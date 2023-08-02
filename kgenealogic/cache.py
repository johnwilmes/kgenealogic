import pandas as pd
import sqlalchemy as sql
import importlib.resources as resources

from kgenealogic.schema import *
import kgenealogic as kg

GENETMAP_PATH='genetmap.csv'
GENETMAP_DTYPE={'chromosome': str, 'position': int, 'cm': float}

def invalidate_cache(engine):
    stmt = (
        sql.update(kgenealogic)
        .where(kgenealogic.c.key=='cache_valid')
        .values(value='N')
    )
    label_sources = sql.update(source).values(has_negative=False)
    with engine.connect() as conn:
        conn.execute(stmt)
        conn.execute(label_sources)
        conn.commit()

def is_cache_valid(engine):
    stmt = (
        sql.select(kgenealogic.c.value=='Y')
        .where(kgenealogic.c.key=='cache_valid')
    )
    with engine.connect() as conn:
        result = conn.execute(stmt)
    return result.scalar_one()

_CACHE_TABLES = [negative, overlap]
def clear_cache(engine):
    invalidate_cache(engine)

    for t in _CACHE_TABLES:
        t.drop(engine)
        t.create(engine)

    drop_gen_segs = sql.delete(segment).where(segment.c.generated)
    with engine.connect() as conn:
        conn.execute(drop_gen_segs)
        conn.commit()

def build_cache(engine):
    from kgenealogic.data import as_internal_segment
    clear_cache(engine)

    valid_neg_targets = (
        sql.select(
            source.c.kit.label("source"),
            triangle.c.kit2.label("target"),
        )
        .join_from(source, triangle, source.c.kit==triangle.c.kit1)
        .where(source.c.has_triangles)
        .distinct()
        .cte()
    )
    t1, t2 = valid_neg_targets.alias(), valid_neg_targets.alias()
    m1, m2 = match.alias(), match.alias()
    s1, s2 = segment.alias(), segment.alias()
    select_overlap = (
        sql.select(
            t1.c.source,
            t1.c.target.label("target1"),
            t2.c.target.label("target2"),
            s1.c.chromosome,
            sql.func.max(s1.c.start, s2.c.start).label("start"),
            sql.func.min(s1.c.end, s2.c.end).label("end"),
        )
        .join_from(t1, t2, sql.and_(
            t1.c.source==t2.c.source,
            t1.c.target!=t2.c.target,
        ))
        .join(m1, sql.and_(
            t1.c.source==m1.c.kit1,
            t1.c.target==m1.c.kit2,
        ))
        .join(m2, sql.and_(
            t2.c.source==m2.c.kit1,
            t2.c.target==m2.c.kit2,
        ))
        .join(s1, m1.c.segment==s1.c.id)
        .join(s2, sql.and_(
            m2.c.segment==s2.c.id,
            s1.c.chromosome==s2.c.chromosome,
            s1.c.start < s2.c.end,
            s2.c.start < s1.c.end,
        ))
    )
    with engine.connect() as conn:
        match_overlap = pd.read_sql(select_overlap, conn)
    match_overlap = as_internal_segment(engine, match_overlap, generated=True)
    match_overlap = match_overlap[["source", "target1", "target2", "segment"]]
    with engine.connect() as conn:
        conn.execute(overlap.insert(), match_overlap.to_dict(orient="records"))
        conn.commit()
        
    select_neg_overlap = (
        sql.select(
            overlap.c.id.label("overlap"),
            s1.c.chromosome,
            s1.c.start,
            s1.c.end,
            s2.c.start.label("tri_start"),
            s2.c.end.label("tri_end"),
        )
        .join_from(overlap, s1, overlap.c.segment==s1.c.id)
        .join(triangle, sql.and_(
            overlap.c.source==triangle.c.kit1,
            overlap.c.target1==triangle.c.kit2,
            overlap.c.target2==triangle.c.kit3,
        ), isouter=True)
        .join(s2, triangle.c.segment==s2.c.id)
        .where(sql.or_(
            s2.c.id.is_(None),
            sql.and_(
                s2.c.chromosome==s1.c.chromosome,
                s2.c.start <= s1.c.end,
                s2.c.end >= s1.c.start,
            ),
        ))
    )
    with engine.connect() as conn:
        neg_overlap = pd.read_sql(select_neg_overlap, conn)

    neg_base = (
        neg_overlap[neg_overlap.tri_start.isnull()]
        .drop(columns=["tri_start", "tri_end"])
        .drop_duplicates()
    )
    neg_remain = []
    neg_group_cols = ["overlap", "chromosome"]
    for _, tri_data in neg_overlap.dropna().groupby(neg_group_cols):
        id_vals = tri_data[neg_group_cols].iloc[0].to_dict()
        start, end = tri_data[['start', 'end']].iloc[0]
        # negative triangulations are the raw triangulation (start to end),
        # excluding all the positive parts (tri_start to tri_end)
        for _, row in tri_data.sort_values('tri_start').iterrows():
            if row.tri_start > start:
                neg_remain.append(dict(
                    start=start,
                    end=int(row.tri_start),
                    **id_vals))
            start = int(row.tri_end)
        if end > start:
             neg_remain.append(dict(
                start=start,
                end=end,
                **id_vals))       
    neg_tri = pd.concat([neg_base, pd.DataFrame(neg_remain)], ignore_index=True).drop_duplicates()
    neg_tri = as_internal_segment(engine, neg_tri, generated=True).drop_duplicates()
    neg_tri = neg_tri[["overlap", "segment"]].rename(columns=dict(segment="neg_segment"))
    with engine.connect() as conn:
        conn.execute(negative.insert(), neg_tri.to_dict(orient="records"))
        conn.commit()

    label_sources = sql.update(source).values(has_negative=True)
    with engine.connect() as conn:
        neg_tri.to_sql("negative", conn, if_exists="append", index=False, chunksize=25000)
        conn.execute(label_sources)
        conn.commit()

    # compute cm lengths for all segments
    null_len_segs_query = sql.select(segment).where(segment.c.length.is_(None))
    with engine.connect() as conn:
        null_len_segs = pd.read_sql(null_len_segs_query, conn)

    genetmap = pd.read_csv(resources.files(kg).joinpath(GENETMAP_PATH), dtype=GENETMAP_DTYPE)
    genetmap = genetmap.sort_values("position")

    seg_len_start = pd.merge_asof(
        null_len_segs.sort_values("start"),
        genetmap,
        left_on="start",
        right_on="position",
        by="chromosome",
        direction="nearest"
    ).set_index("id").cm

    seg_len_end = pd.merge_asof(
        null_len_segs.sort_values("end"),
        genetmap,
        left_on="end",
        right_on="position",
        by="chromosome",
        direction="nearest"
    ).set_index("id").cm

    null_len_segs['length'] = null_len_segs.id.map(seg_len_end-seg_len_start)
    null_len_segs = null_len_segs[['id', 'length']].rename(columns=dict(id="b_id"))

    seg_len_update = (
        segment.update()
        .where(segment.c.id==sql.bindparam("b_id"))
        .values(length=sql.bindparam("length"))
    )
    with engine.connect() as conn:
        conn.execute(seg_len_update, null_len_segs.to_dict(orient="records"))
        conn.commit()

    # all done, set cache as valid
    validate_cache = (
        sql.update(kgenealogic)
        .where(kgenealogic.c.key=='cache_valid')
        .values(value='Y')
    )
    with engine.connect() as conn:
        conn.execute(validate_cache)
        conn.commit()
