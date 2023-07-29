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

_CACHE_TABLES = [negative]
def clear_cache(engine):
    invalidate_cache(engine)

    for t in _CACHE_TABLES:
        t.drop(engine)
        t.create(engine)

def build_cache(engine):
    clear_cache(engine)

    m1, m2 = match.alias(), match.alias()
    s1, s2, st = segment.alias(), segment.alias(), segment.alias()
    select_neg_triangles = (
        sql.select(
            m1.c.kit1.label("source"),
            m1.c.kit2.label("target1"),
            m2.c.kit2.label("target2"),
            s1.c.id.label("segment1"),
            s2.c.id.label("segment2"),
            s1.c.chromosome,
            sql.func.max(s1.c.start, s2.c.start).label("start"),
            sql.func.min(s1.c.end, s2.c.end).label("end"),
            st.c.start.label("start_pos"),
            st.c.end.label("end_pos"),
        )
        .join_from(source, m1, m1.c.kit1==source.c.kit)
        .join(m2, m1.c.kit1==m2.c.kit1)
        .join(s1, m1.c.segment==s1.c.id)
        .join(s2, m2.c.segment==s2.c.id)
        .join(triangle, sql.and_(
            m1.c.kit1==triangle.c.kit1,
            m1.c.kit2==triangle.c.kit2,
            m2.c.kit2==triangle.c.kit3,
        ), isouter=True)
        .join(st, triangle.c.segment==st.c.id, isouter=True)
        .where(s1.c.chromosome==s2.c.chromosome)
        .where(s1.c.start < s2.c.end)
        .where(s2.c.start < s1.c.end)
        .where(sql.or_(
            st.c.id.is_(None),
            sql.and_(
                st.c.chromosome==s1.c.chromosome,
                st.c.start < sql.func.min(s1.c.end, s2.c.end),
                st.c.end > sql.func.max(s1.c.start, s2.c.start),
            ),
        ))
    )
    with engine.connect() as conn:
        neg_candidates = pd.read_sql(select_neg_triangles, conn)
    overlap_segments = neg_candidates[['chromosome', 'start', 'end']].drop_duplicates()

    select_seg_ids = sql.select(segment.c["id", "chromosome", "start", "end"])
    with engine.connect() as conn:
        conn.execute(segment.insert(), overlap_segments.to_dict(orient="records"))
        conn.commit()
        seg_ids = pd.read_sql(select_seg_ids, conn).rename(columns=dict(id="overlap_segment"))

    neg_candidates = neg_candidates.merge(seg_ids, on=["chromosome", "start", "end"])
    neg_group_cols=['source', 'target1', 'target2', 'overlap_segment']
    neg_candidates = neg_candidates.drop_duplicates(subset=neg_group_cols)

    neg_base = (
            neg_candidates[neg_candidates.start_pos.isnull()]
            .drop(columns=["start_pos", "end_pos"])
    )
    neg_remain = []
    neg_id_cols = neg_group_cols + ['segment1', 'segment2', 'chromosome']
    for _, tri_data in neg_candidates.dropna().groupby(neg_group_cols):
        id_vals = tri_data[neg_id_cols].iloc[0].to_dict()
        start, end = tri_data[['start', 'end']].iloc[0]
        # negative triangulations are the raw triangulation (start to end),
        # excluding all the positive parts (start_pos to end_pos)
        for _, row in tri_data.sort_values('start_pos').iterrows():
            if row.start_pos > start:
                neg_remain.append(dict(
                    start=start,
                    end=int(row.start_pos),
                    **id_vals))
            start = int(row.end_pos)
        if end > start:
             neg_remain.append(dict(
                start=start,
                end=end,
                **id_vals))       
    neg_tri = pd.concat([neg_base, pd.DataFrame(neg_remain)]).drop_duplicates()

    neg_segments = neg_tri[['chromosome', 'start', 'end']].drop_duplicates()
    with engine.connect() as conn:
        conn.execute(segment.insert(), neg_segments.to_dict(orient="records"))
        conn.commit()
        seg_ids = pd.read_sql(select_seg_ids, conn).rename(columns=dict(id="neg_segment"))

    neg_tri = neg_tri.merge(seg_ids, on=['chromosome', 'start', 'end'])
    neg_tri = neg_tri[['source', 'target1', 'target2', 'segment1', 'segment2', 'overlap_segment', 'neg_segment']]

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
