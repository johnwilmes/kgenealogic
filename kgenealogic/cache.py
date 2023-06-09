import numpy as np
import scipy as sp
import pandas as pd
import osqp
import sqlalchemy as sql

from kgenealogic.schema import *

SOLVER_PARAMS = dict(
    eps_abs=1e-5,
    eps_rel=1e-5,
    eps_prim_inf=1e-7,
    eps_dual_inf=1e-7,
    scaling=10,
    polish=True,
    max_iter=10000,
    verbose=False,
)
# rounding and numerical errors mean lineqs not exactly solvable
SOLVER_TOLERANCE = 0.2

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

_CACHE_TABLES = [partition, segment_partition, negative]
def clear_cache(engine):
    invalidate_cache(engine)
    delete_imputed = sql.delete(segment).where(segment.c.imputed)
    with engine.connect() as conn:
        conn.execute(delete_imputed)
        conn.commit()

    for t in _CACHE_TABLES:
        t.drop(engine)
        t.create(engine)

def build_partition(segments):
    partition = (
        pd.concat([
            segments[['chromosome', 'start']].rename(columns=dict(start='end')),
            segments[['chromosome', 'end']],
        ])
        .drop_duplicates()
        .sort_values(['chromosome', 'end'])
    )
    partition['start'] = (
        partition
        .groupby('chromosome')
        ['end']
        .shift(1)
    )
    partition = partition.dropna().reset_index(drop=True)
    partition['start'] = partition.start.astype(int)
    partition['id'] = partition.index+1

    return partition

def build_cache(engine, sources=None):
    clear_cache(engine)
    # get segments
    # build partition

    select_segs = sql.select((segment.c.id-1).label("seg_id"), segment.c["chromosome", "start", "end", "length"])
    with engine.connect() as conn:
        all_segs = pd.read_sql(select_segs, conn)
    cells = build_partition(all_segs)
    
    select_seg_part = (
        sql.select(segment.c.id, partition.c.id)
        .join_from(segment, partition, sql.and_(
            segment.c.chromosome==partition.c.chromosome,
            segment.c.start<=partition.c.start,
            segment.c.end>=partition.c.end,
        ))
    )

    populate_segment_partition = (
        segment_partition.insert().from_select(["segment", "partition"], select_seg_part)
    )
    with engine.connect() as conn:
        cells.to_sql("partition", conn, if_exists="append", index=False, chunksize=25000)
        conn.execute(populate_segment_partition)
        conn.commit()

    select_seg_part_size = sql.select(
        (segment_partition.c.segment-1).label("seg_id"),
        (segment_partition.c.partition-1).label("part_id"),
        partition.c.mbp,
    ).join_from(segment_partition, partition)

    chrom_size = (
        sql.select(partition.c.chromosome, sql.func.sum(partition.c.mbp).label("mbp"))
        .group_by(partition.c.chromosome)
        .cte()
    )
    px = partition.alias()
    py = partition.alias()

    select_objective = (
        sql.select(
            (px.c.id-1).label("id_x"),
            (py.c.id-1).label("id_y"),
            sql.case(
                (px.c.id==py.c.id, (px.c.mbp/100)*(chrom_size.c.mbp-px.c.mbp)),
                 else_=(px.c.mbp/100)*py.c.mbp
            ).label("value")
        )
        .join_from(px, py, px.c.chromosome==py.c.chromosome)
        .join(chrom_size, px.c.chromosome==chrom_size.c.chromosome)
    )

    with engine.connect() as conn:
        objective = pd.read_sql(select_objective, conn)
        seg_parts = pd.read_sql(select_seg_part_size, conn)

    # minimize (1/2)x^T P x + q^T x
    # subject to l<= Ax <=u
    # x_i: cM/(fraction of chromosome) for cell i
    n = len(cells)

    # objective: minimize basepair-length weighted variance of x within each chromosome
    P = sp.sparse.csc_matrix((objective.value, (objective.id_x, objective.id_y)))
    q = np.zeros(n)

    # subject to matching overall cM length of known segments, x nonnegative
    A = sp.sparse.vstack([
        sp.sparse.csc_matrix((seg_parts.mbp, (seg_parts.seg_id, seg_parts.part_id))),
        sp.sparse.eye(n, format='csc')
    ], format='csc')

    seg_lengths = np.zeros(all_segs.seg_id.max()+1)
    seg_lengths[all_segs.seg_id] = all_segs.length

    l = np.concatenate([seg_lengths-SOLVER_TOLERANCE, np.zeros(n)])
    u = np.concatenate([seg_lengths+SOLVER_TOLERANCE, np.full((n,), np.Inf)])

    solver = osqp.OSQP()
    solver.setup(P=P, q=q, A=A, l=l, u=u, **SOLVER_PARAMS)
    results = solver.solve()

    cells['rate'] = results.x
    cells['length'] = 1e-6*(cells.end-cells.start)*cells.rate

    set_partition_length = (
        partition.update()
        .where(partition.c.id==sql.bindparam("b_id"))
        .values(length=sql.bindparam("length"))
    )
    with engine.connect() as conn:
        conn.execute(set_partition_length, cells[['id', 'length']].rename(columns=dict(id="b_id")).to_dict(orient="records"))
        conn.commit()

    if sources is not None:
        values = sql.values(sql.column('kit', sql.Integer)).data([(x,) for x in sources])
        source_table = sql.select(values).cte()
        label_sources = (
            sql.update(source)
            .where(source.c.kit==source_table.c.kit)
            .values(has_negative=True)
        )
    else:
        source_table = source
        label_sources = sql.update(source).values(has_negative=True)

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
        .join_from(source_table, m1, m1.c.kit1==source_table.c.kit)
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
    select_imputed_segments = (
        sql.select(
            sql.bindparam("chromosome"),
            sql.bindparam("start"),
            sql.bindparam("end"),
            sql.func.sum(partition.c.length),
            sql.literal(True),
        ).where(partition.c.chromosome==sql.bindparam("chromosome"))
        .where(partition.c.start>=sql.bindparam("start"))
        .where(partition.c.start<=sql.bindparam("end"))
    )

    insert_imputed = (
        segment
        .insert()
        .from_select(["chromosome", "start", "end", "length", "imputed"], select_imputed_segments)
    )
    select_seg_ids = sql.select(segment.c["id", "chromosome", "start", "end"])
    with engine.connect() as conn:
        conn.execute(insert_imputed, overlap_segments.to_dict(orient="records"))
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
        conn.execute(insert_imputed, neg_segments.to_dict(orient="records"))
        conn.commit()
        seg_ids = pd.read_sql(select_seg_ids, conn).rename(columns=dict(id="neg_segment"))

    neg_tri = neg_tri.merge(seg_ids, on=['chromosome', 'start', 'end'])
    neg_tri = neg_tri[['source', 'target1', 'target2', 'segment1', 'segment2', 'overlap_segment', 'neg_segment']]
    validate_cache = (
        sql.update(kgenealogic)
        .where(kgenealogic.c.key=='cache_valid')
        .values(value='Y')
    )

    with engine.connect() as conn:
        neg_tri.to_sql("negative", conn, if_exists="append", index=False, chunksize=25000)
        conn.execute(label_sources)
        conn.execute(validate_cache)
        conn.commit()
