import numpy as np
import scipy as sp
import pandas as pd
import osqp
import importlib.resources

from kgenealogic.data import get_query

SOLVER_PARAMS = dict(
    eps_abs=1e-5,
    eps_rel=1e-5,
    eps_prim_inf=1e-7,
    eps_dual_inf=1e-7,
    scaling=10,
    polish=True,
    max_iter=10000,
)
# rounding and numerical errors mean lineqs not exactly solvable
SOLVER_TOLERANCE = 0.2

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

def estimate_partition_lengths():
    pass


def build_cache(db):
    # get segments
    # build partition

    select_segs = "SELECT id-1 AS seg_id, chromosome, start, end, length FROM segment"
    segments = pd.read_sql(select_segs, db, dtype=dict(seg_id=int, chromosome=str, start=int, end=int, length=float))
    partition = build_partition(segments)

    # reset cache tables of db
    # insert partition into db
    # partition seg/partition relationship
    insert_partition = """
    INSERT INTO partition (id, chromosome, start, end)
    VALUES (:id, :chromosome, :start, :end)
    """
    populate_segment_partition = """
        INSERT INTO segment_partition
        SELECT segment.id, partition.id FROM
        segment JOIN partition
        ON segment.chromosome=partition.chromosome
        AND segment.start<=partition.start
        AND segment.end>=partition.end
    """
    with db:
        db.executescript(get_query("init_cache"))
        db.executemany(insert_partition, partition.to_dict(orient='records'))
        db.execute(populate_segment_partition)

    # solve optimization problem for partition sizes
    select_seg_parts = """
    SELECT segment_partition.segment-1 AS seg_id, segment_partition.partition-1 AS part_id, partition.mbp
    FROM segment_partition JOIN partition on segment_partition.partition=partition.id
    """
    seg_parts = pd.read_sql(select_seg_parts, db)

    select_objective = """
    WITH chrom_size AS (
        SELECT chromosome, SUM(end-start) AS mbp FROM partition GROUP BY chromosome
    )

    SELECT px.id-1 AS id_x, py.id-1 AS id_y,
        CASE WHEN px.id=py.id THEN px.mbp*(chrom_size.mbp-px.mbp)/100 ELSE px.mbp*py.mbp/100 END AS value
    FROM partition AS px
    JOIN partition AS py ON px.chromosome=py.chromosome
    JOIN chrom_size ON px.chromosome=chrom_size.chromosome
    """
    objective = pd.read_sql(select_objective, db)

    # minimize (1/2)x^T P x + q^T x
    # subject to l<= Ax <=u
    # x_i: cM/(fraction of chromosome) for cell i
    n = len(partition)

    # objective: minimize basepair-length weighted variance of x within each chromosome
    P = sp.sparse.csc_matrix((objective.value, (objective.id_x, objective.id_y)))
    q = np.zeros(n)

    # subject to matching overall cM length of known segments, x nonnegative
    A = sp.sparse.vstack([
        sp.sparse.csc_matrix((seg_parts.mbp, (seg_parts.seg_id, seg_parts.part_id))),
        sp.sparse.eye(n, format='csc')
    ], format='csc')

    seg_lengths = np.zeros(segments.seg_id.max()+1)
    seg_lengths[segments.seg_id] = segments.length

    l = np.concatenate([seg_lengths-SOLVER_TOLERANCE, np.zeros(n)])
    u = np.concatenate([seg_lengths+SOLVER_TOLERANCE, np.full((n,), np.Inf)])

    solver = osqp.OSQP()
    solver.setup(P=P, q=q, A=A, l=l, u=u, **SOLVER_PARAMS)
    results = solver.solve()

    partition['rate'] = results.x
    partition['length'] = 1e-6*(partition.end-partition.start)*partition.rate
    set_partition_length = """
    UPDATE partition SET length=:length WHERE id=:id
    """
    with db:
        db.executemany(set_partition_length, partition[['id', 'length']].to_dict(orient='records'))

    select_neg_triangles = """
    SELECT m1.kit1 AS source, m1.kit2 AS target1, m2.kit2 AS target2,
        s1.id AS segment1, s2.id AS segment2,
        s1.chromosome, MAX(s1.start, s2.start) AS start, MIN(s1.end, s2.end) AS end,
        st.start AS start_pos, st.end as end_pos
    FROM match AS m1
    JOIN match AS m2 on m1.kit1=m2.kit1
    JOIN segment AS s1 on m1.segment=s1.id
    JOIN segment AS s2 on m2.segment=s2.id
    LEFT OUTER JOIN triangle ON m1.kit1=triangle.kit1 AND m1.kit2=triangle.kit2 AND m2.kit2=triangle.kit3
    LEFT OUTER JOIN segment AS st ON triangle.segment=st.id
    WHERE s1.chromosome=s2.chromosome AND s1.start < s2.end AND s2.start < s1.end
    AND ((st.id IS NULL) OR (
            st.chromosome=s1.chromosome
            AND st.start < MIN(s1.end, s2.end)
            AND st.end > MAX(s1.start, s2.start)
    ))
    """

    neg_candidates = pd.read_sql(select_neg_triangles, db)
    overlap_segments = neg_candidates[['chromosome', 'start', 'end']].drop_duplicates()
    insert_imputed = """
    INSERT INTO segment (chromosome, start, end, length, imputed) SELECT
    :chromosome, :start, :end, SUM(partition.length), TRUE
    FROM partition WHERE partition.chromosome=:chromosome
        AND partition.start >= :start
        AND partition.end <= :end
    """
    with db:
        db.executemany(insert_imputed, overlap_segments.to_dict(orient='records'))

    seg_ids = pd.read_sql("SELECT id AS overlap_segment, chromosome, start, end FROM segment", db)
    neg_candidates = neg_candidates.merge(seg_ids, on=['chromosome', 'start', 'end'])
    neg_group_cols=['source', 'target1', 'target2', 'overlap_segment']
    neg_candidates = neg_candidates.drop_duplicates(subset=neg_group_cols)

    neg_base = (
            neg_candidates[neg_candidates.start_pos.isnull()]
            .drop(columns=["start_pos", "end_pos"])
    )
    neg_remain = []
    neg_id_cols = neg_group_cols + ['segment1', 'segment2']
    for _, tri_data in neg_candidates.dropna().groupby(neg_group_cols):
        id_vals = neg_candidates[neg_id_cols].to_dict(orient='records')[0]
        start, end = neg_candidates[['start', 'end']]
        # negative triangulations are the raw triangulation (start to end),
        # excluding all the positive parts (start_pos to end_pos)
        for _, row in tri_data.sort_values('start_pos').iterrows():
            if row.start_pos > start:
                neg_tri_remain.append(dict(
                    start=start,
                    end=int(row.start_pos),
                    **id_vals))
            start = int(row.end_pos)
        if end > start:
             neg_tri_remain.append(dict(
                start=start,
                end=end,
                **id_vals))       
    neg_tri = pd.concat([neg_base, pd.DataFrame(neg_remain)]).drop_duplicates()

    neg_segments = neg_tri[['chromosome', 'start', 'end']].drop_duplicates()
    with db:
        db.executemany(insert_imputed, neg_segments.to_dict(orient='records'))
    seg_ids = pd.read_sql("SELECT id AS neg_segment, chromosome, start, end FROM segment", db)
    neg_tri = neg_tri.merge(seg_ids, on=['chromosome', 'start', 'end'])
    neg_tri = neg_tri[['source', 'target1', 'target2', 'segment1', 'segment2', 'overlap_segment', 'neg_segment']]
    insert_neg = """
        INSERT INTO negative
        VALUES (:source, :target1, :target2, :segment1, :segment2, :overlap_segment, :neg_segment)
    """
    with db:
        db.executemany(insert_neg, neg_tri)
        db.execute("UPDATE kgenealogic SET value=TRUE WHERE key='cache_valid'")
