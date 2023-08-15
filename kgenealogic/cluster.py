import numpy as np
import scipy as sp
import pandas as pd
import sqlalchemy as sql
from dataclasses import dataclass, field
import typing
import typer

from kgenealogic.schema import *

# factor for weighting pairwise matches in the graph; triangles have weight 1
PAIRWISE_FACTOR = 0.25

@dataclass
class Seed:
    kit: int
    floating: bool = field(default_factory=bool)
    negative: bool = field(default_factory=bool)

@dataclass
class SeedTree:
    """Data structure for organizing cluster seeds in tree"""
    ahnentafel: int
    values: list[Seed] = field(default_factory=list)
    branches: dict[str, typing.Any] = field(default_factory=dict)

def flatten(tree, seeds=None):
    """Flatten a SeedTree by putting internal kit ids into values list."""
    if seeds is None:
        seeds = {}
    for b in tree.branches.values():
        flatten(b, seeds)
    seeds.update({v.kit: tree.ahnentafel for v in tree.values})
    return seeds

def cluster_data(engine, config):
    """Cluster data according to the config file

    The behavior is documented in the CLI cluster command.

    Args:
        engine: the sqlalchemy engine for the project database
        config: a ClusterConfig object describing the clustering to perform

    Returns:
        a pandas dataframe classifying kits, as described in the CLI cluster command documentation
    """
    with engine.connect() as conn:
        kits = pd.read_sql(sql.select(kit.c['id', 'kitid']), conn)
    kits = kits.set_index('kitid').id

    seeds = set()
    exclude = set()
    for kitid_list, kit_set in ((config.seeds, seeds), (config.exclude, exclude)):
        for kitid in kitid_list:
            if kitid not in kits.index:
                typer.echo(f"Invalid kit id {kitid}", err=True)
                raise typer.Exit(code=1)
            else:
                kit_set.add(int(kits.loc[kitid]))

    if config.include:
        include = set()
        match_include_query = (
            sql.select(sql.distinct(match.c.kit2).label("kit"))
            .join_from(match, segment)
            .where(match.c.kit1==sql.bindparam("kit"))
            .where(segment.c.length > sql.bindparam("length"))
        )
        tri_include_query = (
            sql.select(sql.distinct(triangle.c.kit2).label("kit"))
            .join_from(triangle, segment)
            .where(triangle.c.kit1==sql.bindparam("kit"))
            .where(segment.c.length > sql.bindparam("length"))
        )
        with engine.connect() as conn:
            for k in config.include:
                kitid = k['id']
                if kitid not in kits.index:
                    typer.echo(f"Invliad kit id {kitid}", err=True)
                    raise typer.Exit(code=1)
                else:
                    include.add(int(kits.loc[kitid]))
                kitloc = int(kits.loc[kitid])
                if k.get('matches', None) is not None:
                    query_params = dict(kit=kitloc, length=k['matches'])
                    k_nbr = pd.read_sql(match_include_query, conn, params=query_params)
                    include.update([int(x) for x in k_nbr.kit])
                if k.get('triangles', None) is not None:
                    query_params = dict(kit=kitloc, length=k['triangles'])
                    k_nbr = pd.read_sql(tri_include_query, conn, params=query_params)
                    include.update([int(x) for x in k_nbr.kit])
        kits = kits[kits.isin(include|seeds)]

    kits = kits[~kits.isin(exclude)]

    graph_query = (
        sql.select(
            match.c.kit1,
            match.c.kit2,
            sql.func.sum(segment.c.length).label("weight"),
        )
        .join_from(match, segment, match.c.segment==segment.c.id)
        .where(segment.c.length >= config.min_length)
        .where(match.c.kit1 != match.c.kit2)
        .group_by(match.c['kit1', 'kit2'])
    )
    triangle_query = (
        sql.select(
            triangle.c.kit1,
            triangle.c.kit2,
            sql.func.sum(segment.c.length).label("weight"),
        )
        .join_from(triangle, segment, triangle.c.segment==segment.c.id)
        .where(segment.c.length >= config.min_length)
        .where(triangle.c.kit1 != triangle.c.kit2)
        .where(triangle.c.kit3.not_in(config.exclude))
        .group_by(triangle.c['kit1', 'kit2'])
    )
    trisource_query = (
        sql.select(source.c.kit)
        .where(source.c.triangle.is_not(None))
    )

    with engine.connect() as conn:
        graph = pd.read_sql(graph_query, conn)
        tri_graph = pd.read_sql(triangle_query, conn)
        trisource= set(int(x) for x in pd.read_sql(trisource_query, conn).kit)

    graph['weight'] = PAIRWISE_FACTOR*graph['weight'].astype(float)
    graph = (
        graph
        .set_index(['kit1', 'kit2'])
        .add(tri_graph.set_index(['kit1', 'kit2']), fill_value=0)
        .reset_index()
    )

    autox_query = (
        sql.select(sql.distinct(match.c.kit2).label('kit'))
        .join_from(match, segment, match.c.segment==segment.c.id)
        .where(segment.c.chromosome=='X')
        .where(segment.c.length >= config.min_length)
        .where(match.c.kit1==sql.bindparam('seed'))
    )
    def build_seed_tree(raw, ahnentafel):
        parsed = SeedTree(ahnentafel=ahnentafel)
        autox = set()
        for k in raw.get('kits', []):
            kitlocal = int(kits.loc[k['id']])
            floating = k.get('float')
            if floating is None:
                floating = kitlocal not in trisource
            if k.get('autox'):
                with engine.connect() as conn:
                    k_autox = pd.read_sql(autox_query, conn, params=dict(seed=kitlocal))
                k_autox = set(int(x) for x in k_autox.kit)
                k_autox.difference_update(exclude)
                k_autox.difference_update(seeds)
                seeds.update(k_autox)
                autox.update(k_autox)
            negative = k.get('negative', False)
            parsed.values.append(Seed(kit=kitlocal, floating=floating, negative=negative))

        for a, short, long in ((1, 'M', 'maternal'), (0, 'P', 'paternal')):
            if long in raw:
                parsed.branches[short] = build_seed_tree(raw[long],a+(2*ahnentafel))

        if autox:
            if 'M' not in parsed.branches:
                parsed.branches['M'] = SeedTree(values=[],ahnentafel=1+(2*ahnentafel))
            maternal_values = parsed.branches['M'].values
            for x in autox:
                maternal_values.append(Seed(kit=x, floating=True))
        return parsed
    tree = build_seed_tree(config.tree, 1)

    seed_labels = flatten(tree)

    source_neg = lambda s_kit: get_negative(engine, s_kit, config.min_length)
    clusters = recursive_cluster(kits, tree, graph, source_neg)
    
    clusters['seed'] = clusters.kit.map(seed_labels).astype('Int64')
    kitids = kits.reset_index().set_index('id').kitid
    clusters['kit'] = clusters.kit.map(kitids)

    fixed_cols = ['kit', 'ahnentafel', 'seed']
    clusters = clusters[fixed_cols + [c for c in clusters.columns if c not in fixed_cols]]

    return clusters

def recursive_cluster(kits, tree, graph, source_neg):
    """Recursively partition a list of kits using a given tree structure.

    Args:
        kits: pandas Series of internal kit ids to be be included in the segmentation
        tree: a SeedTree describing the structure of the segmentation and the seeds to be used
        graph: the a pandas DataFrame describing the graph of weighted relationships to be used to
            construct the segmentation, in sparse COO format. Columns 'kit1' and 'kit2' give internal
            kit numbers, and column 'weight' gives the weight of the edge between them. Symmetric,
            so 'kit1' and 'kit2' appear in both orders, if there is an edge between them.
        source_neg: a function mapping internal kit numbers to negative triangles in the same
            format as `graph`

    Returns:
        pandas DataFrame in the same format as the result of cluster_data, but with internal kit
        numbers instead of external kit id strings
    """
    nonfloat = []
    tri_graph = graph
    for s_kit in tree.values:
        if s_kit.negative:
            tri_graph = tri_graph.set_index(['kit1', 'kit2']).add(
                source_neg(s_kit.kit).set_index(['kit1', 'kit2']), fill_value=0
            ).reset_index()
        if not s_kit.floating:
            nonfloat.append(s_kit.kit)
        
    depth = int(np.log2(tree.ahnentafel))
    result = pd.DataFrame(dict(kit=kits))

    kits = kits[~kits.isin(nonfloat)]
    graph = graph[graph.kit1.isin(kits)&graph.kit2.isin(kits)]
    tri_graph = tri_graph[tri_graph.kit1.isin(kits)&tri_graph.kit2.isin(kits)]
    if tree.branches and (len(result)>0) and (len(tri_graph) > 0):
        flat_seeds = []
        for label, tree_branch in tree.branches.items():
            branch_seeds = flatten(tree_branch).keys()
            branch_seeds = pd.DataFrame(dict(kit=branch_seeds))
            branch_seeds['label'] = label
            flat_seeds.append(branch_seeds)
        flat_seeds = pd.concat(flat_seeds, ignore_index=True)

        clusters = get_clusters(tri_graph, flat_seeds).set_index('kit')
        result['label'] = (
            result.kit.map(clusters.label)
            .fillna('')
            .astype(str)
        )
        result['confidence'] = result.kit.map(clusters.confidence)
        result['ahnentafel'] = np.select(
            [result.label=='P', result.label=='M'],
            [2*tree.ahnentafel, 1+(2*tree.ahnentafel)], tree.ahnentafel
        )

        result_branches = []
        for label, kit_branch in result.groupby('label'):
            if label in tree.branches:
                tree_branch = tree.branches[label]
                branch = recursive_cluster(kit_branch.kit, tree_branch, graph, source_neg)
                branch[f"label{depth}"] = branch.kit.map(kit_branch.set_index('kit').label)
                branch[f"confidence{depth}"] = branch.kit.map(kit_branch.set_index('kit').confidence)
            else:
                branch = kit_branch.rename(columns=dict(
                    label=f"label{depth}",
                    confidence=f"confidence{depth}",
                ))
            result_branches.append(branch)
        result = pd.concat(result_branches, ignore_index=True)
    else:
        flat_ahnentafel = flatten(tree)
        result['ahnentafel'] = result.kit.map(flat_ahnentafel).fillna(tree.ahnentafel)

    return result

def get_clusters(graph, seeds, max_rounds=None, fix_seeds=True):
    """Greedily segment each component of graph to respect the labels of seeds."""
    graph = graph.copy()
    labels = (
        graph
        .assign(weight=graph.weight.abs())
        .groupby('kit1')
        .weight
        .sum()
        .fillna(0)
        .reset_index()
        .rename(columns=dict(kit1='kit'))
    )
    labels['label'] = labels.kit.map(seeds.set_index('kit').label)
    labels = labels[labels.weight>0].set_index('kit').sort_index()
    
    max_rounds = max_rounds or 2*len(labels)
    i=0
    while True:
        graph['label'] = graph.kit2.map(labels.label)
        
        # TODO would be nice if seg_sum took confidence-weighted sum, but this is circular. does it converge?
        seg_sum = graph.groupby(['kit1', 'label']).weight.sum().unstack(level=-1).fillna(0)
        for l in ['P', 'M']:
            if l not in seg_sum:
                seg_sum[l] = 0
        labels['paternal'] = seg_sum.P-seg_sum.M
        labels['confidence'] = labels['paternal'].abs()/labels['weight']
        available = np.any([
            labels.label.isnull(),
            (labels.label=='P')&(labels.paternal<0),
            (labels.label=='M')&(labels.paternal>0),
        ], axis=0)
        available = available & (labels.confidence > 0)
        if fix_seeds:
            available = available & (~labels.index.isin(seeds.kit))
        i+=1
        
        if (not np.any(available)) or (i > max_rounds):
            break

        next_kit = labels[available].iloc[labels.confidence[available].argmax()]
        next_label = 'P' if next_kit.paternal>0 else 'M'
        labels.loc[next_kit.name, 'label'] = next_label

    labels['label'] = labels.label.fillna('').astype(str)
    labels = labels.reset_index()[['kit', 'label', 'confidence']]
    isolated_seeds = seeds[~seeds.kit.isin(labels.kit)]
    return pd.concat([labels, isolated_seeds], ignore_index=True)

def get_negative(engine, s_kit, min_length):
    """Get all negative triangles for a kit s_kit of length at least min_length."""
    build_negative(engine, s_kit)

    neg_tri = (
        sql.select(
            overlap.c.target1.label("kit1"),
            overlap.c.target2.label("kit2"),
            sql.func.sum(-segment.c.length).label("weight"),
        )
        .join_from(negative, overlap)
        .join(segment, negative.c.neg_segment==segment.c.id)
        .where(overlap.c.source==s_kit)
        .where(segment.c.length >= min_length)
        .group_by(overlap.c['target1', 'target2'])
    )

    with engine.connect() as conn:
        result = pd.read_sql(neg_tri, conn)
    return result

def build_negative(engine, s_kit):
    """Check if negative triangles are cached for source s_kit, and build if necessary."""
    with engine.connect() as conn:
        status = conn.execute(sql.select(source).where(source.c.kit==s_kit)).mappings().one_or_none()
    if (not status) or (not status['match']) or (not status['triangle']):
        # missing triangles or matches
        return False

    batch = max(status['match'], status['triangle'])
    neg_batch = status['negative'] or 0
    if neg_batch >= batch:
        # up to date
        return True

    update_source = sql.update(source).where(source.c.kit==s_kit).values(negative=batch)
    from kgenealogic.data import as_internal_segment

    valid_neg_targets = (
        sql.select(
            triangle.c.kit2.label("target"),
        )
        .where(triangle.c.kit1==s_kit)
        .distinct()
        .cte()
    )

    valid_matches = (
        sql.select(
            valid_neg_targets.c.target,
            segment.c.chromosome,
            segment.c.start,
            segment.c.end,
        )
        .join_from(valid_neg_targets, match, sql.and_(
            match.c.kit1==s_kit, match.c.kit2==valid_neg_targets.c.target
        ))
        .join(segment, match.c.segment==segment.c.id)
        .cte()
    )

    m1, m2 = valid_matches.alias(), valid_matches.alias()

    select_overlap = (
        sql.select(
            m1.c.target.label("target1"),
            m2.c.target.label("target2"),
            m1.c.chromosome,
            sql.func.max(m1.c.start, m2.c.start).label("start"),
            sql.func.min(m1.c.end, m2.c.end).label("end"),
        )
        .join_from(m1, m2, m1.c.chromosome==m2.c.chromosome)
        .where(m1.c.target!=m2.c.target)
        .where(m1.c.start < m2.c.end)
        .where(m2.c.start < m1.c.end)
    )

    with engine.connect() as conn:
        conn.execute(sql.delete(overlap).where(overlap.c.source==s_kit))
        match_overlap = pd.read_sql(select_overlap, conn)
    if len(match_overlap)==0:
        return True
    match_overlap = as_internal_segment(engine, match_overlap)

    match_overlap = match_overlap[["target1", "target2", "segment"]]
    match_overlap["source"] = s_kit
    with engine.connect() as conn:
        conn.execute(overlap.insert(), match_overlap.to_dict(orient="records"))
        conn.commit()
        
    seg_over, seg_tri = segment.alias(), segment.alias()
    select_neg_overlap = (
        sql.select(
            overlap.c.id.label("overlap"),
            seg_over.c.chromosome,
            seg_over.c.start,
            seg_over.c.end,
            seg_tri.c.start.label("tri_start"),
            seg_tri.c.end.label("tri_end"),
        )
        .join_from(overlap, seg_over, overlap.c.segment==seg_over.c.id)
        .join(triangle, sql.and_(
            overlap.c.source==triangle.c.kit1,
            overlap.c.target1==triangle.c.kit2,
            overlap.c.target2==triangle.c.kit3,
        ), isouter=True)
        .join(seg_tri, triangle.c.segment==seg_tri.c.id, isouter=True)
        .where(sql.or_(
            seg_tri.c.id.is_(None),
            sql.and_(
                seg_tri.c.chromosome==seg_over.c.chromosome,
                seg_tri.c.start <= seg_over.c.end,
                seg_tri.c.end >= seg_over.c.start,
            ),
        ))
        .where(overlap.c.source==s_kit)
    )
    with engine.connect() as conn:
        neg_overlap = pd.read_sql(select_neg_overlap, conn)

    if len(neg_overlap)==0:
        return True
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
    neg_tri = as_internal_segment(engine, neg_tri).drop_duplicates()
    neg_tri = neg_tri[["overlap", "segment"]].rename(columns=dict(segment="neg_segment"))
    if len(neg_tri)==0:
        return True
    with engine.connect() as conn:
        conn.execute(negative.insert(), neg_tri.to_dict(orient="records"))
        conn.execute(update_source)
        conn.commit()

    return True
