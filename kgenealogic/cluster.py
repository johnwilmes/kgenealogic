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

def flatten(tree, values):
    for b in tree.branches.values():
        flatten(b, values)
    values.extend([v.kit for v in tree.values])
    return values

def cluster_data(engine, config):
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

    kits = kits[~kits.index.isin(exclude)]

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
        .where(source.c.has_triangles)
        .where(source.c.has_negative)
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
            floating = k['float']
            if floating is None:
                floating = kitlocal not in trisource
            if k['autox']:
                with engine.connect() as conn:
                    k_autox = pd.read_sql(autox_query, conn, params=dict(seed=kitlocal))
                k_autox = set(int(x) for x in k_autox.kit)
                k_autox.difference_update(exclude)
                k_autox.difference_update(seeds)
                seeds.update(k_autox)
                autox.update(k_autox)
            parsed.values.append(Seed(kit=kitlocal, floating=floating, negative=k['negative']))

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

    source_tri = lambda source: get_triangles(engine, source, config.min_length)
    clusters = recursive_cluster(kits, tree, graph, source_tri)

    kitids = kits.reset_index().set_index('id').kitid
    clusters['kit'] = clusters.kit.map(kitids)

    fixed_cols = ['kit', 'ahnentafel']
    clusters = clusters[fixed_cols + [c for c in clusters.columns if c not in fixed_cols]]

    return clusters

def get_triangles(engine, source, min_length):
    neg_tri = (
        sql.select(
            overlap.c.target1.label("kit1"),
            overlap.c.target2.label("kit2"),
            sql.func.sum(-segment.c.length).label("weight"),
        )
        .join_from(negative, overlap)
        .join(segment, negative.c.neg_segment==segment.c.id)
        .where(overlap.c.source==source)
        .where(segment.c.length >= min_length)
        .group_by(overlap.c['target1', 'target2'])
    )

    with engine.connect() as conn:
        result = pd.read_sql(neg_tri, conn)
    return result

def recursive_cluster(kits, tree, graph, source_tri):
    # graph format: kit1/kit2/weight
    # initial value: all pairwise matches, suitably weighted relative to triangles that will be added?
    nonfloat = []
    tri_graph = graph
    for source in tree.values:
        if source.negative:
            tri_graph = tri_graph.set_index(['kit1', 'kit2']).add(
                source_tri(source.kit).set_index(['kit1', 'kit2']), fill_value=0
            ).reset_index()
        if not source.floating:
            nonfloat.append(source.kit)
        
    depth = int(np.log2(tree.ahnentafel))
    result = pd.DataFrame(dict(kit=kits))

    kits = kits[~kits.isin(nonfloat)]
    graph = graph[graph.kit1.isin(kits)&graph.kit2.isin(kits)]
    tri_graph = tri_graph[tri_graph.kit1.isin(kits)&tri_graph.kit2.isin(kits)]
    if tree.branches and (len(result)>0) and (len(tri_graph) > 0):
        flat_seeds = []
        for label, tree_branch in tree.branches.items():
            branch_seeds = flatten(tree_branch, [])
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
        #result.label.mask(result.kit.isin(flat_seeds.kit), flat_seeds.label, inplace=True)
        result['confidence'] = result.kit.map(clusters.confidence)
        result['ahnentafel'] = np.select(
            [result.label=='P', result.label=='M'],
            [2*tree.ahnentafel, 1+(2*tree.ahnentafel)], tree.ahnentafel
        )

        result_branches = []
        for label, kit_branch in result.groupby('label'):
            if label in tree.branches:
                tree_branch = tree.branches[label]
                branch = recursive_cluster(kit_branch.kit, tree_branch, graph, source_tri)
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
        result['ahnentafel'] = tree.ahnentafel

    return result

def get_adjacency(edges, weights):
    vertex_to_id = (
        pd.concat(edges)
        .drop_duplicates()
        .rename('id')
        .reset_index(drop=True)
    )
    n = len(vertex_to_id)
    id_to_vertex = vertex_to_id.reset_index().set_index('id')['index']
    
    graph = sp.sparse.coo_array((weights, (edges[0].map(id_to_vertex), edges[1].map(id_to_vertex))), shape=(n,n)).toarray()
    return graph, vertex_to_id

def get_clusters(graph, seeds, max_rounds=None, fix_seeds=True):
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
    return labels.reset_index()[['kit', 'label', 'confidence']]
