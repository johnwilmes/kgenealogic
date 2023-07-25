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
class SeedTree:
    """Data structure for organizing cluster seeds in tree"""
    values: list[int] = field(default_factory=list)
    branch: dict[str, typing.Any] = field(default_factory=dict)


def warn_invalid_kitid(kitid):
    typer.echo(f"Invalid kit id {kitid}", err=True)

def build_seed_tree(kits, tree, config):
    for short, long in (('M', 'maternal'), ('P', 'paternal')):
        if long in config:
            child_kits = []
            for k in config[long]['kits']:
                kitid = k['id']
                if kitid not in kits.index:
                    warn_invalid_kitid(kitid)
                else:
                    child_kits.append(kits.loc[kitid])
            child = SeedTree(values=child_kits)
            tree.branch[short] = child
            build_seed_tree(kits, child, config[long])

def flatten(tree, values):
    for b in tree.branch.values():
        flatten(b, values)
    values.extend(tree.values)

def cluster_data(engine, config):
    kits_query = sql.select(kit.c['id', 'kitid'])
    graph_query = (
        sql.select(
            match.c.kit1,
            match.c.kit2,
            sql.func.sum(segment.c.length).label("weight"),
        )
        .join_from(match, segment, match.c.segment==segment.c.id)
        .where(segment.c.length >= config.min_length)
        .where(match.c.kit1 != match.c.kit2)
        .group_by(match.c.kit1, match.c.kit2)
    )

    with engine.connect() as conn:
        kits = pd.read_sql(kits_query, conn).set_index('kitid').id
        graph = pd.read_sql(graph_query, conn)
    kitids = kits.reset_index().set_index('id').kitid

    graph['weight'] = PAIRWISE_FACTOR*graph['weight'].astype(float)

    root_kits = []
    for k in config.tree.get('kits', []):
        kitid = k['id']
        if kitid not in kits.index:
            warn_invalid_kitid(kitid)
        else:
            root_kits.append(kits.loc[kitid])

    seeds = SeedTree(values=root_kits)
    build_seed_tree(kits, seeds, config.tree)

    for kitid in config.exclude:
        if kitid not in kits.index:
            warn_invalid_kitid(kitid)
    kits = kits[~kits.index.isin(config.exclude)]

    flat_seeds = []
    flatten(seeds, flat_seeds)
    unique_seeds = set()
    for s in flat_seeds:
        if s in unique_seeds:
            typer.echo("Duplicated seed {}".format(kitids.loc[s]), err=True)
        else:
            unique_seeds.add(s)

    source_tri = lambda source: get_triangles(engine, int(source), config.min_length)
    clusters = recursive_cluster(kits, seeds, graph, source_tri)

    clusters['kit'] = clusters.kit.map(kitids)

    return clusters

def get_triangles(engine, source, min_length):
    pos_tri = (
        sql.select(
            triangle.c.kit1,
            triangle.c.kit2,
            segment.c.length,
        )
        .join_from(triangle, segment, triangle.c.segment==segment.c.id)
        .where(triangle.c.kit3==source)
        .where(segment.c.length >= min_length)
        .where(triangle.c.kit1!=triangle.c.kit2)
    )
    neg_tri = (
        sql.select(
            negative.c.target1.label("kit1"),
            negative.c.target2.label("kit2"),
            (-segment.c.length).label("length"),
        )
        .join_from(negative, segment, negative.c.neg_segment==segment.c.id)
        .where(negative.c.source==source)
        .where(segment.c.length >= min_length)
        .where(negative.c.target1!=negative.c.target2)
    )
    all_tri = sql.union_all(pos_tri, neg_tri).cte()

    sum_tri = (
        sql.select(
            all_tri.c.kit1,
            all_tri.c.kit2,
            sql.func.sum(all_tri.c.length).label("weight"),
        )
        .group_by(all_tri.c.kit1, all_tri.c.kit2)
    )
    with engine.connect() as conn:
        result = pd.read_sql(sum_tri, conn)
    return result

def recursive_cluster(kits, seeds, graph, source_tri):
    # graph format: kit1/kit2/weight
    # initial value: all pairwise matches, suitably weighted relative to triangles that will be added?
    for source in seeds.values:
        graph = graph.set_index(['kit1', 'kit2']).add(
            source_tri(source).set_index(['kit1', 'kit2']), fill_value=0
        ).reset_index()
        
    kits = kits[~kits.isin(seeds.values)]
    graph = graph[graph.kit1.isin(kits)&graph.kit2.isin(kits)]
    
    paternal_seeds = []
    maternal_seeds = []
    if 'P' in seeds.branch:
        flatten(seeds.branch['P'], paternal_seeds)
    if 'M' in seeds.branch:
        flatten(seeds.branch['M'], maternal_seeds)
    
    paternal_seeds = pd.DataFrame(dict(kit=paternal_seeds))
    paternal_seeds['label'] = 'P'
    maternal_seeds = pd.DataFrame(dict(kit=maternal_seeds))
    maternal_seeds['label'] = 'M'
    flat_seeds = pd.concat([paternal_seeds, maternal_seeds], ignore_index=True)

    clusters = get_clusters(graph, flat_seeds)
    
    for label, branch in clusters.groupby('label'):
        if label in seeds.branch:
            branch = recursive_cluster(branch.kit, seeds.branch[label], graph, source_tri)
            if len(branch) > 0:
                branch_labels = (label + clusters.kit.map(branch.set_index('kit').label))
                clusters.label.where(branch_labels.isnull(), branch_labels, inplace=True)
    return clusters

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

def get_clusters(graph, seeds):
    clusters = (
        graph
        .assign(weight=graph.weight.abs())
        .groupby('kit1')
        .weight
        .sum()
        .reset_index()
        .rename(columns=dict(kit1='kit'))
    )
    clusters['label'] = clusters.kit.map(seeds.set_index('kit').label).fillna('')

    adjacency, vertex_to_kit = get_adjacency((graph.kit1, graph.kit2), graph.weight)
    n_components, components = sp.sparse.csgraph.connected_components(adjacency!=0, directed=False)
    for c_idx in range(n_components):
        c = (components==c_idx)
        c_kits = vertex_to_kit[c].reset_index(drop=True)
        c_clusters = clusters.kit.isin(c_kits)
        c_adj = adjacency[c,:][:,c]

        m_seeds = c_kits.isin(seeds.query("label=='M'").kit)
        p_seeds = c_kits.isin(seeds.query("label=='P'").kit)
        n_m_seeds = m_seeds.sum()
        n_p_seeds = p_seeds.sum()
        
        n = len(c_adj)
        # order vertices according to principal eigenvalue
        _, eig = sp.linalg.eigh(c_adj, subset_by_index=(n-1, n-1))
        eig = eig[:,0].copy()
        
        opt = float('inf')
        for eig0 in (eig, -eig):
            eig0[m_seeds] = float('-inf')
            eig0[p_seeds] = float('inf')
            
            vert_order = eig0.argsort()
            # compute size of cut (sum of crossing weights) for each cut of form vert_order[:i]
            triu = np.triu(c_adj[vert_order, :][:, vert_order])
            cut_size = np.zeros((n+1,))
            cut_size[1:] = triu.sum(axis=1).cumsum()-triu.sum(axis=0).cumsum()
            if n_m_seeds > 0:
                cut_size[:n_m_seeds] = float('inf')
            if n_p_seeds > 0:
                cut_size[-n_p_seeds:] = float('inf')
            val = cut_size.min()
            if val < opt:
                opt = val
                n_M = cut_size.argmin()
                M = vert_order[:n_M]
                P = vert_order[n_M:]
            
        clusters['label'].mask(clusters.kit.isin(c_kits[M]), 'M', inplace=True)
        clusters['label'].mask(clusters.kit.isin(c_kits[P]), 'P', inplace=True)
    
    return clusters[['kit', 'label']]
