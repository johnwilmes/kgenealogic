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

@dataclass
class SeedTree:
    """Data structure for organizing cluster seeds in tree"""
    values: list[Seed] = field(default_factory=list)
    branches: dict[str, typing.Any] = field(default_factory=dict)

def flatten(tree, values):
    for b in tree.branches.values():
        flatten(b, values)
    values.extend([v.kit for v in tree.values])
    return values

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
    trisource_query = (
        sql.select(source.c.kit)
        .where(source.c.has_triangles)
        .where(source.c.has_negative)
    )
    autox_query = (
        sql.select(sql.distinct(match.c.kit2).label('kit'))
        .join_from(match, segment, match.c.segment==segment.c.id)
        .where(segment.c.chromosome=='X')
        .where(segment.c.length >= config.min_length)
        .where(match.c.kit1==sql.bindparam('seed'))
   )

    with engine.connect() as conn:
        kits = pd.read_sql(kits_query, conn).set_index('kitid').id
        trisource= set(int(x) for x in pd.read_sql(trisource_query, conn).kit)
        graph = pd.read_sql(graph_query, conn)

    graph['weight'] = PAIRWISE_FACTOR*graph['weight'].astype(float)

    seeds = set()
    exclude = set()
    for kitid_list, kit_set in ((config.seeds, seeds), (config.exclude, exclude)):
        for kitid in kitid_list:
            if kitid not in kits.index:
                typer.echo(f"Invalid kit id {kitid}", err=True)
                raise typer.Exit(code=1)
            else:
                kit_set.add(int(kits.loc[kitid]))

    kits = kits[~kits.index.isin(exclude)]

    def build_seed_tree(raw):
        parsed = SeedTree()
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
            parsed.values.append(Seed(kit=kitlocal, floating=floating))

        for short, long in (('M', 'maternal'), ('P', 'paternal')):
            if long in raw:
                parsed.branches[short] = build_seed_tree(raw[long])

        if autox:
            maternal_values = parsed.branches.setdefault('M', SeedTree(values=[])).values
            for x in autox:
                maternal_values.append(Seed(kit=x, floating=True))
        return parsed
    tree = build_seed_tree(config.tree)

    source_tri = lambda source: get_triangles(engine, source, config.min_length)
    clusters = recursive_cluster(kits, tree, graph, source_tri)

    kitids = kits.reset_index().set_index('id').kitid
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

def recursive_cluster(kits, tree, graph, source_tri):
    # graph format: kit1/kit2/weight
    # initial value: all pairwise matches, suitably weighted relative to triangles that will be added?
    nonfloat = []
    for source in tree.values:
        if not source.floating:
            graph = graph.set_index(['kit1', 'kit2']).add(
                source_tri(source.kit).set_index(['kit1', 'kit2']), fill_value=0
            ).reset_index()
            nonfloat.append(source.kit)
        
    kits = kits[~kits.isin(nonfloat)]
    graph = graph[graph.kit1.isin(kits)&graph.kit2.isin(kits)]
    
    if tree.branches:
        flat_seeds = []
        for label, branch in tree.branches.items():
            branch_seeds = flatten(branch, [])
            branch_seeds = pd.DataFrame(dict(kit=branch_seeds))
            branch_seeds['label'] = label
            flat_seeds.append(branch_seeds)
        flat_seeds = pd.concat(flat_seeds, ignore_index=True)
        clusters = get_clusters(graph, flat_seeds)

        for label in tree.branches:
            branch = clusters[clusters.label==label]
            if len(branch) > 0:
                branch = recursive_cluster(branch.kit, tree.branches[label], graph, source_tri)
                branch_labels = label + clusters.kit.map(branch.set_index('kit').label)
                clusters.label.where(branch_labels.isnull(), branch_labels, inplace=True)
    else:
        clusters = pd.DataFrame(dict(kit=kits))
        clusters['label'] = ''

    nonfloat = pd.DataFrame(dict(kit=nonfloat))
    nonfloat['label'] = ''

    return pd.concat([clusters, nonfloat], ignore_index=True)

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
