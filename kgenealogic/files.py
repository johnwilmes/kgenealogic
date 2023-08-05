import pandas as pd
import strictyaml as yaml
from collections import defaultdict
from dataclasses import dataclass, field
import typing
import typer

GED_MATCH_COLS = {
    'PrimaryKit': 'kit1',
    'MatchedKit': 'kit2',
    'chr': 'chromosome',
    'B37Start': 'start',
    'B37End': 'end',
    'Segment cM': 'length',
    'MatchedName': 'name',
    'Matched Sex': 'sex',
    'MatchedEmail': 'email',
}
def is_ged_matches(path):
    try:
        csv = pd.read_csv(path, dtype=str)
    except:
        return False
    return set(GED_MATCH_COLS.keys()).issubset(csv.columns)

def read_ged_matches(path):
    dtype = defaultdict(lambda: str)
    dtype['B37Start'] = int
    dtype['B37End'] = int
    matches = pd.read_csv(path, dtype=dtype)
    return matches[GED_MATCH_COLS.keys()].rename(columns=GED_MATCH_COLS)

GED_TRIANGLE_COLS = {
    'Kit1 Number': 'kit2',
    'Kit1 Name': 'name2',
    'Kit1 Email': 'email2',
    'Kit2 Number': 'kit3',
    'Kit2 Name': 'name3',
    'Kit2 Email': 'email3',
    'Chr': 'chromosome',
    'B37 Start': 'start',
    'B37 End': 'end',
    'cM': 'length',
}
def is_ged_triangles(path):
    try:
        csv = pd.read_csv(path, dtype=str)
    except:
        return False
    return set(GED_TRIANGLE_COLS.keys()).issubset(csv.columns)

def read_ged_triangles(path, primary_kit):
    dtype = defaultdict(lambda: str)
    dtype['B37 Start'] = int
    dtype['B37 End'] = int
    tri = pd.read_csv(path, dtype=dtype)
    tri = tri[GED_TRIANGLE_COLS.keys()].rename(columns=GED_TRIANGLE_COLS)
    tri['kit1'] = primary_kit
    return tri

CLUSTER_KIT_SCHEMA = yaml.Str() | yaml.Map({
    'id': yaml.Str(),
    yaml.Optional('autox', default=False): yaml.Bool(),
    yaml.Optional('float', default=None, drop_if_none=False): yaml.EmptyNone() | yaml.Bool(),
    yaml.Optional('negative', default=None, drop_if_none=False): yaml.EmptyNone() | yaml.Bool(),
})
CLUSTER_INCLUDE_SCHEMA = yaml.Str() | yaml.Map({
    'id': yaml.Str(),
    yaml.Optional('matches', default=None, drop_if_none=False): yaml.EmptyNone() | yaml.Float(),
    yaml.Optional('triangles', default=None, drop_if_none=False): yaml.EmptyNone() | yaml.Float(),
})
CLUSTER_TREE_SCHEMA = yaml.Map({
    yaml.Optional('kits'): yaml.EmptyList() | yaml.Seq(CLUSTER_KIT_SCHEMA),
    yaml.Optional('maternal'): yaml.Any(),
    yaml.Optional('paternal'): yaml.Any(),
})
CLUSTER_CONFIG_SCHEMA = yaml.Map({
    yaml.Optional('include'): yaml.EmptyList() | yaml.Seq(CLUSTER_INCLUDE_SCHEMA),
    yaml.Optional('exclude'): yaml.EmptyList() | yaml.UniqueSeq(yaml.Str()),
    yaml.Optional('min_length'): yaml.Float(),
    'tree': CLUSTER_TREE_SCHEMA,
})

@dataclass
class ClusterConfig:
    """Data structure for storing parsed cluster configuration"""
    include: list[dict[str, typing.Any]]
    exclude: list[str]
    min_length: float
    tree: dict[str, typing.Any]
    seeds: set[str]

def validate_config_tree(config, seeds):
    expanded_kits = []
    for k in config.data.get('kits', []):
        if type(k)==str:
            expanded_kits.append(dict(id=k, autox=False, float=None, negative=False))
        else:
            expanded_kits.append(k)
    config['kits']=expanded_kits
    seeds.extend(expanded_kits)

    for branch in ('maternal', 'paternal'):
        if branch in config.data:
            config[branch].revalidate(yaml.EmptyDict() | CLUSTER_TREE_SCHEMA)
            validate_config_tree(config[branch], seeds)

def read_cluster_config(path):
    with open(path) as f:
        config_yaml = f.read()
    parsed = yaml.load(config_yaml, CLUSTER_CONFIG_SCHEMA)
    seeds = []
    validate_config_tree(parsed['tree'], seeds)
    exclude = parsed.data.get('exclude', [])

    include = []
    for k in parsed.data.get('include', []):
        if type(k)==str:
            include.append(dict(id=k, matches=None, triangles=None))
        else:
            include.append(k)

    unique_seeds = set()
    for s in seeds:
        if s['id'] in unique_seeds:
            typer.echo("duplicated seed {}".format(s['id']), err=True)
            raise typer.Exit(code=1)
        else:
            unique_seeds.add(s['id'])
    exclude_seeds = unique_seeds.intersection(exclude)
    if exclude_seeds:
        k = exclude_seeds.pop()
        typer.echo(f"excluded kit {k} is listed as seed", err=True)
        raise typer.Exit(code=1)

    result =  ClusterConfig(
        include=include,
        exclude=exclude,
        min_length=parsed.data.get('min_length', 7.),
        tree=parsed['tree'].data,
        seeds=unique_seeds,
    )
    return result

def write_clusters(clusters, path):
    clusters.to_csv(path, index=False)
