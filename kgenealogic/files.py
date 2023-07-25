import pandas as pd
import strictyaml as yaml
from collections import defaultdict
from dataclasses import dataclass, field
import typing

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

CLUSTER_TREE_SCHEMA = yaml.Map({
    yaml.Optional('kits'): yaml.UniqueSeq(yaml.Str()),
    yaml.Optional('maternal'): yaml.Any(),
    yaml.Optional('paternal'): yaml.Any(),
})
CLUSTER_CONFIG_SCHEMA = yaml.Map({
    yaml.Optional('exclude'): yaml.UniqueSeq(yaml.Str()),
    yaml.Optional('min_length'): yaml.Float(),
    'tree': CLUSTER_TREE_SCHEMA,
})

@dataclass
class ClusterConfig:
    """Data structure for storing parsed cluster configuration"""
    exclude: list[str]
    min_length: float
    tree: dict[str, typing.Any]

def validate_config_tree(config):
    for branch in ('maternal', 'paternal'):
        if branch in config.data:
            config[branch].revalidate(CLUSTER_TREE_SCHEMA)
            validate_config_tree(config[branch])

def read_cluster_config(path):
    with open(path) as f:
        config_yaml = f.read()
    parsed = yaml.load(config_yaml, CLUSTER_CONFIG_SCHEMA)
    validate_config_tree(parsed['tree'])

    result =  ClusterConfig(
        exclude=parsed.data.get('exclude', []),
        min_length=parsed.data.get('min_length', 7.),
        tree=parsed['tree'].data,
    )
    return result

def write_clusters(clusters, path):
    clusters.to_csv(path, index=False)
