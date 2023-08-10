import pandas as pd
import yaml
from collections import defaultdict
from dataclasses import dataclass, field
import typing
import typer

DEFAULT_MIN_LENGTH = 7.

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
    """Check if file at path is a valid GEDmatch pairwise matches file."""
    try:
        csv = pd.read_csv(path, dtype=str)
    except:
        return False
    return set(GED_MATCH_COLS.keys()).issubset(csv.columns)

def read_ged_matches(path):
    """Load GEDmatch pairwise matches file into pandas dataframe."""
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
    """Check if file at path is a valid GEDmatch triangulation file."""
    try:
        csv = pd.read_csv(path, dtype=str)
    except:
        return False
    return set(GED_TRIANGLE_COLS.keys()).issubset(csv.columns)

def read_ged_triangles(path, primary_kit):
    """Load GEDmatch triangulation file into pandas dataframe.

    Args:
        path: path to the file
        primary_kit: the kit from which the triangulations were generated

    Returns:
        A pandas dataframe in the format expected by data.import_triangles.
        Specifically, primary_kit will appear as kit1 in the dataframe
    """
    dtype = defaultdict(lambda: str)
    dtype['B37 Start'] = int
    dtype['B37 End'] = int
    tri = pd.read_csv(path, dtype=dtype)
    tri = tri[GED_TRIANGLE_COLS.keys()].rename(columns=GED_TRIANGLE_COLS)
    tri['kit1'] = primary_kit
    return tri

@dataclass
class ClusterConfig:
    """Data structure for storing parsed cluster configuration"""
    include: list[dict[str, typing.Any]]
    exclude: list[str]
    min_length: float
    tree: dict[str, typing.Any]
    seeds: set[str]

def expect_keys(d, l, error):
    bad_keys = set(d.keys())
    bad_keys.difference_update(l)
    if bad_keys:
        invalid = bad_keys.pop()
        typer.echo(error.format(invalid), err=True)
        raise typer.Exit(code=1)

def parse_config_tree(tree, seeds):
    """Recursively expand and extract all seeds from YAML config tree

    Seeds that appear in the config file as strings are expanded to dictionaries with default
    values for the keys and reinserted into the tree

    Args:
        tree: the dictionary representing the current node of the tree
        seeds: the list of all encountered seeds
    """
    expect_keys(tree, ['kits', 'maternal', 'paternal'], 
                "invalid YAML configuration format: invalid tree key {}")

    expanded_kits = []
    for k in tree.get('kits', []):
        if not k:
            typer.echo(f"invalid (false/missing) tree kits entry", err=True)
            raise typer.Exit(code=1)
        elif type(k)!=dict:
            expanded_kits.append(dict(id=str(k), autox=False, float=None, negative=False))
        else:
            expect_keys(k, ['id', 'autox', 'float', 'negative'],
                        "invalid YAML configuration format: invalid tree kits key {}")
            expanded_kits.append(k)
    tree['kits']=expanded_kits
    seeds.extend(expanded_kits)

    for branch in ('maternal', 'paternal'):
        if branch in tree:
            value = tree[branch]
            if value:
                if type(value)!=dict:
                    typer.echo(f"invalid {branch} branch: {value}", err=True)
                    raise typer.Exit(code=1)
                parse_config_tree(tree[branch], seeds)
            else:
                tree[branch] = {}

def read_cluster_config(path):
    """Parse a clustering config file in YAML format.

    The schema is given above in code, and is documented in the cluster CLI command.

    Args:
        path: path to the cluster config file

    Returns:
        A ClusterConfig object representing the contents of the config file
    """
    with open(path) as f:
        parsed = yaml.safe_load(f)

    if 'min_length' in parsed:
        min_length = parsed['min_length']
        if min_length != 0 and not min_length:
            min_length = DEFAULT_MIN_LENGTH
        else:
            min_length = float(min_length)
    else:
        min_length = DEFAULT_MIN_LENGTH

    tree = dict(parsed.get('tree') or {})
    exclude = list(parsed.get('exclude') or [])

    seeds = []
    parse_config_tree(tree, seeds)
    
    include = []
    include_raw = parsed.get('include', []) or []
    if type(include_raw)!=list:
        typer.echo(f"invalid include value: expected list", err=True)
        raise typer.Exit(code=1)
    for k in include_raw:
        if not k:
            typer.echo(f"invalid (false/missing) include list entry", err=True)
            raise typer.Exit(code=1)
        elif type(k)!=dict:
            include.append(dict(id=str(k), matches=None, triangles=None))
        else:
            expect_keys(k, ['id', 'matches', 'triangles'],
                        "invalid YAML configuration format: invalid include item key {}")
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

    result = ClusterConfig(
        include=include,
        exclude=exclude,
        min_length=min_length,
        tree=tree,
        seeds=unique_seeds,
    )
    return result

def write_clusters(clusters, path):
    """Write out the results of clustering to a CSV file."""
    clusters.to_csv(path, index=False)
