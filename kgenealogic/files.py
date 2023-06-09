import pandas as pd
from collections import defaultdict

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
