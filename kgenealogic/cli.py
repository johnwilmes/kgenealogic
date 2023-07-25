"""CLI for kgenealogic"""


from pathlib import Path
from typing import Optional, List
from typing_extensions import Annotated
import sys

import sqlalchemy as sql
import typer

import kgenealogic
import kgenealogic.files as kg_files

import warnings
warnings.filterwarnings("ignore",
    category=UserWarning,
    message="pandas only supports SQLAlchemy connectable"
)

app = typer.Typer()

DEFAULT_PROJECT = Path("kgenealogic.db")
DEFAULT_OUTPUT = Path("kg_results.csv")
PROJECT_OPTION = Annotated[
    Optional[Path],
    typer.Option(
        "--project",
        "-p",
        file_okay=True,
        dir_okay=False,
        writable=True,
        readable=True,
        resolve_path=True,
        help=f"Project file for storing processed genealogic data",
    ),
]

@app.command()
def init(
    project: PROJECT_OPTION = DEFAULT_PROJECT,
    force:  Annotated[bool, typer.Option(help="Force reinitialization of existing project file")] = False,
):
    """
    Initialize a new project file for storing processed genealogic data.

    If the project file already exists, --force clears and reinitializes the project.
    """
    if project.exists() and not force:
        typer.echo("Project already exists. Use --force to reinitialize.", err=True)
        raise typer.Abort()

    engine = sql.create_engine(f"sqlite:///{project}")
    kgenealogic.initialize(engine)

@app.command()
def add(
    files: List[Path],
    source:  Annotated[
        Optional[str],
        typer.Option(help="For GEDmatch triangulation files, the source kit number that is triangulated")
    ] = None,
    project: PROJECT_OPTION = DEFAULT_PROJECT,
):
    """
    Add data files to a project.

    Recognized file types:

    - GEDmatch segment match files, in CSV format with fields PrimaryKit, MatchedKit, chr,
    B37Start, B37End, Segment cM, MatchedName, Matched Sex, MatchedEmail

    - GEDmatch triangulation files, in CSV format with fields Kit1 Number, Kit1 Name, Kit1 Email,
    Kit2 Number, Kit2 Name, Kit2 Email, Chr, B37 Start, B37 End, cM

    For GEDmatch triangulation files, the source/primary kit number is not given in the data file
    and must be passed separately with the --source option
    """
    engine = sql.create_engine(f"sqlite:///{project}")
    for path in files:
        if kg_files.is_ged_matches(path):
            data = kg_files.read_ged_matches(path)
            kgenealogic.import_matches(engine, data)
        elif kg_files.is_ged_triangles(path):
            if not source:
                typer.echo("Source kit number required for GEDmatch triangulation files, aborting.",
                      err=True)
                raise typer.Exit(code=1)
            data = kg_files.read_ged_triangles(path, source)
            kgenealogic.import_triangles(engine, data)
        else:
            typer.echo(f"Unrecognized file type: {path}", err=True)

@app.command()
def build(
    project: PROJECT_OPTION = DEFAULT_PROJECT,
    force:  Annotated[bool, typer.Option(help="Force re-build, even if already built")] = False
):
    """
    Process all added data files in preparation for clustering.

    - Uses all segments with provided lengths in order to build a model of crossover probability, to
    impute cM lengths of derived segments

    - Finds all "negative" triangles, where two target kits have overlapping matching segments on a
    common source kit, but the segment does not appear as a triangulation for the three kits
    """
    engine = sql.create_engine(f"sqlite:///{project}")
    if not force:
        if kgenealogic.is_cache_valid(engine):
            typer.echo("Project is already built. Use --force to re-build", err=True)
            raise typer.Abort()

    kgenealogic.build_cache(engine)

@app.command()
def build(
    project: PROJECT_OPTION = DEFAULT_PROJECT,
    force:  Annotated[bool, typer.Option(help="Force re-build, even if already built")] = False
):
    """
    Process all added data files in preparation for clustering.

    - Uses all segments with provided lengths in order to build a model of crossover probability, to
    impute cM lengths of derived segments

    - Finds all "negative" triangles, where two target kits have overlapping matching segments on a
    common source kit, but the segment does not appear as a triangulation for the three kits
    """
    engine = sql.create_engine(f"sqlite:///{project}")
    if not force:
        if kgenealogic.is_cache_valid(engine):
            typer.echo("Project is already built. Use --force to re-build", err=True)
            raise typer.Abort()

    kgenealogic.build_cache(engine)

@app.command()
def cluster(
    config: Annotated[Path, typer.Argument(
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help=f"The configuration file describing how to cluster kits",
        )],
    outfile: Annotated[Path, typer.Option(
        "--outfile",
        "-o",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
        help=f"The destination for the output"
        )] = DEFAULT_OUTPUT,
    project: PROJECT_OPTION = DEFAULT_PROJECT,
):
    """
    Cluster kits to predict family tree structure.

    A config file in (strict) YAML format must be provided. An example of such a file is
    distributed with this package. The valid keys for this config file are:

    - exclude: Kit numbers to be excluded are specified as a list. E.g., known duplicate kits can
      be specified here, or individuals that are known to not respect the desired tree structure.
      Optional, defaults to an empty list.

    - min_length: Minimum segment length to be considered, in cM. Optional, defaults to 7cM

    - tree: A tree giving the desired family structure to be produced, along with any kits known to reside at
      each node of the of the tree. The more accurate seeds are available, the more reliable and
      useful the results will be. This root node of the tree is specified with the "tree" key. At
      each node of the tree, you can optionally specify any of:
      1. a "kits" key, associated with a list of kit numbers known to reside at that node of the
      tree. See below for additional options.
      2. a "paternal" key, formatted as a node of the tree, specifying the paternal branch
      2. a "maternal" key, formatted as a node of the tree, specifying the maternal branch
      The tree key is required

    Entries of each "kits" list of the tree can be specified in one of two ways:

    - just the kit number as a string
    - a mapping with the following keys:
      -- id: the kit number as a string (required)
      -- autox: whether to classify kits that match this kit on the X chromosome as necessarily
      maternal (optional, boolean). If the kit corresponds to a male child of the parents
      represented by the paternal and maternal branches of this node of the tree, it makes sense to
      use a value of true. Otherwise omit the key or set it to false (the default)
      -- float: whether to keep the kit at this node of the tree, or attempt to float it into the
      branches (optional, boolean). If float is false, this seed will be kept at this node of the
      tree - this is desired when the kit is known to be a descendant of the parents represented by
      the maternal and paternal branches at this node. If float is true, this seed will be pushed
      out as far as possible into the branches - this is desired when the kit's relation to the
      root is partially but not completely known (e.g., the kit is some relative of the maternal
      grandfather, but nothing more is known). Trianglulations for which the kit is the source are
      NOT used if float is set to true. If float is not set, the default behavior is to treat it as
      false if triangulations are available, and otherwise treat it as true.

    For example:

    exclude: # these kits will be excluded from the clustering
      - T000001
      - T000002
    min_length: 8.5
    tree:
      kits:
        - # me (male)
          id: T000003
          autox: true
          float: false
        - # my sister's son
          id: T000004
          float: false
      paternal:
        kits:
          - T000004 # my father's grand-niece
        paternal:
          kits:
            - T000005 # a distant cousin on my paternal grandfather's side
        maternal:
          # I don't know any relatives of my paternal grandmother, but I'd like any that are
          # discovered to be present in the output
      maternal:
        kits:
          - T000006 # a relative of my mother, but not my father, with no other information

    This will produce an output file separating my maternal and paternal relatives, with the
    paternal relatives further split into paternal grandfather and paternal grandmother relatives.
    The output file will be a CSV file with two columns: kit, specifying the kit number, and label,
    a string consisting of Ms and Ps. The label MMPM means a relative of the root's mother's
    mother's father's mother. Labels may be shorter than the desired tree depth - or even empty
    strings - if there is not enough information to place them deeper in the tree.
    """

    engine = sql.create_engine(f"sqlite:///{project}")
    if not kgenealogic.is_cache_valid(engine):
        typer.echo("Project is not yet built. First run kgenealogic build", err=True)
        raise typer.Abort()
    config = kg_files.read_cluster_config(config)
    clusters = kgenealogic.cluster_data(engine, config)
    kg_files.write_clusters(clusters, outfile)

