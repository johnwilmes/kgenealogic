"""CLI for kgenealogic"""


from pathlib import Path
from typing import Optional, List
from typing_extensions import Annotated
import sys

import sqlalchemy as sql
import typer

import kgenealogic as kg
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
    kg.initialize(engine)

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
            kg.import_matches(engine, data)
        elif kg_files.is_ged_triangles(path):
            if not source:
                typer.echo("Source kit number required for GEDmatch triangulation files, aborting.",
                      err=True)
                raise typer.Exit(code=1)
            data = kg_files.read_ged_triangles(path, source)
            kg.import_triangles(engine, data)
        else:
            typer.echo(f"Unrecognized file type: {path}", err=True)

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

    - include: Kit numbers to be included are specified as a list. Each list entry can be either
      (1) just the kit number to be included, as a string, or (2) a mapping indicating a kit number
      and some of its matches to be included. In the second format, the valid keys are

      --- id: the (base) kit number as a string (required)

      --- matches (optional, floating point): if this is given, it is interpreted as the minimum
      length in cM of a match so that if kit2 matches the base kit anywhere with this length, then
      kit2 will also be included

      --- triangles (optional, floating point): if this is given, it is interpreted as the minimum
      length in cM of a triangle so that if kit2 participates in a triangle the base kit anywhere
      with this length, with any third kit, then kit2 will also be included.

      Seeds are automatically included.

    - exclude: Kit numbers to be excluded are specified as a list. E.g., known duplicate kits can
      be specified here, or individuals that are known to not respect the desired tree structure.
      Optional, defaults to an empty list. Anything listed (directly or by reference) in the
      "include" list and also in the "exclude" list will be excluded.

    - min_length: Minimum segment length to be considered, in cM. Optional, defaults to 7cM

    - tree: A tree giving the desired family structure to be produced, along with any kits known to reside at
      each node of the of the tree. The more accurate seeds are available, the more reliable and
      useful the results will be. This root node of the tree is specified with the "tree" key. At
      each node of the tree, you can optionally specify any of:

      1. a "kits" key, associated with a list of kit numbers known to reside at that node of the
      tree. See below for additional options.

      2. a "paternal" key, formatted as a node of the tree, specifying the paternal branch

      3. a "maternal" key, formatted as a node of the tree, specifying the maternal branch
      The tree key is required

    Entries of each "kits" list of the tree can be specified in one of two ways:

    - just the kit number as a string

    - a mapping with the following keys:

      --- id: the kit number as a string (required)

      --- autox: whether to classify kits that match this kit on the X chromosome as necessarily
      maternal (optional, boolean). If the kit corresponds to a male child of the parents
      represented by the paternal and maternal branches of this node of the tree, it makes sense to
      use a value of true. Otherwise omit the key or set it to false (the default)

      --- negative: whether to use negative triangulations for this node, if available (option,
      boolean, default false)

      --- float: whether to keep the kit at this node of the tree, or attempt to float it into the
      branches (optional, boolean). If float is false, this seed will be kept at this node of the
      tree - this is desired when the kit is known to be a descendant of the parents represented by
      the maternal and paternal branches at this node. If float is true, this seed will be pushed
      out as far as possible into the branches - this is desired when the kit's relation to the
      root is partially but not completely known (e.g., the kit is some relative of the maternal
      grandfather, but nothing more is known). If float is not set, the default behavior is to
      treat it as false if triangulations are available, and otherwise treat it as true.

    An example configuration file is available in the "examples" directory included with the source
    of this package.


    This will produce an output file separating my maternal and paternal relatives, with the
    paternal relatives further split into paternal grandfather and paternal grandmother relatives.
    The output file will be a CSV file with several columns:

    - kit, specifying the kit number

    - ahnentafel, the ahnentafel number of the most distant relative of the root to which the kit
      is determiend to be related

    - seed, the ahnentafel number of the node at which this kit was listed as a seed, if any

    - label<n>, where n varies from 0 to depth of the tree minus 1. This is 'M' if the kit is put
      on the maternal side at that depth, and P if it is put on the paternal side. The list of
      these labels is equivalent to the value in the ahnentafel field

    - confidence<n>, a number between 0 and 1 indicating a rough confidence with which the kit is
      given label<n>. These numbers should NOT be interpreted as probabilities. A value of 1
      indicates that all of the kit's matches have labels consistent with the kit's label at this
      depth.

    """

    engine = sql.create_engine(f"sqlite:///{project}")
    config = kg_files.read_cluster_config(config)
    clusters = kg.cluster_data(engine, config)
    kg_files.write_clusters(clusters, outfile)

