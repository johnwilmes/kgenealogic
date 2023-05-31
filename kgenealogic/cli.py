"""CLI for kgenealogic"""


from pathlib import Path
from typing import Optional, List
from typing_extensions import Annotated
import sys

import pysqlite3
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
        print("Project already exists. Use --force to reinitialize.")
        raise typer.Abort()

    db = pysqlite3.connect(project)
    kgenealogic.initialize(db)

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
    db = pysqlite3.connect(project)
    for path in files:
        if kg_files.is_ged_matches(path):
            data = kg_files.read_ged_matches(path)
            kgenealogic.import_matches(db, data)
        elif kg_files.is_ged_triangles(path):
            if not source:
                print("Source kit number required for GEDmatch triangulation files, aborting.",
                      file=sys.stderr)
                raise typer.Exit(code=1)
            data = kg_files.read_ged_triangles(path, source)
            kgenealogic.import_triangles(db, data)
        else:
            print(f"Unrecognized file type: {path}")

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
    db = pysqlite3.connect(project)
    if not force:
        if kgenealogic.is_cache_valid(db):
            print("Project is already built. Use --force to re-build")
            raise typer.Abort()

    kgenealogic.build_cache(db)

