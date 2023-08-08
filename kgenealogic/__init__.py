"""Base package for kgenealogic"""

__app_name__ = "kgenealogic"
__version__ = "0.1.0"

from kgenealogic.schema import initialize
from kgenealogic.data import import_matches, import_triangles
from kgenealogic.cluster import cluster_data
