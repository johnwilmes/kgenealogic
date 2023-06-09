"""Base package for kgenealogic"""

__app_name__ = "kgenealogic"
__version__ = "0.1.0"

from schema import initialize
from data import import_matches, import_triangles
from cache import is_cache_valid, build_cache
