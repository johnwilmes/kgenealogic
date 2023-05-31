"""Base package for kgenealogic"""

__app_name__ = "kgenealogic"
__version__ = "0.1.0"

from data import initialize, is_cache_valid, import_matches, import_triangles
from cache import build_cache
