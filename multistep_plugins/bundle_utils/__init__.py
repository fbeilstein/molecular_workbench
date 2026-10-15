"""Bundle utilities for Step 5: standalone HTML viewer, high-res QM evaluation, and zip packager."""

from .standalone_template import generate_standalone_html
from .highres_engine import (
    recompute_highres_batch,
    get_highres_grid,
    get_cached_highres_cube,
    _HIGHRES_CUBE_CACHE,
    _HIGHRES_GRID_CACHE,
)
from .packager import (
    generate_xyz,
    build_bundle_manifest,
    create_maxcomp_bundle_zip,
    unpack_molecule_bundle,
    purge_cube_files,
)

__all__ = [
    "generate_standalone_html",
    "recompute_highres_batch",
    "get_highres_grid",
    "get_cached_highres_cube",
    "_HIGHRES_CUBE_CACHE",
    "_HIGHRES_GRID_CACHE",
    "generate_xyz",
    "build_bundle_manifest",
    "create_maxcomp_bundle_zip",
    "unpack_molecule_bundle",
    "purge_cube_files",
]
