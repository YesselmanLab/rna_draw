from importlib.metadata import PackageNotFoundError, version

from .draw import rna_draw

try:
    __version__ = version("rna_draw")
except PackageNotFoundError:  # source tree without install
    __version__ = "0.0.0"

__all__ = ["rna_draw", "__version__"]
