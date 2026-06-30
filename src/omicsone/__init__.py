"""Public package name for the OmicsOne Streamlit application."""

try:
    from importlib.metadata import PackageNotFoundError, version
except ImportError:  # pragma: no cover
    PackageNotFoundError = Exception
    version = None

try:
    __version__ = version("omicsone") if version else "0+unknown"
except PackageNotFoundError:
    __version__ = "0+unknown"

__all__ = ["__version__"]