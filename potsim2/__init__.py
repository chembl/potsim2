from .PotGrid import PotGrid

try:
    from importlib.metadata import version
    __version__ = version("potsim2")
except ImportError:
    __version__ = "unknown"
