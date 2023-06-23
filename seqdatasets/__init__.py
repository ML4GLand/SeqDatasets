try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata

package_name = "seqdatasets"
__version__ = importlib_metadata.version(package_name)

from ._datasets import get_dataset_info, random1000
from ._datasets import jores21, ray13, deBoer20, deAlmeida22
