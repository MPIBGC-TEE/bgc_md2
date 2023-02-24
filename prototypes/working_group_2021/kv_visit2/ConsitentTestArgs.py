import json
from pathlib import Path


class ConsistentTestArgs():
    def __init__(
        """A model (and data) dependent class that expresses
        dependencies between certain parts of the parameterization as
        dictated by a particular choice of estimated and constant
        paramters for parameter estimation """
            self,
            targetPath: Path 
        ):
        self.targetPath = targetPath
    
    @property
    def timelines():
        svs, dvs = msh.get_cached_global_mean_vars(Path(__file__).parent)
        return svs, dvs
