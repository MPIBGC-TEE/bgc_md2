import os
from abc import ABC, abstractmethod
import netCDF4 as nc
from sympy import Symbol
from . import helper as h

def write_driver_data(path, dvs):
    if path.exists():
        os.remove(path)
    dp=path.parent
    if not dp.exists():
        dp.mkdir(parents=True)

    with nc.Dataset(str(path),"w",persist=True) as ds:
        for k,v in dvs._asdict().items():
            dim_name=f"time_{k}"
            ds.createDimension(dim_name,size=len(v)) 
            ds.createVariable(varname=k,datatype=float,dimensions=(
                dim_name,)
            )
            ds.variables[k][:]=v 


class VirtualCachedParameterization(ABC):
    par_dict_path = "param_dict.json"
    func_dict_param_path = "func_param_dict.json"
    func_dict_data_path = "Drivers.nc"
    X_0_data_path = "X_0.json"
    file_names = [
        par_dict_path,
        func_dict_data_path,
        X_0_data_path,
    ]
    # The subclasses have to set some class variables 
    # Drivers= ... (a namedtuple)
    func_dict_param_keys=None # default 
 

    def __init__(
            self,
            parameter_dict,
            drivers,
            X_0_dict,
            func_dict_param_dict=None
        ):
            self.parameter_dict = parameter_dict
            self.drivers=drivers
            self.X_0_dict = X_0_dict
            self.func_dict_param_dict=func_dict_param_dict


    @classmethod
    def from_path(cls,data_path):
        ds=nc.Dataset(
            str(
                data_path.joinpath(
                    cls.func_dict_data_path
                )
            ),
            "r"
        )
        Drivers = cls.Drivers
        dvs = Drivers(*[ds.variables[k][:] for k in Drivers._fields])
        str_dict = h.load_dict_from_json_path(
                data_path.joinpath(cls.par_dict_path)
            )
        par_dict = {Symbol(k): v for k,v in str_dict.items()}
        X_0_dict = h.load_dict_from_json_path((
            data_path.joinpath(cls.X_0_data_path))
        )
        fdpdp = data_path.joinpath(cls.func_dict_param_path)
        func_dict_param_dict = h.load_dict_from_json_path(fdpdp) if fdpdp.exists() else None

        
        return cls(
            parameter_dict=par_dict,
            drivers=dvs,
            X_0_dict=X_0_dict,
            func_dict_param_dict=func_dict_param_dict
        )

    @property
    @abstractmethod
    def func_dict(self):
        pass

    def write(self, data_path):

        h.dump_dict_to_json_path(
            self.parameter_dict,
            data_path.joinpath(self.par_dict_path),
            indent=2
        )


        write_driver_data(
            data_path.joinpath(self.func_dict_data_path),
            self.drivers
        )

        h.dump_dict_to_json_path(
            self.X_0_dict,
            data_path.joinpath(self.X_0_data_path),
            indent=2
        )
        if self.func_dict_param_dict is not None:
            h.dump_dict_to_json_path(
                self.func_dict_param_dict,
                data_path.joinpath(self.func_dict_param_path),
                indent=2
            )


