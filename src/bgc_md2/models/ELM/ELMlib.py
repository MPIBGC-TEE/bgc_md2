from netCDF4 import Dataset
import numpy as np
import pandas as pd

from ModelDataObject import (
    ModelDataObject,
    getStockVariable_from_Density,
    getFluxVariable_from_DensityRate,
)
from ModelStructure import ModelStructure
from SmoothModelRunFromData import SmoothModelRunFromData as SMRFD
from Variable import Variable




################################################################################


def load_model_structure():
    pool_structure = [
        {
            "pool_name": "CWD",
            "stock_var": "CWDC_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
        {
            "pool_name": "Litter 1",
            "stock_var": "LITR1C_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
        {
            "pool_name": "Litter 2",
            "stock_var": "LITR2C_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
        {
            "pool_name": "Litter 3",
            "stock_var": "LITR3C_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
        {
            "pool_name": "Soil 1",
            "stock_var": "SOIL1C_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
        {
            "pool_name": "Soil 2",
            "stock_var": "SOIL2C_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
        {
            "pool_name": "Soil 3",
            "stock_var": "SOIL3C_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
    ]

    external_input_structure = {
        "CWD": [
            "fire_mortality_c_to_cwdc",
            "gap_mortality_c_to_cwdc",
            "harvest_c_to_cwdc",
        ],
        "Litter 1": [
            "gap_mortality_c_to_litr_met_c",
            "harvest_c_to_litr_met_c",
            "m_c_to_litr_met_fire",
            "phenology_c_to_litr_met_c",
        ],
        "Litter 2": [
            "gap_mortality_c_to_litr_cel_c",
            "harvest_c_to_litr_cel_c",
            "m_c_to_litr_cel_fire",
            "phenology_c_to_litr_cel_c",
        ],
        "Litter 3": [
            "gap_mortality_c_to_litr_lig_c",
            "harvest_c_to_litr_lig_c",
            "m_c_to_litr_lig_fire",
            "phenology_c_to_litr_lig_c",
        ],
    }

    horizontal_structure = {
        ("CWD", "Litter 2"): ["CWDC_TO_LITR2C_vr"],
        ("CWD", "Litter 3"): ["CWDC_TO_LITR3C_vr"],
        ("Litter 1", "Soil 1"): ["LITR1C_TO_SOIL1C_vr"],
        ("Litter 2", "Soil 1"): ["LITR2C_TO_SOIL1C_vr"],
        ("Litter 3", "Soil 2"): ["LITR3C_TO_SOIL2C_vr"],
        ("Soil 1", "Soil 2"): ["SOIL1C_TO_SOIL2C_vr"],
        ("Soil 1", "Soil 3"): ["SOIL1C_TO_SOIL3C_vr"],
        ("Soil 2", "Soil 1"): ["SOIL2C_TO_SOIL1C_vr"],
        ("Soil 2", "Soil 3"): ["SOIL2C_TO_SOIL3C_vr"],
        ("Soil 3", "Soil 1"): ["SOIL3C_TO_SOIL1C_vr"],
    }

    vertical_structure = {
        "CWD": {
            "to_below": ["CWD_diffus_down"],
            "from_below": ["CWD_diffus_up"],
            "to_above": [],
            "from_above": [],
        },
        "Litter 1": {
            "to_below": ["LITR1_diffus_down"],
            "from_below": ["LITR1_diffus_up"],
            "to_above": [],
            "from_above": [],
        },
        "Litter 2": {
            "to_below": ["LITR2_diffus_down"],
            "from_below": ["LITR2_diffus_up"],
            "to_above": [],
            "from_above": [],
        },
        "Litter 3": {
            "to_below": ["LITR3_diffus_down"],
            "from_below": ["LITR3_diffus_up"],
            "to_above": [],
            "from_above": [],
        },
        "Soil 1": {
            "to_below": ["SOIL1_diffus_down"],
            "from_below": ["SOIL1_diffus_up"],
            "to_above": [],
            "from_above": [],
        },
        "Soil 2": {
            "to_below": ["SOIL2_diffus_down"],
            "from_below": ["SOIL2_diffus_up"],
            "to_above": [],
            "from_above": [],
        },
        "Soil 3": {
            "to_below": ["SOIL3_diffus_down"],
            "from_below": ["SOIL3_diffus_up"],
            "to_above": [],
            "from_above": [],
        },
    }

    external_output_structure = {
        "CWD": ["CWD_HR_L2_vr", "CWD_HR_L3_vr", "M_CWDC_TO_FIRE_vr"],
        "Litter 1": ["LITR1_HR_vr", "M_LITR1C_TO_FIRE_vr"],
        "Litter 2": ["LITR2_HR_vr", "M_LITR2C_TO_FIRE_vr"],
        "Litter 3": ["LITR3_HR_vr", "M_LITR3C_TO_FIRE_vr"],
        "Soil 1": ["SOIL1_HR_S2_vr", "SOIL1_HR_S3_vr"],
        "Soil 2": ["SOIL2_HR_S1_vr", "SOIL2_HR_S3_vr"],
        "Soil 3": ["SOIL3_HR_vr"],
    }

    model_structure = ModelStructure(
        pool_structure=pool_structure,
        external_input_structure=external_input_structure,
        horizontal_structure=horizontal_structure,
        vertical_structure=vertical_structure,
        external_output_structure=external_output_structure,
    )

    return model_structure


def load_mdo(ds):
    ms = load_model_structure()

    mdo = ModelDataObject(
        model_structure=ms,
        dataset=ds,
        stock_unit=

def load_mdo_12C(parameter_set):
    #    dataset         = parameter_set['dataset']
    ds_filename = parameter_set["ds_filename"]
    stock_unit = parameter_set["stock_unit"]
    cftime_unit_src = parameter_set["cftime_unit_src"]
    calendar_src = parameter_set["calendar_src"]
    nr_layers = parameter_set["nr_layers"]

    dz_var_name = parameter_set.get("dz_var_name", None)
    dz_var_names = parameter_set.get("dz_var_names", dict())
    min_time = parameter_set.get("min_time", 0)
    max_time = parameter_set.get("max_time", None)
    nstep = parameter_set.get("nstep", 1)
    time_shift = parameter_set.get("time_shift", 0)

    mdo_12C = ModelDataObject(
        model_structure=model_structure,
        #        dataset         = dataset,
        ds_filename=ds_filename,
        nstep=nstep,
        stock_unit=stock_unit,
        cftime_unit_src=cftime_unit_src,
        calendar_src=calendar_src,
        min_time=min_time,
        max_time=max_time,
        time_shift=time_shift,
        dz_var_names=dz_var_names,
    )

    return mdo_12C


## stub of 14C model
def load_mdo_14C(parameter_set):
    #    dataset         = parameter_set['dataset']
    ds_filename = parameter_set["ds_filename"]
    stock_unit = parameter_set["stock_unit"]
    cftime_unit_src = parameter_set["cftime_unit_src"]
    calendar_src = parameter_set["calendar_src"]
    nr_layers = parameter_set["nr_layers"]

    dz_var_name = parameter_set.get("dz_var_name", None)
    dz_var_names = parameter_set.get("dz_var_names", dict())
    min_time = parameter_set.get("min_time", 0)
    max_time = parameter_set.get("max_time", None)
    nstep = parameter_set.get("nstep", 1)
    time_shift = parameter_set.get("time_shift", 0)

    pool_structure_14C = [
        {
            "pool_name": "CWD",
            "stock_var": "C14_CWDC_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
        {
            "pool_name": "Litter 1",
            "stock_var": "C14_LITR1C_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
        {
            "pool_name": "Litter 2",
            "stock_var": "C14_LITR2C_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
        {
            "pool_name": "Litter 3",
            "stock_var": "C14_LITR3C_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
        {
            "pool_name": "Soil 1",
            "stock_var": "C14_SOIL1C_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
        {
            "pool_name": "Soil 2",
            "stock_var": "C14_SOIL2C_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
        {
            "pool_name": "Soil 3",
            "stock_var": "C14_SOIL3C_vr",
            "nr_layers": nr_layers,
            "dz_var": dz_var_name,
        },
    ]

    external_input_structure_14C = {
        "CWD": [
            "C14_fire_mortality_c_to_cwdc",
            "C14_gap_mortality_c_to_cwdc",
            "C14_harvest_c_to_cwdc",
        ],
        "Litter 1": [
            "C14_gap_mortality_c_to_litr_met_c",
            "C14_harvest_c_to_litr_met_c",
            "C14_m_c_to_litr_met_fire",
            "C14_phenology_c_to_litr_met_c",
        ],
        "Litter 2": [
            "C14_gap_mortality_c_to_litr_cel_c",
            "C14_harvest_c_to_litr_cel_c",
            "C14_m_c_to_litr_cel_fire",
            "C14_phenology_c_to_litr_cel_c",
        ],
        "Litter 3": [
            "C14_gap_mortality_c_to_litr_lig_c",
            "C14_harvest_c_to_litr_lig_c",
            "C14_m_c_to_litr_lig_fire",
            "C14_phenology_c_to_litr_lig_c",
        ],
    }

    external_output_structure_14C = {
        "CWD": ["C14_CWD_HR_L2_vr", "C14_CWD_HR_L3_vr", "C14_M_CWDC_TO_FIRE_vr"],
        "Litter 1": ["C14_LITR1_HR_vr", "C14_M_LITR1C_TO_FIRE_vr"],
        "Litter 2": ["C14_LITR2_HR_vr", "C14_M_LITR2C_TO_FIRE_vr"],
        "Litter 3": ["C14_LITR3_HR_vr", "C14_M_LITR2C_TO_FIRE_vr"],
        "Soil 1": ["C14_SOIL1_HR_S2_vr", "C14_SOIL1_HR_S3_vr"],
        "Soil 2": ["C14_SOIL2_HR_S1_vr", "C14_SOIL2_HR_S3_vr"],
        "Soil 3": ["C14_SOIL3_HR_vr"],
    }

    model_structure_14C = ModelStructure(
        pool_structure=pool_structure_14C,
        external_input_structure=external_input_structure_14C,
        horizontal_structure=[],
        external_output_structure=external_output_structure_14C,
    )

    mdo_14C = ModelDataObject(
        model_structure=model_structure_14C,
        #        dataset         = dataset,
        ds_filename=ds_filename,
        nstep=nstep,
        stock_unit=stock_unit,
        cftime_unit_src=cftime_unit_src,
        calendar_src=calendar_src,
        min_time=min_time,
        max_time=max_time,
        time_shift=time_shift,
        dz_var_names=dz_var_names,
    )

    return mdo_14C


def save_SMRFD_12C(parameter_set, filename, location):
    mdo = load_mdo_12C(parameter_set)

    smrfd = None
    log = None

    try:
        smrfd = SMRFD.load_from_file(filename)
    except (FileNotFoundError, OSError) as e:
        print(type(e), e, flush=True)

        try:
            smrfd = mdo.create_model_run(
                lat_index=location["lat_index"], lon_index=location["lon_index"]
            )
            if smrfd is not None:
                smrfd.location = location
                smrfd.save_to_file(filename)
                print("saved", filename, flush=True)
            else:
                print("masked data", flush=True)
                log = "masked"
        except (ValueError, np.linalg.LinAlgError) as e:
            print(type(e), e, flush=True)
            log = str(e)

    return (location["cell_nr"], log, filename)


def save_SMRFD_14C(parameter_set, filename_12C, filename_14C, location, xs_14C, us_14C):
    mdo = load_mdo_14C(parameter_set)

    smrfd_14C = None
    log = None

    try:
        smrfd_14C = SMRFD.load_from_file(filename_14C)
    except (FileNotFoundError, OSError) as e:
        print(type(e), e, flush=True)

        try:
            smrfd_12C = SMRFD.load_from_file(filename_12C)
            nr_pools = smrfd_12C.model.nr_pools
            xs_14C_data = np.array(xs_14C).reshape((-1, nr_pools))
            us_14C_data = np.array(us_14C).reshape((-1, nr_pools))

            smrfd_14C = smrfd_12C.to_14C_only(
                xs_14C_data[0], us_14C_data, parameter_set["decay_rate"]
            )
            smrfd_14C.location = location
            smrfd_14C.save_to_file(filename_14C)
            print("saved", filename_14C, flush=True)
        except Exception as e:
            print(type(e), e, flush=True)
            log = str(e)
            raise (e)

    return (location["cell_nr"], log, filename_14C)


def solve_SMRFD(cell_nr, filename):
    soln = None
    log = None

    try:
        smrfd = SMRFD.load_from_file(filename)
        soln = smrfd.solve()
        soln = soln.ravel().tolist()
    except Exception as e:
        print(e, flush=True)
        log = str(e)

    return (cell_nr, log, soln)


def load_model_12C_data(parameter_set, location):
    xs_12C = None
    us_12C = None
    rs_12C = None

    log = None

    lat_index = location["lat_index"]
    lon_index = location["lon_index"]
    try:
        mdo_12C = load_mdo_12C(parameter_set)

        xs_12C = (
            mdo_12C.load_stocks(
                func=getStockVariable_from_Density,
                lat_arr=[lat_index],
                lon_arr=[lon_index],
                data_shift=0,
            )
            .data[:, 0, 0, ...]
            .ravel()
            .tolist()
        )

        us_12C = (
            mdo_12C.load_external_input_fluxes(
                func=getFluxVariable_from_DensityRate,
                lat_arr=[lat_index],
                lon_arr=[lon_index],
                data_shift=1,
            )
            .data[:, 0, 0, ...]
            .ravel()
            .tolist()
        )

        rs_12C = (
            mdo_12C.load_external_output_fluxes(
                func=getFluxVariable_from_DensityRate,
                lat_arr=[lat_index],
                lon_arr=[lon_index],
                data_shift=1,
            )
            .data[:, 0, 0, ...]
            .ravel()
            .tolist()
        )
    except Exception as e:
        print(e, flush=True)
        log = str(e)
        raise (e)

    return (location["cell_nr"], log, xs_12C, us_12C, rs_12C)


def load_model_14C_data(parameter_set, location):
    xs_14C = None
    us_14C = None
    rs_14C = None

    log = None

    lat_index = location["lat_index"]
    lon_index = location["lon_index"]
    try:
        mdo_14C = load_mdo_14C(parameter_set)

        xs_14C = (
            mdo_14C.load_stocks(
                func=getStockVariable_from_Density,
                lat_arr=[lat_index],
                lon_arr=[lon_index],
                data_shift=0,
            )
            .data[:, 0, 0, ...]
            .ravel()
            .tolist()
        )

        us_14C = (
            mdo_14C.load_external_input_fluxes(
                func=getFluxVariable_from_DensityRate,
                lat_arr=[lat_index],
                lon_arr=[lon_index],
                data_shift=1,
            )
            .data[:, 0, 0, ...]
            .ravel()
            .tolist()
        )

        rs_14C = (
            mdo_14C.load_external_output_fluxes(
                func=getFluxVariable_from_DensityRate,
                lat_arr=[lat_index],
                lon_arr=[lon_index],
                data_shift=1,
            )
            .data[:, 0, 0, ...]
            .ravel()
            .tolist()
        )
    except Exception as e:
        print(e, flush=True)
        log = str(e)

    return (location["cell_nr"], log, xs_14C, us_14C, rs_14C)


def load_parameter_set(**keywords):
    ds_filename = keywords.get("ds_filename")
    nstep = keywords.get("nstep", 1)
    min_time = keywords.get("min_time", 0)
    max_time = keywords.get("max_time", None)
    cftime_unit_src = "days since 1901-01-01 00:00:00"
    cftime_unit_src = keywords.get("cftime_unit_src", cftime_unit_src)
    time_shift = keywords.get("time_shift", 0)
    calendar_src = keywords.get("calendar", "noleap")

    parameter_set = {
        "ds_filename": ds_filename,
        "stock_unit": "g/m^2",
        "nr_layers": 10,
        "dz_var_name": "dz",
        "min_time": min_time,  # actually the minimum time INDEX
        "max_time": max_time,  # before time aggregation
        "nstep": nstep,
        "time_shift": time_shift,
        "cftime_unit_src": cftime_unit_src,
        "decay_rate": np.log(2) / 5568.0 / 365.0,  # daily value
        "calendar_src": calendar_src,
    }

    ## depth variable with bounds
    ds_depth = Dataset("../Data/DZSOI.nc", "r")
    dz_var = ds_depth.variables["DZSOI"]
    dz = Variable(
        name=parameter_set["dz_var_name"],
        data=dz_var[: parameter_set["nr_layers"], 0],
        unit=dz_var.units,
    )
    ds_depth.close()
    parameter_set["dz_var_names"] = {"dz": dz}

    #    dataset = Dataset(parameter_set['ds_filename'], 'r')
    #    parameter_set['dataset']   = dataset
    #    parameter_set['calendar_src']  = dataset['time'].calendar

    return parameter_set


################################################################################


if __name__ == "__main__":
    pass
