import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

from bgc_md2.ModelDataObject import (
    ModelDataObject,
    getStockVariable_from_Density,
    getFluxVariable_from_DensityRate,
)
from bgc_md2.ModelStructure import ModelStructure
from CompartmentalSystems.pwc_model_run_fd import (
    PWCModelRunFD,
    PWCModelRunFDError
)
from bgc_md2.Variable import Variable


################################################################################


def load_parameter_set(**keywords):
    nstep = keywords.get("nstep", 1)

    parameter_set = {
        "stock_unit": "g/m^2",
        "nr_layers": 10,
        "nstep": nstep,
    }

    ## depth variable with bounds
    ds_depth = keywords['ds_depth']
    dz_var = ds_depth.variables["DZSOI"]
    dz_var_name = "dz"
    dz = Variable(
        name=dz_var_name,
        data=dz_var[:parameter_set["nr_layers"], 0],
        unit=dz_var.attrs['units'],
    )
    ds_depth.close()
    
    parameter_set["dz_var_name"] = dz_var_name
    parameter_set["dz_var_names"] = {dz_var_name: dz}

    #    dataset = Dataset(parameter_set['ds_filename'], 'r')
    #    parameter_set['dataset']   = dataset
    #    parameter_set['calendar_src']  = dataset['time'].calendar

    return parameter_set


def load_model_structure(nr_layers, dz_var_name):
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


def load_mdo(ds, parameter_set):
    nstep = parameter_set.get("nstep", 1)
    stock_unit = parameter_set["stock_unit"]
    nr_layers = parameter_set["nr_layers"]
    dz_var_name = parameter_set["dz_var_name"]
    dz_var_names = parameter_set["dz_var_names"]

    ms = load_model_structure(nr_layers, dz_var_name)

    time = Variable(
        name="time",
        data=np.arange(len(ds.time)),
        unit="d"
    )

    mdo = ModelDataObject(
        model_structure=ms,
        dataset=ds,
        nstep=nstep,
        stock_unit=stock_unit,
        time=time,
        dz_var_names=dz_var_names,
    )

    return mdo


def compute_ds_pwc_mr_fd(ds, parameter_set):
#    if ds.ens.values.shape == (0,):

    mdo = load_mdo(ds, parameter_set)
    ms = mdo.model_structure
    nr_pools = ms.nr_pools

    error = ''
    try:
        #pwc_mr_fd, err_dict = mdo.create_model_run(errors=True)
        pwc_mr_fd = mdo.create_model_run()
    except PWCModelRunFDError as e:
        error = str(e)

    coords_time = ds.time
    coords_pool = np.arange(nr_pools)

    data_vars = dict()
    data = np.nan * np.ones((nr_pools, ))
    if not error:
        data = pwc_mr_fd.start_values
    data_vars['start_values'] = xr.DataArray(
        data=data,
        coords={'pool': coords_pool},
        dims=['pool'],
        attrs={'units': mdo.stock_unit}
    )

    data = np.nan * np.ones((len(ds.time),))
    if not error:
        data = pwc_mr_fd.times
    data_vars['times'] = xr.DataArray(
        data=data,
        coords={'time': coords_time},
        dims=['time'],
        attrs={'units': mdo.time_agg.unit}
    )

    data = np.nan * np.ones((len(ds.time), nr_pools))
    if not error:
        data[:-1, ...] = pwc_mr_fd.us
    data_vars['us'] = xr.DataArray(
        data=data,
        coords={'time': coords_time, 'pool': coords_pool},
        dims=['time', 'pool'],
        attrs={'units': mdo.stock_unit+'/'+mdo.time_agg.unit}
    )

    data = np.nan * np.ones((len(ds.time), nr_pools, nr_pools))
    if not error:
        data[:-1, ...] = pwc_mr_fd.Bs
    data_vars['Bs'] = xr.DataArray(
        data=data,
        coords={
            'time': coords_time,
            'pool_to': coords_pool,
            'pool_from': coords_pool
        },
        dims=['time', 'pool_to', 'pool_from'],
        attrs={'units': '1/'+mdo.time_agg.unit}
    )

    data = np.array(error, dtype="<U150")
    data_vars['log'] = xr.DataArray(data=data)

    coords = {
        'time': coords_time,
        'pool': coords_pool,
        'pool_to': coords_pool,
        'pool_from': coords_pool
    }
    ds_res = xr.Dataset(
        coords=coords,
        data_vars=data_vars,
    )
    ds_res.close()

    return ds_res


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


def resample_daily_to_monthly(ds):
    ms = load_model_structure(nr_layers=10, dz_var_name='dz')


    data_vars = dict()

    stock_vars = [d['stock_var'] for d in ms.pool_structure]
    for stock_var in stock_vars:
        data_vars[stock_var] = ds[stock_var].resample(
            time='1M',
            keep_attrs=True
        ).first()

    time_resampled = ds.coords['time'].resample(time='1M', label='left').first().data
    print(time_resampled)
    coords_resampled = {
        'time': time_resampled,
        'lat': ds.coords['lat'],
        'lon': ds.coords['lon']
    }

    ds_resampled = xr.Dataset(
#        coords=coords_resampled,
        data_vars=data_vars
    )
    ds_resampled.coords['time'] = time_resampled
    print(type(time_resampled))
    return ds_resampled


################################################################################


if __name__ == "__main__":
    ELMDataDir = "/home/hmetzler/SOIL-R/Manuscripts/Berkeley/2019/Data/"
    runID = "14C_transient_holger_fire.2x2_small"
#    runID = "holger.new3sites.tr_2019_grid"
    fn = runID + ".nc"

    ds = xr.open_dataset(Path(ELMDataDir).joinpath(runID + ".nc"))
    print('A')
    ds_depth = xr.open_dataset(Path(ELMDataDir).joinpath('DZSOI.nc'))
    print('B')
    parameter_set = load_parameter_set(
        nstep       = 365,
        ds_depth    = ds_depth
    )
    print('C')
    ds_single_site = ds.isel(lat=0, lon=0)
    print('D')
    ds_mr_pwc_fd = compute_ds_pwc_mr_fd(ds_single_site, parameter_set)

