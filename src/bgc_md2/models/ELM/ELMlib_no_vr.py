import numpy as np
from pathlib import Path

from bgc_md2.ModelStructure import ModelStructure
from bgc_md2.ModelDataObject import ModelDataObject
from bgc_md2.Variable import Variable


DATA_PATH = Path("../../../../SOIL-R/Manuscripts/Berkeley/2019/Data/Holger_veg_soil/no_vr/")


def load_model_structure_with_vegetation():
    # list of all pools
    pool_structure = [
        # vegetation component
        {"pool_name": "Leaves", "stock_var": "LEAFC_TOT"},
        {"pool_name": "Live stem", "stock_var": "LIVESTEMC_TOT"},
        {"pool_name": "Dead stem", "stock_var": "DEADSTEMC_TOT"},
        {"pool_name": "Fine roots", "stock_var": "FROOTC_TOT"},
        {"pool_name": "Live coarse roots", "stock_var": "LIVECROOTC_TOT"},
        {"pool_name": "Dead coarse roots", "stock_var": "DEADCROOTC_TOT"},      
        # soil component
        {"pool_name": "CWD", "stock_var": "CWDC"},
        {"pool_name": "Litter metabolic", "stock_var": "LITR1C"},
        {"pool_name": "Litter cellulose", "stock_var": "LITR2C"},
        {"pool_name": "Litter lignin", "stock_var": "LITR3C"},
        {"pool_name": "Soil 1", "stock_var": "SOIL1C"},
        {"pool_name": "Soil 2", "stock_var": "SOIL2C"},
        {"pool_name": "Soil 3", "stock_var": "SOIL3C"}
    ]

    # fluxes from outside the model into vegetation
    external_input_structure = {
        # pool_name: list of associated flux variables
        "Leaves": ["PSN_TO_LEAFC_TOT"],
        "Live stem": ["PSN_TO_LIVESTEMC_TOT"],
        "Dead stem": ["PSN_TO_DEADSTEMC_TOT"],
        "Fine roots": ["PSN_TO_FROOTC_TOT"],
        "Live coarse roots": ["PSN_TO_LIVECROOTC_TOT"],
        "Dead coarse roots": ["PSN_TO_DEADCROOTC_TOT"]
    }

    # internal horizontal fluxes
    horizontal_structure = {
        # (pool_from, pool_to): list of associated flux variables

        # vegetation component
        ("Live stem", "Dead stem"): ["LIVESTEMC_TO_DEADSTEMC"],
        ("Live coarse roots", "Dead coarse roots"): ["LIVECROOTC_TO_DEADCROOTC"],

        # vegetation component to soil component
        ("Leaves", "Litter metabolic"): ["LEAFC_TO_LITTER_MET"],
        ("Leaves", "Litter cellulose"): ["LEAFC_TO_LITTER_CEL"],
        ("Leaves", "Litter lignin"): ["LEAFC_TO_LITTER_LIG"],

        ("Live stem", "Litter metabolic"): ["LIVESTEMC_TO_LITTER_MET"],
        ("Live stem", "CWD"): ["LIVESTEMC_TO_CWDC"],

        ("Dead stem", "Litter metabolic"): ["DEADSTEMC_TO_LITTER_MET"],
        ("Dead stem", "CWD"): ["DEADSTEMC_TO_CWDC"],

        ("Fine roots", "Litter metabolic"): ["FROOTC_TO_LITTER_MET"],
        ("Fine roots", "Litter cellulose"): ["FROOTC_TO_LITTER_CEL"],
        ("Fine roots", "Litter lignin"): ["FROOTC_TO_LITTER_LIG"],

        ("Live coarse roots", "Litter metabolic"): ["LIVECROOTC_TO_LITTER_MET"],
        ("Live coarse roots", "CWD"): ["LIVECROOTC_TO_CWDC"],

        ("Dead coarse roots", "Litter metabolic"): ["DEADCROOTC_TO_LITTER_MET"],
        ("Dead coarse roots", "CWD"): ["DEADCROOTC_TO_CWDC"],

        # soil component
        ("CWD", "Litter cellulose"): ["CWDC_TO_LITR2C"],
        ("CWD", "Litter lignin"): ["CWDC_TO_LITR3C"],
        ("Litter metabolic", "Soil 1"): ["LITR1C_TO_SOIL1C"],
        ("Litter cellulose", "Soil 1"): ["LITR2C_TO_SOIL1C"],
        ("Litter lignin", "Soil 2"): ["LITR3C_TO_SOIL2C"],
        ("Soil 1", "Soil 2"): ["SOIL1C_TO_SOIL2C"],
        ("Soil 1", "Soil 3"): ["SOIL1C_TO_SOIL3C"],
        ("Soil 2", "Soil 1"): ["SOIL2C_TO_SOIL1C"],
        ("Soil 2", "Soil 3"): ["SOIL2C_TO_SOIL3C"],
        ("Soil 3", "Soil 1"): ["SOIL3C_TO_SOIL1C"],
    }

    # internal vertical fluxes, soil only
    vertical_structure = {}

    # fluxes leaving the model
    external_output_structure = {
        # pool_name: list of fluxes leaving the model

        # vegetation component
        "Leaves": ["LEAF_MR", "LEAF_GR"],
        "Live stem": ["LIVESTEM_MR", "LIVESTEM_GR"],
        "Dead stem": ["DEADSTEM_GR"],
        "Fine roots": ["FROOT_MR", "FROOT_GR"],
        "Live coarse roots": ["LIVECROOT_MR", "LIVECROOT_GR"],
        "Dead coarse roots": ["DEADCROOT_GR"],

        # soil component
        "CWD": ["CWD_HR_L2", "CWD_HR_L3", "M_CWDC_TO_FIRE"],
        "Litter metabolic": ["LITR1_HR", "M_LITR1C_TO_FIRE"],
        "Litter cellulose": ["LITR2_HR", "M_LITR2C_TO_FIRE"],
        "Litter lignin": ["LITR3_HR", "M_LITR3C_TO_FIRE"],
        "Soil 1": ["SOIL1_HR_S2", "SOIL1_HR_S3"],
        "Soil 2": ["SOIL2_HR_S1", "SOIL2_HR_S3"],
        "Soil 3": ["SOIL3_HR"],
    }

    model_structure = ModelStructure(
        pool_structure=pool_structure,
        external_input_structure=external_input_structure,
        horizontal_structure=horizontal_structure,
        vertical_structure=vertical_structure,
        external_output_structure=external_output_structure,
    )

    return model_structure


def check_data_consistency(ds_single_site, time_step_in_days):
    ms = load_model_structure_with_vegetation()

    # data_test.py says tht for labile 31.0 is constantly right
    time = Variable(
        name="time",
        data=np.arange(len(ds_single_site.time)) * time_step_in_days,
        unit="d"
    )

    mdo = ModelDataObject(
        model_structure=ms,
        dataset=ds_single_site, 
        stock_unit="gC/m2", 
        time=time
    )

    abs_err, rel_err = mdo.check_data_consistency()
    return abs_err, rel_err

