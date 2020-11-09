# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
#spark.stop()
from pyspark.sql import SparkSession
from pyspark.conf import SparkConf 
spark = SparkSession.builder \
.master("local[4]") \
.appName('ELM') \
.config(conf = SparkConf()) \
.getOrCreate()

from pyspark.sql.functions import udf, col
from pyspark.sql.types import StructType, StringType, ArrayType, IntegerType, FloatType, StructField

import importlib
import ELMlib
importlib.reload(ELMlib)

import pandas as pd
import numpy as np

from SmoothModelRunFromData import SmoothModelRunFromData as SMRFD
# -

#spinup_filestem = '14C_spinup_holger.2x2_2_small'
spinup_filestem = '14C_spinup_holger_fire.2x2_small'
df = ELMlib.initialize_spark_df_from_nc('../Data/' + spinup_filestem + '.nc')
df.show(4)

parameter_set = ELMlib.load_parameter_set(
    ds_filename = '../Data/' + spinup_filestem + '.nc',
    time_shift  = -198*365,
    nstep       = 10
)

# +
df2_rdd = df.rdd.map(
    lambda x: ELMlib.load_model_12C_data(
        parameter_set,
        {'cell_nr': x.cell_nr, 'lat_index': x.lat_index, 'lon_index': x.lon_index}
    )
)
schema = StructType([
    StructField('cell_nr', IntegerType(), True),
    StructField('model_12C_log', StringType(), True),
    StructField('xs_12C', ArrayType(FloatType()), True),
    StructField('us_12C', ArrayType(FloatType()), True),
    StructField('rs_12C', ArrayType(FloatType()), True)
])
df2 = spark.createDataFrame(df2_rdd, schema)
df = df.join(df2, on='cell_nr')
print('Columns:', df.columns)
print(
    df.filter(~df.model_12C_log.isNull()).count(),
    "sites could not load the model's 12C stocks, inputs, and outputs."
)

df.show()


# +
df = df.filter(df.model_12C_log.isNull()) \
.drop('model_12C_log')
print('Columns:', df.columns)

df.show()


# +
filestem = '../Data/output/%s/%02d/' % (spinup_filestem, parameter_set['nstep'])
filestem += 'smrfd_12C_%05d.gzip'

#@udf(StructType([StructField('filename', StringType(), True), StructField('log', StringType(), True)]))
#def save_SMRFD_12C_(cell_nr, lat_index, lon_index):
##    return (cell_nr, (lat_index, lon_index))
#    
#    return ELMlib.save_SMRFD_12C(
#            parameter_set,
#            dump_file_stem,
#            location = {'cell_nr': cell_nr, 'lat_index': lat_index, 'lon_index': lon_index}
#    )
#    return (cell_nr, file_name, log)
#
#df = df.withColumn('SMRFD', save_SMRFD_12C_('cell_nr', 'lat_index', 'lon_index')) \
#.withColumn('SMRFD_filename', col('SMRFD').filename) \
#.withColumn('SMRFD_log', col('SMRFD').log) \
#.drop('SMRFD')

df2_rdd = df.rdd.map(
    lambda x: ELMlib.save_SMRFD_12C(
        parameter_set,
        filestem % x.cell_nr,
        {'cell_nr': x.cell_nr, 'lat_index': x.lat_index, 'lon_index': x.lon_index}
    )
)
schema = StructType([
    StructField('cell_nr', IntegerType(), True),
    StructField('smrfd_12C_log', StringType(), True),
    StructField('smrfd_12C_filename', StringType(), True)
])
df2 = spark.createDataFrame(df2_rdd, schema)

df = df.join(df2, on='cell_nr')
print('Columns:', df.columns)
print(df.filter(~df.smrfd_12C_log.isNull()).count(), 'sites could not create a smrfd_12C.')

df.show()

# +
df = df.filter(df.smrfd_12C_log.isNull()) \
.drop('smrfd_12C_log')
print('Columns:', df.columns)

df.show()

# +
df2_rdd = df.rdd.map(lambda x: ELMlib.solve_SMRFD(x.cell_nr, x.smrfd_12C_filename))
schema = StructType([
    StructField('cell_nr', IntegerType(), True),
    StructField('soln_12C_log', StringType(), True),
    StructField('soln_12C', ArrayType(FloatType()), True)
])
df2 = spark.createDataFrame(df2_rdd, schema)
df = df.join(df2, on='cell_nr')
print('Columns:', df.columns)
print(df.filter(~df.soln_12C_log.isNull()).count(), 'sites could not solve the smrfd_12C.')

df.show()


# +
df = df.filter(df.soln_12C_log.isNull()) \
.drop('soln_12C_log')
print('Columns:', df.columns)

df.show()

# +
df2_rdd = df.rdd.map(
    lambda x: ELMlib.load_model_14C_data(
        parameter_set,
        {'cell_nr': x.cell_nr, 'lat_index': x.lat_index, 'lon_index': x.lon_index}
    )
)
schema = StructType([
    StructField('cell_nr', IntegerType(), True),
    StructField('model_14C_log', StringType(), True),
    StructField('xs_14C', ArrayType(FloatType()), True),
    StructField('us_14C', ArrayType(FloatType()), True),
    StructField('rs_14C', ArrayType(FloatType()), True)
])
df2 = spark.createDataFrame(df2_rdd, schema)
df = df.join(df2, on='cell_nr')
print('Columns:', df.columns)
print(
    df.filter(~df.model_14C_log.isNull()).count(),
    "sites could not load the model's 14C stocks, inputs, and outputs."
)

df.show()

# +
df = df.filter(df.model_14C_log.isNull()) \
.drop('model_14C_log')
print('Columns:', df.columns)

df.show()

# +
filestem = '../Data/output/%s/%02d/' % (spinup_filestem, parameter_set['nstep'])
filestem += 'smrfd_14C_%05d.gzip'

df2_rdd = df.rdd.map(
    lambda x: ELMlib.save_SMRFD_14C(
        parameter_set,
        x.smrfd_12C_filename,
        filestem % x.cell_nr,
        {'cell_nr': x.cell_nr, 'lat_index': x.lat_index, 'lon_index': x.lon_index},
        x.xs_14C,
        x.us_14C
    )
)
schema = StructType([
    StructField('cell_nr', IntegerType(), True),
    StructField('smrfd_14C_log', StringType(), True),
    StructField('smrfd_14C_filename', StringType(), True)
])
df2 = spark.createDataFrame(df2_rdd, schema)

df = df.join(df2, on='cell_nr')
print('Columns:', df.columns)
print(df.filter(~df.smrfd_14C_log.isNull()).count(), 'sites could not create a smrfd_14C.')

df.show()


# +
df = df.filter(df.smrfd_14C_log.isNull()) \
.drop('smrfd_14C_log')
print('Columns:', df.columns)

df.show()

# +
df2_rdd = df.rdd.map(lambda x: ELMlib.solve_SMRFD(x.cell_nr, x.smrfd_14C_filename))
schema = StructType([
    StructField('cell_nr', IntegerType(), True),
    StructField('soln_14C_log', StringType(), True),
    StructField('soln_14C', ArrayType(FloatType()), True)
])
df2 = spark.createDataFrame(df2_rdd, schema)
df = df.join(df2, on='cell_nr')
print('Columns:', df.columns)
print(df.filter(~df.soln_14C_log.isNull()).count(), 'sites could not solve the smrfd_14C.')

df.show()

# +
df = df.filter(df.soln_14C_log.isNull()) \
.drop('soln_14C_log')
print('Columns:', df.columns)

df.show()

# +
df.write.json('df', mode='overwrite', compression='gzip')

df.show()
# -

df2 = spark.read.json('df')

df.count()

df2


