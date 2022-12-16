#!/usr/bin/env python3
# vim:set ff=unix expandtab ts=4 sw=4:

# This is the standard pyhon install script.
# Only if it does not work use the shell script ../install_bgc_md.sh in the folder above
# Keep the requirements.freeze clean and to a minimum.

from setuptools import setup, find_packages


def readme():
    with open("README.md") as f:
        return f.read()

packages=find_packages('src'),  # find all packages (multifile modules) recursively
print(packages)

setup(
    name="bgc_md2",
    version="0.1",
    test_suite="example_package.tests",  # http://pythonhosted.org/setuptools/setuptools.html#test
    description="Model Database for Carbon Allocation",
    long_description=readme(),  # avoid duplication
    author="MarkusHolger, Veronika, Thomas,  ",
    author_email="markus.mueller.1.g@gmail.com",
    url="https://github.com/MPIBGC-TEE/bgc_md2",
    packages=find_packages('src'),  # find all packages (multifile modules) recursively
    package_dir={'': 'src'},
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: POSIX :: Linux",
        "Topic :: Education ",
    ],
    # entry_points={
    #'console_scripts': [
    #        'render= bgc_md.reports:render_parse'
    #        ]
    # },
    install_requires=[
        "testinfrastructure",   
        "ComputabilityGraphs",
        "LAPM",
        "CompartmentalSystems",
        "netCDF4",
        "sympy",
        "xarray",
        "networkx",
        "IPython",
        #"pygraphviz",
        "cf-units",
        "frozendict",
        "nbformat", # for testing notebooks
    ]
)
