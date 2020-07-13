#!/usr/bin/env python3
# vim:set ff=unix expandtab ts=4 sw=4:

# This is the standard pyhon install script.
# Only if it does not work use the shell script ../install_bgc_md.sh in the folder above
# Keep the requirements.freeze clean and to a minimum.

from setuptools import setup,find_packages
def readme():
    with open('README.md') as f:
        return f.read()
    
setup(name='bgc_md2',
        version='0.1',
        test_suite="example_package.tests",#http://pythonhosted.org/setuptools/setuptools.html#test
        description='Model Database for Carbon Allocation',
        long_description=readme(),#avoid duplication 
        author='Markus, Thomas, Holger, Veronika ',
        author_email='markus.mueller.1.g@gmail.com',
        url='https://github.com/MPIBGC-TEE/bgc_md2',
        packages=find_packages(), #find all packages (multifile modules) recursively
        #py_modules=['external_module'], # external_module.py is a module living outside this dir as an example how to iclude something not 
        classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: POSIX :: Linux",
        "Topic :: Education "
        ],
        #entry_points={
        #'console_scripts': [
        #        'render= bgc_md.reports:render_parse'
        #        ]
        #},
        install_requires=[
            "CompartmentalSystems"
            ,"testinfrastructure" #also on https://github.com/MPIBGC-TEE/testinfrastructure.githttps://github.com/MPIBGC-TEE/testinfrastructure.git
            ,"concurrencytest"
            #,"requests"
            #,"mendeley "
            #,"PyYAML"
            #,"pandas"
            ,'netCDF4'
            #,'sqlalchemy'
            #,'oslash'
            #,'pypandoc'
            ,"sympy"
            ,"xarray"
            ,"networkx"
            ,"pygraphviz"
            ,'cf-units'
        ],
        include_package_data=True,
        #zip_safe=False
     )

 
 
 
