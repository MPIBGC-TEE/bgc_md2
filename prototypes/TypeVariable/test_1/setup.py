#!/usr/bin/env python3
# vim:set ff=unix expandtab ts=4 sw=4:

from setuptools import setup, find_packages


def readme():
    with open("README.md") as f:
        return f.read()



packages=find_packages('src'),  # find all packages (multifile modules) recursively
print(packages)

setup(
    name="trendy9helpers",
    version="0.1",
    description="Data assimilation for many models with type hints",
    long_description=readme(),  # avoid duplication
    author="Markus",
    author_email="markus.mueller.1.g@gmail.com",
    package_dir={"": "trendy9helpers"},
    #packages=find_packages('src'),  # find all packages (multifile modules) recursively
    #package_dir={'': 'src'},
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
    install_requires=[
        "numpy",
    ] 
)
