#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name="topoflow",
    version="3.6",
    description="d8-based, spatial hydrologic model",
    author="Scott D. Peckham",
    author_email="Scott.Peckham@colorado.edu",
    license="MIT",
    url="https://github.com/peckhams/topoflow36",
    include_package_data=True,
    packages=find_packages("."),
    install_requires=["numpy", "scipy", "h5py", "netCDF4", "cfunits"],
)
