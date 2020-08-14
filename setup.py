#! /usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tcpypi",
    # TODO: Update every time the version updates
    version="1.3",
    author="Daniel M. Gilford, PhD",
    author_email="daniel.gilford@rutgers.edu",
    description="pyPI: Tropical cyclone potential intensity calculations in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dgilford/pyPI",
    download_url="https://github.com/darothen/pyrcel",
    packages=setuptools.find_packages(),
    package_data={"pyrcel": ["data/sample_data.nc"]},
    license="MIT LICENSE",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    # TODO: Update install requirements and corresponding documentation
    install_requires=[
        "xarray==0.15.1",
        "numba==0.48.0",
        "numpy==1.18.1",
    ],
)



