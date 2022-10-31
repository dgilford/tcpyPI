#! /usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tcpypi",
    # TODO: Update every time the version updates
    version="1.3.5",
    author="Daniel M. Gilford, PhD",
    author_email="dgilford@climatecentral.org",
    description="tcpyPI: Tropical cyclone potential intensity calculations in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dgilford/tcpyPI",
    download_url="https://github.com/dgilford/tcpyPI",
    packages=setuptools.find_packages(),
    package_data={"tcpyPI": ["data/sample_data.nc"]},
    license="MIT LICENSE",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    # TODO: Update install requirements and corresponding documentation
    install_requires=[
        "xarray>=0.16.2",
        "numba>=0.51.2",
        "numpy>=1.19.5",
    ],
)
