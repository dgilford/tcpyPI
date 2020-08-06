import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyPI",
    version="1.3",
    author="Daniel M. Gilford",
    author_email="daniel.gilford@rutgers.edu",
    description="pyPI: Tropical cyclone potential intensity calculations in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dgilford/pyPI",
    packages=setuptools.find_packages(),
    license="MIT LICENSE",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5',
)