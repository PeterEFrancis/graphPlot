from setuptools import setup
from os import path

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="graphPlot",
    packages=["graphPlot"],
    version="1.0.0",
    author="Peter E. Francis",
    author_email="francisp336@gmail.com",
    description="Plot (di)graphs using a timestep simulation of charged particles and springs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/francisp336/graphPlot",
    classifiers=[],
    python_requires=">=3.0.0",
    install_requires=["matplotlib", "numpy"],
    keywords="Graph",
)
