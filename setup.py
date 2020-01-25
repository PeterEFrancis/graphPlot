from setuptools import setup
from os import path

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="graph_plot",
    packages=["graph_plot"],
    version="0.0.1",
    author="Peter Francis",
    author_email="franpe02@gettysburg.edu",
    description="For plotting Graphs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/francisp336/graph",
    classifiers=[],
    python_requires=">=3.0.0",
    install_requires=["matplotlib", "numpy"],
    keywords="Graph",
)
