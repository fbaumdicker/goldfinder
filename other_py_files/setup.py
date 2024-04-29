from setuptools import setup

setup(
    name="goldfinder",
    version="0.1",
    description="a useful module",
    author="lord",
    author_email="christian-resl@gmx.at",
    packages=['goldfinder'],
    install_requires=["ete3", "tqdm"]
)