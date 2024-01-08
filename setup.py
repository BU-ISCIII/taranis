#!/usr/bin/env python

from setuptools import setup, find_packages

version = "3.0.0"

with open("README.md") as f:
    readme = f.read()

with open("requirements.txt") as f:
    required = f.read().splitlines()

setup(
    name="taranis",
    version=version,
    description="Tools for gene-by-gene allele calling analysis",
    long_description=readme,
    long_description_content_type="text/markdown",
    keywords=[
        "buisciii",
        "bioinformatics",
        "pipeline",
        "sequencing",
        "NGS",
        "next generation sequencing",
    ],
    author="Sara Monzon",
    author_email="smonzon@isciii.es",
    url="https://github.com/BU-ISCIII/taranis",
    license="GNU GENERAL PUBLIC LICENSE v.3",
    entry_points={"console_scripts": ["taranis=taranis.__main__:run_taranis"]},
    install_requires=required,
    packages=find_packages(exclude=("docs")),
    include_package_data=True,
    zip_safe=False,
)
