"""
ncOrtho - Targeted ortholog search for miRNAs
Copyright (C) 2021 Felix Langschied

ncOrtho is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ncOrtho is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ncOrtho.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
from setuptools import setup, find_packages


here = os.path.abspath(os.path.dirname(__file__))
readme_path = os.path.join(here, "README.md")
if os.path.isfile(readme_path):
    with open(readme_path, encoding="utf-8") as fh:
        long_description = fh.read()
else:
    long_description = ""

setup(
    name="ncOrtho",
    version="1.0.0",
    python_requires=">=3.10.0",
    description="Targeted ortholog search for miRNAs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Felix Langschied",
    author_email="langschied@bio.uni-frankfurt.de",
    url="https://github.com/felixlangschied/ncortho",
    packages=find_packages(),
    package_data={"": ["*"]},
    install_requires=[
        "PyYAML>=5.0",
        "biopython>=1.79",
        "pyfaidx>=0.6.0",
        "pyfiglet",
        "tqdm>=4.0",
    ],
    entry_points={
        "console_scripts": [
            "ncCreate = ncOrtho.coreset.coreset:main",
            "ncSearch = ncOrtho.ncortho:main",
            "ncAnalyze = ncOrtho.analysis.ncAnalyze:main",
            "ncCheck = ncOrtho.check.core_synteny:main",
        ],
    },
    license="GPL-3.0",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)