from setuptools import find_packages, setup

setup(name = 'seqdatasets', packages = find_packages())

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

requirements = []

setup(
    name="seqdatasets",
    version="0.1.2",
    author="Adam Klie",
    author_email="aklie@ucsd.edu",
    description="Datasets for benchmarking, testing and developing in EUGENe",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/adamklie/SeqDatasets",
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)
