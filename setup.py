import setuptools


INSTALL_REQUIRES=[
        "biopython>=1.7", 
        "pandas>=0.25",
        "joblib>=0.13",
        "scikit-learn>=0.21"
        ]

#TEST_REQUIRES = [
#    # testing and coverage
#    "pytest",
#    "coverage",
#    "pytest-cov",
#    # to be able to run `python setup.py checkdocs`
#    "collective.checkdocs",
#    "pygments",
#]

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("bacphlip/__init__.py", "r") as f:
    init = f.readlines()
for line in init:
    if "__version__" in line:
        __version__ = line.split('=')[-1]

setuptools.setup(
    name="bacphlip",
    version=__version__,
    author="Adam J Hockenberry",
    author_email="adam.hockenberry@utexas.edu",
    description="A Random Forest classifier to predict bacteriophage lifestyle",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=INSTALL_REQUIRES,
    url="https://github.com/adamhockenberry/bacphlip-py",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
    python_requires='>=3.6',
)


"""
setup(
    download_url="https://github.com/clauswilke/PeptideBuilder/releases",
    platforms="Tested on Mac OS X and Windows 10",
    packages=["PeptideBuilder"],
    extras_require={"test": TEST_REQUIRES + INSTALL_REQUIRES,},
)

"""
