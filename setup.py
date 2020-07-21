import setuptools


INSTALL_REQUIRES=[
        "biopython>=1.7", 
        "pandas>=0.25",
        "joblib>=0.13",
        "scikit-learn==0.23.1"
        ]

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("bacphlip/__init__.py", "r") as f:
    init = f.readlines()
    for line in init:
        if "__version__" in line:
            __version__ = line.split('"')[-2]
            #__version__ = line.split('=')[-1].strip('"')

setuptools.setup(
    name="bacphlip",
    packages=['bacphlip'],
    entry_points={ 'console_scripts': ['Package = bacphlip.__main__:main', 'bacphlip=bacphlip.command_line:main']}, 
    version=__version__,
    author="Adam J Hockenberry",
    author_email="adam.hockenberry@utexas.edu",
    description="A Random Forest classifier to predict bacteriophage lifestyle",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=INSTALL_REQUIRES,
    include_package_data=True,
    url="https://github.com/adamhockenberry/bacphlip",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
    python_requires='>=3.6',
)
