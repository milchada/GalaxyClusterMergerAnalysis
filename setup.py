import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gcmerge",
    version="0.1",
    author="Urmila Chadayammuri",
    author_email="uchadaya@gmail.com",
    description="Multi-wavelength comparisons of galaxy cluster simulations and observations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/milchada/gcmerge",
    packages=setuptools.find_packages(),
    classifiers=["Intended Audience :: Science/Research",
                 "License :: OSI Approved :: BSD License",
                 "Operating System :: MacOS :: MacOS X",
                 "Operating System :: POSIX :: AIX",
                 "Operating System :: POSIX :: Linux",
                 "Programming Language :: C",
                 "Programming Language :: Python :: 2",
                 "Programming Language :: Python :: 2.7",
                 "Programming Language :: Python :: 3",
                 "Programming Language :: Python :: 3.4",
                 "Programming Language :: Python :: 3.5",
                 "Topic :: Scientific/Engineering :: Astronomy",
                 "Topic :: Scientific/Engineering :: Physics",
                 "Topic :: Scientific/Engineering :: Visualization"],
)
