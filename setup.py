from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup


__version__ = "0.1.0"


ext_modules = [
    Pybind11Extension(
        "potsim2.skinlib",
        ["potsim2/src/skin.cpp"],
        # Example: passing in the version to the compiled code
        define_macros=[("VERSION_INFO", __version__)],
    ),
]

setup(
    name="potsim2",
    version=__version__,
    author="Eloy Felix",
    author_email="eloyfelix@gmail.com",
    url="https://github.com/chembl/potsim2",
    license="MIT",
    packages=["potsim2"],
    description="",
    long_description=open("README.md", encoding="utf-8").read(),
    ext_modules=ext_modules,
    install_requires=["pybind11>=2.4", "MDAnalysis~=2.0.0", "GridDataFormats~=0.6.0"],
    setup_requires=["pybind11>=2.4"],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
