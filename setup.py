from setuptools import setup, find_packages
import versioneer

with open("README.md", "r") as fh:
    long_description = fh.read()
with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="rna_draw",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    packages=["rna_draw"],
    py_modules=[
        "rna_draw/colorer",
        "rna_draw/data",
        "rna_draw/draw",
        "rna_draw/parameters",
        "rna_draw/render_rna",
    ],
    include_package_data=True,
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    entry_points={
        "console_scripts": [
            "rna_draw = rna_draw.draw:main",
        ]
    },
    python_requires=">=3.6",
)
