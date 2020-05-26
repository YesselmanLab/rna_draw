from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='rna_draw',
    version='0.0.01',
    author='Joe Yesselman',
    author_email='jyesselm@unl.edu',
    packages=["rna_draw"],
    py_modules=["rna_draw/colorer", "rna_draw/data", "rna_draw/draw",
                "rna_draw/parameters", "rna_draw/render_rna", "rna_draw/svg"],
    include_package_data=True,
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'rna_draw = rna_draw.draw:main',
        ]
    }

)