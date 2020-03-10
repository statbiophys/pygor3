import setuptools
setuptools.setup(
    name='pygor3',
    version='0.0.1',
    description='Python package to manipulate and run IGoR data files',
    license="GNU GPLv3",
    python_requires='>=3.5',
    install_requires = [
        "pandas",
        "holoviz",
        "numpy",
        "xarray",
        "beautifulsoup4",
        "biopython"
    ],
    packages=setuptools.find_packages(),
    entry_points= {
        'console_scripts' : [
            'igor-bs_export=pygor3.bs_export:main',
            'igor-pgen_sequences=pygor3.pgen_sequences:main',
            'igor-model_export=pygor3.model_export:main',
        ],
    }
)