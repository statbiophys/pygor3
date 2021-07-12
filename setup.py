import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

data_files_to_include = [('', ['README.md', 'LICENSE'])]


setuptools.setup(
    name='pygor3',
    url="https://github.com/alfaceor/pygor3",
    author="Carlos Olivares",
    author_email="carlos.olivares@phys.ens.fr",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # version='0.0.5.dev1',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    description='Python package to manipulate and run IGoR data files',
    license="GNU GPLv3",
    python_requires='>=3.5',
    install_requires = [
        "pandas",
        #"holoviz",
        "numpy",
        "xarray",
        "beautifulsoup4",
        "biopython",
        "Click",
        "networkx",
        "requests",
        "matplotlib",
        "airr",
        "appdirs"
    ],
    packages=setuptools.find_packages(),
    package_data = {
            'pygor3': ['demo/data/IgL/*.txt', "config.json"],
            # 'pygor3.demo.data': ['*'],
            # 'pygor3.demo.data.IgL': ['*.txt']
            },
    data_files = data_files_to_include,
    include_package_data=True,
    entry_points= {
        'console_scripts' : [
            'pygor=pygor3.scripts.cli:cli',
            'pygor3-bs_pairwise_prob=pygor3.scripts.bs_pairwise_prob:main',
            'pygor3-observable=pygor3.scripts.observable:main',
            #'pygor3-scriptTest=pygor3.scripts.scriptTest:main',
            #'pygor3-scriptTest02=pygor3.scripts.scriptTest02:main',
            'pygor3-stas_from_bs=pygor3.scripts.stas_from_bs:main',
            'pygor3-load_database=pygor3.scripts.load_database:main',
            'pygor3-naive_align=pygor3.scripts.naive_align:main',
            'pygor3-all_gene_alignments_id=pygor3.scripts.all_gene_alignments_id:main',
            'pygor3-naive_align_id=pygor3.scripts.naive_align_id:main',
            'pygor3-align_export=pygor3.scripts.align_export:main',
            'pygor3-bs_export=pygor3.scripts.bs_export:main',
            #'igor-pgen_sequences=pygor3.scripts.pgen_sequences:main',
            #'igor-infer_sequences=pygor3.scripts.infer_sequences:main',
            #'igor-model_plot=pygor3.scripts.model_plot:main',
            #'igor-model_export=pygor3.scripts.model_export:main',
            'pygor3-model_export=pygor3.scripts.model_export:main',
            'pygor3-plot_marginals=pygor3.scripts.plot_marginals:main',
            #'igor-model_create=pygor3.scripts.model_create:main',
            #'igor-get_imgt_data=pygor3.scripts.get_imgt_data:main',
            #'pygor3-download_imgt_gene_templates=pygor3.scripts.download_imgt_gene_templates:main',
            #'pygor3-short_names_imgt_templates=pygor3.scripts.short_names_imgt_templates:main'
            #'igor-model_import=pygor3.model_import:main',
            #'igor-get_imgt_data=pygor3.get_imgt_data:main'
        ],
    }
)
