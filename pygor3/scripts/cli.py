import click

from pygor3.IgorIO import IgorTask

pass_igortask = click.make_pass_decorator(IgorTask, ensure=True)

@click.group()
# @click.pass_context
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse") # FIXME:
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA") # FIXME:
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')  # , default='./ref_genome')
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated') #, required=True)
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")
@pass_igortask
def cli(igortask, igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch, igor_fln_db):

    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.igor_species = igor_species
    igortask.igor_chain = igor_chain
    igor_model_parms = igor_model[0]
    igor_model_marginals = igor_model[1]
    igortask.igor_model_parms_file = igor_model[0]
    igortask.igor_model_marginals_file = igor_model[1]
    igortask.igor_model_dir_path = igor_model_path
    igortask.igor_wd = igor_wd
    igortask.igor_batchname = igor_batch
    igortask.igor_fln_db = igor_fln_db

    Q_species_chain = (not (igor_species is None) and not (igor_chain is None))
    Q_model_files = (not (igor_model_parms is None) and not (igor_model_marginals is None))
    if Q_species_chain:
        igortask.igor_species = igor_species
        igortask.igor_chain = igor_chain
    elif Q_model_files:
        igortask.igor_model_parms_file = igor_model_parms
        igortask.igor_model_marginals_file = igor_model_marginals
    elif igor_fln_db is not None:
        igortask.create_db()
    else:
        print("WARNING: No model provided!")
    click.echo("--------------------------------")



########### IGoR's run commands ###########
@click.group("igor-read-seqs")
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
@pass_igortask
def run_read_seqs(igortask, igor_read_seqs):
    """IGoR's call to read_seqs"""
    igortask.run_read_seqs(igor_read_seqs=igor_read_seqs)

@click.group("igor-align")
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
@pass_igortask
def run_align(igortask, igor_read_seqs):
    """IGoR's call to aligns"""
    igortask.run_align(igor_read_seqs=igor_read_seqs)

@click.command("igor-infer")
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
@click.option("-o", "--output-prefix", "output_fln_prefix", default=None, help="Output database file.")
@pass_igortask
def run_infer(igortask, igor_read_seqs, output_fln_prefix):
    """IGoR's call to infer model from input sequences and model"""
    print("===== Running inference =====")
    igortask.update_batch_filenames()
    igortask.update_model_filenames()
    print(igortask)

    igortask.load_IgorRefGenome()
    igortask.load_IgorModel()
    print(igortask)
    output = igortask.run_infer(igor_read_seqs=igor_read_seqs)
    print(output)

    print("===== Saving files in database : =====")
    igortask.create_db()
    igortask.load_db_from_indexed_sequences()
    igortask.load_db_from_indexed_cdr3()
    igortask.load_db_from_genomes()
    igortask.load_db_from_alignments()
    igortask.load_IgorModel_from_infer_files()
    igortask.load_db_from_models()

    # if inference succesfull add files to db

    if output_fln_prefix is not None:
        import os
        output_fln_db = output_fln_prefix+".db"
        output_fln_parms = output_fln_prefix + "_parms.txt"
        output_fln_marginals = output_fln_prefix + "_marginals.txt"

        os.rename(igortask.igor_fln_db, output_fln_db)
        igortask.mdl.write_model(output_fln_parms, output_fln_marginals)
        # copy files
        print("Database file : ", output_fln_prefix)
    else:
        print("Database file : ", igortask.igor_fln_db)
        # mv


@click.command("igor-evaluate")
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
@click.option("-o", "--output-db", "output_db", default=None, help="Output database file.")
@pass_igortask
def run_evaluate(igortask, igor_read_seqs, output_db):
    """IGoR's call to evaluate input sequences"""
    click.echo("Running IGoR evaluation process...")
    igortask.update_batch_filenames()
    igortask.update_model_filenames()

    igortask.load_IgorRefGenome()
    igortask.load_IgorModel()
    # batchname_model/
    # batchname_model/models
    # batchname_model/ref_genome
    output = igortask.run_evaluate(igor_read_seqs=igor_read_seqs)
    print(output)

    print("===== Saving files in database : =====")
    igortask.create_db()
    igortask.load_db_from_indexed_sequences()
    igortask.load_db_from_indexed_cdr3()
    igortask.load_db_from_genomes()
    igortask.load_db_from_alignments()
    igortask.load_IgorModel()
    igortask.load_db_from_models()
    igortask.load_db_from_bestscenarios()
    igortask.load_db_from_pgen()

@click.command("igor-scenarios")
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
@click.option("-o", "--output-db", "output_db", default=None, help="Output database file.")
@pass_igortask
def run_get_scenarios(igortask, igor_read_seqs, output_db):
    """IGoR's call to get best scenarios."""
    click.echo("Running IGoR scenarios process...")
    igortask.run_evaluate(igor_read_seqs=igor_read_seqs)

@click.command("igor-pgen")
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None,
              help="Input sequences in FASTA, TXT or CSV formats.")
@click.option("-o", "--output-db", "output_db", default=None, help="Output database file.")
@pass_igortask
def run_get_pgen(igortask, igor_read_seqs, output_db):
    """IGoR's call to calculate pgen of input sequences"""
    click.echo("Get IGoR pgen process...")
    igortask.run_evaluate(igor_read_seqs=igor_read_seqs)

@click.command("igor-generate")
@pass_igortask
def run_generate(igortask):
    """IGoR's call to generate sequences"""
    igortask.run_generate()

# run_read_seqs.add_command(get_ref_genome)
cli.add_command(run_read_seqs)
cli.add_command(run_align)
cli.add_command(run_infer)
cli.add_command(run_evaluate)
cli.add_command(run_get_scenarios)
cli.add_command(run_get_pgen)
cli.add_command(run_generate)


########### IMGT commands ###########
# @cli.command()
@click.group(invoke_without_command=True)
@click.option("--info", "info", help="List species and chain avialable in imgt website.", is_flag=True, default=False)
@pass_igortask
def imgt(igortask, info): #igor_species, igor_chain, igor_model):
    """Download genomic information from imgt website --info for more details"""

    import pygor3.imgt as p3imgt
    print(p3imgt.imgt_params['url.home'])
    species_list = p3imgt.get_species_list()

    if info:
        click.echo("Downloading data from ... ")
        print("List of IMGT available species:")
        print("\n".join(species_list))
        print("For more details access:")
        print(p3imgt.imgt_params['url.genelist'])

@click.command("imgt-get-genomes")
@click.option("--info", "info", help="List species and chain avialable in imgt website.", is_flag=True, default=False)
@click.option("-t", "--recombination_type", "rec_type", type=click.Choice(['VJ', 'VDJ']),
              help="Igor recombination type.")
@click.option("--imgt-species", "imgt_species", #type=click.Choice(get_species_list()),
              help="IMGT species name for name specifications run imgt --info.")
@click.option("--imgt-chain", "imgt_chain", #type=click.Choice(['VJ', 'VDJ']),
              help="IMGT chain name e.g. TRA, TRB.")
@pass_igortask
def get_ref_genome(igortask, info, rec_type, imgt_species, imgt_chain):
    """
    Get genomes from imgt website of specifing species and chain in imgt format.
    """

    import pygor3.imgt as p3imgt
    print(p3imgt.imgt_params['url.home'])
    species_list = p3imgt.get_species_list()

    if info:
        click.echo("Downloading data from ... ")
        print("List of IMGT available species:")
        print("\n".join(species_list))
        print("For more details access:")
        print(p3imgt.imgt_params['url.genelist'])
    else:
        print("get_ref_genome")
        # print(igortask.to_dict())
        try:
            import pygor3.imgt as p3imgt
            if (rec_type == 'VDJ'):
                if igortask.igor_path_ref_genome is None:
                    p3imgt.download_ref_genome_VDJ(imgt_species, imgt_chain)
                else:
                    p3imgt.download_ref_genome_VDJ(imgt_species, imgt_chain, modelspath=None)
            elif (rec_type == 'VJ'):
                if igortask.igor_path_ref_genome is None:
                    p3imgt.download_ref_genome_VJ(imgt_species, imgt_chain)
                else:
                    p3imgt.download_ref_genome_VJ(imgt_species, imgt_chain, modelspath=None)
            else:
                click.echo("ERROR: Type " + str(rec_type) + " not valid, please choose between VDJ or VJ.")
        except Exception as e:
            print("ERROR: get_ref_genome ")
            print(e)

# imgt.add_command(get_ref_genome)
# cli.add_command(imgt)
cli.add_command(get_ref_genome)


########### model commands ###########
# @click.group() #invoke_without_command=True)
# @pass_igortask
# def model(igortask):
#     """Manipulations of models"""
#     pass

@click.command("model-export")
@click.option("--from-txt", "fln_from_txt", nargs=2, default=[None, None],
              help="Export Igor's model from txt files model_parms.txt and model_marginals.txt.")
@click.option("--from-db", "fln_from_db", nargs=1, default=None,
              help="Export Igor's model from database file.")
@click.option("--to-txt", "fln_to_txt", nargs=2, default=[None, None],
              help="Output filename of Igor recombination model to  <model_parms.txt> <model_marginals.txt>.")
@click.option("--to-db", "fln_to_db", nargs=1, default=None,
              help="Output filename of Igor recombination model to  <model.db>.")
@pass_igortask
def model_export(igortask, fln_from_txt, fln_from_db, fln_to_txt, fln_to_db):
    """
    Export IGoR's models from txt (model_parms.txt, model_marginals.txt) files to db viceversa
    """
    b_from_db = False
    b_from_txt = False
    b_to_db = False
    b_to_txt = False
    if (fln_from_db is not None) : b_from_db = True
    if (fln_to_db is not None): b_to_db = True
    if (fln_from_txt[0] is not None ): b_from_txt = True
    if (fln_to_txt[0] is not None ): b_to_txt = True

    # print(fln_from_txt, fln_from_db, fln_to_txt, fln_to_db)
    # print(b_from_db, b_from_txt, b_to_db, b_to_txt)

    if b_from_db and b_from_txt:
        print("ERROR: --from-txt and --from-db are excludent")
        return 0

    if b_to_db and b_to_txt :
        print("ERROR: --to-txt and --to-db are excludent")
        return 0

    import pygor3 as p3
    mdl_from = p3.IgorModel()
    mdl_to = p3.IgorModel()
    if b_from_txt:
        try:
            mdl_from = p3.IgorModel(fln_from_txt[0], fln_from_txt[1])
            print(mdl_from)
        except Exception as e:
            print("ERROR: ", fln_from_txt)
            print(e)
    elif b_from_db:
        try:
            igor_db = p3.IgorSqliteDB.create_db(fln_from_db)
            mdl_from = igor_db.get_IgorModel()
            print(mdl_from)
        except Exception as e:
            print("ERROR: ", fln_from_db)
            print(e)
    else:
        print("A file is need it!")
        return 0

    if b_to_txt:
        try:
            if fln_to_txt[1] is not None:
                mdl_to = mdl_from
                mdl_to.parms.write_model_parms(fln_to_txt[0])
                mdl_to.marginals.write_model_marginals(fln_to_txt[1])
            else:
                mdl_to = mdl_from
                mdl_to.parms.write_model_parms(fln_to_txt[0])
                # TODO: MAKE A UNIFORM MARGINALS DISTRIBUTION
                mdl_to.marginals.initialize_uniform_event_from_model_parms(mdl_to.parms)
                fln_to_txt[1] = fln_to_txt[0]+"_marginals.txt"
                mdl_to.marginals.write_model_marginals(fln_to_txt[1])
                print("WARNING: No marginals file was specified, a uniform distribution was set and save in file ", fln_to_txt[1])
            print(mdl_to)
        except Exception as e:
            print("ERROR: ", fln_to_txt)
            print(e)
    elif b_to_db:
        try:
            mdl_to = mdl_from
            igor_db = p3.IgorSqliteDB.create_db(fln_to_db)
            # mdl_to = igor_db.get_IgorModel()
            igor_db.delete_IgorModel_Tables()
            igor_db.load_IgorModel(mdl_to)
        except Exception as e:
            print("ERROR: ", fln_to_db)
            print(e)


@click.command("model-create")
@click.option("-t", "--recombination_type", "rec_type", type=click.Choice(['VJ', 'VDJ']),
              help="Igor recombination type.")
# @click.option("-o", "--output-prefix", "fln_output_prefix", default=None, help="Prefix for models files.")
@pass_igortask
def model_create(igortask, rec_type, fln_output_prefix):
    """Make a new VJ or VDJ model with uniform probability distribution"""
    # click.echo( "igor_model_dir_path: "+igortask.igor_model_dir_path )
    import pygor3 as p3
    if rec_type == 'VJ':
        # load genomics
        # igortask.igor_model_dir_path
        igortask.update_model_filenames(model_path=igortask.igor_model_dir_path)
        igortask.load_IgorRefGenome()
        igortask.make_model_default_VJ_from_genomes()
        igortask.mdl.parms.write_model_parms(igortask.igor_model_parms_file)
        igortask.mdl.marginals.write_model_marginals(igortask.igor_model_marginals_file, igortask.mdl.parms)
    elif rec_type == 'VDJ':
        igortask.update_model_filenames(model_path=igortask.igor_model_dir_path)
        # print(igortask.to_dict())
        igortask.load_IgorRefGenome()

        igortask.make_model_default_VDJ_from_genomes()
        print(igortask.igor_model_dir_path)
        # igortask.igor_model_parms_file
        igortask.mdl.parms.write_model_parms(igortask.igor_model_parms_file)
        igortask.mdl.marginals.write_model_marginals(igortask.igor_model_marginals_file, igortask.mdl.parms)
    else:
        print("Model type ", rec_type)

    print("igortask.igor_model_dir_path: ", igortask.igor_model_dir_path)


@click.command("model-plot")
@click.option("-o", "--output-prefix", "fln_output_prefix", default=None, help="Prefix to pdf files with model plots.")
@pass_igortask
def model_plot(igortask, fln_output_prefix):
    """Plot real marginals of the bayesian network events """
    if igortask.igor_fln_db is not None:
        igortask.create_db(igortask.igor_fln_db)
        igortask.load_mdl_from_db()
        print("Model loaded from ", igortask.igor_fln_db)
    else:
        igortask.update_model_filenames()
        igortask.load_IgorModel()

    igortask.mdl.export_plot_Pmarginals(fln_output_prefix)


cli.add_command(model_export)
cli.add_command(model_create)
cli.add_command(model_plot)
# cli.add_command(model)


# pygor3-cli [GENERAL_OPTIONS] database export all
########### database commands ###########
@click.command("db-cp") #invoke_without_command=True)
@click.option("-o", "--output-db", "fln_output_db", default=None, help="output attached database.")

@click.option("--igor-reads", "b_igor_reads", is_flag=True,
              help='Copy sequences reads to output-db.')
@click.option("--igor-genomes", "b_igor_genomes", is_flag=True,
              help='Copy V, (D) and J genetic data to database')
@click.option("--igor-alignments", "b_igor_alignments", is_flag=True,
              help='Copy all available alignments tables to database.')
@click.option("--igor-model", "b_igor_model", is_flag=True,
              help='Copy IGoR model to database file.')
@click.option("--igor-scenarios", "b_igor_scenarios", is_flag=True,
              help='Copy scenarios table to database')
@click.option("--igor-pgen", "b_igor_pgen", is_flag=True,
              help='Copy pgen table to database')
@pass_igortask
def database_copy(igortask, fln_output_db, b_igor_reads,
                  b_igor_genomes,
                  b_igor_alignments,
                  b_igor_model,
                  b_igor_scenarios, b_igor_pgen):
    """Testing function before commit - export database"""
    # Create a database
    from pygor3.IgorSqliteDB import IgorSqliteDB
    # OPEN CONEXION TO output_db
    output_db = IgorSqliteDB.create_db(fln_output_db)

    # Get list of tables in igortask.igor_db
    tablename_ctsql_dict = igortask.igor_db.get_dict_of_Igortablename_sql()
    print("Tables in source database : ", tablename_ctsql_dict.keys())
    if b_igor_reads and igortask.igor_db.Q_sequences_in_db():
        try:
            # Ask to sqlite_master for the way to create it in new database
            tablename_to_copy = 'IgorIndexedSeq'
            # Create table in database destiny.
            output_db.execute_query(tablename_ctsql_dict[tablename_to_copy])
            # Copy table
            fln_source_db = igortask.igor_db.flnIgorDB
            output_db.copytable_from_source(tablename_to_copy, fln_source_db)
        except Exception as e:
            print("ERROR: igor-reads")
            print(e)

    if b_igor_genomes and igortask.igor_db.Q_ref_genome_in_db():
        try:
            # Ask to sqlite_master for the way to create it in new database
            # igortask.igor_db.
            from pygor3.IgorSQL import sql_tablename_patterns_dict
            for table_pattern in sql_tablename_patterns_dict['ref_genome']:
                tablename_list = igortask.igor_db.get_list_of_tables_with_name(table_pattern)
                for tablename_to_copy in tablename_list:
                    # tablename_to_copy = 'IgorIndexedSeq'
                    # Create table in database destiny.
                    output_db.execute_query(tablename_ctsql_dict[tablename_to_copy])
                    # Copy table
                    fln_source_db = igortask.igor_db.flnIgorDB
                    output_db.copytable_from_source(tablename_to_copy, fln_source_db)
        except Exception as e:
            print("ERROR: igor-genomes")
            print(e)

    if b_igor_alignments and igortask.igor_db.Q_align_in_db():
        try:
            # Ask to sqlite_master for the way to create it in new database
            # igortask.igor_db.
            from pygor3.IgorSQL import sql_tablename_patterns_dict
            for table_pattern in sql_tablename_patterns_dict['align']:
                tablename_list = igortask.igor_db.get_list_of_tables_with_name(table_pattern)
                for tablename_to_copy in tablename_list:
                    # tablename_to_copy = 'IgorIndexedSeq'
                    # Create table in database destiny.
                    output_db.execute_query(tablename_ctsql_dict[tablename_to_copy])
                    # Copy table
                    fln_source_db = igortask.igor_db.flnIgorDB
                    output_db.copytable_from_source(tablename_to_copy, fln_source_db)
        except Exception as e:
            print("ERROR: igor-genomes")
            print(e)

    if b_igor_model and igortask.igor_db.Q_model_in_db():
        try:
            # Ask to sqlite_master for the way to create it in new database
            # igortask.igor_db.
            from pygor3.IgorSQL import sql_tablename_patterns_dict
            for table_pattern in sql_tablename_patterns_dict['model']:
                tablename_list = igortask.igor_db.get_list_of_tables_with_name(table_pattern)
                for tablename_to_copy in tablename_list:
                    # tablename_to_copy = 'IgorIndexedSeq'
                    # Create table in database destiny.
                    output_db.execute_query(tablename_ctsql_dict[tablename_to_copy])
                    # Copy table
                    fln_source_db = igortask.igor_db.flnIgorDB
                    output_db.copytable_from_source(tablename_to_copy, fln_source_db)
        except Exception as e:
            print("ERROR: igor-model")
            print(e)

        # b_igor_scenarios

    if b_igor_pgen and igortask.igor_db.Q_IgorPgen_in_db():
        print("COPY PGEN")
        try:
            # Ask to sqlite_master for the way to create it in new database
            # igortask.igor_db.
            from pygor3.IgorSQL import sql_tablename_patterns_dict
            for table_pattern in sql_tablename_patterns_dict['pgen']:
                tablename_list = igortask.igor_db.get_list_of_tables_with_name(table_pattern)
                for tablename_to_copy in tablename_list:
                    # tablename_to_copy = 'IgorIndexedSeq'
                    # Create table in database destiny.
                    output_db.execute_query(tablename_ctsql_dict[tablename_to_copy])
                    # Copy table
                    fln_source_db = igortask.igor_db.flnIgorDB
                    output_db.copytable_from_source(tablename_to_copy, fln_source_db)
        except Exception as e:
            print("ERROR: igor-pgen")
            print(e)

        #

    if b_igor_scenarios and igortask.igor_db.Q_IgorBestScenarios_in_db():
        try:
            # Ask to sqlite_master for the way to create it in new database
            # igortask.igor_db.
            from pygor3.IgorSQL import sql_tablename_patterns_dict
            for table_pattern in sql_tablename_patterns_dict['scenarios']:
                tablename_list = igortask.igor_db.get_list_of_tables_with_name(table_pattern)
                for tablename_to_copy in tablename_list:
                    # tablename_to_copy = 'IgorIndexedSeq'
                    # Create table in database destiny.
                    output_db.execute_query(tablename_ctsql_dict[tablename_to_copy])
                    # Copy table
                    fln_source_db = igortask.igor_db.flnIgorDB
                    output_db.copytable_from_source(tablename_to_copy, fln_source_db)
        except Exception as e:
            print("ERROR: igor-scenarios")
            print(e)

    print("Tables in destiny database: ", fln_output_db)
    output_db.list_from_db()




@click.command("db-ls")
@pass_igortask
def database_ls(igortask):
    """List tables in database by groups and show  number of records."""
    igortask.create_db()
    igortask.igor_db.list_from_db()

    # Q_seqs = igortask.igor_db.Q_output_in_db()
    # print(Q_seqs)

    # igortask.igor_db.write_IgorIndexedSeq_to_CSV("jojo.csv")

@click.command("db-attach")
@click.option("--from-db", "fln_from_db", default=None, help="Database copy source filename.")
@click.option("--from-batch", "fln_from_db", default=None, help="Database copy source filename.")
@click.option("--from-genome-dir", "fln_from_db", default=None, help="Database copy source filename.")
@click.option("--from-model-path", "fln_from_db", default=None, help="Database copy source filename.")

@click.option("--igor-model-dir", "igor_model_db", default=None,
              metavar="<model.db>",
              help='IGoR model database file.')
@click.option("--igor-model-parms", "igor_model_parms", default=None,
              metavar="<model_parms.txt>",
              help='IGoR model parms (or params) file.')
@click.option("--igor-model-marginals", "igor_model_marginals", default=None,
              metavar="<model_marginals.txt>",
              help='IGoR model marginals file.')
@click.option("--scenarios", "scenarios", default=None,
              help='If --from-db no need to add filename')
@click.option("--pgen", "pgen", default=None,
              help='If --from-db no need to add filename')
@click.option("--genomes", "genomes", default=None,
              help='Copy V, (D) and J genetic data to database')
@click.option("--genomesV", "genomesV", default=None,
              help='Copy just V genomes tables to database.')
@click.option("--genomesD", "genomesD", default=None,
              help='Copy just D genomes tables to database.')
@click.option("--genomesJ", "genomesJ", default=None,
              help='Copy just J genomes tables to database.')
@click.option("--genomesCDR3", "genomesCDR3", default=None,
              help='Copy just CDR3 anchors tables to database.')
@click.option("--alignments", "alignments", default=None,
              help='Copy all available alignments tables to database.')
@click.option("--alignmentsV", "alignmentsV", default=None,
              help='Copy V alignments tables to database.')
@click.option("--alignmentsD", "alignmentsD", default=None,
              help='Copy D alignments tables to database.')
@click.option("--alignmentsJ", "alignmentsJ", default=None,
              help='Copy J alignments tables to database.')
@click.option("--alignmentsCDR3", "alignmentsCDR3", default=None,
              help='Copy indexed cdr3 table to database.')
# @click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
# @click.option("-o", "--output-prefix", "output_prefix", default=None, help="Filename output prefix.")
@pass_igortask
def database_attach(igortask, fln_from_db):
    """
    Attach tables to database.
    """
    igortask.create_db()
    from pygor3 import IgorSqliteDB
    other_igor_db = IgorSqliteDB.create_db(fln_from_db)
    dicto = other_igor_db.get_dict_of_Igortablename_sql()
    print(dicto.keys())
    import pygor3 as p3
    print("sql_tablename_patterns_dict : ", p3.sql_tablename_patterns_dict)
    other_igor_db.close_db()
    # igortask.igor_db.attach_table_from_db("ddd")
    # igortask.igor_db.write_IgorIndexedSeq_to_CSV("jojo.csv")



@click.command("db-rm")
@click.option("--igor-reads", "b_igor_reads", is_flag=True,
              help='Delete sequences reads in database.')
@click.option("--igor-model", "b_igor_model", is_flag=True,
              help='IGoR model database file.')
@click.option("--igor-model-parms", "b_igor_model_parms", is_flag=True,
              help='IGoR model parms (or params) file.')
@click.option("--igor-model-marginals", "b_igor_model_marginals", is_flag=True,
              help='IGoR model marginals file.')

@click.option("--igor-scenarios", "b_igor_scenarios", is_flag=True,
              help='If --from-db no need to add filename')
@click.option("--igor-pgen", "b_igor_pgen", is_flag=True,
              help='If --from-db no need to add filename')
@click.option("--igor-genomes", "b_igor_genomes", is_flag=True,
              help='Delete V, (D) and J genetic data.')
@click.option("--igor-genomesV", "b_igor_genomesV", is_flag=True,
              help='Delete just V genomes table.')
@click.option("--igor-genomesD", "b_igor_genomesD", is_flag=True,
              help='Delete just D genomes table.')
@click.option("--igor-genomesJ", "b_igor_genomesJ", is_flag=True,
              help='Delete just J genomes table.')
@click.option("--igor-genomesCDR3", "b_igor_genomesCDR3", is_flag=True,
              help='Delete just CDR3 anchors table.')
@click.option("--igor-alignments", "b_igor_alignments", is_flag=True,
              help='Delete all available alignments table.')
@click.option("--igor-alignmentsV", "b_igor_alignmentsV", is_flag=True,
              help='Delete V alignments table.')
@click.option("--igor-alignmentsD", "b_igor_alignmentsD", is_flag=True,
              help='Delete D alignments table.')
@click.option("--igor-alignmentsJ", "b_igor_alignmentsJ", is_flag=True,
              help='Delete J alignments table.')
@click.option("--igor-alignmentsCDR3", "b_igor_alignmentsCDR3", is_flag=True,
              help='Delete indexed cdr3 table.')
# @click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
# @click.option("-o", "--output-prefix", "output_prefix", default=None, help="Filename output prefix.")
@pass_igortask
def database_rm(igortask, b_igor_reads, b_igor_model,
                b_igor_model_parms, b_igor_model_marginals,
                b_igor_scenarios, b_igor_pgen,
                b_igor_genomes,
                b_igor_genomesV, b_igor_genomesD, b_igor_genomesJ, b_igor_genomesCDR3,
                b_igor_alignments,
                b_igor_alignmentsV, b_igor_alignmentsD, b_igor_alignmentsJ,
                b_igor_alignmentsCDR3):
    """
    Delete tables in database by groups.
    """
    igortask.create_db()
    # from pygor3.IgorSqliteDB import *
    if b_igor_reads:
        try:
            igortask.igor_db.delete_IgorIndexedSeq_Tables()
        except Exception as e:
            print("ERROR: igor-reads")
            print(e)

    if b_igor_genomes:
        try:
            igortask.igor_db.delete_IgorGeneTemplate_Tables()
            igortask.igor_db.delete_IgorGeneAnchors_Tables()
        except Exception as e:
            print("ERROR: igor-genomes")
            print(e)


    if b_igor_alignments:
        try:
            igortask.igor_db.delete_IgorAlignments_Tables()
            igortask.igor_db.delete_IgorIndexedCDR3_Tables()
        except Exception as e:
            print("ERROR: igor-alignments")
            print(e)

    if b_igor_model:
        try:
            igortask.igor_db.delete_IgorModel_Tables()
        except Exception as e:
            print("ERROR: igor-model")
            print(e)

    if b_igor_pgen:
        try:
            igortask.igor_db.delete_IgorPgen_Tables()
        except Exception as e:
            print("ERROR: igor-pgen")
            print(e)
    if b_igor_scenarios:
        try:
            igortask.igor_db.delete_IgorBestScenarios_Tables()
        except Exception as e:
            print("ERROR: igor-scenarios")
            print(e)


@click.command("db-export") #invoke_without_command=True)
# @click.option("--igor-all", "b_igor_all", is_flag=True,
#               help='Export all available data in db with prefix.')
# @click.option("--igor-reads", "b_igor_reads", is_flag=True,
#               help='Delete sequences reads in database.')
# @click.option("--igor-model", "b_igor_model", is_flag=True,
#               help='IGoR model database file.')
# @click.option("--igor-model-parms", "b_igor_model_parms", is_flag=True,
#               help='IGoR model parms (or params) file.')
# @click.option("--igor-model-marginals", "b_igor_model_marginals", is_flag=True,
#               help='IGoR model marginals file.')
#
# @click.option("--igor-scenarios", "b_igor_scenarios", is_flag=True,
#               help='If --from-db no need to add filename')
# @click.option("--igor-pgen", "b_igor_pgen", is_flag=True,
#               help='If --from-db no need to add filename')
# @click.option("--igor-genomes", "b_igor_genomes", is_flag=True,
#               help='Copy V, (D) and J genetic data to database')
# @click.option("--igor-genomesV", "b_igor_genomesV", is_flag=True,
#               help='Copy just V genomes tables to database.')
# @click.option("--igor-genomesD", "b_igor_genomesD", is_flag=True,
#               help='Copy just D genomes tables to database.')
# @click.option("--igor-genomesJ", "b_igor_genomesJ", is_flag=True,
#               help='Copy just J genomes tables to database.')
# @click.option("--igor-genomesCDR3", "b_igor_genomesCDR3", is_flag=True,
#               help='Copy just CDR3 anchors tables to database.')
# @click.option("--igor-alignments", "b_igor_alignments", is_flag=True,
#               help='Copy all available alignments tables to database.')
# @click.option("--igor-alignmentsV", "b_igor_alignmentsV", is_flag=True,
#               help='Copy V alignments tables to database.')
# @click.option("--igor-alignmentsD", "b_igor_alignmentsD", is_flag=True,
#               help='Copy D alignments tables to database.')
# @click.option("--igor-alignmentsJ", "b_igor_alignmentsJ", is_flag=True,
#               help='Copy J alignments tables to database.')
# @click.option("--igor-alignmentsCDR3", "b_igor_alignmentsCDR3", is_flag=True,
#               help='Copy indexed cdr3 table to database.')
@pass_igortask
def database_export(igortask):
    """Export database model in igor formatted files"""
    # fln_prefix = "tmp_export_"
    print("batchname: ", igortask.igor_batchname)
    if igortask.igor_batchname is None:
        print("ERROR: batch option required.")
        return 0

    # TODO: IgorTask should help me to pass from batch to db and from db to batch and model_directory.
    igor_fln_db = igortask.igor_fln_db
    igortask.update_batch_filenames()
    igortask.update_model_filenames(model_path=igortask.igor_batchname)

    igortask.create_db(igor_fln_db)
    ii_db = igortask.igor_db
    print("ii_db : ", ii_db, type(ii_db), type(igortask.igor_db))
    # import pygor3 as p3
    # ii_db = p3.IgorSqliteDB()
    ii_db.write_IgorIndexedSeq_to_CSV(igortask.igor_fln_indexed_sequences)

    # b_igor_genomes
    print(igortask.fln_genomicVs)
    igortask.genomes.update_fln_names()

    if ii_db.Q_ref_genome_in_db_by_gene("V"):
        igortask.fln_genomicVs = igortask.genomes.fln_genomicVs
        ii_db.write_IgorGeneTemplate_to_fasta("V", igortask.fln_genomicVs)
    if ii_db.Q_ref_genome_in_db_by_gene("J"):
        igortask.fln_genomicJs = igortask.genomes.fln_genomicJs
        ii_db.write_IgorGeneTemplate_to_fasta("J", igortask.fln_genomicJs)
    if ii_db.Q_ref_genome_in_db_by_gene("D"):
        igortask.fln_genomicDs = igortask.genomes.fln_genomicDs
        ii_db.write_IgorGeneTemplate_to_fasta("D", igortask.fln_genomicDs)

    if ii_db.Q_align_in_db():
        # b_igor_alignments
        if ii_db.Q_align_in_db_by_gene("V"):
            ii_db.write_IgorAlignments_to_CSV("V", igortask.igor_fln_align_V_alignments)
        if ii_db.Q_align_in_db_by_gene("J"):
            ii_db.write_IgorAlignments_to_CSV("J", igortask.igor_fln_align_J_alignments)
        if ii_db.Q_align_in_db_by_gene("D"):
            ii_db.write_IgorAlignments_to_CSV("D", igortask.igor_fln_align_D_alignments)
        try:
            ii_db.write_IgorIndexedCDR3_to_CSV(igortask.igor_fln_indexed_CDR3)
        except Exception as e:
            print("WARNING: No indexed CDR3 files found", igortask.igor_fln_indexed_CDR3)
            print(e)
            pass

    # b_igor_model
    if ii_db.Q_model_in_db():
        print("MODEL : ",igortask.igor_model_parms_file, igortask.igor_model_marginals_file)
        ii_db.write_IgorModel_to_TXT(igortask.igor_model_parms_file, igortask.igor_model_marginals_file)
    # b_igor_model_parms
    # b_igor_model_marginals

    # b_igor_pgen
    if ii_db.Q_IgorPgen_in_db():
        ii_db.write_IgorPgen_to_CSV(igortask.igor_fln_output_pgen)
    # b_igor_scenarios
    if ii_db.Q_IgorBestScenarios_in_db():
        ii_db.write_IgorBestScenarios_to_CSV(igortask.igor_fln_output_scenarios)

    igortask.export_to_igorfiles()

    """
    igortask.load_db_from_indexed_cdr3()

    igortask.igor_db.write_IgorIndexedSeq_to_CSV(fln_prefix+"_indexed_seqs.csv")
    strGene = "J"
    igortask.igor_db.write_IgorGeneTemplate_to_fasta(strGene, fln_prefix+"_genomics"+strGene+".fasta")
    igortask.igor_db.write_IgorGeneAnchors_to_CSV(strGene, fln_prefix+"_genomics"+strGene+"_anchors.csv")
    # TODO: FINISH IT

    igortask.igor_db.write_IgorAlignments_to_CSV(strGene, fln_prefix+"_"+strGene+"_aligns.csv")
    # print(igortask.igor_db.get_columns_type_of_tables('IgorBestScenarios'))
    igortask.igor_db.write_IgorBestScenarios_to_CSV(fln_prefix+'_scenarios.csv')
    igortask.igor_db.write_IgorPgen_to_CSV(fln_prefix+'_pgen.csv')
    """


@click.command("db-naive-align") #invoke_without_command=True)
@click.option("-o", "--output-csv", "output_csv", default=None, help="Filename of output file.")
@pass_igortask
def database_naive_align(igortask, output_csv):
    """
    Get a naive alignment (no scenarios) from igor alignments
    """
    import pygor3 as p3

    def generate_csv_line(indexed_sequence: p3.IgorIndexedSequence, indexed_cdr3_record: list,
                          list_vdj_alignments: dict, sep=';',
                          header_list=['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score', 'd_score',
                                       'j_score']):
        csv_line = ""
        indexed_sequence.sequence = indexed_sequence.sequence.lower()
        fields_dict = dict()
        fields_dict['sequence_id'] = str(indexed_sequence.seq_index)
        fields_dict['sequence'] = indexed_sequence.sequence
        fields_dict['v_call'] = list_vdj_alignments['V'].strGene_name
        fields_dict['d_call'] = list_vdj_alignments['D'].strGene_name
        fields_dict['j_call'] = list_vdj_alignments['J'].strGene_name
        fields_dict['v_score'] = str(list_vdj_alignments['V'].score)
        fields_dict['d_score'] = str(list_vdj_alignments['D'].score)
        fields_dict['j_score'] = str(list_vdj_alignments['J'].score)
        # FIXME: CREATE A BETTER WAY TO DO THIS
        fields_dict['v_anchor'] = ""
        fields_dict['j_anchor'] = ""
        fields_dict['junction'] = ""
        fields_dict['junction_aa'] = ""
        try:
            fields_dict['v_anchor'] = str(indexed_cdr3_record[1])  # ""
            fields_dict['j_anchor'] = str(indexed_cdr3_record[2])  # ""
            fields_dict['junction'] = str(indexed_cdr3_record[3])  # ""
            fields_dict['junction_aa'] = str(indexed_cdr3_record[4])  # ""
        except Exception as e:
            print("No junction for sequence : " + fields_dict['sequence_id'])
            print(e)
            pass

        # align = p3.IgorAlignment_data()

        for field in header_list:
            csv_line = csv_line + fields_dict[field] + sep
        return csv_line

    #### STARTS HERE
    db = p3.IgorSqliteDB()
    # db.flnIgorDB = task.igor_wd+"/"+task.igor_batchname+".db"
    db.flnIgorDB = igortask.igor_db.flnIgorDB #args.database
    db.connect_db()

    # Make a loop over all sequences
    if output_csv is None:
        fln_output = db.flnIgorDB.split(".db")[0] + "_na.csv"
    else:
        fln_output = output_csv

    ofile = open(fln_output, 'w')
    seq_index_list = db.execute_select_query("SELECT seq_index FROM IgorIndexedSeq;")
    seq_index_list = map(lambda x: x[0], seq_index_list)
    header_list = ['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score', 'd_score', 'j_score', 'v_anchor',
                   'j_anchor', 'junction', 'junction_aa']
    sep = ";"
    str_header = sep.join(header_list)
    ofile.write(str_header + '\n')
    for seq_index in seq_index_list:
        try:
            # seq_index = args.seq_index
            indexed_sequence = db.get_IgorIndexedSeq_By_seq_index(seq_index)
            indexed_sequence.offset = 0

            indexed_cdr3_record = db.fetch_IgorIndexedCDR3_By_seq_index(seq_index)
            print(indexed_cdr3_record)

            best_v_align_data = db.get_best_IgorAlignment_data_By_seq_index('V', indexed_sequence.seq_index)
            best_d_align_data = db.get_best_IgorAlignment_data_By_seq_index('D', indexed_sequence.seq_index)
            best_j_align_data = db.get_best_IgorAlignment_data_By_seq_index('J', indexed_sequence.seq_index)

            vdj_naive_alignment = {'V': best_v_align_data,
                                   'D': best_d_align_data,
                                   'J': best_j_align_data}

            v_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
            # print('V', len(v_align_data_list), [ ii.score for ii in v_align_data_list])
            try:
                d_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('D', indexed_sequence.seq_index)
            except Exception as e:
                print(e)
                print("No D sequences alignments found!")
                pass
            # print('D', len(d_align_data_list), [ ii.score for ii in d_align_data_list])
            j_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('J', indexed_sequence.seq_index)
            # print('J', len(j_align_data_list), [ ii.score for ii in j_align_data_list])

            # 1. Choose the highest score then check if this one is the desire range.
            # if there is an overlap
            # calculate score without overlap. If overlap
            # if hightest score
            for i, d_align_data in enumerate(d_align_data_list):
                # Check if D is btwn V and J position
                if (best_v_align_data.offset_3_p <= d_align_data.offset_5_p) and (
                        d_align_data.offset_3_p <= best_j_align_data.offset_5_p):
                    # vdj_naive_alignment['D'+str(i)] = d_align_data
                    vdj_naive_alignment['D'] = d_align_data
                    break

            # ofile = open(task.igor_batchname+'__'+str(indexed_sequence.seq_index)+'_na.fasta', 'w')
            # str_fasta = generate_str_fasta(indexed_sequence, vdj_naive_alignment)
            # ofile.write(str_fasta)

            str_csv_line = generate_csv_line(indexed_sequence, indexed_cdr3_record, vdj_naive_alignment, sep=sep,
                                             header_list=header_list)
            ofile.write(str_csv_line + "\n")
        except Exception as e:
            print("ERROR : with sequence id : " + str(seq_index))
            print(e)
            pass
    ofile.close()


# database.add_command(database_create)
cli.add_command(database_copy)
cli.add_command(database_ls)
cli.add_command(database_attach)
cli.add_command(database_rm)
cli.add_command(database_export)
cli.add_command(database_naive_align)



if __name__ == '__main__':
    cli()


# pygor3-cli igor-infer -i sequences.txt -o model
# pygor3-cli igor-generate -human -TRB -o seqs -n 1000
# pygor3-cli olga-compute -ppost --VJmodel model/ -i
# pygor3-cli igor-model -i model/ -o csv_files
# pygor3-cli igor-model -i model/ -o png_files
