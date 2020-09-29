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
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./')
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
    else:
        print("WARNING: No model provided!")
    click.echo("--------------------------------")


    pass


########### IGoR's run commands ###########
@click.group("igor_read_seqs")
@pass_igortask
def run_read_seqs(igortask):
    """IGoR's call to read_seqs"""
    igortask.run_read_seqs()

@click.group("igor_align")
@pass_igortask
def run_align(igortask):
    """IGoR's call to aligns"""
    igortask.run_align()

@click.command("igor_infer")
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
@click.option("-o", "--output-db", "output_db", default=None, help="Output database file.")
@pass_igortask
def run_infer(igortask, igor_read_seqs, output_db):
    """IGoR's call to infer"""
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

    # if inference succesfull add files to db

    if output_db is not None:
        import os
        os.rename(igortask.igor_fln_db, output_db)
        print("Database file : ", output_db)
    else:
        print("Database file : ", igortask.igor_fln_db)
        # mv


@click.group("igor_evaluate")
@pass_igortask
def run_evaluate(igortask):
    """IGoR's call to infer"""
    igortask.run_evaluate()

@click.group("igor_generate")
@pass_igortask
def run_generate(igortask):
    """IGoR's call to infer"""
    igortask.run_generate()

# run_read_seqs.add_command(get_ref_genome)
cli.add_command(run_read_seqs)
cli.add_command(run_align)
cli.add_command(run_infer)
cli.add_command(run_evaluate)
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

@click.command()
@click.option("-t", "--recombination_type", "rec_type", type=click.Choice(['VJ', 'VDJ']),
              help="Igor recombination type.")
@click.option("--imgt-species", "imgt_species", #type=click.Choice(get_species_list()),
              help="IMGT species name for name specifications run imgt --info.")
@click.option("--imgt-chain", "imgt_chain", #type=click.Choice(['VJ', 'VDJ']),
              help="IMGT chain name e.g. TRA, TRB.")
@pass_igortask
def get_ref_genome(igortask, rec_type, imgt_species, imgt_chain):
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

imgt.add_command(get_ref_genome)
cli.add_command(imgt)

########### model commands ###########
# @click.group() #invoke_without_command=True)
# @pass_igortask
# def model(igortask):
#     """Manipulations of models"""
#     pass

@click.command("model-export")
@pass_igortask
def model_export(igortask):
    pass

@click.command("model-create")
@click.option("-t", "--recombination_type", "rec_type", type=click.Choice(['VJ', 'VDJ']),
              help="Igor recombination type.")
@pass_igortask
def model_create(igortask, rec_type):
    """Make a new default model VJ or VDJ with uniform probability distribution"""
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

    pass

@click.command("model-plot")
@click.option("-o", "--output-prefix", "fln_output_prefix", default=None, help="Prefix to pdf files with model plots.")
@pass_igortask
def model_plot(igortask, fln_output_prefix):
    if igortask.igor_fln_db is not None:

        igortask.create_db()
        igortask.load_mdl_from_db()
        print("Model loaded from ", igortask.igor_fln_db)
    else:
        igortask.update_model_filenames()
        igortask.load_IgorModel()


    import matplotlib.pyplot as plt

    for event_nickname in igortask.mdl.Pmarginal.keys():
        fig, ax = plt.subplots()
        # df = task.mdl.Pmarginal[event_nickname].to_dataframe(name=event_nickname)
        igortask.mdl.plot_Event_Marginal(event_nickname, ax=ax)
        fig.tight_layout()
        flnOutput = fln_output_prefix + "_" + event_nickname + ".pdf"
        fig.savefig(flnOutput)
        print("***** Marginal plot of ", event_nickname, " in ", flnOutput)
    pass

cli.add_command(model_export)
cli.add_command(model_create)
cli.add_command(model_plot)
# cli.add_command(model)


# pygor3-cli [GENERAL_OPTIONS] database export all
########### database commands ###########
@click.group() #invoke_without_command=True)
@pass_igortask
def database(igortask):
    """Manipulations of models"""
    pass

@click.command("export")
@pass_igortask
def database_export(igortask):
    pass

database.add_command(database_export)
# database.add_command(database_create)
cli.add_command(database)



#
# @cli.command()
# @click.option("-o", "--output", "model_output", default=None, help="Model csv file output prefix.")
# @pass_igortask
# def export_model(igortask, model_output): #igor_species, igor_chain, igor_model):
#     """Export model in csv format with a pdf file of marginals plot"""
#     click.echo('Model utilities.')
#     igortask.load_IgorModel()
#     igortask.mdl.export_csv(model_output)
#
#
# #click.echo(IgorTask.to_dict())
# @cli.command()
# @click.option("-t", "--recombination_type", "rec_type", type=click.Choice(['VJ', 'VDJ']),
#               help="Igor recombination type.")
# @pass_igortask
# def new_model(igortask, rec_type):
#     """Creates new VJ or VDJ model parms from set_path_ref_genome
#     and save it in  set_model_path """
#     click.echo("New model ")
#
#
# @cli.command()
# @click.option("-t", "--recombination_type", "rec_type", type=click.Choice(['VJ', 'VDJ']),
#               help="Igor recombination type from default VJ or VDJ.")
# @click.option("--info", "info", help="List species and chain avialable in imgt website.", is_flag=True, default=False)
# @pass_igortask
# def download_imgt_genomics(igortask, rec_type, info):
#     """Download new VJ or VDJ from imgt website model parms from set_path_ref_genome
#     and save it in  set_model_path """
#     click.echo("Downloading data from ... ")
#     import pygor3.imgt as imgt
#     print(imgt.imgt_params['url.home'])
#     if info:
#         species_list = imgt.get_species_list()
#
#         print("List of IMGT available species from :")
#         print(imgt.imgt_params['url.genelist'])
#         print(species_list)
#     else:
#
#         if (rec_type == 'VDJ'):
#             if igortask.igor_path_ref_genome is None:
#                 imgt.download_ref_genome_VDJ(igortask.igor_species, igortask.igor_chain)
#             else:
#                 imgt.download_ref_genome_VDJ(igortask.igor_species, igortask.igor_chain, modelspath=None)
#         elif(rec_type == 'VJ'):
#             if igortask.igor_path_ref_genome is None:
#                 imgt.download_ref_genome_VJ(igortask.igor_species, igortask.igor_chain)
#             else:
#                 imgt.download_ref_genome_VJ(igortask.igor_species, igortask.igor_chain, modelspath=None)
#         else:
#             click.echo("ERROR: Type "+str(rec_type)+" not valid, please choose between VDJ or VJ.")
#
# @cli.command()
# @pass_igortask
# def load_database(igortask): #task):
#     click.echo('Initialized the database')
#     # click.echo(task.to_dict())
#
# @cli.command()
# # @pass_igortask
# def initdb(): #task):
#     click.echo('Initialized the database')
#     # click.echo(task.to_dict())
#
# @cli.command()
# def dropdb():
#     click.echo('Dropped the database')

# cli.add_command(dropdb)
# cli.add_command(export_model)

if __name__ == '__main__':
    cli()


# pygor3-cli igor-infer -i sequences.txt -o model
# pygor3-cli igor-generate -human -TRB -o seqs -n 1000
# pygor3-cli olga-compute -ppost --VJmodel model/ -i
# pygor3-cli igor-model -i model/ -o csv_files
# pygor3-cli igor-model -i model/ -o png_files


# @click.command()
# @click.option("--V", "Vgene", default=True, help="Download gene V")
# def main(Vgene):
#     from pygor3.imgt import get_species_list
#     species_list = get_species_list()
#     print(species_list, Vgene)