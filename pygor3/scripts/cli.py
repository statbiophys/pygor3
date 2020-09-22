import click

from pygor3.IgorIO import IgorTask

pass_igortask = click.make_pass_decorator(IgorTask, ensure=True)


@click.group()
# @click.pass_context
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse") # FIXME:
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA") # FIXME:
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2, help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path",
                    help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome",
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')  # , default='./ref_genome')
@click.option("-w", "--set_wd", "igor_wd", help="Path where files gonna be created.", default='./')
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch",
                    help='Batchname to identify run. If not set random name is generated') #, required=True)
@click.option("-D", "--set_database", "igor_db", help="Igor database created with database script.")
@pass_igortask
def cli(igortask, igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome,igor_wd, igor_batch,igor_db):
    igortask.igor_species = igor_species
    igortask.igor_chain = igor_chain
    igor_model_parms = igor_model[0]
    igor_model_marginals = igor_model[0]
    igortask.igor_model_parms_file = igor_model[0]
    igortask.igor_model_marginals_file = igor_model[1]

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
    pass


@cli.command()
@click.option("-o", "--output", "model_output", default=None, help="Model csv file output prefix.")
# @click.option("-c", "--igor_chain", "igor_chain", default=None, help="Species in IGoR's format")
# @click.option("-m", "--igor_model", "igor_model", default=["model_parms.txt", "model_marginals.txt"], help='IGoR model_params.txt')
# @click.argument("outfilename", type=click.File('w'))
@pass_igortask
def export_model(igortask, model_output): #igor_species, igor_chain, igor_model):
    """Export model in csv format with a pdf file of marginals plot"""
    click.echo('Model utilities.')
    igortask.load_IgorModel()
    igortask.mdl.export_csv(model_output)


#click.echo(IgorTask.to_dict())
@cli.command()
@click.option("-t", "--recombination_type", "rec_type", help="Igor recombination type from default VJ or VDJ.")
@pass_igortask
def new_model(igortask, rec_type):
    """Creates new VJ or VDJ model parms from set_path_ref_genome
    and save it in  set_model_path """
    click.echo("New model ")


@cli.command()
@click.option("-t", "--recombination_type", "rec_type", help="Igor recombination type from default VJ or VDJ.")
@click.option("--info", "info", help="List species and chain avialable in imgt website.")
@pass_igortask
def download_imgt_genomics(igortask, rec_type, species_list):
    """Download new VJ or VDJ from imgt website model parms from set_path_ref_genome
    and save it in  set_model_path """
    click.echo("Download from imgt ")

@cli.command()
# @pass_igortask
def initdb(): #task):
    click.echo('Initialized the database')
    # click.echo(task.to_dict())

@cli.command()
def dropdb():
    click.echo('Dropped the database')

cli.add_command(initdb)
cli.add_command(dropdb)
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