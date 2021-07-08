import click

import os

# def get_env_vars(ctx, args, incomplete):
#     return [k for k in os.environ.keys() if incomplete in k]

class RegisterReaderOption(click.Option):
    """ Mark this option as getting a _set option """
    register_reader = True

class RegisterWriterOption(click.Option):
    """ Fix the help for the _set suffix """
    def get_help_record(self, ctx):
        help = super(RegisterWriterOption, self).get_help_record(ctx)
        return (help[0].replace('_set ', '='),) + help[1:]

class RegisterWriterCommand(click.Command):
    def parse_args(self, ctx, args):
        """ Translate any opt= to opt_set= as needed """
        options = [o for o in ctx.command.params
                   if getattr(o, 'register_reader', None)]
        prefixes = {p for p in sum([o.opts for o in options], [])
                    if p.startswith('--')}
        for i, a in enumerate(args):
            a = a.split('=')
            if a[0] in prefixes and len(a) > 1:
                a[0] += '_set'
                args[i] = '='.join(a)

        return super(RegisterWriterCommand, self).parse_args(ctx, args)


# from pygor3.IgorIO import IgorTask
# pass_igortask = click.make_pass_decorator(IgorTask, ensure=True)
from pygor3 import __version__
@click.group()
@click.version_option(version=__version__)
def cli():
    # igortask, igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch, igor_fln_db,
    #     fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors):
    click.echo("--------------------------------")

#######################################################################
import functools
def common_options(f):
    options = [
        click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse"),
        click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA"),
        click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
                      metavar="<model_parms.txt> <model_marginals.txt>",
                      help='IGoR model_params.txt and model_marginals.txt filenames.'),
        click.option("-M", "--set_model_path", "igor_model_path", default=None,
                      metavar="<model_directory_path>",
                      help='IGoR model directory path, this path include ref_genomes and model_parms'),
        click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                            help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv'),
        click.option("--Vgene", "fln_genomicVs", default=None,
                      help="Genes template for V gene"),
        click.option("--Dgene", "fln_genomicDs", default=None,
                      help="Gimport functoolsenes template for D gene"),
        click.option("--Jgene", "fln_genomicJs", default=None,
                      help="Genes template for J gene"),
        click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
                      help="V gene anchors filename"),
        click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
                      help="V gene anchors filename"),
        click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True),
        click.option("-b", "--set_batch", "igor_batch", default=None,
                            help='Sets batchname to identify run. If not set random name is generated'),
        click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

    ]
    return functools.reduce(lambda x, opt: opt(x), options, f)



########### IGoR's run commands ###########
@click.group("igor-read-seqs")
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse")
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA")
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')
@click.option("--Vgene", "fln_genomicVs", default=None,
              help="Filename of genes template for V gene")
@click.option("--Dgene", "fln_genomicDs", default=None,
              help="Filename of genes template for D gene")
@click.option("--Jgene", "fln_genomicJs", default=None,
              help="Filename of genes template for J gene")
@click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
              metavar="<V_gene_CDR3_anchors.csv>",
              help="Anchors filename of V gene")
@click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
              metavar="<J_gene_CDR3_anchors.csv>",
              help="Anchors filename of J gene")
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated')
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
def run_read_seqs(igor_read_seqs, igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch, igor_fln_db,
        fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors):
    """IGoR's call to read_seqs"""
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

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
    igortask.run_read_seqs(igor_read_seqs=igor_read_seqs)







@click.group("igor-align")
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse") # FIXME:
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA") # FIXME:
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')
@click.option("--Vgene", "fln_genomicVs", default=None,
              help="Filename of genes template for V gene")
@click.option("--Dgene", "fln_genomicDs", default=None,
              help="Filename of genes template for D gene")
@click.option("--Jgene", "fln_genomicJs", default=None,
              help="Filename of genes template for J gene")
@click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
              metavar="<V_gene_CDR3_anchors.csv>",
              help="Anchors filename of V gene")
@click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
              metavar="<J_gene_CDR3_anchors.csv>",
              help="Anchors filename of J gene")
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated') #, required=True)
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
def run_align(igor_read_seqs,
              igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch,
              igor_fln_db,
              fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors):
    """IGoR's call to aligns"""
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

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

    igortask.run_align(igor_read_seqs=igor_read_seqs)










@click.command("igor-infer")
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse")
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA")
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')
@click.option("--Vgene", "fln_genomicVs", default=None,
              metavar="<genomicVs.fasta>",
              help="Filename of genes template for V gene")
@click.option("--Dgene", "fln_genomicDs", default=None,
              metavar="<genomicDs.fasta>",
              help="Filename of genes template for D gene")
@click.option("--Jgene", "fln_genomicJs", default=None,
              metavar="<genomicJs.fasta>",
              help="Filename of genes template for J gene")
@click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
              metavar="<V_gene_CDR3_anchors.csv>",
              help="Anchors filename of V gene")
@click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
              metavar="<J_gene_CDR3_anchors.csv>",
              help="Anchors filename of J gene")
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated')
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")
############## NO COMMON options ##############
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
@click.option("-o", "--output-prefix", "output_fln_prefix", default=None, help="Output database file.")
def run_infer(igor_read_seqs, output_fln_prefix,
    igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch, igor_fln_db,
    fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors):
    """IGoR's call to infer model from input sequences and model"""

    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

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
    elif igor_model_path is not None:
        igortask.igor_model_dir_path = igor_model_path
    elif igor_fln_db is not None:
        igortask.create_db()
    else:
        print("WARNING: No model provided!")


    print("===== Running inference =====")
    igortask.update_batch_filenames()
    igortask.update_model_filenames()
    print(igortask)

    igortask.load_IgorRefGenome()
    igortask.load_IgorModel()
    print(igortask)
    output = igortask._run_infer(igor_read_seqs=igor_read_seqs)
    igortask._run_clean_batch_mdldata()
    # print(output)

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
        igortask.mdl.plot_Bayes_network(filename=output_fln_prefix+"_BN.pdf")
        igortask.mdl.export_plot_Pmarginals(output_fln_prefix+"_MP")

        os.rename(igortask.igor_fln_db, output_fln_db)
        igortask.mdl.write_model(output_fln_parms, output_fln_marginals)
        # copy files
        print("Database file : ", output_fln_prefix)
        igortask.run_clean_batch()
    else:
        # igortask.mdl.write_model(output_fln_parms, output_fln_marginals)
        print("Database file : ", igortask.igor_fln_db)
        base_fln_output = igortask.igor_fln_db.split(".db")[0]
        output_fln_prefix = base_fln_output
        # output_fln_db = output_fln_prefix + ".db"
        output_fln_parms = output_fln_prefix + "_parms.txt"
        output_fln_marginals = output_fln_prefix + "_marginals.txt"

        if len(igortask.mdl.V_anchors) > 0:
            output_fln_V_gene_CDR3_anchors = output_fln_prefix + "_V_gene_CDR3_anchors.csv"
        else:
            output_fln_V_gene_CDR3_anchors = None

        if len(igortask.mdl.J_anchors) > 0:
            output_fln_J_gene_CDR3_anchors = output_fln_prefix + "_J_gene_CDR3_anchors.csv"
        else:
            output_fln_J_gene_CDR3_anchors = None

        igortask.mdl.plot_Bayes_network(filename=output_fln_prefix+"_BN.pdf")
        igortask.mdl.export_plot_Pmarginals(output_fln_prefix + "_RM")
        igortask.mdl.write_model(output_fln_parms, output_fln_marginals,
                                 fln_V_gene_CDR3_anchors=output_fln_V_gene_CDR3_anchors,
                                 fln_J_gene_CDR3_anchors=output_fln_J_gene_CDR3_anchors)
        igortask.run_clean_batch()
        # mv









@click.command("igor-evaluate")
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse")
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA")
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')
@click.option("--Vgene", "fln_genomicVs", default=None,
              help="Filename of genes template for V gene")
@click.option("--Dgene", "fln_genomicDs", default=None,
              help="Filename of genes template for D gene")
@click.option("--Jgene", "fln_genomicJs", default=None,
              help="Filename of genes template for J gene")
@click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
              metavar="<V_gene_CDR3_anchors.csv>",
              help="Anchors filename of V gene")
@click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
              metavar="<J_gene_CDR3_anchors.csv>",
              help="Anchors filename of J gene")
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated')
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
@click.option("-o", "--output-prefix", "output_fln_prefix", default=None, help="Output prefix for database file a scenarios file.")
def run_evaluate(igor_read_seqs, output_fln_prefix,
            igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch, igor_fln_db,
            fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors):
    """IGoR's call to evaluate input sequences"""
    ########################
    from pygor3 import IgorTask
    b_clean_igortask_input=False
    if not (igor_fln_db is None):
        igortask_input = IgorTask()
        igortask_input.create_db(igor_fln_db=igor_fln_db)
        igortask_input.igor_db.list_from_db()
        igortask_input.update_batch_filenames()
        path_mdl_data = igortask_input.igor_batchname + "_mdldata"
        igortask_input.update_model_filenames(igor_model_dir_path=path_mdl_data)
        igortask_input.update_ref_genome()

        igor_model_path = path_mdl_data
        print(igortask_input.to_dict())
        igortask_input.db_export_to_igorfiles()
        igortask_input.igor_db.close_db()
        b_clean_igortask_input = True

    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

    igortask.igor_species = igor_species
    igortask.igor_chain = igor_chain
    igor_model_parms = igor_model[0]
    igor_model_marginals = igor_model[1]
    igortask.igor_model_parms_file = igor_model[0]
    igortask.igor_model_marginals_file = igor_model[1]
    igortask.igor_model_dir_path = igor_model_path
    if igor_wd is not None:
        igortask.igor_wd = igor_wd
    if igor_batch is not None:
        igortask.igor_batchname = igor_batch
    igortask.igor_fln_db = igor_fln_db



    Q_species_chain = (not (igor_species is None) and not (igor_chain is None))
    Q_model_files = (not (igor_model_parms is None) and not (igor_model_marginals is None))
    if Q_species_chain:
        igortask.igor_species = igor_species
        igortask.igor_chain = igor_chain
        # load models from default
        igortask.load_IgorModel()
    elif Q_model_files:
        igortask.igor_model_parms_file = igor_model_parms
        igortask.igor_model_marginals_file = igor_model_marginals
        igortask.load_IgorModel()
    elif igor_fln_db is not None:
        igortask.create_db(igor_fln_db=igor_fln_db)
        igortask.load_mdl_from_db()
    else:
        print("ERROR: No model provided!")
        return 0

    click.echo("Running IGoR evaluation process...")
    # TODO:  with the loaded model generate the mdl_datadir
    igortask.write_mdldata_dir()
    igortask.load_IgorRefGenome()


    igortask.update_batch_filenames()
    # igortask.update_model_filenames()
    # igortask.update_ref_genome()

    # igortask.load_IgorRefGenome()

    import json
    print(json.dumps(igortask.to_dict()))

    # igortask.load_IgorModel()
    # batchname_model/
    # batchname_model/models
    # batchname_model/ref_genome
    try:
        output = igortask._run_evaluate(igor_read_seqs=igor_read_seqs)
    except Exception as e:
        raise e
    print(output)
    if b_clean_igortask_input:
        igortask_input.run_clean_batch()


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

    # Use the database to get bestscenarios


    if output_fln_prefix is not None:
        import os
        output_fln_db = output_fln_prefix + ".db"
        output_fln_parms = output_fln_prefix + "_parms.txt"
        output_fln_marginals = output_fln_prefix + "_marginals.txt"
        output_fln_scenarios = output_fln_prefix+"_scenarios.csv"
        output_fln_pgen = output_fln_prefix + "_pgen.csv"
        output_fln_airr = output_fln_prefix + ".tsv"
        # igortask.mdl.plot_Bayes_network(filename=output_fln_prefix + "_BN.pdf")
        # igortask.mdl.export_plot_Pmarginals(output_fln_prefix + "_RM")
        os.rename(igortask.igor_fln_db, output_fln_db)
        os.rename(igortask.igor_fln_output_scenarios, output_fln_scenarios)
        os.rename(igortask.igor_fln_output_pgen, output_fln_pgen)

        igortask.igor_db.connect_db(output_fln_db)
        # copy files
        print("Database file : ", output_fln_db)
        print("airr rearrangement : ", output_fln_airr)
        igortask.igor_db.export_IgorBestScenarios_to_AIRR(output_fln_airr)
        igortask.run_clean_batch()
        igortask._run_clean_batch_mdldata()
    else:
        import os
        # igortask.mdl.write_model(output_fln_parms, output_fln_marginals)
        print("Database file : ", igortask.igor_fln_db)
        base_fln_output = igortask.igor_fln_db.split(".db")[0]
        output_fln_prefix = base_fln_output
        # output_fln_db = output_fln_prefix + ".db"
        output_fln_parms = output_fln_prefix + "_parms.txt"
        output_fln_marginals = output_fln_prefix + "_marginals.txt"
        output_fln_scenarios = output_fln_prefix + "_scenarios.csv"
        output_fln_pgen = output_fln_prefix + "_pgen.csv"
        output_fln_airr = output_fln_prefix + ".tsv"

        os.rename(igortask.igor_fln_output_scenarios, output_fln_scenarios)
        print("scenarios file: ", output_fln_scenarios)
        os.rename(igortask.igor_fln_output_pgen, output_fln_pgen)
        print("pgen file: ", output_fln_pgen)
        print("airr rearrangement : ", output_fln_airr)
        igortask.igor_db.export_IgorBestScenarios_to_AIRR(output_fln_airr)
        igortask.run_clean_batch()
        igortask._run_clean_batch_mdldata()
        # igortask.mdl.plot_Bayes_network(filename=output_fln_prefix + "_BN.pdf")
        # igortask.mdl.export_plot_Pmarginals(output_fln_prefix + "_RM")
        # igortask.mdl.write_model(output_fln_parms, output_fln_marginals)






@click.command("igor-scenarios")
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse")
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA")
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')
@click.option("--Vgene", "fln_genomicVs", default=None,
              help="Filename of genes template for V gene")
@click.option("--Dgene", "fln_genomicDs", default=None,
              help="Filename of genes template for D gene")
@click.option("--Jgene", "fln_genomicJs", default=None,
              help="Filename of genes template for J gene")
@click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
              metavar="<V_gene_CDR3_anchors.csv>",
              help="Anchors filename of V gene")
@click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
              metavar="<J_gene_CDR3_anchors.csv>",
              help="Anchors filename of J gene")
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated')
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None, help="Input sequences in FASTA, TXT or CSV formats.")
@click.option("-o", "--output-prefix", "output_fln_prefix", default=None, help="Output prefix for database file a scenarios file.")
def run_get_scenarios(igor_read_seqs, output_fln_prefix,
                      igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch,
                      igor_fln_db,
                      fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors):
    """IGoR's call to get best scenarios."""
    ########################
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

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

    ########################
    igortask.update_batch_filenames()
    igortask.update_model_filenames()
    print(igortask)

    igortask.load_IgorRefGenome()
    igortask.load_IgorModel()

    click.echo("Running IGoR scenarios process...")
    # igortask.run_evaluate(igor_read_seqs=igor_read_seqs)
    igortask.run_scenarios(igor_read_seqs=igor_read_seqs)
    # cp igortask.igor_fln_output_scenarios

    if output_fln_prefix is not None:
        import os
        output_fln_db = output_fln_prefix + ".db"
        output_fln_parms = output_fln_prefix + "_parms.txt"
        output_fln_marginals = output_fln_prefix + "_marginals.txt"
        # igortask.mdl.plot_Bayes_network(filename=output_fln_prefix + "_BN.pdf")
        # igortask.mdl.export_plot_Pmarginals(output_fln_prefix + "_RM")

        os.rename(igortask.igor_fln_db, output_fln_db)
        os.rename(igortask.igor_fln_output_scenarios, output_fln_prefix+"_scenarios.csv")
        # igortask.mdl.write_model(output_fln_parms, output_fln_marginals)
        # copy files
        print("Database file : ", output_fln_prefix)
    else:
        import os
        # igortask.mdl.write_model(output_fln_parms, output_fln_marginals)
        print("Database file : ", igortask.igor_fln_db)
        base_fln_output = igortask.igor_fln_db.split(".db")[0]
        output_fln_prefix = base_fln_output
        # output_fln_db = output_fln_prefix + ".db"
        output_fln_parms = output_fln_prefix + "_parms.txt"
        output_fln_marginals = output_fln_prefix + "_marginals.txt"
        os.rename(igortask.igor_fln_output_scenarios, output_fln_prefix + "_scenarios.csv")
        print("scenarios file: ", output_fln_prefix + "_scenarios.csv")
        # igortask.mdl.plot_Bayes_network(filename=output_fln_prefix + "_BN.pdf")
        # igortask.mdl.export_plot_Pmarginals(output_fln_prefix + "_RM")
        # igortask.mdl.write_model(output_fln_parms, output_fln_marginals)


@click.command("igor-pgen")
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse")
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA")
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')
@click.option("--Vgene", "fln_genomicVs", default=None,
              help="Filename of genes template for V gene")
@click.option("--Dgene", "fln_genomicDs", default=None,
              help="Filename of genes template for D gene")
@click.option("--Jgene", "fln_genomicJs", default=None,
              help="Filename of genes template for J gene")
@click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
              metavar="<V_gene_CDR3_anchors.csv>",
              help="Anchors filename of V gene")
@click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
              metavar="<J_gene_CDR3_anchors.csv>",
              help="Anchors filename of J gene")
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated')
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
@click.option("-i", "--input-sequences", "igor_read_seqs", default=None,
              help="Input sequences in FASTA, TXT or CSV formats.")
@click.option("-o", "--output-db", "output_db", default=None, help="Output database file.")
def run_get_pgen(igor_read_seqs, output_db,
                 igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch,
                 igor_fln_db,
                 fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors
                 ):
    """IGoR's call to calculate pgen of input sequences"""
    click.echo("Get IGoR pgen process...")
    # TODO: TO RUN THE PGEN I NEED :
    # - INPUT SEQUENCES
    # - ALL THE INFORMATION ABOUT ALIGNMENTS AND MODELS

    ########################
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

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

    ########################
    igortask._run_evaluate(igor_read_seqs=igor_read_seqs)












@click.command("igor-generate")
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse")
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA")
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')
@click.option("--Vgene", "fln_genomicVs", default=None,
              help="Filename of genes template for V gene")
@click.option("--Dgene", "fln_genomicDs", default=None,
              help="Filename of genes template for D gene")
@click.option("--Jgene", "fln_genomicJs", default=None,
              help="Filename of genes template for J gene")
@click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
              metavar="<V_gene_CDR3_anchors.csv>",
              help="Anchors filename of V gene")
@click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
              metavar="<J_gene_CDR3_anchors.csv>",
              help="Anchors filename of J gene")
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated')
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
@click.option("-N","N", default=None, help="Number of sequences to generate.")
@click.option("-o", "--output-prefix", "fln_output_prefix", default=None, help="Prefix for models files.")
def run_generate(N,
                 igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch,
                 igor_fln_db,
                 fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors,
                 fln_output_prefix):
    """IGoR's call to generate sequences"""
    ########################
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

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
        igortask.load_IgorModel()
    elif Q_model_files:
        igortask.igor_model_parms_file = igor_model_parms
        igortask.igor_model_marginals_file = igor_model_marginals
        igortask.load_IgorModel()
    elif igor_fln_db is not None:
        igortask.create_db(igor_fln_db)
        igortask.load_mdl_from_db()
    else:
        print("WARNING: No model provided!")

    ########################
    igortask.update_batch_filenames()
    igortask.update_model_filenames()

    print(igortask.to_dict())
    igortask._run_generate(N_seqs=N)
    igortask._run_clean_batch_mdldata()
    import os
    if fln_output_prefix is None:
        igortask.igor_fln_generated_realizations_werr = None
        igortask.igor_fln_generated_seqs_werr = None
    else:
        igortask.igor_fln_generated_realizations_werr = None
        igortask.igor_fln_generated_seqs_werr = igortask.igor_wd + "/" + igortask.igor_batchname + "_generated/generated_seqs_werr.csv"
        igortask.igor_fln_generated_realizations_werr = igortask.igor_wd + "/" + igortask.igor_batchname + "_generated/generated_realizations_werr.csv"
        output_generated_sequences = fln_output_prefix+"_sequences.csv"
        output_generated_realizations = fln_output_prefix + "_realizations.csv"
        # output_generated_sequences_airr = fln_output_prefix + "_sequences.tsv"
        # TODO: EXPORT IN AIRR REARRANGEMENT
        os.rename(igortask.igor_fln_generated_seqs_werr, output_generated_sequences)
        os.rename(igortask.igor_fln_generated_realizations_werr, output_generated_realizations)
        igortask.run_clean_batch()







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
# @click.group(invoke_without_command=True)
# @click.option("--info", "info", help="List species and chain avialable in imgt website.", is_flag=True, default=False)
# @pass_igortask
# def imgt(igortask, info): #igor_species, igor_chain, igor_model):
#     """Download genomic information from imgt website --info for more details"""
#
#     import pygor3.imgt as p3imgt
#     print(p3imgt.imgt_params['url.home'])
#     species_list = p3imgt.get_species_list()
#
#     if info:
#         click.echo("Downloading data from ... ")
#         print("List of IMGT available species:")
#         print("\n".join(species_list))
#         print("For more details access:")
#         print(p3imgt.imgt_params['url.genelist'])
#





@click.command("imgt-get-genomes")
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse")
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA")
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')
@click.option("--Vgene", "fln_genomicVs", default=None,
              help="Filename of genes template for V gene")
@click.option("--Dgene", "fln_genomicDs", default=None,
              help="Filename of genes template for D gene")
@click.option("--Jgene", "fln_genomicJs", default=None,
              help="Filename of genes template for J gene")
@click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
              metavar="<V_gene_CDR3_anchors.csv>",
              help="Anchors filename of V gene")
@click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
              metavar="<J_gene_CDR3_anchors.csv>",
              help="Anchors filename of J gene")
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated')
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
@click.option("--info", "info", help="List species and chain avialable in imgt website.", is_flag=True, default=False)
@click.option("-t", "--recombination_type", "rec_type", type=click.Choice(['VJ', 'VDJ']),
              help="Igor recombination type.")
@click.option("--imgt-species", "imgt_species", #type=click.Choice(get_species_list()),
              help="IMGT species name for name specifications run imgt --info.")
@click.option("--imgt-chain", "imgt_chain", #type=click.Choice(['VJ', 'VDJ']),
              help="IMGT chain name e.g. TRA, TRB.")
def get_ref_genome(info, rec_type, imgt_species, imgt_chain,
                   igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch,
                   igor_fln_db,
                   fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors):
    """
    Get genomes from imgt website of specifing species and chain in imgt format.
    """

    ########################
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

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
        # print("WARNING: No model provided!")
        pass

    ########################

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




cli.add_command(get_ref_genome)







########### model commands ###########
# @click.group() #invoke_without_command=True)
# @pass_igortask
# def model(igortask):
#     """Manipulations of models"""
#     pass

@click.command("model-export")
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse")
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA")
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')
@click.option("--Vgene", "fln_genomicVs", default=None,
              help="Filename of genes template for V gene")
@click.option("--Dgene", "fln_genomicDs", default=None,
              help="Filename of genes template for D gene")
@click.option("--Jgene", "fln_genomicJs", default=None,
              help="Filename of genes template for J gene")
@click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
              metavar="<V_gene_CDR3_anchors.csv>",
              help="Anchors filename of V gene")
@click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
              metavar="<J_gene_CDR3_anchors.csv>",
              help="Anchors filename of J gene")
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated')
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
@click.option("--from-txt", "fln_from_txt", nargs=2, default=[None, None],
              help="Export Igor's model from txt files model_parms.txt and model_marginals.txt.")
@click.option("--from-db", "fln_from_db", nargs=1, default=None,
              help="Export Igor's model from database file.")
@click.option("--to-txt", "fln_to_txt", nargs=2, default=[None, None],
              help="Output filename of Igor recombination model to  <model_parms.txt> <model_marginals.txt>.")
@click.option("--to-db", "fln_to_db", nargs=1, default=None,
              help="Output filename of Igor recombination model to  <model.db>.")
def model_export(fln_from_txt, fln_from_db, fln_to_txt, fln_to_db,
                 igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch,
                 igor_fln_db,
                 fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors):
    """
    Export IGoR's models from txt (model_parms.txt, model_marginals.txt) files to db viceversa
    """
    ########################
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

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
    ########################

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
                mdl_to.write_model(fln_to_txt[0], fln_to_txt[1])
                # mdl_to.parms.write_model_parms(fln_to_txt[0])
                # mdl_to.marginals.write_model_marginals(fln_to_txt[1])
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
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse")
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA")
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')
@click.option("--Vgene", "fln_genomicVs", default=None,
              help="Filename of genes template for V gene")
@click.option("--Dgene", "fln_genomicDs", default=None,
              help="Filename of genes template for D gene")
@click.option("--Jgene", "fln_genomicJs", default=None,
              help="Filename of genes template for J gene")
@click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
              metavar="<V_gene_CDR3_anchors.csv>",
              help="Anchors filename of V gene")
@click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
              metavar="<J_gene_CDR3_anchors.csv>",
              help="Anchors filename of J gene")
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated')
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
@click.option("-t", "--recombination_type", "rec_type", type=click.Choice(['VJ', 'VDJ']),
              help="Igor recombination type.", required=True)
@click.option("-o", "--output-prefix", "fln_output_prefix", default=None, help="Prefix for models files.")
def model_create(rec_type, fln_output_prefix,
                 igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch,
                 igor_fln_db,
                 fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors):
    """Make a new VJ or VDJ model with uniform probability distribution"""
    # click.echo( "igor_model_dir_path: "+igortask.igor_model_dir_path )
    ########################
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

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
    elif igor_model_path is not None:
        igortask.igor_model_dir_path = igor_model_path
        igortask.update_model_filenames(igor_model_dir_path=igortask.igor_model_dir_path)
    else:
        print("WARNING: No model provided!")
    ########################
    import pygor3 as p3
    if rec_type == 'VJ':
        # load genomics
        # igortask.igor_model_dir_path
        igortask.update_model_filenames(igor_model_dir_path=igortask.igor_model_dir_path)
        igortask.load_IgorRefGenome()
        igortask.make_model_default_VJ_from_genomes_dir()
    elif rec_type == 'VDJ':
        igortask.update_model_filenames(igor_model_dir_path=igortask.igor_model_dir_path)
        # print(igortask.to_dict())
        igortask.load_IgorRefGenome()
        igortask.make_model_default_VDJ_from_genomes_dir()
        # print(igortask.igor_model_dir_path)
    else:
        print("Model type ", rec_type)

    print("igortask.igor_model_dir_path: ", igortask.igor_model_dir_path)
    if fln_output_prefix is None:
        # print(igortask.igor_model_dir_path)
        igortask.mdl.parms.write_model_parms(igortask.igor_model_parms_file)
        igortask.mdl.marginals.write_model_marginals(igortask.igor_model_marginals_file, igortask.mdl.parms)
    else:
        igortask.igor_model_parms_file = fln_output_prefix + "_parms.txt"
        igortask.igor_model_marginals_file = fln_output_prefix + "_marginals.txt"
        igortask.mdl.parms.write_model_parms(igortask.igor_model_parms_file)
        igortask.mdl.marginals.write_model_marginals(igortask.igor_model_marginals_file, igortask.mdl.parms)








@click.command("model-plot")
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse")
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA")
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')
@click.option("--Vgene", "fln_genomicVs", default=None,
              help="Filename of genes template for V gene")
@click.option("--Dgene", "fln_genomicDs", default=None,
              help="Filename of genes template for D gene")
@click.option("--Jgene", "fln_genomicJs", default=None,
              help="Filename of genes template for J gene")
@click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
              metavar="<V_gene_CDR3_anchors.csv>",
              help="Anchors filename of V gene")
@click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
              metavar="<J_gene_CDR3_anchors.csv>",
              help="Anchors filename of J gene")
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated')
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
@click.option("-o", "--output-prefix", "fln_output_prefix", default=None, help="Prefix to pdf files with model plots.")
def model_plot(fln_output_prefix,
               igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch,
               igor_fln_db,
               fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors):
    """Plot real marginals of the bayesian network events """
    ########################
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

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
    elif igor_model_path is not None:
        igortask.igor_model_dir_path = igor_model_path
        igortask.update_model_filenames(igor_model_dir_path=igortask.igor_model_dir_path)
    else:
        print("WARNING: No model provided!")
    ########################

    if igortask.igor_fln_db is not None:
        igortask.create_db(igortask.igor_fln_db)
        igortask.load_mdl_from_db()
        print("Model loaded from ", igortask.igor_fln_db)
    else:
        try:
            igortask.update_model_filenames()
        except Exception as e:
            print("WARNING: update_model_filenames ", e)

        igortask.load_IgorModel()

    igortask.mdl.export_plot_events(fln_output_prefix+"_CP")
    igortask.mdl.export_plot_Pmarginals(fln_output_prefix+"_MP")


cli.add_command(model_export)
cli.add_command(model_create)
cli.add_command(model_plot)







# pygor3-cli [GENERAL_OPTIONS] database export all
########### database commands ###########
@click.command("db-import") #invoke_without_command=True)
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse")
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA")
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')
@click.option("-g", "--set_path_ref_genome", "igor_path_ref_genome", default=None,
                    help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv')
@click.option("--Vgene", "fln_genomicVs", default=None,
              help="Filename of genes template for V gene")
@click.option("--Dgene", "fln_genomicDs", default=None,
              help="Filename of genes template for D gene")
@click.option("--Jgene", "fln_genomicJs", default=None,
              help="Filename of genes template for J gene")
@click.option("--Vanchors", "fln_V_gene_CDR3_anchors", default=None,
              metavar="<V_gene_CDR3_anchors.csv>",
              help="Anchors filename of V gene")
@click.option("--Janchors", "fln_J_gene_CDR3_anchors", default=None,
              metavar="<J_gene_CDR3_anchors.csv>",
              help="Anchors filename of J gene")
@click.option("-w", "--set_wd", "igor_wd", help="Sets the working directory to path", default='./', show_default=True)
# To load all data files use batchname
@click.option("-b", "--set_batch", "igor_batch", default=None,
                    help='Sets batchname to identify run. If not set random name is generated')
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
@click.option("-o", "--output-db", "fln_output_db", default=None, help="output attached database.")
def database_import(fln_output_db,
                    igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd, igor_batch,
                    igor_fln_db,
                    fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors):
    """
    Import igor files to database
    """
    ########################
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

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
    ########################

    # load a task object to get all data in one structure
    igortask.update_batch_filenames()
    print(igortask.igor_fln_db)
    igortask.create_db(igor_fln_db=fln_output_db)

    try:
        igortask.load_db_from_indexed_sequences()
    except Exception as e:
        print("WARNING: indexed_sequences not found.", e)

    try:
        igortask.load_db_from_indexed_cdr3()
    except Exception as e:
        print("WARNING: indexed_sequences not found.", e)


    import pygor3 as p3
    if (igortask.igor_species is not None) and (igortask.igor_chain is not None):
        print("species : ", igortask.igor_species, " chain: ", igortask.igor_chain)
        try:
            igortask.run_datadir()
            igor_model_path = igortask.igor_models_root_path + igortask.igor_species + "/" \
                              + p3.igor_option_path_dict[igortask.igor_chain] + "/"
            # igor_model_path = task.igor_models_root_path + task.igor_species + "/" \
            #                 + p3.igor_option_path_dict[task.igor_chain] + "/"
            igortask.igor_model_dir_path = igor_model_path
        except Exception as e:
            print("WARNING: Default model not found!", e)

    # IF MODEL_PATH PROVIDED THEN
    # task.igor_model_dir_path = args.model_path
    if (igortask.igor_model_dir_path is not None):
        try:
            igortask.update_model_filenames(igortask.igor_model_dir_path)
            # task.igor_model_parms_file = task.igor_model_dir_path+"/models/model_parms.txt"
            # task.igor_model_marginals_file = task.igor_model_dir_path+"/models/model_marginals.txt"
            # task.igor_path_ref_genome = task.igor_model_dir_path+"/ref_genome/"
            # task.load_IgorModel()
            # args.path_ref_genome = igortask.igor_path_ref_genome
            igortask.load_IgorModel()
            igortask.load_db_from_models()
        except Exception as e:
            print("Couldn't load models to database")
            print(e)

    # IF GENOME PATH PROVIDED
    if igortask.igor_path_ref_genome is not None:
        try:
            # igortask.igor_path_ref_genome = args.path_ref_genome  # "/home/alfaceor/Dropbox/PosDoc/IGoR/dev/MyGithub/pygor3/demo/thi/genomics_repseqio_F"
            igortask.load_IgorRefGenome()
            # print(task.igor_fln_indexed_sequences)
            # print(task.genomes.df_genomicVs) # is where all this data is collected
            # print(task.genomes.df_genomicDs)
            # print(task.genomes.df_genomicJs)

            igortask.load_db_from_genomes()
        except Exception as e:
            print("Couldn't load genome templates to database")
            print("ERROR: ", e)

    try:
        igortask.load_db_from_alignments()
    except Exception as e:
        print("WARNING: No alignments found.", e)

    try:
        igortask.load_db_from_bestscenarios()
    except Exception as e:
        print("WARNING: No alignments found.", e)

    try:
        igortask.load_db_from_pgen()
    except Exception as e:
        print("WARNING: No alignments found.", e)











@click.command("db-cp") #invoke_without_command=True)
############## COMMON options ##############
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
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
def database_copy(fln_output_db, b_igor_reads,
                  b_igor_genomes,
                  b_igor_alignments,
                  b_igor_model,
                  b_igor_scenarios, b_igor_pgen,
                  igor_fln_db):
    """Testing function before commit - export database"""
    ########################
    from pygor3 import IgorTask
    igortask = IgorTask()

    igortask.igor_fln_db = igor_fln_db
    igortask.create_db()

    ########################

    # Create a database
    from pygor3.IgorSqliteDB import IgorSqliteDB
    # OPEN CONEXION TO output_db
    output_db = IgorSqliteDB.create_db(fln_output_db)

    # Get list of tables in igortask.igor_db
    tablename_ctsql_dict = igortask.igor_db.get_dict_of_Igortablename_sql()
    print("**** Tables in source database : ", igor_fln_db) #, tablename_ctsql_dict.keys())
    igortask.igor_db.list_from_db()
    if b_igor_reads and igortask.igor_db.Q_sequences_in_db():
        try:
            # Ask to sqlite_master for the way to create it in new database
            tablename_to_copy = 'IgorIndexedSeq'
            # Create table in database destiny.
            output_db.execute_query(tablename_ctsql_dict[tablename_to_copy])
            # Copy table
            fln_source_db = igortask.igor_db.fln_db
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
                    fln_source_db = igortask.igor_db.fln_db
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
                    fln_source_db = igortask.igor_db.fln_db
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
                    fln_source_db = igortask.igor_db.fln_db
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
                    fln_source_db = igortask.igor_db.fln_db
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
                    fln_source_db = igortask.igor_db.fln_db
                    output_db.copytable_from_source(tablename_to_copy, fln_source_db)
        except Exception as e:
            print("ERROR: igor-scenarios")
            print(e)

    print("**** Tables in destiny database: ", fln_output_db)
    output_db.list_from_db()















@click.command("db-ls")
############## COMMON options ##############
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")
############## NO COMMON options ##############
def database_ls(igor_fln_db):
    """List tables in database by groups and show  number of records."""
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_fln_db = igor_fln_db
    igortask.create_db()
    igortask.igor_db.list_from_db()

    # Q_seqs = igortask.igor_db.Q_output_in_db()
    # print(Q_seqs)

    # igortask.igor_db.write_IgorIndexedSeq_to_CSV("jojo.csv")











@click.command("db-attach")
############## COMMON options ##############
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
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
def database_attach(fln_from_db,
                    igor_fln_db,):
    """
    Attach tables to database.
    """

    ########################
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_path_ref_genome = igor_path_ref_genome
    igortask.fln_genomicVs = fln_genomicVs
    igortask.fln_genomicDs = fln_genomicDs
    igortask.fln_genomicJs = fln_genomicJs
    igortask.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
    igortask.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

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
    ########################

    igortask.create_db()
    from pygor3 import IgorSqliteDB
    other_igor_db = IgorSqliteDB.create_db(fln_from_db)
    dicto = other_igor_db.get_dict_of_Igortablename_sql()
    print(dicto.keys())
    # import pygor3 as p3
    import pygor3 as p3
    print("sql_tablename_patterns_dict : ", p3.sql_tablename_patterns_dict)
    other_igor_db.close_db()
    # igortask.igor_db.attach_table_from_db("ddd")
    # igortask.igor_db.write_IgorIndexedSeq_to_CSV("jojo.csv")












@click.command("db-rm")
############## COMMON options ##############
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
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
def database_rm(b_igor_reads, b_igor_model,
                b_igor_model_parms, b_igor_model_marginals,
                b_igor_scenarios, b_igor_pgen,
                b_igor_genomes,
                b_igor_genomesV, b_igor_genomesD, b_igor_genomesJ, b_igor_genomesCDR3,
                b_igor_alignments,
                b_igor_alignmentsV, b_igor_alignmentsD, b_igor_alignmentsJ,
                b_igor_alignmentsCDR3,
                igor_fln_db):
    """
    Delete tables in database by groups.
    """
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_fln_db = igor_fln_db
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















@click.command("db-export", cls=RegisterWriterCommand)
############## COMMON options ##############
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")
@click.option("-b", "--set_batch", "igor_batch", default=None, type=str, required=True,
                    help='Sets batchname to identify run. If not set random name is generated')
############## NO COMMON options ##############
@click.option("--airr-rearrangement", "b_airr_rearrangement", cls=RegisterReaderOption,
              is_flag=True,
              default=False, help='Export IGoR scenarios and pgen in tsv airr format file')
@click.option("--airr-rearrangement_set", "fln_airr_rearrangement", cls=RegisterWriterOption,
              default=None,
              help='Export to airr file.')
@click.option("--igor-all", "b_igor_all", is_flag=True,
              default=False,
              help='Export all available data in db with prefix.')
@click.option("--igor-reads", "b_igor_reads", cls=RegisterReaderOption,
              default=False, is_flag=True,
              help='Export sequences reads to batch.')
@click.option("--igor-reads_set", "fln_igor_reads", cls=RegisterWriterOption,
              default=None,
              help='Export sequences reads to filename.')
@click.option("--igor-genomes", "b_igor_genomes", cls=RegisterReaderOption,
              default=False, is_flag=True,
              help='Export genomes to batch.')
@click.option("--igor-genomes_set", "fln_igor_genomes", cls=RegisterWriterOption,
              default=None,
              help='Export genomes to filename prefix.')

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
@click.option("--igor-genomesCDR3", "b_igor_genomesCDR3", cls=RegisterReaderOption,
              default=False, is_flag=True,
              help='Export CDR3 anchors to batch.')
@click.option("--igor-genomesCDR3_set", "fln_igor_genomesCDR3", cls=RegisterWriterOption,
              default=None,
              help='Export CDR3 anchors to filename prefix.')
# @click.option("--igor-alignments", "b_igor_alignments", is_flag=True,
#               help='Copy all available alignments tables to database.')
@click.option("--igor-alignments", "b_igor_alignments", cls=RegisterReaderOption,
              default=False, is_flag=True,
              help='Export alignments to batch.')
@click.option("--igor-alignments_set", "fln_igor_alignments", cls=RegisterWriterOption,
              default=None,
              help='Export alignments to filename prefix.')
# @click.option("--igor-alignmentsV", "b_igor_alignmentsV", is_flag=True,
#               help='Copy V alignments tables to database.')
# @click.option("--igor-alignmentsD", "b_igor_alignmentsD", is_flag=True,
#               help='Copy D alignments tables to database.')
# @click.option("--igor-alignmentsJ", "b_igor_alignmentsJ", is_flag=True,
#               help='Copy J alignments tables to database.')
# @click.option("--igor-alignmentsCDR3", "b_igor_alignmentsCDR3", is_flag=True,
#               help='Copy indexed cdr3 table to database.')
@click.option("--igor-alignmentsCDR3", "b_igor_alignmentsCDR3", cls=RegisterReaderOption,
              default=False, is_flag=True,
              help='Export CDR3 segments from best alignments to batch.')
@click.option("--igor-alignmentsCDR3_set", "fln_igor_alignmentsCDR3", cls=RegisterWriterOption,
              default=None,
              help='Export segments from best alignments to filename.')
# @click.option("--igor-model", "b_igor_model", is_flag=True,
#               help='IGoR model database file.')
@click.option("--igor-model", "b_igor_model", cls=RegisterReaderOption,
              default=False, is_flag=True,
              help="Export IGoR's model to batch.")
@click.option("--igor-model_set", "fln_igor_model", cls=RegisterWriterOption,
              default=None,
              help="Export IGoR's model to filename prefix.")
# @click.option("--igor-model-parms", "b_igor_model_parms", is_flag=True,
#               help='IGoR model parms (or params) file.')
# @click.option("--igor-model-marginals", "b_igor_model_marginals", is_flag=True,
#               help='IGoR model marginals file.')
# @click.option("--igor-scenarios", "b_igor_scenarios", is_flag=True,
#               help='If --from-db no need to add filename')
@click.option("--igor-scenarios", "b_igor_scenarios", cls=RegisterReaderOption,
              default=False, is_flag=True,
              help="Export IGoR's scenarios to batch.")
@click.option("--igor-scenarios_set", "fln_igor_scenarios", cls=RegisterWriterOption,
              default=None,
              help="Export IGoR's scenarios to filename prefix.")
# @click.option("--igor-pgen", "b_igor_pgen", is_flag=True,
#               help='If --from-db no need to add filename')
@click.option("--igor-pgen", "b_igor_pgen", cls=RegisterReaderOption,
              default=False, is_flag=True,
              help="Export IGoR's pgen to batch.")
@click.option("--igor-pgen_set", "fln_igor_pgen", cls=RegisterWriterOption,
              default=None,
              help="Export IGoR's pgen to filename prefix.")
def database_export(b_airr_rearrangement, fln_airr_rearrangement,
                    b_igor_all, b_igor_reads, fln_igor_reads,
                    b_igor_genomes, fln_igor_genomes, b_igor_genomesCDR3, fln_igor_genomesCDR3,
                    b_igor_alignments, fln_igor_alignments,
                    b_igor_alignmentsCDR3, fln_igor_alignmentsCDR3,
                    b_igor_model, fln_igor_model,
                    b_igor_scenarios, fln_igor_scenarios,
                    b_igor_pgen, fln_igor_pgen,
                    # igor_species, igor_chain, igor_model, igor_model_path, igor_path_ref_genome, igor_wd,
                    igor_batch,
                    # fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors,
                    igor_fln_db):
    """Export database model in igor formatted files"""
    # fln_prefix = "tmp_export_"
    ########################
    print("db-export")
    from pygor3 import IgorTask
    igortask = IgorTask(igor_fln_indexed_sequences=fln_igor_reads, igor_batchname=igor_batch)

    # Q_species_chain = (not (igor_species is None) and not (igor_chain is None))
    # Q_model_files = (not (igor_model_parms is None) and not (igor_model_marginals is None))
    # if Q_species_chain:
    #     igortask.igor_species = igor_species
    #     igortask.igor_chain = igor_chain
    # elif Q_model_files:
    #     igortask.igor_model_parms_file = igor_model_parms
    #     igortask.igor_model_marginals_file = igor_model_marginals
    # elif igor_fln_db is not None:
    #     igortask.create_db()
    # else:
    #     print("WARNING: No model provided!")
    ########################
    # print("Export from database")
    if b_igor_all:
        if igor_batch is None:
            print("ERROR: Batchname is needed. Use -b <batchname> or --set_batch <batchname>")
            return 0

    if fln_airr_rearrangement is not None:
        b_airr_rearrangement = True

    if fln_igor_reads is not None:
        b_igor_reads = True

    if fln_igor_genomes is not None:
        b_igor_genomes = True

    if fln_igor_genomesCDR3 is not None:
        b_igor_genomesCDR3 = True

    if fln_igor_alignments is not None:
        b_igor_alignments = True

    if fln_igor_alignmentsCDR3 is not None:
        b_igor_alignmentsCDR3 = True

    if fln_igor_model is not None:
        b_igor_model = True

    if fln_igor_scenarios is not None:
        b_igor_scenarios = True

    if fln_igor_pgen is not None:
        b_igor_pgen = True

    print("batchname: ", igortask.igor_batchname)
    print("b_igor_reads, fln_igor_reads", b_igor_reads, fln_igor_reads)
    print("b_igor_genomes, fln_igor_genomes", b_igor_genomes, fln_igor_genomes)

    if igortask.igor_batchname is None:
        print("ERROR: batch option required.")
        return 0

    # TODO: IgorTask should help me to pass from batch to db and from db to batch and model_directory.
    # igor_fln_db = igortask.igor_fln_db
    igortask.update_batch_filenames()
    path_mdl_data = igortask.igor_batchname + "_mdldata"
    igortask.update_model_filenames(igor_model_dir_path=path_mdl_data)


    igortask.create_db(igor_fln_db)
    ii_db = igortask.igor_db
    ii_db.list_from_db()
    print("igor_fln_db: " + igor_fln_db, "ii_db : ",
          ii_db, type(ii_db), type(igortask.igor_db))

    try:
        # import pygor3 as p3
        # ii_db = p3.IgorSqliteDB()
        igortask.db_export_IgorIndexedSeq()
        ii_db.write_IgorIndexedSeq_to_CSV(igortask.igor_fln_indexed_sequences)
    except Exception as e:
        print(igortask)
        raise e

    ################################################3
    try:
        igortask.update_ref_genome()
    except Exception as e:
        print(igortask)
        raise e

    ################################################3

    # igortask.genomes.update_fln_names(path_ref_genome=igortask.igor_path_ref_genome)

    if fln_igor_genomes is not None:
        b_igor_genomes = True
        igortask.fln_genomicVs = fln_igor_genomes + "Vs.fasta"
        igortask.fln_genomicJs = fln_igor_genomes + "Js.fasta"
        igortask.fln_genomicDs = fln_igor_genomes + "Ds.fasta"

    if b_igor_all or b_igor_genomes :
        if ii_db.Q_ref_genome_in_db_by_gene("V"):
            igortask.fln_genomicVs = igortask.genomes.fln_genomicVs
            ii_db.write_IgorGeneTemplate_to_fasta("V", igortask.fln_genomicVs)
            try:
                ii_db.write_IgorGeneAnchors_to_CSV("V", igortask.fln_V_gene_CDR3_anchors)
            except Exception as e:
                pass
        if ii_db.Q_ref_genome_in_db_by_gene("J"):
            igortask.fln_genomicJs = igortask.genomes.fln_genomicJs
            ii_db.write_IgorGeneTemplate_to_fasta("J", igortask.fln_genomicJs)
            try:
                ii_db.write_IgorGeneAnchors_to_CSV("J", igortask.fln_J_gene_CDR3_anchors)
            except Exception as e:
                pass
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

    # igortask.export_to_igorfiles()

    if fln_airr_rearrangement is not None:
        b_airr_rearrangement = True
    else:
        fln_prefix = igortask.igor_fln_db.split(".db")[0]
        fln_airr_rearrangement =  fln_prefix + ".tsv"


    if b_airr_rearrangement:
        igortask.igor_db.export_IgorBestScenarios_to_AIRR(fln_airr_rearrangement)

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
############## COMMON options ##############
############## COMMON options ##############
@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")

############## NO COMMON options ##############
@click.option("-o", "--output-filename", "output_fln", default=None, help="Filename of output file csv or fasta if --seq_index is use.")
@click.option("-s", "--seq-index", "seq_index", default=None, type=int, help="Sequence id (seq_index) in db.", required=False)
def database_naive_align(output_fln, seq_index, igor_fln_db):
    """
    Get a naive alignment (no scenarios) from igor alignments
    """
    from pygor3 import IgorTask
    igortask = IgorTask()
    igortask.igor_fln_db = igor_fln_db

    import pygor3 as p3

    if seq_index is None:
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
        db.fln_db = igortask.igor_db.fln_db #args.database
        db.connect_db()

        # Make a loop over all sequences
        if output_fln is None:
            fln_output = db.fln_db.split(".db")[0] + "_na.csv"
        else:
            fln_output = output_fln

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
    else:
        ###############################################
        def generate_str_fasta(indexed_sequence, list_vdj_alignments: dict):
            """ Given an Sequence index and the corresponding alignments vj/ vdj
            return a string with considering only offset"""

            indexed_sequence.sequence = indexed_sequence.sequence.lower()
            # add mismatches in sequence.
            s = list(indexed_sequence.sequence)
            for key_align in list_vdj_alignments.keys():
                for pos_mis in list_vdj_alignments[key_align].mismatches:
                    s[pos_mis] = s[pos_mis].upper()
            indexed_sequence.sequence = "".join(s)

            str_fasta = ""
            min_offset_key = min(list_vdj_alignments.keys(), key=lambda x: list_vdj_alignments[x].offset)  # .offset
            min_offset = list_vdj_alignments[min_offset_key].offset
            min_offset = min(indexed_sequence.offset, min_offset)

            delta_offset = indexed_sequence.offset - min_offset
            str_prefix = '-' * (delta_offset)
            str_fasta_sequence = str_prefix + indexed_sequence.sequence
            print(str_fasta_sequence)
            str_fasta = str_fasta + "> " + str(indexed_sequence.seq_index) + "\n"
            str_fasta = str_fasta + str_fasta_sequence + "\n"
            for key in list_vdj_alignments.keys():
                list_vdj_alignments[key].strGene_seq = list_vdj_alignments[key].strGene_seq.lower()
                delta_offset = list_vdj_alignments[key].offset - min_offset
                str_prefix = '-' * (delta_offset)
                str_fasta_sequence = str_prefix + list_vdj_alignments[key].strGene_seq
                print(str_fasta_sequence)
                str_fasta = str_fasta + "> " + key + ", " + list_vdj_alignments[key].strGene_name + "\n"
                str_fasta = str_fasta + str_fasta_sequence + "\n"
                offset_5_p = list_vdj_alignments[key].offset_5_p - min_offset
                offset_3_p = list_vdj_alignments[key].offset_3_p - min_offset
                print("delta_offset : ", delta_offset)
                print("offset_5_p : ", list_vdj_alignments[key].offset_5_p, offset_5_p)
                print("offset_3_p : ", list_vdj_alignments[key].offset_3_p, offset_3_p)
                str_prefix_2 = '-' * (offset_5_p + 1)
                str_fasta_sequence2 = str_prefix_2 + str_fasta_sequence[offset_5_p + 1:offset_3_p + 1]

                str_fasta = str_fasta + "> " + list_vdj_alignments[key].strGene_name + ", score : " + str(
                    list_vdj_alignments[key].score) + "\n"
                str_fasta = str_fasta + str_fasta_sequence2 + "\n"

                # TODO ADD MISMATCHES
                align = list_vdj_alignments[key]
                # align mismatches are in indexed sequence reference I need to convert it to gene reference given the alignment
                # given the align.offset
                # pos_in_gene  = pos_in_seq - align.offset
                # pos_in_gene = cdr3 - align.offset

            return str_fasta

        def generate_csv_line(indexed_sequence: p3.IgorIndexedSequence, list_vdj_alignments: dict, sep=';',
                              header_list=['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score',
                                           'd_score', 'j_score', 'junction']):
            csv_line = ""
            indexed_sequence.sequence = indexed_sequence.sequence.lower()
            fields_dict = dict()
            fields_dict['sequence_id'] = indexed_sequence.seq_index
            fields_dict['sequence'] = indexed_sequence.sequence
            fields_dict['v_call'] = list_vdj_alignments['V'].strGene_name
            fields_dict['d_call'] = list_vdj_alignments['D'].strGene_name
            fields_dict['j_call'] = list_vdj_alignments['J'].strGene_name
            fields_dict['v_score'] = list_vdj_alignments['V'].score
            fields_dict['d_score'] = list_vdj_alignments['D'].score
            fields_dict['j_score'] = list_vdj_alignments['J'].score

            align = p3.IgorAlignment_data()

            for field in header_list:
                csv_line = csv_line + fields_dict[field] + sep
            return csv_line

            # mixcr :
            # targetSequences	targetQualities	allVHitsWithScore
            # allDHitsWithScore	allJHitsWithScore	allCHitsWithScore	allVAlignments	allDAlignments
            # allJAlignments	allCAlignments	nSeqFR1	minQualFR1	nSeqCDR1	minQualCDR1	nSeqFR2	minQualFR2	nSeqCDR2	minQualCDR2	nSeqFR3	minQualFR3	nSeqCDR3	minQualCDR3	nSeqFR4	minQualFR4	aaSeqFR1	aaSeqCDR1	aaSeqFR2	aaSeqCDR2	aaSeqFR3	aaSeqCDR3	aaSeqFR4	refPoints

            # AIRR :
            # sequence_id (req)
            # sequence (req)
            # rev_comp (req) //false by default # FIXME: for now
            # productive (req) // false by default # FIXME: a naive alignmnet of productive is possible.
            # locus (chain type) // eg. TCR IG
            # v_call
            # d_call
            # d2_call // FIXME: not gonna use this by the moment.
            # j_call
            # sequence alignment // FIXME: Aligned portion of query sequence. IMGT-gaps check the best way to do this
            # germline_alignment // FIXME: Assembled, aligned, full-length inferred germline sequence spanning the same region as the sequence_alignment field (typically the V(D)J region) and including the same set of corrections and spacers (if any)
            # junction
            # junction_aa // Our CDR3
            # np1 (nts btwn V and D or btwn V and J  // vd_ins  or vj_ins
            # np2 (nts btwn first D and J gene or btwn first D and second D // dj_ins or d1d2_ins
            # np3 (nts btwn second D and J gene
            # v_score
            # v_cigar
            # d_score
            # d_cigar
            # j_score
            # j_cigar
            # junction_length

        def generate_str_fasta_simple(indexed_sequence, list_vdj_alignments: dict):
            """ Given an Sequence index and the corresponding alignments vj/ vdj
            return a string with considering only offset"""

            str_fasta = ""
            str_fasta = str_fasta + "> " + str(indexed_sequence.seq_index) + "\n"
            str_fasta = str_fasta + indexed_sequence.sequence + "\n"
            for key in list_vdj_alignments.keys():
                str_fasta = str_fasta + "> " + key + ", " + list_vdj_alignments[key].strGene_name + "\n"
                str_fasta = str_fasta + list_vdj_alignments[key].strGene_seq + "\n"

            # print(list_vdj_alignments['V'])
            # print(list_vdj_alignments['D        '])
            return str_fasta

        ############################################
        #### naive align with id seq_index
        db = p3.IgorSqliteDB()

        db.fln_db = igortask.igor_fln_db
        db.connect_db()

        #seq_index = args.seq_index
        try:
            indexed_sequence = db.get_IgorIndexedSeq_By_seq_index(seq_index)
        except Exception as e:
            print("ERROR: seq_index ", seq_index, " not found in db ", igortask.igor_fln_db)
            raise
        indexed_sequence.offset = 0

        best_v_align_data = db.get_best_IgorAlignment_data_By_seq_index('V', indexed_sequence.seq_index)
        best_j_align_data = db.get_best_IgorAlignment_data_By_seq_index('J', indexed_sequence.seq_index)

        try:
            best_d_align_data = db.get_best_IgorAlignment_data_By_seq_index('D', indexed_sequence.seq_index)
            vdj_naive_alignment = {'V': best_v_align_data,
                                   'D': best_d_align_data,
                                   'J': best_j_align_data}
            v_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
            print('V', len(v_align_data_list), [ii.score for ii in v_align_data_list])
            d_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('D', indexed_sequence.seq_index)
            print('D', len(d_align_data_list), [ii.score for ii in d_align_data_list])
            j_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('J', indexed_sequence.seq_index)
            print('J', len(j_align_data_list), [ii.score for ii in j_align_data_list])
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

        except Exception as e:
            print(e)
            print("No d gene alignments found!")
            vdj_naive_alignment = {'V': best_v_align_data,
                                   'J': best_j_align_data}
            v_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
            print('V', len(v_align_data_list), [ii.score for ii in v_align_data_list])
            j_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('J', indexed_sequence.seq_index)
            print('J', len(j_align_data_list), [ii.score for ii in j_align_data_list])
            pass

        if output_fln is None:
            batchname = db.fln_db.split(".db")[0]
            fln_output = batchname + '__' + str(indexed_sequence.seq_index) + '_na.fasta'
            # fln_output = args.database.split(".db")[0]+"_na.csv"
        else:
            fln_output = output_fln

        ofile = open(fln_output, 'w')
        str_fasta = generate_str_fasta(indexed_sequence, vdj_naive_alignment)
        ofile.write(str_fasta)
        ofile.close()




# TODO;
def database_plot_pgen():
    """Plot Pgen distribution"""
    import pygor3 as p3
    import numpy as np

    igortask = p3.IgorTask()
    igortask.igor_db.fetch_IgorPgen()

    igortask.get_pgen_pd()

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    ###### BEGIN PLOT DECORATION VARIABLES
    font = {'family': 'normal',
            'weight': 'bold',
            'size': 18}
    plt.rc('font', **font)
    plt.rc('text', usetex=True)
    ###### END PLOT DECORATION VARIABLES
    fig, ax = plt.subplots()
    ax.set_xlabel("CDR3nt length")
    ax.set_ylabel("counts")

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-d", "--density", dest="density")
    (options, args) = parser.parse_args()

    for flnIGoRCDR3 in args:
        # flnIGoRCDR3 = "Barb_indexed_CDR3.csv"
        try:
            lblHist = flnIGoRCDR3.split("_indexed_CDR3.csv")[0]
        except:
            lblHist = "No Title"
            pass
        dataCDR3 = pd.read_csv(flnIGoRCDR3, sep=';')  # , sep=";")
        flnIGoRCDR3_hist = flnIGoRCDR3.split(".csv")[0] + "_hist" + ".csv"
        flnIGoRCDR3_hist

        dataCDR3LenArr = dataCDR3['CDR3nt'].dropna().map(len).values
        hist, bin_edges = np.histogram(dataCDR3LenArr, bins=range(dataCDR3LenArr.min(), dataCDR3LenArr.max()))
        xbins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        print
        bin_edges[:-1].dtype
        print
        type(bin_edges[1:])
        xbins
        np.savetxt(flnIGoRCDR3_hist, zip(bin_edges[:-1], bin_edges[1:], hist), fmt='%d')
        ax.plot(xbins, hist, marker='o', label=lblHist)

    ax.legend()
    fig.tight_layout()
    plt.show()


# TODO:
def database_plot_CDR3_len():
    """Plot CDR3 length distribution from alignments"""
    pass



# database.add_command(database_create)
cli.add_command(database_import)
cli.add_command(database_copy)
cli.add_command(database_ls)
cli.add_command(database_attach)
cli.add_command(database_rm)
cli.add_command(database_export)
cli.add_command(database_naive_align)
# cli.add_command(database_plot_pgen)
# cli.add_command(database_plot_CDR3_len)



@click.command("test") #invoke_without_command=True)
############## COMMON options ##############
@click.option("-s", "--set_igor_species", "igor_species", default=None, help="Species in IGoR's format: human, mouse")
@click.option("-c", "--set_igor_chain", "igor_chain", default=None, help="Chain in IGoR's format, e.g. alpha, beta, TRB, TRA")
@click.option("-m", "--set_igor_model", "igor_model", default=[None, None], nargs=2,
              metavar="<model_parms.txt> <model_marginals.txt>",
              help='IGoR model_params.txt and model_marginals.txt filenames.')

@click.option("-M", "--set_model_path", "igor_model_path", default=None,
              metavar="<model_directory_path>",
              help='IGoR model directory path, this path include ref_genomes and model_parms')

@click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")
############## COMMON options ##############
@click.option("--igor-scenarios", "igor_fln_output_scenarios", default=None, help="igor generated scenarios filename.", required=True)
@click.option("--airr-scenarios", "airr_fln_output_scenarios", default=None, help="AIRR formatted scenarios filename.", required=True)
############## NO COMMON options ##############
def test(igor_species, igor_chain, igor_model, igor_model_path, igor_fln_db,
         igor_fln_output_scenarios, airr_fln_output_scenarios):
    """
    Get a naive alignment (no scenarios) from igor alignments
    """
    from pygor3 import IgorTask

    if not (igor_fln_db is None):
        igortask_input = IgorTask(igor_fln_db=igor_fln_db)


    igortask = IgorTask()
    igor_model_parms = igor_model[0]
    igor_model_marginals = igor_model[1]
    igortask.igor_model_parms_file = igor_model[0]
    igortask.igor_model_marginals_file = igor_model[1]
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
        igortask.create_db(igor_fln_db=igor_fln_db)
    else:
        print("WARNING: No model provided!")

    # igortask.load_IgorModel()
    igortask.load_mdl_from_db()

    igortask.igor_db.export_IgorBestScenarios_to_AIRR("nose.tsv")



@click.command("demo-get-data") #invoke_without_command=True)
def demo_get_data():
    """
    Copy demo directory in current directory
    """
    from importlib import resources
    import shutil
    import os

    try:
        with resources.path("pygor3", "demo") as path:
            out_path = os.getcwd() + "/demo"
            print("Copy data from : ", path)
            print("to: ", out_path)
            shutil.copytree(path, out_path)
    except Exception as e:
        print("ERROR: ", e)


    """
    # V_gene_list
    vregion, vdregion, dregion, djregion, jregion = igortask.mdl.construct_sequence_VDJ_from_realization_dict(realization_dict)
    with open("teste.fasta", "w") as ofile:
        print("scenario: ", scenario.to_dict())
        str_fasta_offset = ""
        ofile.write("> " + str(indexed_seq.seq_index) + " : \n")
        str_gap = (0 - v_best_aln.offset)*'-'
        ofile.write(str_gap+indexed_seq.sequence + "\n")

        ofile.write("> vregion : \n")
        ofile.write(str_fasta_offset + vregion + "\n")
        str_fasta_offset = str_fasta_offset + "-"*len(vregion)

        ofile.write("> vdregion : \n")
        ofile.write(str_fasta_offset + vdregion + "\n")
        str_fasta_offset = str_fasta_offset + "-" * len(vdregion)

        ofile.write("> dregion : \n")
        ofile.write(str_fasta_offset + dregion + "\n")
        str_fasta_offset = str_fasta_offset + "-" * len(dregion)

        ofile.write("> djregion : \n")
        ofile.write(str_fasta_offset + djregion + "\n")
        str_fasta_offset = str_fasta_offset + "-" * len(djregion)

        ofile.write("> jregion : \n")
        ofile.write(str_fasta_offset + jregion + "\n")
        str_fasta_offset = str_fasta_offset + "-" * len(jregion)
    
    """


    # igortask.mdl.get_AIRR_rearragement_dict_from_scenario()
    # gene_seq = scenario['v_choice'].value
    # gene_3_del = scenario['v_3_del'].value
    #
    # aaa = igortask.mdl.generate_sequence_construction_list()
    # print(aaa)

    # print([event.nickname for event in aaa])


    # igortask.mdl.get_realizations_dict_from_scenario_dict(scenario_realization_dict)
    # print(event_realization)

    #
    # igortask.igor_fln_output_scenarios = igor_fln_output_scenarios
    # mdl = igortask.mdl
    # sep=';'
    # with open(igor_fln_output_scenarios, 'r') as ofile:
    #     header_line = ofile.readline().replace('\n', '')
    #     header_list = header_line.split(sep)
    #
    #     scenario_ordered_list = list()
    #     for ii, colname in enumerate(header_list):
    #         if colname == 'seq_index':
    #             scenario_ordered_list.append(colname)
    #         elif colname == 'scenario_rank':
    #             scenario_ordered_list.append(colname)
    #         elif colname == 'scenario_proba_cond_seq':
    #             scenario_ordered_list.append(colname)
    #         elif colname == 'Errors':
    #             scenario_ordered_list.append('Mismatches')
    #         else:
    #             scenario_ordered_list.append(mdl.parms.dictNameNickname[colname])
    #
    #
    #
    #     str_line = ofile.readline()
    #     print(str_line)
    #     p3.IgorScenario.load_FromLineBestScenario()



    # import pygor3 as p3
    #
    # if seq_index is None:
    #     def generate_csv_line(indexed_sequence: p3.IgorIndexedSequence, indexed_cdr3_record: list,
    #                           list_vdj_alignments: dict, sep=';',
    #                           header_list=['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score', 'd_score',
    #                                        'j_score']):
    #         csv_line = ""
    #         indexed_sequence.sequence = indexed_sequence.sequence.lower()
    #         fields_dict = dict()
    #         fields_dict['sequence_id'] = str(indexed_sequence.seq_index)
    #         fields_dict['sequence'] = indexed_sequence.sequence
    #         fields_dict['v_call'] = list_vdj_alignments['V'].strGene_name
    #         fields_dict['d_call'] = list_vdj_alignments['D'].strGene_name
    #         fields_dict['j_call'] = list_vdj_alignments['J'].strGene_name
    #         fields_dict['v_score'] = str(list_vdj_alignments['V'].score)
    #         fields_dict['d_score'] = str(list_vdj_alignments['D'].score)
    #         fields_dict['j_score'] = str(list_vdj_alignments['J'].score)
    #         # FIXME: CREATE A BETTER WAY TO DO THIS
    #         fields_dict['v_anchor'] = ""
    #         fields_dict['j_anchor'] = ""
    #         fields_dict['junction'] = ""
    #         fields_dict['junction_aa'] = ""
    #         try:
    #             fields_dict['v_anchor'] = str(indexed_cdr3_record[1])  # ""
    #             fields_dict['j_anchor'] = str(indexed_cdr3_record[2])  # ""
    #             fields_dict['junction'] = str(indexed_cdr3_record[3])  # ""
    #             fields_dict['junction_aa'] = str(indexed_cdr3_record[4])  # ""
    #         except Exception as e:
    #             print("No junction for sequence : " + fields_dict['sequence_id'])
    #             print(e)
    #             pass
    #
    #         # align = p3.IgorAlignment_data()
    #
    #         for field in header_list:
    #             csv_line = csv_line + fields_dict[field] + sep
    #         return csv_line
    #
    #     #### STARTS HERE
    #     db = p3.IgorSqliteDB()
    #     # db.flnIgorDB = task.igor_wd+"/"+task.igor_batchname+".db"
    #     db.flnIgorDB = igortask.igor_db.flnIgorDB #args.database
    #     db.connect_db()
    #
    #     # Make a loop over all sequences
    #     if output_fln is None:
    #         fln_output = db.flnIgorDB.split(".db")[0] + "_na.csv"
    #     else:
    #         fln_output = output_fln
    #
    #     ofile = open(fln_output, 'w')
    #     seq_index_list = db.execute_select_query("SELECT seq_index FROM IgorIndexedSeq;")
    #     seq_index_list = map(lambda x: x[0], seq_index_list)
    #     header_list = ['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score', 'd_score', 'j_score', 'v_anchor',
    #                    'j_anchor', 'junction', 'junction_aa']
    #     sep = ";"
    #     str_header = sep.join(header_list)
    #     ofile.write(str_header + '\n')
    #     for seq_index in seq_index_list:
    #         try:
    #             # seq_index = args.seq_index
    #             indexed_sequence = db.get_IgorIndexedSeq_By_seq_index(seq_index)
    #             indexed_sequence.offset = 0
    #
    #             indexed_cdr3_record = db.fetch_IgorIndexedCDR3_By_seq_index(seq_index)
    #             print(indexed_cdr3_record)
    #
    #             best_v_align_data = db.get_best_IgorAlignment_data_By_seq_index('V', indexed_sequence.seq_index)
    #             best_d_align_data = db.get_best_IgorAlignment_data_By_seq_index('D', indexed_sequence.seq_index)
    #             best_j_align_data = db.get_best_IgorAlignment_data_By_seq_index('J', indexed_sequence.seq_index)
    #
    #             vdj_naive_alignment = {'V': best_v_align_data,
    #                                    'D': best_d_align_data,
    #                                    'J': best_j_align_data}
    #
    #             v_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
    #             # print('V', len(v_align_data_list), [ ii.score for ii in v_align_data_list])
    #             try:
    #                 d_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('D', indexed_sequence.seq_index)
    #             except Exception as e:
    #                 print(e)
    #                 print("No D sequences alignments found!")
    #                 pass
    #             # print('D', len(d_align_data_list), [ ii.score for ii in d_align_data_list])
    #             j_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('J', indexed_sequence.seq_index)
    #             # print('J', len(j_align_data_list), [ ii.score for ii in j_align_data_list])
    #
    #             # 1. Choose the highest score then check if this one is the desire range.
    #             # if there is an overlap
    #             # calculate score without overlap. If overlap
    #             # if hightest score
    #             for i, d_align_data in enumerate(d_align_data_list):
    #                 # Check if D is btwn V and J position
    #                 if (best_v_align_data.offset_3_p <= d_align_data.offset_5_p) and (
    #                         d_align_data.offset_3_p <= best_j_align_data.offset_5_p):
    #                     # vdj_naive_alignment['D'+str(i)] = d_align_data
    #                     vdj_naive_alignment['D'] = d_align_data
    #                     break
    #
    #             # ofile = open(task.igor_batchname+'__'+str(indexed_sequence.seq_index)+'_na.fasta', 'w')
    #             # str_fasta = generate_str_fasta(indexed_sequence, vdj_naive_alignment)
    #             # ofile.write(str_fasta)
    #
    #             str_csv_line = generate_csv_line(indexed_sequence, indexed_cdr3_record, vdj_naive_alignment, sep=sep,
    #                                              header_list=header_list)
    #             ofile.write(str_csv_line + "\n")
    #         except Exception as e:
    #             print("ERROR : with sequence id : " + str(seq_index))
    #             print(e)
    #             pass
    #     ofile.close()
    # else:
    #     ###############################################
    #     def generate_str_fasta(indexed_sequence, list_vdj_alignments: dict):
    #         """ Given an Sequence index and the corresponding alignments vj/ vdj
    #         return a string with considering only offset"""
    #
    #         indexed_sequence.sequence = indexed_sequence.sequence.lower()
    #         # add mismatches in sequence.
    #         s = list(indexed_sequence.sequence)
    #         for key_align in list_vdj_alignments.keys():
    #             for pos_mis in list_vdj_alignments[key_align].mismatches:
    #                 s[pos_mis] = s[pos_mis].upper()
    #         indexed_sequence.sequence = "".join(s)
    #
    #         str_fasta = ""
    #         min_offset_key = min(list_vdj_alignments.keys(), key=lambda x: list_vdj_alignments[x].offset)  # .offset
    #         min_offset = list_vdj_alignments[min_offset_key].offset
    #         min_offset = min(indexed_sequence.offset, min_offset)
    #
    #         delta_offset = indexed_sequence.offset - min_offset
    #         str_prefix = '-' * (delta_offset)
    #         str_fasta_sequence = str_prefix + indexed_sequence.sequence
    #         print(str_fasta_sequence)
    #         str_fasta = str_fasta + "> " + str(indexed_sequence.seq_index) + "\n"
    #         str_fasta = str_fasta + str_fasta_sequence + "\n"
    #         for key in list_vdj_alignments.keys():
    #             list_vdj_alignments[key].strGene_seq = list_vdj_alignments[key].strGene_seq.lower()
    #             delta_offset = list_vdj_alignments[key].offset - min_offset
    #             str_prefix = '-' * (delta_offset)
    #             str_fasta_sequence = str_prefix + list_vdj_alignments[key].strGene_seq
    #             print(str_fasta_sequence)
    #             str_fasta = str_fasta + "> " + key + ", " + list_vdj_alignments[key].strGene_name + "\n"
    #             str_fasta = str_fasta + str_fasta_sequence + "\n"
    #             offset_5_p = list_vdj_alignments[key].offset_5_p - min_offset
    #             offset_3_p = list_vdj_alignments[key].offset_3_p - min_offset
    #             print("delta_offset : ", delta_offset)
    #             print("offset_5_p : ", list_vdj_alignments[key].offset_5_p, offset_5_p)
    #             print("offset_3_p : ", list_vdj_alignments[key].offset_3_p, offset_3_p)
    #             str_prefix_2 = '-' * (offset_5_p + 1)
    #             str_fasta_sequence2 = str_prefix_2 + str_fasta_sequence[offset_5_p + 1:offset_3_p + 1]
    #
    #             str_fasta = str_fasta + "> " + list_vdj_alignments[key].strGene_name + ", score : " + str(
    #                 list_vdj_alignments[key].score) + "\n"
    #             str_fasta = str_fasta + str_fasta_sequence2 + "\n"
    #
    #             # TODO ADD MISMATCHES
    #             align = list_vdj_alignments[key]
    #             # align mismatches are in indexed sequence reference I need to convert it to gene reference given the alignment
    #             # given the align.offset
    #             # pos_in_gene  = pos_in_seq - align.offset
    #             # pos_in_gene = cdr3 - align.offset
    #
    #         return str_fasta
    #
    #     def generate_csv_line(indexed_sequence: p3.IgorIndexedSequence, list_vdj_alignments: dict, sep=';',
    #                           header_list=['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score',
    #                                        'd_score', 'j_score', 'junction']):
    #         csv_line = ""
    #         indexed_sequence.sequence = indexed_sequence.sequence.lower()
    #         fields_dict = dict()
    #         fields_dict['sequence_id'] = indexed_sequence.seq_index
    #         fields_dict['sequence'] = indexed_sequence.sequence
    #         fields_dict['v_call'] = list_vdj_alignments['V'].strGene_name
    #         fields_dict['d_call'] = list_vdj_alignments['D'].strGene_name
    #         fields_dict['j_call'] = list_vdj_alignments['J'].strGene_name
    #         fields_dict['v_score'] = list_vdj_alignments['V'].score
    #         fields_dict['d_score'] = list_vdj_alignments['D'].score
    #         fields_dict['j_score'] = list_vdj_alignments['J'].score
    #
    #         align = p3.IgorAlignment_data()
    #
    #         for field in header_list:
    #             csv_line = csv_line + fields_dict[field] + sep
    #         return csv_line
    #
    #         # mixcr :
    #         # targetSequences	targetQualities	allVHitsWithScore
    #         # allDHitsWithScore	allJHitsWithScore	allCHitsWithScore	allVAlignments	allDAlignments
    #         # allJAlignments	allCAlignments	nSeqFR1	minQualFR1	nSeqCDR1	minQualCDR1	nSeqFR2	minQualFR2	nSeqCDR2	minQualCDR2	nSeqFR3	minQualFR3	nSeqCDR3	minQualCDR3	nSeqFR4	minQualFR4	aaSeqFR1	aaSeqCDR1	aaSeqFR2	aaSeqCDR2	aaSeqFR3	aaSeqCDR3	aaSeqFR4	refPoints
    #
    #         # AIRR :
    #         # sequence_id (req)
    #         # sequence (req)
    #         # rev_comp (req) //false by default # FIXME: for now
    #         # productive (req) // false by default # FIXME: a naive alignmnet of productive is possible.
    #         # locus (chain type) // eg. TCR IG
    #         # v_call
    #         # d_call
    #         # d2_call // FIXME: not gonna use this by the moment.
    #         # j_call
    #         # sequence alignment // FIXME: Aligned portion of query sequence. IMGT-gaps check the best way to do this
    #         # germline_alignment // FIXME: Assembled, aligned, full-length inferred germline sequence spanning the same region as the sequence_alignment field (typically the V(D)J region) and including the same set of corrections and spacers (if any)
    #         # junction
    #         # junction_aa // Our CDR3
    #         # np1 (nts btwn V and D or btwn V and J  // vd_ins  or vj_ins
    #         # np2 (nts btwn first D and J gene or btwn first D and second D // dj_ins or d1d2_ins
    #         # np3 (nts btwn second D and J gene
    #         # v_score
    #         # v_cigar
    #         # d_score
    #         # d_cigar
    #         # j_score
    #         # j_cigar
    #         # junction_length
    #
    #     def generate_str_fasta_simple(indexed_sequence, list_vdj_alignments: dict):
    #         """ Given an Sequence index and the corresponding alignments vj/ vdj
    #         return a string with considering only offset"""
    #
    #         str_fasta = ""
    #         str_fasta = str_fasta + "> " + str(indexed_sequence.seq_index) + "\n"
    #         str_fasta = str_fasta + indexed_sequence.sequence + "\n"
    #         for key in list_vdj_alignments.keys():
    #             str_fasta = str_fasta + "> " + key + ", " + list_vdj_alignments[key].strGene_name + "\n"
    #             str_fasta = str_fasta + list_vdj_alignments[key].strGene_seq + "\n"
    #
    #         # print(list_vdj_alignments['V'])
    #         # print(list_vdj_alignments['D        '])
    #         return str_fasta
    #
    #     ############################################
    #     #### naive align with id seq_index
    #     db = p3.IgorSqliteDB()
    #
    #     db.flnIgorDB = igortask.igor_fln_db
    #     db.connect_db()
    #
    #     #seq_index = args.seq_index
    #     try:
    #         indexed_sequence = db.get_IgorIndexedSeq_By_seq_index(seq_index)
    #     except Exception as e:
    #         print("ERROR: seq_index ", seq_index, " not found in db ", igortask.igor_fln_db)
    #         raise
    #     indexed_sequence.offset = 0
    #
    #     best_v_align_data = db.get_best_IgorAlignment_data_By_seq_index('V', indexed_sequence.seq_index)
    #     best_j_align_data = db.get_best_IgorAlignment_data_By_seq_index('J', indexed_sequence.seq_index)
    #
    #     try:
    #         best_d_align_data = db.get_best_IgorAlignment_data_By_seq_index('D', indexed_sequence.seq_index)
    #         vdj_naive_alignment = {'V': best_v_align_data,
    #                                'D': best_d_align_data,
    #                                'J': best_j_align_data}
    #         v_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
    #         print('V', len(v_align_data_list), [ii.score for ii in v_align_data_list])
    #         d_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('D', indexed_sequence.seq_index)
    #         print('D', len(d_align_data_list), [ii.score for ii in d_align_data_list])
    #         j_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('J', indexed_sequence.seq_index)
    #         print('J', len(j_align_data_list), [ii.score for ii in j_align_data_list])
    #         # 1. Choose the highest score then check if this one is the desire range.
    #         # if there is an overlap
    #         # calculate score without overlap. If overlap
    #         # if hightest score
    #         for i, d_align_data in enumerate(d_align_data_list):
    #             # Check if D is btwn V and J position
    #             if (best_v_align_data.offset_3_p <= d_align_data.offset_5_p) and (
    #                     d_align_data.offset_3_p <= best_j_align_data.offset_5_p):
    #                 # vdj_naive_alignment['D'+str(i)] = d_align_data
    #                 vdj_naive_alignment['D'] = d_align_data
    #                 break
    #
    #     except Exception as e:
    #         print(e)
    #         print("No d gene alignments found!")
    #         vdj_naive_alignment = {'V': best_v_align_data,
    #                                'J': best_j_align_data}
    #         v_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
    #         print('V', len(v_align_data_list), [ii.score for ii in v_align_data_list])
    #         j_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('J', indexed_sequence.seq_index)
    #         print('J', len(j_align_data_list), [ii.score for ii in j_align_data_list])
    #         pass
    #
    #     if output_fln is None:
    #         batchname = db.flnIgorDB.split(".db")[0]
    #         fln_output = batchname + '__' + str(indexed_sequence.seq_index) + '_na.fasta'
    #         # fln_output = args.database.split(".db")[0]+"_na.csv"
    #     else:
    #         fln_output = output_fln
    #
    #     ofile = open(fln_output, 'w')
    #     str_fasta = generate_str_fasta(indexed_sequence, vdj_naive_alignment)
    #     ofile.write(str_fasta)
    #     ofile.close()
    #
    #
    #

cli.add_command(demo_get_data)
# cli.add_command(test)


if __name__ == '__main__':
    cli()


#
# @click.command("test") #invoke_without_command=True)
# ############## COMMON options ##############
# ############## COMMON options ##############
# @click.option("-D", "--set_database", "igor_fln_db", default=None, help="Igor database created with database script.")
#
# ############## NO COMMON options ##############
# @click.option("-o", "--output-filename", "output_fln", default=None, help="Filename of output file csv or fasta if --seq_index is use.")
# @click.option("-s", "--seq-index", "seq_index", default=None, type=int, help="Sequence id (seq_index) in db.", required=False)
# def test(output_fln, seq_index, igor_fln_db):
#     """
#     Get a naive alignment (no scenarios) from igor alignments
#     """
#     from pygor3 import IgorTask
#     igortask = IgorTask()
#     igortask.igor_fln_db = igor_fln_db
#
#     import pygor3 as p3
#
#     if seq_index is None:
#         def generate_csv_line(indexed_sequence: p3.IgorIndexedSequence, indexed_cdr3_record: list,
#                               list_vdj_alignments: dict, sep=';',
#                               header_list=['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score', 'd_score',
#                                            'j_score']):
#             csv_line = ""
#             indexed_sequence.sequence = indexed_sequence.sequence.lower()
#             fields_dict = dict()
#             fields_dict['sequence_id'] = str(indexed_sequence.seq_index)
#             fields_dict['sequence'] = indexed_sequence.sequence
#             fields_dict['v_call'] = list_vdj_alignments['V'].strGene_name
#             fields_dict['d_call'] = list_vdj_alignments['D'].strGene_name
#             fields_dict['j_call'] = list_vdj_alignments['J'].strGene_name
#             fields_dict['v_score'] = str(list_vdj_alignments['V'].score)
#             fields_dict['d_score'] = str(list_vdj_alignments['D'].score)
#             fields_dict['j_score'] = str(list_vdj_alignments['J'].score)
#             # FIXME: CREATE A BETTER WAY TO DO THIS
#             fields_dict['v_anchor'] = ""
#             fields_dict['j_anchor'] = ""
#             fields_dict['junction'] = ""
#             fields_dict['junction_aa'] = ""
#             try:
#                 fields_dict['v_anchor'] = str(indexed_cdr3_record[1])  # ""
#                 fields_dict['j_anchor'] = str(indexed_cdr3_record[2])  # ""
#                 fields_dict['junction'] = str(indexed_cdr3_record[3])  # ""
#                 fields_dict['junction_aa'] = str(indexed_cdr3_record[4])  # ""
#             except Exception as e:
#                 print("No junction for sequence : " + fields_dict['sequence_id'])
#                 print(e)
#                 pass
#
#             # align = p3.IgorAlignment_data()
#
#             for field in header_list:
#                 csv_line = csv_line + fields_dict[field] + sep
#             return csv_line
#
#         #### STARTS HERE
#         db = p3.IgorSqliteDB()
#         # db.flnIgorDB = task.igor_wd+"/"+task.igor_batchname+".db"
#         db.flnIgorDB = igortask.igor_db.flnIgorDB #args.database
#         db.connect_db()
#
#         # Make a loop over all sequences
#         if output_fln is None:
#             fln_output = db.flnIgorDB.split(".db")[0] + "_na.csv"
#         else:
#             fln_output = output_fln
#
#         ofile = open(fln_output, 'w')
#         seq_index_list = db.execute_select_query("SELECT seq_index FROM IgorIndexedSeq;")
#         seq_index_list = map(lambda x: x[0], seq_index_list)
#         header_list = ['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score', 'd_score', 'j_score', 'v_anchor',
#                        'j_anchor', 'junction', 'junction_aa']
#         sep = ";"
#         str_header = sep.join(header_list)
#         ofile.write(str_header + '\n')
#         for seq_index in seq_index_list:
#             try:
#                 # seq_index = args.seq_index
#                 indexed_sequence = db.get_IgorIndexedSeq_By_seq_index(seq_index)
#                 indexed_sequence.offset = 0
#
#                 indexed_cdr3_record = db.fetch_IgorIndexedCDR3_By_seq_index(seq_index)
#                 print(indexed_cdr3_record)
#
#                 best_v_align_data = db.get_best_IgorAlignment_data_By_seq_index('V', indexed_sequence.seq_index)
#                 best_d_align_data = db.get_best_IgorAlignment_data_By_seq_index('D', indexed_sequence.seq_index)
#                 best_j_align_data = db.get_best_IgorAlignment_data_By_seq_index('J', indexed_sequence.seq_index)
#
#                 vdj_naive_alignment = {'V': best_v_align_data,
#                                        'D': best_d_align_data,
#                                        'J': best_j_align_data}
#
#                 v_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
#                 # print('V', len(v_align_data_list), [ ii.score for ii in v_align_data_list])
#                 try:
#                     d_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('D', indexed_sequence.seq_index)
#                 except Exception as e:
#                     print(e)
#                     print("No D sequences alignments found!")
#                     pass
#                 # print('D', len(d_align_data_list), [ ii.score for ii in d_align_data_list])
#                 j_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('J', indexed_sequence.seq_index)
#                 # print('J', len(j_align_data_list), [ ii.score for ii in j_align_data_list])
#
#                 # 1. Choose the highest score then check if this one is the desire range.
#                 # if there is an overlap
#                 # calculate score without overlap. If overlap
#                 # if hightest score
#                 for i, d_align_data in enumerate(d_align_data_list):
#                     # Check if D is btwn V and J position
#                     if (best_v_align_data.offset_3_p <= d_align_data.offset_5_p) and (
#                             d_align_data.offset_3_p <= best_j_align_data.offset_5_p):
#                         # vdj_naive_alignment['D'+str(i)] = d_align_data
#                         vdj_naive_alignment['D'] = d_align_data
#                         break
#
#                 # ofile = open(task.igor_batchname+'__'+str(indexed_sequence.seq_index)+'_na.fasta', 'w')
#                 # str_fasta = generate_str_fasta(indexed_sequence, vdj_naive_alignment)
#                 # ofile.write(str_fasta)
#
#                 str_csv_line = generate_csv_line(indexed_sequence, indexed_cdr3_record, vdj_naive_alignment, sep=sep,
#                                                  header_list=header_list)
#                 ofile.write(str_csv_line + "\n")
#             except Exception as e:
#                 print("ERROR : with sequence id : " + str(seq_index))
#                 print(e)
#                 pass
#         ofile.close()
#     else:
#         ###############################################
#         def generate_str_fasta(indexed_sequence, list_vdj_alignments: dict):
#             """ Given an Sequence index and the corresponding alignments vj/ vdj
#             return a string with considering only offset"""
#
#             indexed_sequence.sequence = indexed_sequence.sequence.lower()
#             # add mismatches in sequence.
#             s = list(indexed_sequence.sequence)
#             for key_align in list_vdj_alignments.keys():
#                 for pos_mis in list_vdj_alignments[key_align].mismatches:
#                     s[pos_mis] = s[pos_mis].upper()
#             indexed_sequence.sequence = "".join(s)
#
#             str_fasta = ""
#             min_offset_key = min(list_vdj_alignments.keys(), key=lambda x: list_vdj_alignments[x].offset)  # .offset
#             min_offset = list_vdj_alignments[min_offset_key].offset
#             min_offset = min(indexed_sequence.offset, min_offset)
#
#             delta_offset = indexed_sequence.offset - min_offset
#             str_prefix = '-' * (delta_offset)
#             str_fasta_sequence = str_prefix + indexed_sequence.sequence
#             print(str_fasta_sequence)
#             str_fasta = str_fasta + "> " + str(indexed_sequence.seq_index) + "\n"
#             str_fasta = str_fasta + str_fasta_sequence + "\n"
#             for key in list_vdj_alignments.keys():
#                 list_vdj_alignments[key].strGene_seq = list_vdj_alignments[key].strGene_seq.lower()
#                 delta_offset = list_vdj_alignments[key].offset - min_offset
#                 str_prefix = '-' * (delta_offset)
#                 str_fasta_sequence = str_prefix + list_vdj_alignments[key].strGene_seq
#                 print(str_fasta_sequence)
#                 str_fasta = str_fasta + "> " + key + ", " + list_vdj_alignments[key].strGene_name + "\n"
#                 str_fasta = str_fasta + str_fasta_sequence + "\n"
#                 offset_5_p = list_vdj_alignments[key].offset_5_p - min_offset
#                 offset_3_p = list_vdj_alignments[key].offset_3_p - min_offset
#                 print("delta_offset : ", delta_offset)
#                 print("offset_5_p : ", list_vdj_alignments[key].offset_5_p, offset_5_p)
#                 print("offset_3_p : ", list_vdj_alignments[key].offset_3_p, offset_3_p)
#                 str_prefix_2 = '-' * (offset_5_p + 1)
#                 str_fasta_sequence2 = str_prefix_2 + str_fasta_sequence[offset_5_p + 1:offset_3_p + 1]
#
#                 str_fasta = str_fasta + "> " + list_vdj_alignments[key].strGene_name + ", score : " + str(
#                     list_vdj_alignments[key].score) + "\n"
#                 str_fasta = str_fasta + str_fasta_sequence2 + "\n"
#
#                 # TODO ADD MISMATCHES
#                 align = list_vdj_alignments[key]
#                 # align mismatches are in indexed sequence reference I need to convert it to gene reference given the alignment
#                 # given the align.offset
#                 # pos_in_gene  = pos_in_seq - align.offset
#                 # pos_in_gene = cdr3 - align.offset
#
#             return str_fasta
#
#         def generate_csv_line(indexed_sequence: p3.IgorIndexedSequence, list_vdj_alignments: dict, sep=';',
#                               header_list=['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score',
#                                            'd_score', 'j_score', 'junction']):
#             csv_line = ""
#             indexed_sequence.sequence = indexed_sequence.sequence.lower()
#             fields_dict = dict()
#             fields_dict['sequence_id'] = indexed_sequence.seq_index
#             fields_dict['sequence'] = indexed_sequence.sequence
#             fields_dict['v_call'] = list_vdj_alignments['V'].strGene_name
#             fields_dict['d_call'] = list_vdj_alignments['D'].strGene_name
#             fields_dict['j_call'] = list_vdj_alignments['J'].strGene_name
#             fields_dict['v_score'] = list_vdj_alignments['V'].score
#             fields_dict['d_score'] = list_vdj_alignments['D'].score
#             fields_dict['j_score'] = list_vdj_alignments['J'].score
#
#             align = p3.IgorAlignment_data()
#
#             for field in header_list:
#                 csv_line = csv_line + fields_dict[field] + sep
#             return csv_line
#
#             # mixcr :
#             # targetSequences	targetQualities	allVHitsWithScore
#             # allDHitsWithScore	allJHitsWithScore	allCHitsWithScore	allVAlignments	allDAlignments
#             # allJAlignments	allCAlignments	nSeqFR1	minQualFR1	nSeqCDR1	minQualCDR1	nSeqFR2	minQualFR2	nSeqCDR2	minQualCDR2	nSeqFR3	minQualFR3	nSeqCDR3	minQualCDR3	nSeqFR4	minQualFR4	aaSeqFR1	aaSeqCDR1	aaSeqFR2	aaSeqCDR2	aaSeqFR3	aaSeqCDR3	aaSeqFR4	refPoints
#
#             # AIRR :
#             # sequence_id (req)
#             # sequence (req)
#             # rev_comp (req) //false by default # FIXME: for now
#             # productive (req) // false by default # FIXME: a naive alignmnet of productive is possible.
#             # locus (chain type) // eg. TCR IG
#             # v_call
#             # d_call
#             # d2_call // FIXME: not gonna use this by the moment.
#             # j_call
#             # sequence alignment // FIXME: Aligned portion of query sequence. IMGT-gaps check the best way to do this
#             # germline_alignment // FIXME: Assembled, aligned, full-length inferred germline sequence spanning the same region as the sequence_alignment field (typically the V(D)J region) and including the same set of corrections and spacers (if any)
#             # junction
#             # junction_aa // Our CDR3
#             # np1 (nts btwn V and D or btwn V and J  // vd_ins  or vj_ins
#             # np2 (nts btwn first D and J gene or btwn first D and second D // dj_ins or d1d2_ins
#             # np3 (nts btwn second D and J gene
#             # v_score
#             # v_cigar
#             # d_score
#             # d_cigar
#             # j_score
#             # j_cigar
#             # junction_length
#
#         def generate_str_fasta_simple(indexed_sequence, list_vdj_alignments: dict):
#             """ Given an Sequence index and the corresponding alignments vj/ vdj
#             return a string with considering only offset"""
#
#             str_fasta = ""
#             str_fasta = str_fasta + "> " + str(indexed_sequence.seq_index) + "\n"
#             str_fasta = str_fasta + indexed_sequence.sequence + "\n"
#             for key in list_vdj_alignments.keys():
#                 str_fasta = str_fasta + "> " + key + ", " + list_vdj_alignments[key].strGene_name + "\n"
#                 str_fasta = str_fasta + list_vdj_alignments[key].strGene_seq + "\n"
#
#             # print(list_vdj_alignments['V'])
#             # print(list_vdj_alignments['D        '])
#             return str_fasta
#
#         ############################################
#         #### naive align with id seq_index
#         db = p3.IgorSqliteDB()
#
#         db.flnIgorDB = igortask.igor_fln_db
#         db.connect_db()
#
#         #seq_index = args.seq_index
#         try:
#             indexed_sequence = db.get_IgorIndexedSeq_By_seq_index(seq_index)
#         except Exception as e:
#             print("ERROR: seq_index ", seq_index, " not found in db ", igortask.igor_fln_db)
#             raise
#         indexed_sequence.offset = 0
#
#         best_v_align_data = db.get_best_IgorAlignment_data_By_seq_index('V', indexed_sequence.seq_index)
#         best_j_align_data = db.get_best_IgorAlignment_data_By_seq_index('J', indexed_sequence.seq_index)
#
#         try:
#             best_d_align_data = db.get_best_IgorAlignment_data_By_seq_index('D', indexed_sequence.seq_index)
#             vdj_naive_alignment = {'V': best_v_align_data,
#                                    'D': best_d_align_data,
#                                    'J': best_j_align_data}
#             v_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
#             print('V', len(v_align_data_list), [ii.score for ii in v_align_data_list])
#             d_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('D', indexed_sequence.seq_index)
#             print('D', len(d_align_data_list), [ii.score for ii in d_align_data_list])
#             j_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('J', indexed_sequence.seq_index)
#             print('J', len(j_align_data_list), [ii.score for ii in j_align_data_list])
#             # 1. Choose the highest score then check if this one is the desire range.
#             # if there is an overlap
#             # calculate score without overlap. If overlap
#             # if hightest score
#             for i, d_align_data in enumerate(d_align_data_list):
#                 # Check if D is btwn V and J position
#                 if (best_v_align_data.offset_3_p <= d_align_data.offset_5_p) and (
#                         d_align_data.offset_3_p <= best_j_align_data.offset_5_p):
#                     # vdj_naive_alignment['D'+str(i)] = d_align_data
#                     vdj_naive_alignment['D'] = d_align_data
#                     break
#
#         except Exception as e:
#             print(e)
#             print("No d gene alignments found!")
#             vdj_naive_alignment = {'V': best_v_align_data,
#                                    'J': best_j_align_data}
#             v_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
#             print('V', len(v_align_data_list), [ii.score for ii in v_align_data_list])
#             j_align_data_list = db.get_IgorAlignment_data_list_By_seq_index('J', indexed_sequence.seq_index)
#             print('J', len(j_align_data_list), [ii.score for ii in j_align_data_list])
#             pass
#
#         if output_fln is None:
#             batchname = db.flnIgorDB.split(".db")[0]
#             fln_output = batchname + '__' + str(indexed_sequence.seq_index) + '_na.fasta'
#             # fln_output = args.database.split(".db")[0]+"_na.csv"
#         else:
#             fln_output = output_fln
#
#         ofile = open(fln_output, 'w')
#         str_fasta = generate_str_fasta(indexed_sequence, vdj_naive_alignment)
#         ofile.write(str_fasta)
#         ofile.close()


# pygor3-cli igor-infer -i sequences.txt -o model
# pygor3-cli igor-generate -human -TRB -o seqs -n 1000
# pygor3-cli olga-compute -ppost --VJmodel model/ -i
# pygor3-cli igor-model -i model/ -o csv_files
# pygor3-cli igor-model -i model/ -o png_files
