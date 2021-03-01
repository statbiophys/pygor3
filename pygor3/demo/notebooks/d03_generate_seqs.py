#!/usr/bin/env python3
import pygor3 as p3
import numpy as np

# all sequences
# run preporo
# task.make_model_default_VDJ_from_genomes_dir()
# task.make_model_default_VDJ_from_dataframes()
# task.make_model_default_VDJ_from_fasta(fln_genomicVs:Union[None,str,Path]=None,
#                                         fln_genomicDs, fln_genomicJs)
# task = p3.IgorTask(igor_read_seqs='')
# df = task.run_generate(N_seqs=10, igor_species="human", igor_chain="tcr_beta")
# task.run_infer('input_file')

df_genes_dict = p3.imgt.download_ref_genome_VDJ('Mus+musculus', 'TRB')
df_genes_dict['J'].columns
records = p3.imgt.get_gene_template('Mus+musculus', 'TRBV')
for record in records:
    print(record)

print(df_genes_dict['V'])

genomes0 = p3.IgorRefGenome.load_VDJ_from_IMGT_website('Mus+musculus', 'TRB')

"""
# Edit genes in any form you want

# Only use genes with
ddf_V = df_genes_dict['V'][df_genes_dict['V']['gfunction'].notna()]
ddf_V = ddf_V.reset_index().drop(columns='index')
ddf_V
ddf_V['gfunction']

ddf_V2 = ddf_V[ddf_V['gfunction'].str.contains('F')]
ddf_V2 = ddf_V2.reset_index().drop(columns='index')
ddf_V2

ddf_V3 = ddf_V[ ~ddf_V['gfunction'].str.contains('F')]
ddf_V3 = ddf_V3.reset_index().drop(columns='index')
ddf_V3

notna_series = df_genes_dict['V']['gfunction'].notna()
not_gfunction = ~(df_genes_dict['V'][notna_series]['gfunction']).str.contains('F')
df_genes_dict['V'][df_genes_dict['V']['gfunction'][not_gfunction]]
df_genes_dict['V'][notna_series][not_gfunction]
notna_series

df_genes_dict['V'].dropna()
# mdl = p3.IgorModel.make_default_VDJ(df_genes_dict['V'], df_genes_dict['D'], df_genes_dict['J'])
mdl = p3.IgorModel.make_default_from_Dataframe_dict(df_genes_dict)


mdl_parms = p3.IgorModel_Parms.make_default_VDJ(df_genes_dict['V'], df_genes_dict['D'], df_genes_dict['J'])
for event in mdl_parms.Event_list:
    # print(event.to_dict())
    if event.event_type == 'GeneChoice':
        if event.seq_type == 'V_gene':
            print("Save V genes in fasta file and anchors")
        elif event.seq_type == 'D_gene':
            print("Save D gene in fasta file")
        elif event.seq_type == 'J_gene':
            print("Save J genes in fasta file and anchors")
        else:
            print(event.seq_type, " Not recognized")

        print(event.event_type, event.seq_type)

event = mdl_parms.get_Event('j_choice')
event.to_dict()
realizations = event.realizations
index = 7
tmp_list = [ realization for realization in realizations if realization.id == index]
tmp_list[0]

ref_genome_path = 'models/Mus+musculus/TRB/ref_genome/'
genomes = p3.IgorRefGenome.load_from_path(ref_genome_path)
anchors = p3.IgorAnchors(path_ref_genome=ref_genome_path)
anchors.flnJanchors
anchors.df_Vanchors
anchors.df_Vanchors.loc['TRBV31*02']

fln_db = '/home/olivares/testing_pygor/demo/new_IgL_naive.db'
# p3.IgorRefGenome.load

anchors.update_default_filenames()

"""
