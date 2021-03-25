#!/usr/bin/python3
from __future__ import print_function
import os
# import IgorModel
flag_verbose = False
import requests
from bs4 import BeautifulSoup
from Bio import SeqIO

from .utils import *

try:
    from StringIO import StringIO # Python 2
except ImportError:
    from io import StringIO # Python 3

imgt_params={
    'url.home' : "http://www.imgt.org",
    'url.genelist' : "http://www.imgt.org/download/GENE-DB/IMGTGENEDB-GeneList",
    'url.genedb' : "http://www.imgt.org/genedb/GENElect?",
    'url.release' : "http://www.imgt.org/download/GENE-DB/RELEASE",
    'url.readme' : "http://www.imgt.org/download/GENE-DB/README.txt",
    'Genesymbol' : [],
    'data.chainlist' : ['TRA', 'TRB', 'BCR']
}

def get_species_list():
    """
    Returns list of all species available in the IMGT database.
    :return : species list to make queries to IMGT GENEDB.
    """
    response = requests.get(imgt_params['url.genelist'])
    tmp_species_list = set()
    for linea in response.text.split('\n')[1:]:
        tmp_species_list.add(linea.split(';')[0])
    species_list = [specie.replace(" ", "+") for specie in tmp_species_list]
    return species_list

def get_genedb_query72(specie: str, gene: str, imgt_genedb=imgt_params['url.genedb']):
    """
    Returns imgt link to download genomic template according to imgt database.
    :param specie: imgt specie name check get_species_list() to get the names of species in IMGT DB.
    :param gene: imgt nomenclature for gene.
    :return : string of the requested link
    """
    strQuery = "query=7.2+"
    url_gene = imgt_genedb + strQuery + gene + "&species=" + specie
    return url_gene

def get_genedb_query81_imgtlabel(specie: str, gene: str, imgtlabel:str, imgt_genedb=imgt_params['url.genedb']):
    """
    Returns imgt link to download genomic template according to imgt database with the corresponding IMGTlabel
    :return : string of the requested link
    """
    strQuery = "query=8.1+"
    url_gene = imgt_genedb + strQuery + gene + "&species=" + specie+"&IMGTlabel="+imgtlabel
    return url_gene

def extract_fasta_from_url(url):
    """ Extract fasta information from the typical IMGT format
    """
    print(url)
    # Get content from imgt url
    response = requests.get(url)

    soup = BeautifulSoup(response.text, 'html.parser')
    result = soup.find_all('pre')[1].text
    return result

#def genSeqRecords(url):
def get_records_list(url):
    """
    :return : a list of Sequence records.
    """
    rawfasta = extract_fasta_from_url(url)
    # FIXME: Is better to preserve the generator but ... you know maybe later XD
    return list(SeqIO.parse(StringIO(rawfasta), "fasta"))

def save_records2fasta(records, filename:str):
    # flnGenomicJs = DIR_REF_GENOME + "genomicJs.fasta"
    with open(filename, "w") as ofile:
        SeqIO.write(records, ofile, "fasta")

    return filename

def load_records_from_fasta(filename:str):
    # flnGenomicJs = DIR_REF_GENOME + "genomicJs.fasta"
    return list(SeqIO.parse(filename, "fasta"))


def get_gene_template(specie: str, gene: str, modelspath=None, filename=None, imgt_genedb=imgt_params['url.genedb']):
    """
    Download genomic template according to specie and gene from imgt database

    :param species: Species name in IMGT nomenclature.
    :type species: str
    :param gene: IMGT gene nomenclature like TRAV, IGHJ, etc.
    :type species: str
    :return : string of the requested link
    """
    try:
        url_query = get_genedb_query72(specie, gene, imgt_genedb=imgt_genedb)
        print(specie, gene, url_query)

        # 2. load data from IMGT url
        records = get_records_list(url_query)
        records = map(lambda x: x.upper(), records)
        return records
    except Exception as e:
        raise e

def download_gene_template(specie: str, gene: str, ref_genes_path=None, modelspath=None, filename=None, imgt_genedb=imgt_params['url.genedb']):
    """
    Create a file in IGoR default format
    and returns output filename.
    :param specie: IMGT specie nomenclature use "+" instead of a space " " Mus+musculus
    :type species: str
    :param gene: IMGT gene nomenclature like TRAV, IGHJ, etc
    :type species: str
    :param modelspath: root paths for all models.
    :type modelspath: str, optional
    :param filename: gene template fasta filename.
    :type filename: str, optional
    :param imgt_genedb: Url of IMGT GeneDB web application.
    :type imgt_genedb: str
    """
    # 3. Write it in a fasta file.
    records = get_gene_template(specie, gene, imgt_genedb=imgt_genedb)

    if gene[-1] in ['V', 'D', 'J']:
        if ref_genes_path is None:
            ref_genes_path = ""

            if modelspath is None:
                modelspath = "models"

            gene_dir = gene[:-1]
            make_igor_directories(gene_dir, specie, modelspath=modelspath)
            ref_genes_path = modelspath + "/" + specie + "/" + gene_dir + "/ref_genome/"
        else:
            # create directory
            ref_genes_path = ref_genes_path + "/"
            os.system("mkdir -p " + ref_genes_path)

        if filename is None :
            filename = ref_genes_path + "genomic" + gene[-1] + "s__imgt.fasta"
            #filename = modelspath + "genomic__" + specie + "_" + gene + ".fasta"

    #flnGenomicJs = DIR_REF_GENOME + "genomicJs.fasta"
    flnGenome = save_records2fasta(records, filename)
    fln_dict = dict()
    fln_dict["fln"+gene[-1]+"Genome"] = flnGenome
    fln_dict['ref_genes_path'] = ref_genes_path
    fln_dict['modelspath'] = modelspath
    fln_dict['records'] = records
    return fln_dict

def genKey(seqDescription):
    """
    :param seqDescritption : is the header description in fasta file
    :return key : a generated key to create a dictionary
    :return startPos : the start position of the sequence based on ACCesion number sequence.
    """
    fields = seqDescription.split("|")
    key    = "|".join(fields[:4])
    return key

def getStartPos(seqDescription):
    """
    Return the start position of sequence based on the fasta file description
    """
    fields = seqDescription.split("|")
    #startPosLIGMDB = map(int, fields[5].split(".."))[0]
    startPosLIGMDB = int ( fields[5].split("..")[0] )
    return startPosLIGMDB

def getFunction(seqDescription):
    """
    Return Functionality
    """
    fields = seqDescription.split("|")
    strFunction = fields[3]
    return strFunction

def genAnchDict(url):
    recordsAnch = get_records_list(url)
    seqsDictAnch = {}
    for rec in recordsAnch:
        seqsDictAnch[genKey(rec.description)] = getStartPos(rec.description)
    #print("genAnchDict", seqsDictAnch)
    return seqsDictAnch

def get_gene_anchors(specie: str, gene: str, imgtlabel: str, modelspath=None, filename=None, imgt_genedb=imgt_params['url.genedb']):
    """
    Download genomic anchors for V or J according to specie and gene from imgt database
    :return : string of the requested link
    """
    try:
        #QUERY_V_2CYS = "query=8.1+" + gene + "&species=" + specie + "&IMGTlabel=2nd-CYS"
        url_query = get_genedb_query81_imgtlabel(specie, gene, imgtlabel, imgt_genedb=imgt_genedb)
        print(specie, gene, url_query)

        # 2. load data from IMGT url
        records = get_records_list(url_query)
        records = map(lambda x: x.upper(), records)
        return records
    except Exception as e:
        raise e

def download_genes_anchors(specie: str, chain: str, flnVGenome, flnJGenome, ref_genes_path=None, modelspath=None, imgt_genedb=imgt_params['url.genedb']):
    """
    Download gene anchors from IMGT website in a IGoR structure directory (modelspath+"/"+specie+"/"+chain)
    """

    V_dict = download_Vgene_anchors(specie, chain, flnVGenome, ref_genes_path=ref_genes_path, modelspath=modelspath, imgt_genedb=imgt_genedb)
    J_dict = download_Jgene_anchors(specie, chain, flnJGenome, ref_genes_path=ref_genes_path, modelspath=modelspath, imgt_genedb=imgt_genedb)
    Anchors_dict = {**V_dict, **J_dict}
    return Anchors_dict

def download_Vgene_anchors(specie: str, chain: str, flnVGenome, ref_genes_path=None, modelspath=None, imgt_genedb=imgt_params['url.genedb'], sep=';'):
    #records = get_gene_template(specie, gene, imgt_genedb=imgt_genedb)
    if ref_genes_path is None:
        ref_genes_path = ""
        if modelspath is None:
            modelspath = "models"
        make_igor_directories(chain, specie, modelspath=modelspath)
        ref_genes_path = modelspath + "/" + specie + "/" + chain + "/ref_genome/"
    else:
        ref_genes_path = ref_genes_path + "/"


    # generate dictionaries with the anchors
    # 2nd-CYS
    urlV_2CYS = get_genedb_query81_imgtlabel(specie, chain + "V", imgtlabel="2nd-CYS", imgt_genedb=imgt_genedb)
    dictV_2CYS = genAnchDict(urlV_2CYS)

    flnAnchors = ref_genes_path + "V_gene_CDR3_anchors__imgt.csv"
    CSVDELIM = sep
    ofileAnch = open(flnAnchors, "w")
    ofileAnch.write("gene"+CSVDELIM+"anchor_index"+CSVDELIM+"gfunction" + "\n")


    Vrecords = SeqIO.parse(flnVGenome, "fasta")
    for rec in Vrecords:
        key = genKey(rec.description)
        # Initiate to avoid bad designation.
        posV_2CYS = -1
        if key in dictV_2CYS.keys():
            posV_2CYS = dictV_2CYS[key] - getStartPos(rec.description)
            posAnch = posV_2CYS
            seqFunc = getFunction(rec.description)
            ofileAnch.write(rec.description + CSVDELIM + str(posAnch) +  CSVDELIM+str(seqFunc)+"\n")
        else:
            if flag_verbose:
                print("No anchor is found for : " + rec.description)

    ofileAnch.close()
    fln_dict = dict()
    fln_dict['flnVGenome'] = flnVGenome
    fln_dict['flnVAnchors'] = flnAnchors
    fln_dict['ref_genes_path'] = ref_genes_path
    fln_dict['modelspath'] = modelspath
    return fln_dict

def download_Vgene_anchors_bk(specie: str, chain: str, modelspath=None, imgt_genedb=imgt_params['url.genedb']):
    #records = get_gene_template(specie, gene, imgt_genedb=imgt_genedb)
    if modelspath is None:
        modelspath = "models"
    make_igor_directories(chain, specie, modelspath=modelspath)
    ref_genes_path = modelspath + "/" + specie + "/" + chain + "/ref_genome/"

    # generate dictionaries with the anchors
    # 2nd-CYS
    urlV_2CYS = get_genedb_query81_imgtlabel(specie, chain + "V", imgtlabel="2nd-CYS", imgt_genedb=imgt_genedb)
    dictV_2CYS = genAnchDict(urlV_2CYS)

    flnAnchors = ref_genes_path + "V_gene_CDR3_anchors__imgt.csv"
    return dictV_2CYS
    # ofileAnch = open(flnAnchors, "w")
    # CSVDELIM = ";"
    # ofileAnch.write("gene"+CSVDELIM+"anchor_index"+CSVDELIM+"function" + "\n")

    # Vrecords = SeqIO.parse(flnVGenome, "fasta")
    # for rec in Vrecords:
    #     key = genKey(rec.description)
    #     # Initiate to avoid bad designation.
    #     posV_2CYS = -1
    #     if key in dictV_2CYS.keys():
    #         posV_2CYS = dictV_2CYS[key] - getStartPos(rec.description)
    #         posAnch = posV_2CYS
    #         seqFunc = getFunction(rec.description)
    #         ofileAnch.write(rec.description + CSVDELIM + str(posAnch) + CSVDELIM+str(seqFunc)+"\n")
    #     else:
    #         print("No anchor is found for : " + rec.description)
    #
    # ofileAnch.close()
    # fln_dict = dict()
    # fln_dict['flnVGenome'] = flnVGenome
    # fln_dict['flnVAnchors'] = flnAnchors
    # fln_dict['ref_genes_path'] = ref_genes_path
    # fln_dict['modelspath'] = modelspath
    # return fln_dict


def download_Jgene_anchors(specie: str, chain: str, flnJGenome, ref_genes_path=None, modelspath=None, imgt_genedb=imgt_params['url.genedb']):
    # records = get_gene_template(specie, gene, imgt_genedb=imgt_genedb)
    if ref_genes_path is None:
        ref_genes_path = ""
        if modelspath is None:
            modelspath = "models"
            make_igor_directories(chain, specie, modelspath=modelspath)
        ref_genes_path = modelspath + "/" + specie + "/" + chain + "/ref_genome/"
    else:
        ref_genes_path = ref_genes_path + "/"


    # J-PHE
    urlJ_PHE = get_genedb_query81_imgtlabel(specie, chain + "J", imgtlabel="J-PHE", imgt_genedb=imgt_genedb)
    dictJ_PHE = genAnchDict(urlJ_PHE)

    # J-TRP
    urlJ_TRP = get_genedb_query81_imgtlabel(specie, chain + "J", imgtlabel="J-TRP", imgt_genedb=imgt_genedb)
    dictJ_TRP = genAnchDict(urlJ_TRP)

    # write the files
    flnAnchors = ref_genes_path + "J_gene_CDR3_anchors__imgt.csv"
    ofileAnch = open(flnAnchors, "w")
    # ofileAnch.write("gene;anchor_index" + "\n")
    ofileAnch.write("gene;anchor_index;function" + "\n")
    CSVDELIM = ";"
    Jrecords = SeqIO.parse(flnJGenome, "fasta")
    for rec in Jrecords:
        key = genKey(rec.description)
        # Initiate to avoid bad designation.
        posJ_PHE = -1
        posJ_TRP = -1
        if key in dictJ_PHE.keys():
            posJ_PHE = dictJ_PHE[key] - getStartPos(rec.description)
            posAnch = posJ_PHE
            seqFunc = getFunction(rec.description)
            ofileAnch.write(rec.description + CSVDELIM + str(posAnch) + CSVDELIM + str(seqFunc) + "\n")
            # ofileAnch.write(rec.description + CSVDELIM + str(posAnch) + "\n")

        # FIXME: this should be elif, but I want to check If the system is consistent
        elif key in dictJ_TRP.keys():
            posJ_TRP = dictJ_TRP[key] - getStartPos(rec.description)
            posAnch = posJ_TRP
            seqFunc = getFunction(rec.description)
            ofileAnch.write(rec.description + CSVDELIM + str(posAnch) + CSVDELIM + str(seqFunc) + "\n")
            # ofileAnch.write(rec.description + CSVDELIM + str(posAnch) + "\n")
        else:
            if flag_verbose:
                print("No anchor is found for : " + rec.description)

    ofileAnch.close()
    fln_dict = dict()
    fln_dict['flnJGenome'] = flnJGenome
    fln_dict['flnJAnchors'] = flnAnchors
    fln_dict['ref_genes_path'] = ref_genes_path
    fln_dict['modelspath'] = modelspath
    return fln_dict

############ Donwload all genetic information if VDJ or VJ
def download_ref_genome_VDJ(species: str, chain: str, **kwargs):
    """
    Return a dictionary with genomics dataframes and also save files in a IGoR directory structure.
    :param species: IMGT species name
    :param chain: IMGT chain receptor name
    :return : dictionary of pandas DataFrame with
    """
    # flag_verbose = False
    dictVGenome = download_gene_template(species, chain + 'V', **kwargs)
    dictDGenome = download_gene_template(species, chain + 'D', **kwargs)
    dictJGenome = download_gene_template(species, chain + 'J', **kwargs)

    flnVGenome = dictVGenome["flnVGenome"]
    flnDGenome = dictDGenome["flnDGenome"]
    flnJGenome = dictJGenome["flnJGenome"]

    # FIXME: ONCE THE GENE TEMPLATES ARE DOWNLOADED CHANGE THE NAME TO
    # write anchors
    Anchors_dict = download_genes_anchors(species, chain, flnVGenome, flnJGenome, **kwargs)

    # TODO: filter sequences for OLGA compatibility
    strGene = chain + "V"
    urlV_2CYS = get_genedb_query81_imgtlabel(species, strGene, imgtlabel="2nd-CYS")
    dict2CYS = genAnchDict(urlV_2CYS)
    list2CYS = list(dict2CYS.keys())
    # print("------------> list2CYS : ", list2CYS)

    v_genomes_list = load_records_from_fasta(flnVGenome)
    # print("------------D v_genomes_list: ", v_genomes_list)
    v_genomes_trim_list = list()
    for v_genome in v_genomes_list:
        if genKey(v_genome.description) in list2CYS:
        # if v_genome in list2CYS:
            v_genomes_trim_list.append(v_genome)
        else:
            if flag_verbose:
                print(v_genome.description)
    save_records2fasta(v_genomes_trim_list, flnVGenome + "_trim")

    # J-PHE
    urlJ_PHE = get_genedb_query81_imgtlabel(species, species + "J", imgtlabel="J-PHE")
    dictJ_PHE = genAnchDict(urlJ_PHE)
    listJ_PHE = list(dictJ_PHE.keys())

    # J-TRP
    urlJ_TRP = get_genedb_query81_imgtlabel(species, species + "J", imgtlabel="J-TRP")
    dictJ_TRP = genAnchDict(urlJ_TRP)
    listJ_TRP = list(dictJ_TRP.keys())

    j_genomes_list = load_records_from_fasta(flnJGenome)
    j_genomes_trim_list = list()
    for j_genome in j_genomes_list:
        if genKey(j_genome.description) in listJ_PHE:
        #if j_genome in listJ_PHE:
            j_genomes_trim_list.append(j_genome)
        elif genKey(j_genome.description) in listJ_TRP:
        # elif j_genome in listJ_TRP:
            j_genomes_trim_list.append(j_genome)
        else:
            if flag_verbose:
                print(j_genome.description)

    save_records2fasta(j_genomes_trim_list, flnJGenome + "_trim")
    if flag_verbose:
        print("----------------------")
        print("Genomic VDJ templates in files: ")
        print(flnVGenome, flnDGenome, flnJGenome)
    dict_V = gen_short_names(flnVGenome, flnAnchors=Anchors_dict["flnVAnchors"])
    dict_J = gen_short_names(flnJGenome, flnAnchors=Anchors_dict["flnJAnchors"])
    dict_D = gen_short_names(flnDGenome)

    import os.path
    V_path = os.path.dirname(dict_V["genomics"])
    J_path = os.path.dirname(dict_J["genomics"])
    D_path = os.path.dirname(dict_D["genomics"])
    import shutil
    shutil.copy(dict_V["genomics"], V_path + "/genomicVs.fasta")
    shutil.copy(dict_J["genomics"], J_path + "/genomicJs.fasta")
    shutil.copy(dict_D["genomics"], D_path + "/genomicDs.fasta")
    shutil.copy(dict_V["anchors"], V_path + "/V_gene_CDR3_anchors.csv")
    shutil.copy(dict_J["anchors"], J_path + "/J_gene_CDR3_anchors.csv")

    from .utils import get_dataframe_from_fasta_and_csv_anchors
    df_genomics_dict = dict()
    df_genomics_dict['V'] = get_dataframe_from_fasta_and_csv_anchors(dict_V["genomics"], dict_V["anchors"])
    df_genomics_dict['D'] = get_dataframe_from_fasta_and_csv_anchors(dict_D["genomics"])
    df_genomics_dict['J'] = get_dataframe_from_fasta_and_csv_anchors(dict_J["genomics"], dict_J["anchors"])
    return df_genomics_dict

    # print(kwargs)

def download_ref_genome_VJ(species: str, chain: str, **kwargs):
    """Download gene templates and anchors from IMGT and creates
    different files with original and short names.

    :param species: Species name in IMGT nomenclature,
    :type species: str
    :param chain: Chain name in IMGT nomenclature,
    :type chain: str
    :param modelspath: root paths for all models.
    :type modelspath: str, optional
    :param filename: gene template fasta filename.
    :type filename: str, optional
    :param imgt_genedb: Url of IMGT GeneDB web application.
    :type imgt_genedb: str
    """
    dictVGenome = download_gene_template(species, chain + 'V', **kwargs)
    dictJGenome = download_gene_template(species, chain + 'J', **kwargs)

    flnVGenome = dictVGenome["flnVGenome"]
    flnJGenome = dictJGenome["flnJGenome"]

    # FIXME: ONCE THE GENE TEMPLATES ARE DOWNLOADED CHANGE THE NAME TO
    # write anchors
    Anchors_dict = download_genes_anchors(species, chain, flnVGenome, flnJGenome, **kwargs)

    # TODO: filter sequences for OLGA compatibility
    strGene = chain + "V"
    urlV_2CYS = get_genedb_query81_imgtlabel(species, strGene, imgtlabel="2nd-CYS")
    dict2CYS = genAnchDict(urlV_2CYS)
    list2CYS = list(dict2CYS.keys())
    # print("------------> list2CYS : ", list2CYS)

    v_genomes_list = load_records_from_fasta(flnVGenome)
    # print("------------D v_genomes_list: ", v_genomes_list)
    v_genomes_trim_list = list()
    for v_genome in v_genomes_list:
        if genKey(v_genome.description) in list2CYS:
        # if v_genome in list2CYS:
            v_genomes_trim_list.append(v_genome)
        else:
            print(v_genome)
    save_records2fasta(v_genomes_trim_list, flnVGenome + "_trim")

    # J-PHE
    urlJ_PHE = get_genedb_query81_imgtlabel(species, species + "J", imgtlabel="J-PHE")
    dictJ_PHE = genAnchDict(urlJ_PHE)
    listJ_PHE = list(dictJ_PHE.keys())

    # J-TRP
    urlJ_TRP = get_genedb_query81_imgtlabel(species, species + "J", imgtlabel="J-TRP")
    dictJ_TRP = genAnchDict(urlJ_TRP)
    listJ_TRP = list(dictJ_TRP.keys())

    j_genomes_list = load_records_from_fasta(flnJGenome)
    j_genomes_trim_list = list()
    for j_genome in j_genomes_list:
        if genKey(j_genome.description) in listJ_PHE:
        #if j_genome in listJ_PHE:
            j_genomes_trim_list.append(j_genome)
        elif genKey(j_genome.description) in listJ_TRP:
        # elif j_genome in listJ_TRP:
            j_genomes_trim_list.append(j_genome)
        else:
            print(j_genome)

    save_records2fasta(j_genomes_trim_list, flnJGenome + "_trim")

    print("----------------------")
    print("Genomic VJ templates in files: ")
    print(flnVGenome, flnJGenome)
    dict_V = gen_short_names(flnVGenome, flnAnchors=Anchors_dict["flnVAnchors"])
    dict_J = gen_short_names(flnJGenome, flnAnchors=Anchors_dict["flnJAnchors"])

    import os.path
    V_path = os.path.dirname(dict_V["genomics"])
    J_path = os.path.dirname(dict_J["genomics"])
    import shutil
    shutil.copy(dict_V["genomics"], V_path + "/genomicVs.fasta")
    shutil.copy(dict_J["genomics"], J_path + "/genomicJs.fasta")
    shutil.copy(dict_V["anchors"], V_path + "/V_gene_CDR3_anchors.csv")
    shutil.copy(dict_J["anchors"], J_path + "/J_gene_CDR3_anchors.csv")

    from .utils import get_dataframe_from_fasta_and_csv_anchors
    df_genomics_dict = dict()
    df_genomics_dict['V'] = get_dataframe_from_fasta_and_csv_anchors(dict_V["genomics"], dict_V["anchors"])
    df_genomics_dict['J'] = get_dataframe_from_fasta_and_csv_anchors(dict_J["genomics"], dict_J["anchors"])
    return df_genomics_dict


# FIXME: WORKING ON THIS
def download_ref_genome(species: str, chain: str, **kwargs):
    """
    FIXME: IN DEV
    TODO: CHECK THE IMGT STANDARDS TO DEFINE IG AND TR
    """
    dictVGenome = download_gene_template(species, chain + 'V', **kwargs)
    dictJGenome = download_gene_template(species, chain + 'J', **kwargs)
    pass

def get_dict_from_imgt_description(str_description):
    """
    Return an OderedDict with following the fields
    The IMGT FASTA header of nucleotide IMGT  reference sequences contains 15 fields separated by '|':
    1. IMGT/LIGM-DB accession number(s)
    2. IMGT gene and allele name
    3. species
    4. IMGT allele functionality
    5. exon(s), region name(s), or extracted label(s)
    6. start and end positions in the IMGT/LIGM-DB accession number(s)
    7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
    8. codon start, or 'NR' (not relevant) for non coding labels
    9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
    10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
    11. +n, -n, and/or nS: number of added, deleted, and/or substituted nucleotides to correct sequencing errors, or 'not corrected' if non corrected sequencing errors
    12. number of amino acids (AA): this field indicates that the sequence is in amino acids
    13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
    14. partial (if it is)
    15. reverse complementary (if it is)
    """

    description_field_list = ["accession_number",  # 1
             "allele_name",  # 2
             "species",  # 3
             "function",  # 4
             "exon_region_name_label",  # 5
             "start_end",  # 6
             "n_nt",  # 7
             "codon_start",  # 8
             "n_nt_added_5p", # 9
             "n_nt_added_removed_3p", # 10
             "n_nt_added_deleted_substituted",  # 11
             "n_aa",  # 12
             "n_char_plus_gaps",  # 13
             "partial",  # 14
             "reverse_complementary"  # 15
             ]
    from collections import OrderedDict
    dicto = OrderedDict()
    description_list = str_description.split("|")
    print(len(description_field_list), description_field_list)
    print(len(description_list), description_list)
    dicto["description"] = str_description
    for ii, elem in enumerate(description_field_list):
        print(ii, elem, description_list[ii])
        dicto[elem] = description_list[ii]

    return dicto



def make_VDJ_model():
    from .IgorIO import IgorModel_Parms
    mdl_parms = IgorModel_Parms

def gen_short_names(flnGenome, flnAnchors=None):
    import pandas as pd
    from .IgorIO import genLabel
    # flnGenome = options.gene_template
    # flnAnchors = options.gene_anchors

    records = load_records_from_fasta(flnGenome)
    new_records = list()
    for record in records:
        record.description = genLabel(record.description)
        record.id = record.description
        new_records.append(record)
    flnGenome_short = flnGenome + "_short"
    save_records2fasta(new_records, flnGenome_short)

    gene_dict = dict()
    gene_dict["genomics"] = flnGenome_short

    try:
        if flnAnchors is not None:
            df = pd.read_csv(flnAnchors, sep=';')
            #print(df['gene'])
            df['gene'] = df['gene'].map(genLabel)
            # print(df['gene'])
            flnAnchors_short = flnAnchors + "_short"
            df.to_csv(flnAnchors_short, sep=';', index=False)
            gene_dict["anchors"] = flnAnchors_short
    except Exception as e:
        print(e)


    return gene_dict

