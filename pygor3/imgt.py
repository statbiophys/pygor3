#!/usr/bin/python3
from __future__ import print_function
import os
# import IgorModel

import requests
from bs4 import BeautifulSoup
from Bio import SeqIO

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

# def extract_fasta_from_url(url):
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

def makeDirectories(gene: str, specie: str, modelspath=None):
    if modelspath is None:
        modelspath = "models"

    os.system("mkdir -p " + modelspath)
    os.system("mkdir -p " + modelspath + "/" + specie)
    os.system("mkdir -p " + modelspath + "/" + specie + "/" + gene)
    os.system("mkdir -p " + modelspath + "/" + specie + "/" + gene + "/ref_genome")
    os.system("mkdir -p " + modelspath + "/" + specie + "/" + gene + "/ref_genome")
    os.system("mkdir -p " + modelspath + "/" + specie + "/" + gene + "/models")

def get_gene_template(specie: str, gene: str, modelspath=None, filename=None, imgt_genedb=imgt_params['url.genedb']):
    """
    Download genomic template according to specie and gene from imgt database
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

def download_gene_template(specie: str, gene: str, modelspath=None, filename=None, imgt_genedb=imgt_params['url.genedb']):
    # 3. Write it in a fasta file.
    records = get_gene_template(specie, gene, imgt_genedb=imgt_genedb)
    if modelspath is None:
        modelspath = "models"

    if gene[-1] in ['V', 'D', 'J']:
        if filename is None :
            makeDirectories(gene, specie, modelspath=modelspath)
            modelspath = modelspath + "/"+specie+"/"+gene+"/ref_genome/"
            filename = modelspath + "genomic" + gene[-1] + "s.fasta"
            #filename = modelspath + "genomic__" + specie + "_" + gene + ".fasta"
    #flnGenomicJs = DIR_REF_GENOME + "genomicJs.fasta"
    save_records2fasta(records, filename)

def write_D_TemplateFiles(gene, specie):
    records = get_records_list(urlDgene)

def write_V_TemplateFiles(gene, specie, ):
    DIR_REF_GENOME = specie + "/" + gene + "/ref_genome/"
    gene = gene + "V"
    print("DIR_REF_GENOME = ", DIR_REF_GENOME)
    # 1. Get the Alleles/Gene sequences from the IMGT website (IMGT/GENE-DB)
    URL_IMGT_GENEDB = "http://www.imgt.org/genedb/GENElect?"
    QUERY_GENE = "query=7.2+" + gene + "&species=" + specie
    QUERY_V_2CYS = "query=8.1+" + gene + "&species=" + specie + "&IMGTlabel=2nd-CYS"
    urlVgene = URL_IMGT_GENEDB + QUERY_GENE
    urlV_2CYS = URL_IMGT_GENEDB + QUERY_V_2CYS

    urlVgene = get_genedb_query72(specie, gene, imgt_genedb=imgt_params['url.genedb'])

    # 2. load data from IMGT url
    records = get_records_list(urlVgene)
    #records = genSeqRecords(urlVgene)
    records = map(lambda x: x.upper(), records)

    # 3. Write it in a fasta file.
    flnGenomicVs = DIR_REF_GENOME + "genomicVs.fasta"
    with open(flnGenomicVs, "w") as ofile:
        SeqIO.write(records, ofile, "fasta")

    # generate dictionaries with the anchors
    dictV_2CYS = genAnchDict(urlV_2CYS)
    # print( len(dictV_2CYS), len(records) )
    # write the files
    flnAnchors = DIR_REF_GENOME + "V_gene_CDR3_anchors.csv"
    ofileAnch = open(flnAnchors, "w")
    ofileAnch.write("gene;anchor_index" + "\n")
    CSVDELIM = ";"
    for rec in records:
        key = genKey(rec.description)
        # Initiate to avoid bad designation.
        posV_2CYS = -1
        if key in dictV_2CYS.keys():
            posV_2CYS = dictV_2CYS[key] - getStartPos(rec.description)
            posAnch = posV_2CYS
            ofileAnch.write(rec.description + CSVDELIM + str(posAnch) + "\n")

        else:
            print("No anchor is found for : " + rec.description)

    ofileAnch.close()

def download_gene_anchors(specie: str, gene: str, imgt_genedb=imgt_params['url.genedb'], modelpath=".", filename=None):
    url_query = get_genedb_query81_imgtlabel(specie, gene, imgt_genedb=imgt_genedb)

def write_V_TemplateFiles(gene, specie ):
    DIR_REF_GENOME = specie + "/" + gene + "/ref_genome/"
    gene = gene + "V"
    print("DIR_REF_GENOME = ", DIR_REF_GENOME)
    # 1. Get the Alleles/Gene sequences from the IMGT website (IMGT/GENE-DB)
    URL_IMGT_GENEDB = "http://www.imgt.org/genedb/GENElect?"
    QUERY_GENE = "query=7.2+" + gene + "&species=" + specie
    QUERY_V_2CYS = "query=8.1+" + gene + "&species=" + specie + "&IMGTlabel=2nd-CYS"
    urlVgene = URL_IMGT_GENEDB + QUERY_GENE
    urlV_2CYS = URL_IMGT_GENEDB + QUERY_V_2CYS

    # 2. load data from IMGT url
    records = genSeqRecords(urlVgene)
    records = map(lambda x: x.upper(), records)

    # 3. Write it in a fasta file.
    flnGenomicVs = DIR_REF_GENOME + "genomicVs.fasta"
    with open(flnGenomicVs, "w") as ofile:
        SeqIO.write(records, ofile, "fasta")

    # generate dictionaries with the anchors
    dictV_2CYS = genAnchDict(urlV_2CYS)
    # print( len(dictV_2CYS), len(records) )
    # write the files
    flnAnchors = DIR_REF_GENOME + "V_gene_CDR3_anchors.csv"
    ofileAnch = open(flnAnchors, "w")
    ofileAnch.write("gene;anchor_index" + "\n")
    CSVDELIM = ";"
    for rec in records:
        key = genKey(rec.description)
        # Initiate to avoid bad designation.
        posV_2CYS = -1
        if key in dictV_2CYS.keys():
            posV_2CYS = dictV_2CYS[key] - getStartPos(rec.description)
            posAnch = posV_2CYS
            ofileAnch.write(rec.description + CSVDELIM + str(posAnch) + "\n")

        else:
            print("No anchor is found for : " + rec.description)

    ofileAnch.close()


"""
    DIR_REF_GENOME = specie+"/"+gene+"/ref_genome/"
    gene = gene + "V"
    print ("DIR_REF_GENOME = ", DIR_REF_GENOME)
    # 1. Get the Alleles/Gene sequences from the IMGT website (IMGT/GENE-DB)
    URL_IMGT_GENEDB = "http://www.imgt.org/genedb/GENElect?"
    QUERY_GENE      = "query=7.2+"+gene+"&species="+specie
    QUERY_V_2CYS    = "query=8.1+"+gene+"&species="+specie+"&IMGTlabel=2nd-CYS"
    urlVgene        = URL_IMGT_GENEDB + QUERY_GENE
    urlV_2CYS       = URL_IMGT_GENEDB + QUERY_V_2CYS
    
    
    
    
    
    
    # 3. Write it in a fasta file.
    flnGenomicJs = DIR_REF_GENOME+"genomicJs.fasta"
    with open(flnGenomicJs, "w") as ofile:
        SeqIO.write(records, ofile, "fasta")
    
    # generate dictionaries with the anchors
    dictJ_PHE = genAnchDict(urlJ_PHE)
    dictJ_TRP = genAnchDict(urlJ_TRP)
    # write the files 
    flnAnchors = DIR_REF_GENOME+"J_gene_CDR3_anchors.csv"
    ofileAnch  = open(flnAnchors, "w")
    ofileAnch.write("gene;anchor_index"+"\n")
    CSVDELIM = ";"
    for rec in records:
        key = genKey(rec.description)
        # Initiate to avoid bad designation.
        posJ_PHE = -1
        posJ_TRP = -1
        if key in dictJ_PHE.keys():
            posJ_PHE =  dictJ_PHE[key] - getStartPos(rec.description) 
            posAnch  = posJ_PHE
            ofileAnch.write(rec.description+CSVDELIM+str(posAnch)+"\n")
    
        # FIXME: this should be elif, but I want to check If the system is consistent
        elif key in dictJ_TRP.keys():
            posJ_TRP = dictJ_TRP[key] -  getStartPos(rec.description) 
            posAnch  = posJ_TRP
            ofileAnch.write(rec.description+CSVDELIM+str(posAnch)+"\n")
        else:
            print("No anchor is found for : "+ rec.description)
    
    ofileAnch.close()
    
"""