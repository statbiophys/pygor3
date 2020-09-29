#!/usr/bin/env python3
# modify genomic templates and anchors to short_names
import pygor3 as p3
from Bio import SeqIO
import pandas as pd

def main():
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog [options] ")
    parser.add_option("-t", "--gene_template",   dest="gene_template",   help="fasta file with gene templates")
    parser.add_option("-a", "--gene_anchors",   dest="gene_anchors",   help="csv file with gene anchors")

    (options, args) = parser.parse_args()

    flnGenome = options.gene_template
    flnAnchors = options.gene_anchors

    records = p3.imgt.load_records_from_fasta(flnGenome)
    new_records = list()
    for record in records:
        record.description = p3.genLabel(record.description)
        record.id = record.description
        new_records.append(record)
    p3.imgt.save_records2fasta(new_records, flnGenome+"_short")

    try:
        df = pd.read_csv(flnAnchors, sep=';')
        print(df['gene'])
        df['gene'] = df['gene'].map(p3.genLabel)
        print(df['gene'])
        df.to_csv(flnAnchors+"_short", sep=';', index=False)
    except Exception as e:
        print(e)


if __name__ == "__main__":
    main()
