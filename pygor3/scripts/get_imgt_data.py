#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 18:58:54 2019

@author: alfaceor
"""
import pygor3.imgt as imgt

def main():
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog [options] ")
    parser.add_option("-t", "--type",   dest="type",   help="VJ or VDJ")
    parser.add_option("-g", "--gene",   dest="gene",   help="Type of gene like TRB")
    parser.add_option("-s", "--specie", dest="specie", help="Type of specie")

    (options, args) = parser.parse_args()

    #options.type   = "VDJ"
    #options.gene   = "TRB"
    #options.specie = "Mus+musculus"


    if options.type == "VDJ":
        if options.gene == None:
            print("Gene option is missing: --gene TRB")
        else:
            if options.specie == None:
                print("Specie option is missing: --specie Mus+musculus")
            else:
                imgt.make_igor_directories(options.gene, options.specie)
                imgt.write_V_TemplateFiles(options.gene, options.specie)
                imgt.write_D_TemplateFiles(options.gene, options.specie)
                imgt.write_J_TemplateFiles(options.gene, options.specie)
                # Now write the model. by using the template files.
                imgt.writeModelVDJ_Parms(options.gene, options.specie)
                imgt.writeModelVDJ_Marginals(options.gene, options.specie)

    elif options.type == "VJ":
        if options.gene == None:
            print("Gene option is missing: --gene TRB")
        else:
            if options.specie == None:
                print("Specie option is missing: --specie Mus+musculus")
            else:
                imgt.make_igor_directories(options.gene, options.specie)
                imgt.write_V_TemplateFiles(options.gene, options.specie)
                imgt.write_J_TemplateFiles(options.gene, options.specie)
                imgt.writeModelVJ_parms(options.gene, options.specie)
                imgt.writeModelVJ_Marginals(options.gene, options.specie)
    else:
        print("Too bad you need to specify a type : --type VDJ")


if __name__ == "__main__":
    main()