# pre-processing script for sequence alignment through pyir, wrapper of Igblast
# Giulio Isacchini, Francesco Camaglia 

import os, sys
import subprocess
import numpy as np
import pandas as pd

class pre_process_gear( pygor_class, WDpath, igdata=None, verbose=False ) :
    
        '''
        It creates an IgBlast alignment file ???.csv.gz in the WDpath using the wrapper pyir.
        
        Parameters :
        ------------
        
        ???
        pygor_class :
        
        {
            reads_data_frame : pd.DataFrame
                Pandas DatFrame containing indexed reads.

            receptor : str, choice    
               Choose the receptor, "Ig" for immunoglobulin or "TCR" for T cell receptor (default).

            specie : str, choice
                Choose the specie, "mouse" or "human" (default).
        }
            
        WDpath : path/to/dir
            Sepcify the working directory path.

        igdata : str, optional
            Path to your custom IGDATA directory. Default is None

        verbose : bool
            Wheter to provide all passages informations.   
            
        '''
    
    # file where to store things (without extension)
    batchname = f"{WDPATH}/" ???

    if options.verbose : print ( 'Loading indexed sequences...' )
    self._dataFrame_to_fasta( self, batchname )

    if options.verbose : print ( 'Aligning sequences...' )
    blastname = self._align_seqs( self, batchname, igdata )

    if options.verbose : print ( '> Processing IgBLAST output sequences...' )
    self._process_seqs( blastname )

    if options.verbose : print( '> Completed.\n' )

        
    # >>>>>>>>>>>>>>>>>
    #  SAVE TO FASTA  #
    # >>>>>>>>>>>>>>>>>

    def _dataFrame_to_fasta( self, batchname ) :
        '''
        It saves the reads to a fasta file named after the batchname.
        Indeces are taken from the dataframe.
        '''

        path_to_fasta = f"{batchname}.fasta"
        with open( path_to_fasta, "w" ) as fw:
            for indx, seq in zip( reads_data_frame.index, reads_data_frame.values ):
                fw.write( ">{}\n".format( indx ) )
                fw.write( f"{seq}\n" )
    ###


    # >>>>>>>>>>>>>
    #  ALIGN SEQS #
    ## >>>>>>>>>>>>

    def _align_seqs( self, batchname, igdata=None ):
        '''
        Call pyir wrap of IgBLAST to perform sequence alignment.
        ''' 

        filein = f"{batchname}.fasta"
        fileout = f"{batchname}-blasted"

        # define the pre-installed pyir command for alignment
        ???
        pyir_cmd = f"pyir {filein} -o {fileout} --outfmt tsv  -r {self.receptor} -s {self.specie}"
        if igdata is not None : pyir_cmd = f"{pyir_cmd} --igdata {igdata}"

        # Run the terminal instruction of pyir within Python
        ???
        subprocess.Popen( pyir_cmd, shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE ).communicate( )[1]
        os.remove( filein )       

        return fileout
    ###


    # >>>>>>>>>>>>>>>>
    #  PROCESS SEQS  #
    # >>>>>>>>>>>>>>>>

    def _process_seqs( blastname ):
        '''
        It takes PyIR output and translates into a csv working file.
        '''
        
        # open temporary igblast alignment file
        aligned = pd.read_csv( f"{blastname}.tsv.gz", sep="\t", compression="gzip", dtype=str ) 

        # choose which igblast output to consider
        keep = [ 'sequence_id', 'junction_aa', 'v_call', 'j_call', 'junction',
                'locus', 'stop_codon', 'vj_in_frame', 'productive', 'rev_comp', 'sequence' ]
        
        processed = aligned[ keep ]
        processed.to_csv( f'{blastname}.csv.gz' , sep=",", index=False, compression="gzip" )

        # remove temporary igblast alignment file
        os.remove( f"{blastname}.tsv.gz" )                                           
    ### 
###