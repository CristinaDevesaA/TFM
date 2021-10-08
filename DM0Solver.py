#!/usr/bin/env python
# coding: utf-8

# Module metadata variables
__author__ = "Cristina Amparo Devesa Arbiol"
__credits__ = ["Cristina Amparo Devesa Arbiol", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "cristinaamparo.devesa@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"

# Import modules
from os import remove
import sys
from optparse import OptionParser
import configparser
import pandas as pd
import numpy as np
import argparse
import os
import logging
from pathlib import Path

###################
# Local functions #
###################

def readInfile(infile):
    '''    
    Read input file to dataframe.
    '''
    
    df = pd.read_csv(infile, sep="\t", float_precision='high',low_memory=False)
    return df

def DM0Solver(Theo_mh,Exp_mh,seq,Error,dic_DM0,SequenceMassmod,sequence):
    
    """
    DM0Solver function returns the sequence and label, that match with the smallest error
    according to the dic_DM0 (dictionary with all possible options).
    Input:
         Theo_mh             > Theoretical mh  
         Exp_mh              > Experimental mh
         seq                 > Sequence with DM
         Error               > Error in ppms stablished by the user
         dic_DM0             > Dictionary with the options of DM0List, in which label is the key and the corresponding mass the value.

         
         
    The function returns:
        DM0Sequence          > Sequence with the DM that best fits
        DM0Label             > Label option that best fits
    
    
    """
    
    
    # These variables are assigned for the cases in which there is no match.
    DM0Sequence = sequence
    DM0Label = ""
    DM0Label_ppm = ""


    
  
                    
    # The dictionary option that obtains the smallest error is saved.
    minimun_DiffPPM = Error  # At the beginning the minimum difference (ppms) is considered as the Error  
    for Label in dic_DM0:
        DiffPPM = abs(((Theo_mh+float(dic_DM0[Label])-Exp_mh)*1000000)/(Theo_mh+float(dic_DM0[Label])))
        if DiffPPM <= Error:
            if DiffPPM < minimun_DiffPPM:
                minimun_DiffPPM = DiffPPM
                DM0Sequence = seq+"_"+str(SequenceMassmod) 
                DM0Label = Label
                DM0Label_ppm = minimun_DiffPPM

   
    return DM0Sequence,DM0Label,DM0Label_ppm

def applySolver(row, Seq_column_name, Theo_mh_column_name, Exp_mh_column_name, Error, dic_DM0,
                DM0Sequence_output_column_name, DM0Label_output_column_name, DM0Label_ppm_output_column_name):
    seq = row[Seq_column_name]
    SequenceMassmod = seq[seq.find("[")+1:]
    SequenceMassmod = SequenceMassmod[:SequenceMassmod.find("]")]  # Mass modification of the sequence is obtained.
    if row[Seq_column_name].find("_") != -1:  # If  scan has already been corrected DM0Seqeuence remains Seq_column_name
  	    seq = seq[:seq.find("_")]
        
    else:

        seq = seq[:seq.find("[")]+seq[seq.find("]")+1:]  # Clean sequence is obtained.# if  scan has not yet been corrected DM0Sover function is executed

    DM0Sequence,DM0Label,DM0Label_ppm = DM0Solver(row[Theo_mh_column_name],row[Exp_mh_column_name],seq,Error,dic_DM0,SequenceMassmod,row[Seq_column_name])
   
        
    # New columns are completed  
    row[DM0Sequence_output_column_name] = DM0Sequence
    row[DM0Label_output_column_name] = DM0Label.upper()
    row[DM0Label_ppm_output_column_name] = DM0Label_ppm
    return row

##################
# Main functions #
##################

def main(solverconfig,infile):
    
    """
    Reading configuration file and processing file

    """

    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(solverconfig) # Reading configuration file
    
    logging.info("Reading DM0Solver configuration file")
    

    Exp_mh_column_name = config["DM0Solver_Parameters"].get("exp_mh_column_name") # Experimental mh  column name
    Theo_mh_column_name = config["DM0Solver_Parameters"].get("theo_mh_column_name")# Theoretical mh column name 
    Seq_column_name = config["DM0Solver_Parameters"].get("Sequence_column_name") # Sequence column name

    DM0Sequence_output_column_name =  config["DM0Solver_Parameters"].get("DM0Sequence_output_column_name") # Column name of the output where the chosen sequence is annotated
    DM0Label_output_column_name =  config["DM0Solver_Parameters"].get("DM0Label_output_column_name") # Column name of the output where the chosen label is annotated
    DM0Label_ppm_output_column_name =  config["DM0Solver_Parameters"].get("DM0Label_ppm_output_column_name") # Column name of the output where the calculated error in ppm is annotated
    output_file_suffix = config["DM0Solver_Parameters"].get("output_file_suffix") # Chosen suffix for output file 
    Error = config["DM0Solver_Parameters"].getfloat("Relative_Error_ppm") # Relative error (ppm)


 
    # A dictionary that save the labels and their corresponding masses is created
    dic_DM0 = {}   
    for option, value in config["DM0Solver_DM0List"].items():
        dic_DM0[option.strip("\n")] = value.strip("\n")
            
            
        
    # Input file is read as data frame and column names are added to the header
    df = readInfile(infile)
    try: 
        df.drop([DM0Sequence_output_column_name,DM0Label_output_column_name,DM0Label_ppm_output_column_name],axis=1)
    except:
        pass
    df[DM0Sequence_output_column_name] = ""
    df[DM0Label_output_column_name] = ""
    df[DM0Label_ppm_output_column_name] = np.nan


    #cont = 0
    logging.info("Processing input file")
    df = df.apply(lambda x: applySolver(x, Seq_column_name, Theo_mh_column_name, Exp_mh_column_name, Error, dic_DM0,
                                        DM0Sequence_output_column_name, DM0Label_output_column_name, DM0Label_ppm_output_column_name), axis = 1)

    # write outputfile
    logging.info("Writing output file")

    outfilename = infile1[:-4]+output_file_suffix
    outfile= outfilename+".txt"
    try:
        remove("outfile")
    except:
        None
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')

    logging.info('end script')



if __name__ == '__main__':
    


    # parse arguments
    parser = argparse.ArgumentParser(
        description='DM0Solver',
        epilog='''
        Example:
            python DM0Solver.py
        ''')
      
    # Default DM0Solver configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom configuration file')
    
    # These will overwrite the config if specified
    parser.add_argument('-r', '--relativeerror', help='Maximum allowable relative error (ppm)')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
   
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.relativeerror is not None:
        config.set('DM0Solver_Parameters', 'relative_error_ppm', str(args.relativeerror))
        config.set('Logging', 'create_ini', '1')
   
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/config/Solver.ini', 'w') as newconfig:
            config.write(newconfig)
        
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + '_DM0Solved_log.txt'
    log_file_debug = outfile = args.infile[:-4] + '_DM0Solved_log_debug.txt'
    if args.verbose:
        logging.basicConfig(level = logging.DEBUG,
                            format = '%(asctime)s - %(levelname)s - %(message)s',
                            datefmt = '%m/%d/%Y %I:%M:%S %p',
                            handlers = [logging.FileHandler(log_file_debug),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(level = logging.INFO,
                            format = '%(asctime)s - %(levelname)s - %(message)s',
                            datefmt = '%m/%d/%Y %I:%M:%S %p',
                            handlers = [logging.FileHandler(log_file),
                                      logging.StreamHandler()])

    infile1 = args.infile

    #start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    

    main(args.config, infile1)
