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

#imports
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





def all_labels_joiner(all_labels_joined):
    
    """
    all_labels_joiner returs a clean string with unique labels  separated by comas
    """
    all_labels_joined = pd.unique(all_labels_joined)
    all_labels_joined = [i for i in all_labels_joined if i ]
    all_labels_joined = ";".join(all_labels_joined)
    
    return  str(all_labels_joined)


##################
# Main functions #
##################

def main(solverconfig,infile):
    
    """
    Reading configuration file and processing file
    """
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(solverconfig) # Reading configuration file
    
    Output_colum_name = config["Joiner_Parameters"].get("Output_column_name")
    Output_file_suffix = config["Joiner_Parameters"].get("Output_file_suffix")
    
    logging.info("Reading Joiner configuration file")
    column_names=[]
    for option,value in config["Joiner_Columns"].items():
        column_names.append(value.strip("\n"))

    logging.info("Processing input file")
    
    df = readInfile(infile)
    df.fillna('', inplace  = True)
    df[Output_colum_name] = ""

    cont = 0
    for index, row in df.iterrows():
       
        joined_labels = []
        for item in column_names:
            item = row[item].strip(" ")
            joined_labels.append(item)



        All_labels_joined = all_labels_joiner(joined_labels)
   

        df.loc[cont,Output_colum_name] = All_labels_joined
        cont = cont+1

   

    logging.info("Writing output file")
    outfilename = infile[:-4]+Output_file_suffix
    outfile = outfilename + ".txt"
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')






if __name__ == '__main__':
    
    try:
        remove('Solver.ini')
    except:
        None

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Joiner',
        epilog='''
        Example:
            python Joiner.py
        ''')
      
    # Default DM0Solver configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    
    # These will overwrite the config if specified
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
   
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)

    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + '_Joiner_log.txt'
    log_file_debug = outfile = args.infile[:-4] + '_Joiner_log_debug.txt'
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
    main(args.config,args.infile)