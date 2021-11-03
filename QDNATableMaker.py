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
import pandas as pd 
from optparse import OptionParser
import configparser
import argparse
import os
import logging
import sys

###################
# Local functions #
###################
def readInfile(infile):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile, sep="\t", float_precision='high',low_memory=False)

    return df



##################
# Main functions #
##################
def main(infile,file):
        
    """
    Reading configuration file
    """

    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(file)
    logging.info("Reading QDNATableMaker configuration file")  
        
    output_name = config["QDNATableMaker_Parameters"].get("Output_file_name") # Sequence with DM column name
    
    logging.info("Processing input file")
    df1 = readInfile(infile) 

    df1.drop(['p', 'pdm',"pFreq","m","l","M","L","ScanFreq","Theo_mh","pd","f","qf","qfFreq"], axis = 'columns', inplace=True) # The columns that are not desired are deleted

    df1_without_duplicates = df1.drop_duplicates(subset=['qdna']) # qdna duplicates are deleted

    df1_correct_order = df1_without_duplicates.reindex(columns=['q','qdna','qFreq','d','a','n','b','e','first_b','first_n','k','qna','qk','qdnaFreq','qnaFreq','qkFreq','A','N','qdNA','qNA'])# new colum order is configured
    
    logging.info("Writing output file")
    df1_correct_order.to_csv(output_name+".txt", index=False, sep='\t', encoding='utf-8')





if __name__ == '__main__':
    
    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='QDNATableMaker',
        epilog='''
        Example:
            python QDNATbleMaker.py
        ''')
      
    # default DM0Solver configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='Input file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    
    args = parser.parse_args()
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + 'QDNATable_log.txt'
    log_file_debug = outfile = args.infile[:-4] + 'QDNATable_log_debug.txt'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file_debug),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file),
                                      logging.StreamHandler()])

    
    

    
    
   
    #start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    

    main (args.infile,args.config)

