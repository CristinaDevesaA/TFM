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
import re
import sys
from pathlib import Path
from optparse import OptionParser
import configparser
import argparse
import os
import logging
from pathlib import Path




###################
# Local functions #
###################
def filterapply (infile,GlobalThres, LocalThres, PeakThres, GlobalFDR_colum_name,PeakFDR_colum_name,LocalFDR_colum_name,Label_colum_name,Decoys_naming):
    
    """
    Filterapply function returns an output file containing only the lines that meet determined 
    conditions.  
    """
    for line in open(infile):
        line = line.strip("\n")
        file = open(line,"r")
        newfile = open((line.split("/")[-1].strip(".txt")) + "_FDRfiltered.txt","w")
        cont = 0
        for line1 in file:
            if cont == 0:
                newfile.write(line1)
                header = line1.strip("\n").split("\t")
            elif cont >= 1:
                fields = line1.strip("\n").split("\t")
  
                if float(fields[header.index(GlobalFDR_colum_name)]) < GlobalThres and (float(fields[header.index(LocalFDR_colum_name)]) < LocalThres or float(fields[header.index(PeakFDR_colum_name)]) < PeakThres) and str(fields[header.index(Label_colum_name)]).upper() != Decoys_naming.upper():        
                    newfile.write(line1)
            cont = cont+1
             


##################
# Main functions #
##################

def main(file,infile):
    
    """
    Reading configuration file
    """
    config = configparser.ConfigParser()
    config = configparser.ConfigParser(inline_comment_prefixes='#')

    config.read(file)
    logging.info("Reading  configuration file")
    
    GlobalThres = config["Parameters"].getfloat("globalthres")
    LocalThres = config["Parameters"].getfloat("localthres")
    PeakThres = config["Parameters"].getfloat("peakthres")
    GlobalFDR_colum_name = config["Parameters"].get("GlobalFDR_colum_name")
    LocalFDR_colum_name = config["Parameters"].get("LocalFDR_colum_name")
    PeakFDR_colum_name = config["Parameters"].get("PeakFDR_colum_name") 
    Label_colum_name = config["Parameters"].get("Label_colum_name") 
    Decoys_naming = config["Parameters"].get("Decoys_naming")
    
   
    
    logging.info('Porcessing input file')
    filterapply (infile,GlobalThres, LocalThres, PeakThres, GlobalFDR_colum_name,PeakFDR_colum_name,LocalFDR_colum_name,Label_colum_name,Decoys_naming)

    logging.info('end script')




if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='TableFilterer',
        epilog='''
        Example:
            python FDRFilterer.py
        ''')
      
    # default TableFilterer configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/FDRFilterer.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='Input file with paths of the file(s) wanted to be filtered')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
 
    args = parser.parse_args()
   
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
   
        
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + '_log.txt'
    log_file_debug = outfile = args.infile[:-4] + '_log_debug.txt'
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

    infile=args.infile
    #start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))


    main(args.config,args.infile)
  


