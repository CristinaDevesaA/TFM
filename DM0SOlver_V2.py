
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





def proxim(apex_file,dic_DM0,Error): 
    dic_selected_DMapex ={}
    for Label in dic_DM0:
        selected_Dmapex = ""
        DMlist= dic_DM0[Label]
        
        minimun_DiffAbs = Error  # At the beginning the minimum difference (ppms) is considered as the Error
        for line in open(apex_file,"r"): 
            DMapex =float(line.strip("\n").strip(" "))
            DiffAbs = abs(DMlist-DMapex)
            if DiffAbs <= Error:
                if DiffAbs < minimun_DiffAbs:
                    minimun_DiffAbs = DiffAbs
                    selected_DMapex = DMapex
                    
        if selected_DMapex!="":
            if round(selected_DMapex,6) in dic_selected_DMapex.keys():
                if dic_selected_DMapex[round(selected_DMapex,6)][1]>minimun_DiffAbs:
                    dic_selected_DMapex[round(selected_DMapex,6)]=Label,minimun_DiffAbs
            else:

                dic_selected_DMapex[round(selected_DMapex,6)]=Label,minimun_DiffAbs
    return dic_selected_DMapex
    



###################
# Local functions #
###################

def readInfile(infile,DM_column_name):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile, sep="\t",                               
    float_precision='high',low_memory=False,dtype={DM_column_name:str})
    df[DM_column_name].astype("float64").dtypes
    return df

def DM0Solver(dic_apex_selected,DM,seq,Error,dic_DM0,SequenceMassmod,sequence,PeakAssignation,PeakNaming):
    
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
    DM0Label_error = ""

    if PeakAssignation ==PeakNaming:
        DM2 =round(float(DM),6)

        find="NO"
        for key in dic_apex_selected.keys():
            if DM2 ==key:
 
                DM0Sequence = seq+"_"+str(SequenceMassmod)
                DM0Label = dic_apex_selected[DM2][0]
                DM0Label_error = dic_apex_selected[DM2][1]
        



   
    return DM0Sequence,DM0Label,DM0Label_error






def applySolver(row, Seq_column_name, DM_column_name,dic_apex_selected, Error, dic_DM0,
                DM0Sequence_output_column_name, DM0Label_output_column_name, DM0Label_error_output_column_name,PeakAssignation_Column_name,PeakNaming):
    seq = row[Seq_column_name]
    SequenceMassmod = seq[seq.find("[")+1:]
    SequenceMassmod = SequenceMassmod[:SequenceMassmod.find("]")]  # Mass modification of the sequence is obtained.

    if row[Seq_column_name].find("_") != -1:  # If  scan has already been corrected DM0Seqeuence remains Seq_column_name
        seq = seq[:seq.find("_")]
        DM0Sequence=row[Seq_column_name]
        DM0Label=""
        DM0Label_error =""
        
    else:

        seq = seq[:seq.find("[")]+seq[seq.find("]")+1:]  # Clean sequence is obtained.# if  scan has not yet been corrected DM0Sover function is executed

        DM0Sequence,DM0Label,DM0Label_error = DM0Solver(dic_apex_selected,row[DM_column_name],seq,Error,dic_DM0,SequenceMassmod,row[Seq_column_name],row[PeakAssignation_Column_name],PeakNaming)
   
        
    # New columns are completed  
    row[DM0Sequence_output_column_name] = DM0Sequence
    row[DM0Label_output_column_name] = DM0Label.upper()
    row[DM0Label_error_output_column_name] = DM0Label_error
    return row





##################
# Main functions #
##################

def main(solverconfig,infile,apex_file):
    
    """
    Reading configuration file and processing file

    """

    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(solverconfig) # Reading configuration file
    
    logging.info("Reading DM0Solver configuration file")
    

    DM_column_name = config["DM0Solver_Parameters"].get("DM_column_name") # Experimental mh  column name
    Seq_column_name = config["DM0Solver_Parameters"].get("Sequence_column_name") # Sequence column name

    DM0Sequence_output_column_name =  config["DM0Solver_Parameters"].get("DM0Sequence_output_column_name") # Column name of the output where the chosen sequence is annotated
    DM0Label_output_column_name =  config["DM0Solver_Parameters"].get("DM0Label_output_column_name") # Column name of the output where the chosen label is annotated
    DM0Label_error_output_column_name =  config["DM0Solver_Parameters"].get("DM0Label_error_output_column_name") # Column name of the output where the calculated error in ppm is annotated
    output_file_suffix = config["DM0Solver_Parameters"].get("output_file_suffix") # Chosen suffix for output file 
    Error = config["DM0Solver_Parameters"].getfloat("Absolute_Error") # Relative error (ppm)
    PeakAssignation_Column_name = config["DM0Solver_Parameters"].get("PeakAssignation_Column_name") # Relative error (ppm)
    PeakNaming = config["DM0Solver_Parameters"].get("PeakNaming") # Relative error (ppm)    


 
    # A dictionary that save the labels and their corresponding masses is created
    dic_DM0 = {}   
    for option, value in config["DM0Solver_DM0List"].items():
        dic_DM0[option.strip("\n")] = float(value.strip("\n"))
            

    dic_apex_selected =  proxim(apex_file,dic_DM0,Error)


    # Input file is read as data frame and column names are added to the header
    df = readInfile(infile,DM_column_name)
    try: 
        df.drop([DM0Sequence_output_column_name,DM0Label_output_column_name,DM0Label_output_column_name],axis=1)
    except:
        pass
    df[DM0Sequence_output_column_name] = ""
    df[DM0Label_output_column_name] = ""
    df[DM0Label_error_output_column_name] = np.nan


    #cont = 0
    logging.info("Processing input file")
    df = df.apply(lambda x: applySolver(x, Seq_column_name, DM_column_name,dic_apex_selected, Error, dic_DM0,
                                        DM0Sequence_output_column_name, DM0Label_output_column_name, DM0Label_error_output_column_name,PeakAssignation_Column_name,PeakNaming), axis = 1)

    # write outputfile
    logging.info("Writing output file")

    outfilename = infile[:-4]+output_file_suffix
    outfile= outfilename+".txt"
    try:
        remove("outfile")
    except:
        None
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')

    logging.info('end script')







def truncate(num, n):
    integer = int(num * (10**n))/(10**n)
    return float(integer)



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
    parser.add_argument('-a', '--apexfile', default=defaultconfig, help='Path toapex file')    
    
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
    

    main(args.config, infile1,args.apexfile)

