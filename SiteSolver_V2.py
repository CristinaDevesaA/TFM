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

import pandas as pd
import numpy as np
from optparse import OptionParser
import configparser
import argparse
import os
import logging
from pathlib import Path
import sys
import re





def readInfile(infile,Dm_column_name):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile, sep="\t", float_precision='high',low_memory=False,dtype={Dm_column_name:str})
    df[Dm_column_name].astype("float64").dtypes
    return df





def breakUp1(string,DM):
    """
    breakUp1 function extract the mass modifications, first and second half of te sequence (being the separation 
    the mass modification), their length, and the residue (amino acid found on the left attached to the mass
    modification). 
    """
    start = string.find("[")
    end = string.find("]")
    massMod = round(float(DM),6)
    DM=float(string[(start)+1: (end)])
    secondHalf = string[end+1:]

    
    if start != 0:
        fisrtHalf = string[0:start]
        fisrtHalf = re.sub('\(\-\D+\)', '', fisrtHalf)
        residue = fisrtHalf[-1]
    elif start == 0: 
        residue = "NA"
        fisrtHalf = "0"
    elif end == (len(string)-1):
        secondHalf = "0"
 
    
    lenseconfHalf=len(secondHalf)
    lenfirstHalf=len(fisrtHalf)
    
    return(fisrtHalf, secondHalf, massMod, residue, lenfirstHalf, lenseconfHalf)





def applySS(row,Error,Theo_column_name,seq_column_name,x,cal_Dm_mh_column_name,SiteSequence_column_name,SiteCorrection_column_name,SiteDM_column_name,SiteDMError_ppm_column_name,dicPL):

    minimum = Error
    final = ""
    SiteSequence=row[seq_column_name]
    SiteCorrection=""
    finalkey=""
    SiteDM=""
    SiteDMError_ppm=""
    if row[seq_column_name].find("_")==-1:
        for key in dicPL:
            DiffPPM = abs(((float(key)-float(row[cal_Dm_mh_column_name]))*1000000)/(row[Theo_column_name]+float(key)))
            if DiffPPM <= minimum:
                minimum = DiffPPM
                finalkey = key
        
            
        if finalkey!="" and row[seq_column_name][row[seq_column_name].find("[")-1]  not in dicPL[finalkey].split(","):

            firstHalf, secondHalf, massMod, residue, lenfirstHalf, lenseconfHalf= breakUp1(row[seq_column_name],row[cal_Dm_mh_column_name])
            for i in range(x): # All positions are traversed starting from the positions closest to the  original DM site
                if final == "": # If a site has not yet been found

                    try:
                        # It is verified that in that position there exists and amino acid 
                        aasecond = secondHalf[i] # First amino acid starting from the second part
                        aafirst = firstHalf[-(i+2)] # Last amino acid starting from the first part
                        if aasecond in dicPL[finalkey].split(",") and aafirst in dicPL[finalkey].split(","):
                            if str(dicPL[finalkey]).find(aasecond)<str(dicPL[finalkey]).find(aafirst):
                                secondHalf2 = list(secondHalf)
                                secondHalf2[i] = secondHalf2[i]+"["+str(massMod)+"]"
                                secondHalf2 = "".join(secondHalf2)
                                final = firstHalf+secondHalf2
                            elif str(dicPL[finalkey]).find(aasecond)>str(dicPL[finalkey]).find(aafirst):
                                firstHalf2 = list(firstHalf)
                                firstHalf2[-(i+2)] = firstHalf2[-(i+2)]+"["+str(massMod)+"]"                                   
                                firstHalf2 = "".join(firstHalf2)                                   
                                final = firstHalf2+secondHalf
                                
                        else:
                            aasecond = secondHalf[-55] 
                    except:                
                        try: 
                            aasecond = secondHalf[i] # First amino acid starting from the second part
                            
                            if aasecond in dicPL[finalkey].split(","):      
                                secondHalf2 = list(secondHalf)
                                secondHalf2[i] = secondHalf2[i]+"["+str(massMod)+"]"
                                secondHalf2 = "".join(secondHalf2)
                                final = firstHalf+secondHalf2
                            else:
                                aasecond = secondHalf[-55]                     
                        except:                            
                            try:

                                aafirst = firstHalf[-(i+2)] # Last amino acid starting from the first part
                                if aafirst in  dicPL[finalkey].split(","): 
                                    firstHalf2 = list(firstHalf)
                                    firstHalf2[-(i+2)] = firstHalf2[-(i+2)]+"["+str(massMod)+"]"                                   
                                    firstHalf2 = "".join(firstHalf2)                                   
                                    final = firstHalf2+secondHalf
                                    
                                else:
                                    aasecond = secondHalf[-55] 
                            except:
                                continue
            

            if final!="":
                SiteSequence=final
                SiteCorrection=row[seq_column_name][row[seq_column_name].find("[")-1] +">"+ final[final.find("[")-1]
                SiteDM=float(finalkey)
                SiteDMError_ppm=DiffPPM
                      

    row[SiteSequence_column_name]=SiteSequence
    row[SiteCorrection_column_name] = SiteCorrection
    row[SiteDM_column_name] = SiteDM
    row[SiteDMError_ppm_column_name] = SiteDMError_ppm
    return row




##################
# Main functions #
##################


def main(file,infile,PrimaryList):
    
    """
    Reading configuration file
    """


    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(file)

    logging.info("Reading SiteSolver configuration file") 
    
    Error = config["SiteSolver_Parameters"].getfloat("Relative_Error_ppm") # Relative error (ppm)
    Theo_column_name = config["SiteSolver_Parameters"].get("theo_mh_column_name")# Theoretical mh column name
    seq_column_name = config["SiteSolver_Parameters"].get("Sequence_column_name") # Sequence with DM column name
    Output_file_suffix = config["SiteSolver_Parameters"].get("Output_file_suffix") # Chosen suffix for output file
    x = config["SiteSolver_Parameters"].getint("x") # Parameter that indicates the extension (left and right) of the amino acids wanted to be analyzed 
    cal_Dm_mh_column_name =  config["SiteSolver_Parameters"].get("cal_Dm_mh_column_name")# Calibrated DM MH name
    
    SiteSequence_column_name = config["SiteSolver_Parameters"].get("SiteSequence_column_name") # Column name of the output where the sequence is annotated 
    SiteCorrection_column_name = config["SiteSolver_Parameters"].get("SiteCorrection_column_name")# Column name of the output where correction site is annotated # Column name of the output where correction site is annotated 
    SiteDM_column_name = config["SiteSolver_Parameters"].get("SiteDM_column_name") # Column name of the output where selected DM is annotated 
    SiteDMError_ppm_column_name = config["SiteSolver_Parameters"].get("SiteDMError_ppm_column_name") # Column name of the output where the error of the selected DM is annotated 
     
    
    logging.info("Reading input files")   
    dfPL = readInfile(PrimaryList,"DM")
    dicPL = dfPL.set_index('DM').to_dict()['Residue']
    df =readInfile(infile,cal_Dm_mh_column_name)

    df[SiteSequence_column_name]=""
    df[SiteCorrection_column_name] = ""
    df[SiteDM_column_name] = ""
    df[SiteDMError_ppm_column_name] = ""

    logging.info("Processing input file")
    df = df.apply(lambda y: applySS(y,Error,Theo_column_name,seq_column_name,x,cal_Dm_mh_column_name,SiteSequence_column_name,SiteCorrection_column_name,SiteDM_column_name,SiteDMError_ppm_column_name,dicPL), axis = 1)
   
    # write outputfile
    logging.info("Writing output file")

    outfilename = infile[:-4]
    outfile= outfilename+Output_file_suffix+".txt"
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')


    logging.info('end script')



if __name__ == '__main__':



    # parse arguments
    parser = argparse.ArgumentParser(
        description='SiteSolver',
        epilog='''
        Example:
            python SiteSolverr.py
        ''')
      
    # default DM0Solver configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-pl', '--SiteList', required=True, help='Path to sitelist file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    
    # these will overwrite the config if specified
    parser.add_argument('-r', '--relerror', help='Maximum allowable relative error (ppm)')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
   
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.relerror is not None:
        config.set('SiteSolver_Parameters', 'Relative_Error_ppm', str(args.relerror))
        config.set('Logging', 'create_ini', '1')
   
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/config/Solver.ini', 'w') as newconfig:
            config.write(newconfig)
        
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + 'SiteSolved_log.txt'
    log_file_debug = outfile = args.infile[:-4] + 'SiteListSolved_log_debug.txt'
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

        
    main(args.config, args.infile, args.SiteList)

