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

# Imports
import pandas as pd
from os import remove
import numpy as np
from pandas import ExcelWriter
from optparse import OptionParser
import configparser
import argparse
import os
import logging
from pathlib import Path
import sys
import re
import operator




def readInfile(infile,Dm_column_name):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile, sep="\t", float_precision='high',low_memory=False,dtype={Dm_column_name:str})
    df[Dm_column_name].astype("float64").dtypes
    return df

def DMSelection (DMPDM,Theo_mh,Error,userlist_unique): 
    """
    Detecta la DM más cercana de la PDMTable a la introducida por el usuraio  
    """
    # DM from PSL that gives rise to the smallest error is selected 
    minimum = 2*10**10
    finalkey =""
    for DMuser in userlist_unique:
        DiffPPM =abs(((float(DMuser)-float(DMPDM))*1000000)/(float(Theo_mh)+float(DMuser)))

        if DiffPPM <= minimum and DiffPPM<=Error:
            minimum = DiffPPM
            finalkey = DMuser
        
    return finalkey




def applyGM (row1,dfUserFile,Error,Theo_mh_column_name,DM_PDM_column_name,userlist_unique,group_name):
    cond = 0 
    e = "YES"
    final_group_name = " "
    try: 
        row1[Theo_mh_column_name]
    except: 
        e = "NO"
    if e != "NO":

        initial_list = list(dfUserFile)
        number_of_conditions = len(initial_list)-1
        listconditions = initial_list[:-1]
        DM = "NO"
        if DM_PDM_column_name in listconditions:
            finalkey = DMSelection(row1[DM_PDM_column_name],row1[Theo_mh_column_name],Error,userlist_unique)
            listconditions.remove(DM_PDM_column_name)
            DM = "YES"
  
        for index, row in dfUserFile.iterrows():
            cond = 0 
            if DM =="YES" : 

                if row[DM_PDM_column_name] == finalkey: 
                    cond = cond+1
                    for e in listconditions: 
                        if row1[e]==row[e]:
                            cond = cond+1

            elif DM == "NO": 
                for e in listconditions: 
                    if row1[e]==row[e]:
                        cond = cond+1


            if cond == number_of_conditions: 
                final_group_name=row[group_name]
        row1[group_name]=final_group_name

        return row1
                        
        

##################
# Main functions #
##################

def main(configurationfile,PDMTable,UserFile,):
    """
    Reading configuration file
    """
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(configurationfile)

    logging.info("Reading configuration file")
    Error = config["GroupMaker_Parameters"].getfloat("Relative_Error_ppm") # Relative error is saved     
    Theo_mh_column_name = config["GroupMaker_Parameters"].get("theo_mh_column_name") # Theoretical mh name
    DM_PDM_column_name = config["GroupMaker_Parameters"].get("DM_PDM_column_name") # DM of the PDMTable file column name
    output_file_suffix = config["GroupMaker_Parameters"].get("output_file_suffix") # DM of the PDMTable file column name
    
    dfPDM = readInfile(PDMTable,DM_PDM_column_name)
    dfUserFile =readInfile(UserFile,DM_PDM_column_name)
    group_name =(list(dfUserFile.columns.values)[-1]) #Cogemos el nombre del grupo que siempre estará en último lugar en el fichero
    userlist_unique = list(dict.fromkeys(dfUserFile[DM_PDM_column_name].tolist()))
    dfPDM[group_name]=""
    logging.info("Processing input file")
    

    dfPDM = dfPDM.apply(lambda y:applyGM(y,dfUserFile,float(Error),Theo_mh_column_name,DM_PDM_column_name,userlist_unique,group_name), axis = 1)


    # write outputfile
    logging.info("Writing output file")

    outfilename = PDMTable[:-4]
    outfile= outfilename+output_file_suffix+".txt"
    dfPDM.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
    logging.info('end script')


if __name__ == '__main__':


    # parse arguments
    parser = argparse.ArgumentParser(
        description='TrunkSolver',
        epilog='''
        Example:
            python TrunkSolver.py
        ''')
      
    # default TrunkSolver configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    defaultconfigM = os.path.join(os.path.dirname(__file__), "config/MassMod.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='Path to input file(PDMTable)')
    parser.add_argument('-u', '--UserFile', required=True, help='Path to input UserFile')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')

    
    # these will overwrite the config if specified
    parser.add_argument('-r', '--relerror', help='Maximum allowable relative error (ppm)')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
   
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.relerror is not None:
        config.set('TrunkSolver_Parameters', 'relative_error', str(args.relerror))
        config.set('Logging', 'create_ini', '1')
   
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/config/Solver.ini', 'w') as newconfig:
            config.write(newconfig)
        
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + 'TrunkSolver_log.txt'
    log_file_debug = outfile = args.infile[:-4] + 'TrunkSolver_log_debug.txt'
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

    PDMTable=args.infile
    UserFile=args.UserFile
    
    #start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))

    


    main(args.config,PDMTable,UserFile)
