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
def readInfile(infile,Dm_colum_name):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile, sep="\t", float_precision='high',low_memory=False,dtype={Dm_colum_name:str}, na_values=[" "])
    df[Dm_colum_name] = df[Dm_colum_name].astype(float)

    return df


def StickerSolver( Theo_mh,seq,DM_user_selection,Error,dfCondFile,listConditions,listLabels,listInfile):
    
    """
    StickerSolver function returns the sequence and label, that match with the smallest error
    according to the dfCondFile and according to listConditions.
    Input:
         Theo_mh             > Theoretical mh  
         seq                 > Sequence with DM
         Error               > Error in ppms stablished by the user
         dfCondFile          > Dataframe that contains all conditions
         listConditios       > Conditions columns names list
         listLabels          > Labels list
         listInfile          > List of infile conditions columns names

         
         
    The function returns:
        StickerLabel             > Label option that best fits
        StickerLable_2           > Label description if it exist
        StickerLable_ppm         > Error (ppm) of label selection
        
    
    """

    
    #These variables are assigned for the cases in which there is no match.

    StickerLabel = ""
    StickerLabel_ppm = ""
    StickerLabel_description = ""
    selected_Label2 = ""
    
    seq = seq[:seq.find("[")]+seq[seq.find("]")+1:]  # Clean sequence is obtained.
    
    dic ={} # A dictionary is created with the aim of save all possible options              
    minimun = 2*10**10 
    n = 0
    matches = 0
    for index, row in dfCondFile.iterrows():
        if pd.isnull(row[0]):
            
            final_label = row[listLabels[0]]
        else:
            DiffAbs = abs(row[0]-DM_user_selection)
            if DiffAbs <= minimun: 
                selected_conditions = []
                minimun = DiffAbs
                selected_DM = row[0]
                selected_Label = row[listLabels[0]] 
                if len(listLabels) > 1: # If there is more than one label
                    selected_Label2 = row[listLabels[1]] 
                for element in listConditions:
                    selected_conditions.append(row[element])

                DiffPPM = abs(((selected_DM-DM_user_selection)*1000000)/(Theo_mh+selected_DM))
                if DiffPPM  < Error: # If is lower than the error this label optio is saved 
                    n = n+1
                    dic[n] = selected_conditions,selected_DM,selected_Label,selected_Label2,DiffPPM

    for key in dic:
        if StickerLabel == "":

            c = 0
            matches = 0
            for element in listConditions: # Each conditions is checked 
                if listInfile[c] == dic[key][0][c] or  pd.isnull(dic[key][0][c]):
                    matches = matches+1
                c = c+1

            if len(listConditions)==matches:# If all the conditions are met
                StickerLabel = dic[key][2]
                StickerLabel_ppm = dic[key][4]
                StickerLabel_description = dic[key][3]

   
    if StickerLabel == "":
        StickerLabel = final_label
    return StickerLabel,StickerLabel_ppm, StickerLabel_description

##################
# Main functions #
##################

def main(solverconfig,infile,CondFile):
    
    """
    Reading configuration file and processing file
    """
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(solverconfig) # Reading configuration file
    
    logging.info("Reading Sticker configuration file")
    
    Error = config["Sticker_Parameters"].getfloat("Relative_Error_ppm") # Relative error (ppm)
    Theo_column_name = config["Sticker_Parameters"].get("theo_mh_column_name")# Theoretical mh column name
    Seq_column_name = config["Sticker_Parameters"].get("Sequence_column_name") # Sequence with DM column
    Selected_DM_column_name = config["Sticker_Parameters"].get("Selected_DM_column_name")# DM column
    StickerLabel_User_output_column_name = config["Sticker_Parameters"].get("StickerLabel_User_output_column_name") # Column name of the output where the chosen label  is annotated  
    StickerLabel_ppm_User_output_column_name =  config["Sticker_Parameters"].get("StickerLabel_ppm_User_output_column_name")# Column name of the output where the calculated error in ppm for the selected label  is annotated 
    output_file_suffix = config["Sticker_Parameters"].get("output_file_suffix") # Chosen suffix for output file
    
   
    dfCondFile=readInfile(CondFile,"d") # User ConditionFile
    headerdfCondFile = list(dfCondFile)
    listConditions= [] # List to save all conditions 
    listLabels= [] # List to save all labels
    for element in headerdfCondFile: 
        if element[0]=="C":
            listConditions.append(element)
            
        elif element[0]=="L": 
            listLabels.append(element)
         
    
    # Input file is read as data frame and column names are added to the header
    dfInfile = readInfile(infile,Selected_DM_column_name)
    dfInfile[StickerLabel_User_output_column_name] = ""
    dfInfile[StickerLabel_ppm_User_output_column_name] = np.nan
    
    if len(listLabels)>1: # If there is more than one label
        dfInfile.loc[StickerLabel_User_output_column_name+"_2"] = np.nan

    cont = 0
    logging.info("Processing input file")
    for index, row in dfInfile.iterrows():
        if pd.notnull(row["p"]):

            listInfile = []
            for element in listConditions: 
                listInfile.append(row[element[element.find("_")+1:]])

     
            StickerLabel,StickerLabel_ppm, StickerLabel_description= StickerSolver(row[Theo_column_name],row[Seq_column_name],row[Selected_DM_column_name],Error,dfCondFile,listConditions,listLabels,listInfile)

            # New columns are completed  
            dfInfile.loc[cont,StickerLabel_User_output_column_name] = StickerLabel
            dfInfile.loc[cont,StickerLabel_ppm_User_output_column_name] = StickerLabel_ppm
            if len(listLabels)>1:# If there is more than one label
                dfInfile.loc[cont,StickerLabel_User_output_column_name+"_2"] = StickerLabel_description
        cont=cont+1
    
    # writting outputfile
    logging.info("Writing output file")

    outfilename = infile[:-4]+output_file_suffix+".txt"

    
    dfInfile.to_csv(outfilename, index=False, sep='\t', encoding='utf-8')


    logging.info('end script')
   
   


if __name__ == '__main__':
    
    try:
        remove('Solver.ini')
    except:
        None

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Sticker',
        epilog='''
        Example:
            python Sticker.py
        ''')
      
    # Default DM0Solver configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    parser.add_argument('-i', '--infile', required=True, help='File wanted to be labelled')
    parser.add_argument('-d', '--conditionsfile', required=True, help='Path to Unimod  file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    
    # These will overwrite the config if specified
    parser.add_argument('-r', '--relerror', help='Maximum allowable relative error (ppm)')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
   
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    
    if args.relerror is not None:
        config.set('DM0Solver_Parameters', 'Relative_Error', str(args.relerror))
        config.set('Logging', 'create_ini', '1')
   
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/Solver.ini', 'w') as newconfig:
            config.write(newconfig)
        
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + '_Sticker_log.txt'
    log_file_debug = outfile = args.infile[:-4] + '_Sticker_log_debug.txt'
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



    #start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    
    # Configuration files are read   
    try:
        open('Solver.ini',"r")
        solverini ='Solver.ini'
        logging.info("Modified Solver configuration file is going to be use")


    except:
        open("config/Solver.ini","r")
        solverini = "config/Solver.ini"
    main(solverini, args.infile,args.conditionsfile)