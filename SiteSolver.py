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

###################
# Local functions #
###################


def readInfile(infile,Dm_column_name):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile, sep="\t", float_precision='high',low_memory=False,dtype={Dm_column_name:str})
    df[Dm_column_name].astype("float64").dtypes
    return df



def readSiteList(SiteListfile ,SiteList_column_name,SiteCorrection_label):
    """
    readSiteList function extracts all the information of SiteList file
    """
 
    # A dicctionary which save all SiteList content is created
    dicc_SiteList = {}
    dicc_SiteList_list_label = {}
    filecount = 0

    #SiteList file is read
    SiteListfile = open(SiteListfile,"r")
    c = 0
    for line in SiteListfile:
        c = c+1
        if c == 1:
            header = line.strip("\n").split("\t")     
        else:
            fields = line.strip("\n").split("\t")
            DM = float(fields[0].replace(",","."))
            aa = fields[1]
            
            # N-terminad la C-terimnal, is changed to ^ and - respectively 
            if aa == "NT":
                aa = "^"
            elif aa == "CT":
                aa = "-"
            freq = float(fields[header.index(SiteList_column_name)].replace(",","."))

            if freq > 0 : 
                if DM not in dicc_SiteList.keys():
                    dicc_SiteList[DM] = {}
                    dicc_SiteList[DM][aa] = freq
                else:
                    dicc_SiteList[DM][aa] = freq
                    
                if DM not in dicc_SiteList_list_label.keys():
                    dicc_SiteList_list_label[DM] = {}
                    dicc_SiteList_list_label[DM][aa] = SiteCorrection_label
                else:    
                    dicc_SiteList_list_label[DM][aa] = SiteCorrection_label
 
    
    return dicc_SiteList, dicc_SiteList_list_label



def readPeptideSiteList(PeptideSiteList,MinScanFreq,MaxNSite):
    """
    Extracting all the information of PeptideSiteList file
    """

    df_PSL = pd.read_csv(PeptideSiteList, sep="\t", float_precision='high',low_memory=False) # PepetideSiteList is read as dataframe
    dicc_PSL = {}   # A dicctionary which save all PeptideSiteList content is created
    for index, row in df_PSL.iterrows():
        if row["pdm"].find("_") == -1:
            if row["d"] not in dicc_PSL.keys():
                dicc_PSL[row["d"]] = {}
            if row["p"] not in dicc_PSL[row["d"]].keys():
                dicc_PSL[row["d"]][row["p"]] = {}
            if row["m"] not in dicc_PSL[row["d"]][row["p"]].keys() and row["ScanFreq"] >= MinScanFreq : # Those that do not meet MinScanFreq will be discarded
                dicc_PSL[row["d"]][row["p"]][int(row["m"])] = {}
                if  row["a"] not in dicc_PSL[row["d"]][row["p"]][int(row["m"])].keys() :
                    dicc_PSL[row["d"]][row["p"]][int(row["m"])][row["a"]] = row["ScanFreq"]
    # Only those that meet MaxNSite condition are saved
    dicc_del = {}
    for d in dicc_PSL.keys():
        dicc_del[d] = {}
        for p in dicc_PSL[d].keys():
            dicc_del[d][p]={}
            subdicc= {}
            lista = []
            for m in dicc_PSL[d][p].keys():
                for a in dicc_PSL[d][p][m].keys():
                    subdicc[a+"-"+str(m)]=dicc_PSL[d][p][m][a]

            sort = sorted((subdicc).items(), key=operator.itemgetter(1),reverse= True)
            
            dicc_del[d][p]=sort[:MaxNSite]


    dicc_PSL2 = {}
    for d in dicc_PSL.keys():
        dicc_PSL2[d] = {}
        for p in dicc_PSL[d].keys():
            dicc_PSL2[d][p] = {}
            for position in  dicc_PSL[d][p].keys():
                dicc_PSL2[d][p][position] = {}
                for aa in dicc_PSL[d][p][position].keys():
                    for item in dicc_del[d][p]:
                        if aa == item[0].split("-")[0] and  float(position) ==  float(item[0].split("-")[1]) and dicc_PSL[d][p][position][aa] == item[1]:
                            dicc_PSL2[d][p][position][aa] = {}
                            dicc_PSL2[d][p][position][aa] = dicc_PSL[d][p][position][aa]


    return dicc_PSL2




def PSL_analysis(aasecond,aafirst,firstHalf,secondHalf,dicc_PSL,key1PSL,dicc_SiteList,key1):
     
    """
    PSL_analysis function will return the amino acid with the highest frequency according to the user's input files.
    If the frequencies are equal it will take into account the order of appearance in the input files.
    """
    if dicc_PSL[key1PSL][firstHalf+secondHalf][aafirst][0] > dicc_PSL[key1PSL][firstHalf+secondHalf][aasecond][0]:
        best = "aafirst"
        add = ""
    elif dicc_PSL[key1PSL][firstHalf+secondHalf][aafirst][0] < dicc_PSL[key1PSL][firstHalf+secondHalf][aasecond][0]:
        best = "aasecond"
        add = ""
    elif dicc_PSL[key1PSL][firstHalf+secondHalf][aafirst][0] == dicc_PSL[key1PSL][firstHalf+secondHalf][aasecond][0]:
        best = "tie"
        add = aasecond+"/"+aafirst
    
    if best == "tie": 
        lista = list(dicc_SiteList[key1].keys())
        if lista.index(aafirst)>lista.index(aasecond):
            best = "aasecond"
            add = ""
        elif lista.index(aafirst)<lista.index(aasecond)  :
            best = "aafirst"
            add = ""
    return add, best 


def breakUp1(string,cal_DM):
    """
    breakUp1 function extract the mass modifications, first and second half of te sequence (being the separation 
    the mass modification), their length, and the residue (amino acid found on the left attached to the mass
    modification). 
    """
    start = string.find("[")
    end = string.find("]")
    massMod = float(cal_DM)
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
    
    return(fisrtHalf, secondHalf, massMod, residue, lenfirstHalf, lenseconfHalf,DM)



def SimpleSiteSolver(string,cal_dm_mh,dicc_SiteList,dicc_PSL,x,Theo_mh,Error,DM_user_selection):
    """
    SimpleSiteSolver function will check the best option of the site based on PPDMTable and SiteList and will return,
    if there is a better option than the initial one, the sequence with the DM located in the best place.
    """
 
    firstHalf, secondHalf, massMod, residue, lenfirstHalf, lenseconfHalf,DM=breakUp1(string,cal_dm_mh)
 
    # Initial parameters 
    x = x*2
    final = ""
    SiteSolverOptions = " " 
    lastI = ""
    minimum = 2*10**10
    minimumPSL = 2*10**10
    key1 = ""
    key1PSL = ""
    finalkey = ""
    finalkeyPSL = ""
    
    # DM from SL that gives rise to the smallest error is selected 
    for key in dicc_SiteList:
        a = ""
        DiffPPM = abs(((float(key)-DM_user_selection)*1000000)/(Theo_mh+float(key)))
        if DiffPPM <= minimum:
            minimum = DiffPPM
            key1 = key
            finalkey = key1

    # DM from PSL that gives rise to the smallest error is selected 
    for keyPSL in dicc_PSL:
        a = ""
        DiffPPM_PSL =abs(((float(keyPSL)-DM_user_selection)*1000000)/(Theo_mh+float(keyPSL)))
        if DiffPPM_PSL <= minimumPSL:
            minimumPSL = DiffPPM_PSL
            key1PSL = keyPSL
            finalkeyPSL = key1PSL
   
    try:
        a = dicc_SiteList[key1][residue] # If this residue (amino acid in zero position) is an allowed amino acid (appears in SiteList dictionary)

        aPLS = dicc_PSL[key1PSL][firstHalf+secondHalf][len(firstHalf)][residue] # If this residue (amino acid in zero position) is an allowed amino acid (appears in PeptideSiteList dictionary)
        
        final = " " # the site will not change, so the final variable will be empty since no ca will be applied
    except:
   
        if key1!="" and key1PSL!="": # If any of the DMs have  crossed the threshold in both lists (SL, PSL)
            if secondHalf=="" and firstHalf=="":
               
                if "-" in  dicc_SiteList[key1].keys() and   firstHalf[-1] not in dicc_SiteList[key1].keys():
                    searchsecondHalf =  "_"
                else:
                    searchsecondHalf= secondHalf
              
             
                if "^" in  dicc_SiteList[key1].keys() and   secondHalf[0] not in dicc_SiteList[key1].keys():
                    searchfirstHalf =  "^"
                else:
                   searchfirstHalf= firstHalf
             
            elif secondHalf=="":
                if "-" in  dicc_SiteList[key1].keys() and   firstHalf[-1] not in dicc_SiteList[key1].keys():
                    searchsecondHalf =  "_"
                else:
                        searchsecondHalf= secondHalf
                if "^" in  dicc_SiteList[key1].keys() and  firstHalf[0] not in dicc_SiteList[key1].keys():
                    searchfirstHalf = firstHalf.replace(firstHalf[0],"^")
                else:
                    searchfirstHalf=firstHalf
            elif firstHalf=="":
                if "^" in  dicc_SiteList[key1].keys() and   secondHalf[0] not in dicc_SiteList[key1].keys():
                    searchfirstHalf =  "^"
                else:
                   searchfirstHalf= firstHalf
                if "-" in  dicc_SiteList[key1].keys() and   secondHalf[-1] not in dicc_SiteList[key1].keys():
                    searchsecondHalf =  secondHalf.replace( secondHalf[-1],"-")
                else:
                    searchsecondHalf= secondHalf

               
            else: 
                if "^" in  dicc_SiteList[key1].keys() and  firstHalf[0] not in dicc_SiteList[key1].keys():
                    searchfirstHalf = firstHalf.replace(firstHalf[0],"^")
                else:
                    searchfirstHalf=firstHalf
                if "-" in  dicc_SiteList[key1].keys() and   secondHalf[-1] not in dicc_SiteList[key1].keys():
                    searchsecondHalf =  secondHalf.replace( secondHalf[-1],"-")
                else:
                    searchsecondHalf= secondHalf


            for i in range(x): # All positions are traversed starting from the positions closest to the  original DM site
                if final == "": # If a site has not yet been found
                    try:
                        # It is verified that in that position there exists and amino acid 
                        aasecond = searchsecondHalf[i] # First amino acid starting from the second part
                        aafirst = searchfirstHalf[-(i+2)] # Last amino acid starting from the first part
                        
                        # It is checked if the amino acid and it position are allowed for that DM
                        a = dicc_SiteList[float(key1)][aasecond]
                        b = dicc_SiteList[float(key1)][aafirst]
                        aPSL = dicc_PSL[key1PSL][firstHalf+secondHalf][len(firstHalf)+i+1][secondHalf[i]]
                        bPSL = dicc_PSL[key1PSL][firstHalf+secondHalf][len(firstHalf)-i-1][firstHalf[-(i+2)]]
                        
                        add,best = PSL_analysis(aasecond,aafirst,firstHalf,secondHalf,dicc_PSL,key1PSL,dicc_SiteList,key1) # if both parts contain, in that position, an allowed amino acid PSL_analysis will be employed
                        if best != None:# If some result is obtain from PSL_anlysis function it will be determined the best option
                            if best == "aasecond":
                                secondHalf2 = list(secondHalf)
                                secondHalf2[i] = secondHalf2[i]+"["+str(DM)+"]"
                                secondHalf2 = "".join(secondHalf2)
                                final = firstHalf+secondHalf2
                                finalkey=key1
                                lastI=i
                            elif best == "tie":
                                final = " "
                                finalkey=key1
                                SiteSolverOptions = add
                                lastI=i
                            elif best == "aafirst":
                                firstHalf2 = list(firstHalf)
                                firstHalf2[-(i+2)] = firstHalf2[-(i+2)]+"["+str(DM)+"]"

                                firstHalf2 = "".join(firstHalf2)
                                final = firstHalf2+secondHalf
                                lastI=i

                    except: # If both parts do not have an allowed amino acid It will be checked each one , separately
                        try: 
                            aasecond = searchsecondHalf[i] # First amino acid starting from the second part
                            
                            # It is checked if the amino acid and it position are allowed for that DM 
                            a = dicc_SiteList[float(key1)][aasecond] 
                            aPSL = dicc_PSL[key1PSL][firstHalf+secondHalf][len(firstHalf)+i+1][secondHalf[i]]
        
                            secondHalf2 = list(secondHalf)
                            secondHalf2[i] = secondHalf2[i]+"["+str(DM)+"]"
                            secondHalf2 = "".join(secondHalf2)
                            final = firstHalf+secondHalf2
                            lastI=i

                        except:
                            try:
                                aafirst = searchfirstHalf[-(i+2)] # Last amino acid starting from the first part
                                
                                # It is checked if the aminoacid and it position are allowed for that DM 
                                a = dicc_SiteList[float(key1)][aafirst]
                                bPSL = dicc_PSL[key1PSL][firstHalf+secondHalf][len(firstHalf)-i-1][firstHalf[-(i+2)]]
                  
                                firstHalf2 = list(firstHalf)
                                firstHalf2[-(i+2)] = firstHalf2[-(i+2)]+"["+str(DM)+"]"
                                firstHalf2 = "".join(firstHalf2)
                                final = firstHalf2+secondHalf
                                lastI=i

                            except:
                                continue
                
                
    return final,SiteSolverOptions,lastI,finalkey,minimum

##################
# Main functions #
##################


def main(file,infile1,PrimaryList,peptidesitelist,SecondaryList):
    
    """
    Reading configuration file
    """


    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(file)

    logging.info("Reading SiteSolver configuration file") 
    
    Error = config["SiteSolver_Parameters"].getfloat("Relative_Error_ppm") # Relative error (ppm)
    Theo_column_name = config["SiteSolver_Parameters"].get("theo_mh_column_name")# Theoretical mh column name
    seq_column_name = config["SiteSolver_Parameters"].get("Sequence_column_name") # Sequence with DM column name
    PeakNaming = config["SiteSolver_Parameters"].get("PeakNaming") # Parameter that indicates how peaks are named
    PeakAssignation_column_name = config["SiteSolver_Parameters"].get("PeakAssignation_column_name") # Name of column that contains peak assignation
    Output_file_suffix = config["SiteSolver_Parameters"].get("Output_file_suffix") # Chosen suffix for output file
    x = config["SiteSolver_Parameters"].getint("x") # Parameter that indicates the extension (left and right) of the amino acids wanted to be analyzed 
    cal_Dm_mh_column_name =  config["SiteSolver_Parameters"].get("cal_Dm_mh_column_name")# Calibrated DM MH name
    MinScanFreq = config["SiteSolver_Parameters"].getint("MinScanFreq") # Parameter that indicates frequency thershold for PSL
    MaxNSite = config["SiteSolver_Parameters"].getint("MaxNSite") # Parameter that indicates number of amino acid allowed when  PSL is analyzed
     
    SiteSequence_column_name = config["SiteSolver_Parameters"].get("SiteSequence_column_name") # Column name of the output where the sequence is annotated 
    SiteCorrection_column_name = config["SiteSolver_Parameters"].get("SiteCorrection_column_name")# Column name of the output where correction site is annotated # Column name of the output where correction site is annotated 
    SiteOption_column_name = config["SiteSolver_Parameters"].get("SiteOption_column_name") # Column name of the output where possible options are annotated 
    SiteDM_column_name = config["SiteSolver_Parameters"].get("SiteDM_column_name") # Column name of the output where selected DM is annotated 
    SiteDMError_ppm_column_name = config["SiteSolver_Parameters"].get("SiteDMError_ppm_column_name") # Column name of the output where the error of the selected DM is annotated 

       
    PrimaryList_column_name = config["SiteSolver_Parameters"].get("PrimaryList_column_name") # Column from Primary List wanted to be used for SiteSolver
    SecondaryList_column_name = config["SiteSolver_Parameters"].get("SecondaryList_column_name")# Column from Secondary List wanted to be used for SiteSolver
    SiteCorrection_PrimaryList_label = config["SiteSolver_Parameters"].get("SiteCorrection_PrimaryList_label") # Site List selected as a  Primary list
    SiteCorrection_SecondaryList_label = config["SiteSolver_Parameters"].get("SiteCorrection_SecondaryList_label")# Site List selected as a  Secondary list
  

    logging.info("Reading PrimaryList")
    dicc_PrimaryList, dicc_PrimaryList_label = readSiteList(PrimaryList,PrimaryList_column_name, SiteCorrection_PrimaryList_label)
    logging.info("Reading SecondaryLis")
    dicc_SecondaryList, dicc_SecondaryList_label = readSiteList(SecondaryList,SecondaryList_column_name, SiteCorrection_SecondaryList_label)
    logging.info("Reading PeptideSiteList")
    dicc_PSL = readPeptideSiteList(peptidesitelist,MinScanFreq,MaxNSite)

    logging.info("Reading input file")
    df=readInfile(infile1,cal_Dm_mh_column_name)
    df[SiteSequence_column_name]=""
    df[SiteCorrection_column_name] = ""
    df[SiteOption_column_name] = ""
    df[SiteDM_column_name] = ""
    df[SiteDMError_ppm_column_name] = ""

    logging.info("Processing input file") 
    cont=0
    for row in df.iterrows():

        
        # If the scan has been already analyzed
        if df.loc[cont,seq_column_name].find("_")!=-1 or df.loc[cont,PeakAssignation_column_name].upper()!= PeakNaming.upper():
            final= df.loc[cont,seq_column_name]
            SiteSolverOptions = " "
            corrected="NO"
            correction  = ""
            DMSelection = ""
            Error_ppm_selection = ""
           
   

        else :
            # SimpleSiteSolver function is applied using Primary List
            
            final,SiteSolverOptions,position,DMSelected,Selection_Error=SimpleSiteSolver(df.loc[cont,seq_column_name],df.loc[cont,cal_Dm_mh_column_name],dicc_PrimaryList,dicc_PSL,x,df.loc[cont,Theo_column_name],Error,float(df.loc[cont,cal_Dm_mh_column_name]))
            DMSelected1 = DMSelected
            Selection_Error1=Selection_Error
            ListCorrection = "Primary"
            
            # If there is no matching DM at primary list Secondary List will be analyzed
            if Selection_Error1>Error:
                ListCorrection="Secondary"
                final,SiteSolverOptions,position,DMSelected,Selection_Error=SimpleSiteSolver(df.loc[cont,seq_column_name],df.loc[cont,cal_Dm_mh_column_name],dicc_SecondaryList,dicc_PSL,x,df.loc[cont,Theo_column_name],Error,float(df.loc[cont,cal_Dm_mh_column_name]))
                DMSelected2 = DMSelected
                Selection_Error2=Selection_Error
            
       
               
                    
            if final=="" or final==" ": # If no reassignation has been found the DM that has the lowest error will be selected
                corrected="NO"
                final= df.loc[cont,seq_column_name]
                correction = ""
                try:
                    if Selection_Error2 > Selection_Error1:
                        DMSelection = DMSelected1
                        Error_ppm_selection = Selection_Error1
                    else:
                        DMSelection = DMSelected2
                        Error_ppm_selection = Selection_Error2
                except:
                    DMSelection = DMSelected1
                    Error_ppm_selection = Selection_Error1

            else:# If a  resignation has been found NT and CT are annotated and correction variable is completed
                
                initial_position = df.loc[cont,seq_column_name][df.loc[cont,seq_column_name].find("[")-1]
                last_position = final[final.find("[")-1]
                DMSelection = DMSelected
                Error_ppm_selection= Selection_Error
                corrected="YES"
                correction = initial_position+"<"+last_position+";"+str(position+1)
 
                if ListCorrection == "Secondary":
                    dicc_SiteList_list_label=dicc_SecondaryList_label

                elif ListCorrection == "Primary":
                    dicc_SiteList_list_label=dicc_PrimaryList_label

                
                if final[-1]=="]":
                    final_F = final[final.find("[")-1]
                else:
                    final_F = final[-1]
                if final[0] =="[":
                    final_P = final[final.find("]")+1]
                        
                else:
                    final_P = final[0]
  
                try:
                    correction = dicc_SiteList_list_label[DMSelection][last_position]+correction 
                except:
                    if final_F==last_position and final_P==last_position:
                        try:
                            correction = dicc_SiteList_list_label[DMSelection]["-"]+correction
                        except:
                            correction = dicc_SiteList_list_label[DMSelection]["^"]+correction
                    elif final_F==last_position:
                        correction = dicc_SiteList_list_label[DMSelection]["-"]+correction
                    elif final_P==last_position:
                        correction = dicc_SiteList_list_label[DMSelection]["^"]+correction
                
                
                if Error_ppm_selection > Error:
                    correction = ""
                    final = df.loc[cont,seq_column_name]
                    SiteSolverOptions = ""
                    if Selection_Error2 > Selection_Error1:
                        DMSelection = DMSelected1
                        Error_ppm_selection = Selection_Error1
                    else:
                        DMSelection = DMSelected2
                        Error_ppm_selection = Selection_Error2
              
        # All columnas are anotated
        df.loc[cont, SiteSequence_column_name] = final
        df.loc[cont,SiteCorrection_column_name] = correction
        df.loc[cont,SiteOption_column_name] = SiteSolverOptions
        df.loc[cont,SiteDM_column_name] =  DMSelection
        df.loc[cont,SiteDMError_ppm_column_name] =  Error_ppm_selection



        cont=cont+1
    logging.info("Writing output file") 
    name = infile1[:-4]+Output_file_suffix
    try:
        remove(name+".txt")
    except:
        None

    
    df.to_csv(name+'.txt', index=False, sep='\t', encoding='utf-8')



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
    parser.add_argument('-pl', '--Primarylist', required=True, help='Path to sitelist file')
    parser.add_argument('-p', '--peptidesitelist', required=True, help='Path to peptidesitelist file')
    parser.add_argument('-sl', '--Secondarylist', required=True, help='Path to  Unimod sitelist file')
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

    infile1 = args.infile
    

    
    
    #start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))

        
    main(args.config, infile1, args.Primarylist, args.peptidesitelist, args.Secondarylist)
