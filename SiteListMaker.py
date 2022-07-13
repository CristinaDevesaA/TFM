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
import matplotlib.pyplot as plt
import operator
from os import rmdir
from shutil import rmtree
import seaborn as sns





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


def calculations(AA0,summ,x):
    """
    Caclculations functions returns the average off all positions and the frequency of position zero. 
    """
    
    average = summ/((x*2)+1) # Average of all positions
    freq = AA0-average # Background correction
    return average, freq


def SiteListMaker(seq,dic_SL,x,MassMod):
    
    """
    SisteListMaker function complets the dictionary in which the  DM, aminoacid and the repetitions for
    each postions are saved. It is necessary to know the sequence that is being analyzed,dictionary that is being completed (dic_SL), extension of the positions wanted to be
    anlalyzed (x).
    """
    
    listaj = [0]*((x*2)+1) # A list is created with as many zeros as the parameter "x" specifies

    # If the sequence does not have the requered extension (left and right), missing positions are completed      
    # with points with the aim of not taking them into account
    while len(seq[0:seq.find("[")]) < x+1:
        seq = "."+seq
    while len(seq[seq.find("]")+1:])< x:
        seq = seq+"." 
        
    seq = seq[seq.find("[")-x-1:seq.find("]")+x+1]
    seq=seq[0:seq.find("[")]+ seq[seq.find("]")+1:] # Compleated and cleaned sequence

    
    cont2 = 0
    for aa in seq:     
        if aa == ".": # If the aminoacid does not exist, the function will continue
            cont2 = cont2+1
            continue         
        if aa not in dic_SL[MassMod].keys(): # If is the first time that the aminoacid is being analyzed with this DM, it is added to aminoacid list
            dic_SL[MassMod][aa] = [0]*((x*2)+1)
                
        # All this information is saved in addictionary with a list         
        dic_SL[MassMod][aa][cont2] = dic_SL[MassMod][aa][cont2]+1
        listaj[cont2] = listaj[cont2]+1
        cont2 = cont2+1

    return dic_SL

def tablesMaker(dic,dic_NM,x):
    dic_SL=dic
    listdf=[]
    listdf_corr=[]
    listdf_corr_p0=[]
    x=5
    position_list = []
    for i in range((x*2)+1):
        position_list.append("P"+str((x-i)*-1))
    for DM in dic_SL:

        dfA_P= pd.DataFrame(dic_SL[DM]) # Para el cáclculo de A, para cada aa. Total de ese aa para esa DM y Para el cáclulo de P, para cada posición. Total posicción cero para esa DM. Y apra el calculo de Ptotal, suma de todos los P. 
        dfFreq_transpose = dfA_P.transpose() #Hacmeos la traspuesta para tener posición por columna y aa por fila
        dfFreq_transpose_DM_DMfreq= dfA_P.transpose() #Hacmos lo mismo para tener un anueva dfd donde  añaduir la Dm y su DMFreq y que no INETRfiera en klos claculos esa columamnn
        dfA_P["Total_pos"]=list(dfA_P.sum(axis=1)) #dfA_P con  el total por fila, es decir por aa
        dfA_P.loc[len(dfA_P.index)+1]=list(dfA_P.sum()) #dfA_P con el total por columna, es decir por posición
        dfA_P_transpose = dfA_P.transpose()


        Ptotal = dfA_P_transpose.iloc[len(dfA_P_transpose.index)-1,len(dfA_P_transpose.columns)-1]
        FreqDM = 0
        df_corr=dfFreq_transpose
        for col in dfFreq_transpose:
            controw = 0 
            for e in dfFreq_transpose[col]:

                A = dfA_P_transpose.iloc[controw,len(dfA_P_transpose.columns)-1]
                P = dfA_P_transpose.iloc[len(dfA_P_transpose.index)-1,col]
                new_value = e-((P*A)/Ptotal) #Claculo del esperado-observadp
                FreqDM = FreqDM+e
                df_corr.iloc[controw,col]=new_value
                controw=controw+1


        df_corr_p0 =df_corr.iloc[:,[5]] 
        df_p0 = dfFreq_transpose_DM_DMfreq.iloc[:,[5]]
        df_p0
        suma = list(df_p0.sum())
        Frequency_DM= int(suma[0])
        df_corr_p0_trans = df_corr_p0.transpose()

        df_corr["DMFreq"] = 0
        df_corr["DMFreq"] = suma*len(df_corr.index)
        df_corr["DM"] = 0
        df_corr["DM"] = DM
        df_corr['aa'] = df_corr.index
        df_corr_p0_trans ["DMFreq"] = 0
        df_corr_p0_trans ["DMFreq"] = suma*len(df_corr_p0_trans.index)
        df_corr_p0_trans ["DM"] = 0
        df_corr_p0_trans ["DM"] = DM
        dfFreq_transpose_DM_DMfreq['aa'] = dfFreq_transpose_DM_DMfreq.index
        dfFreq_transpose_DM_DMfreq["DM"] = DM
        dfFreq_transpose_DM_DMfreq["DMFreq"] = 0
        dfFreq_transpose_DM_DMfreq["DMFreq"] =  suma*len(dfFreq_transpose_DM_DMfreq.index)
        listdf_corr.append(df_corr)
        listdf_corr_p0.append(df_corr_p0_trans)
        listdf.append(dfFreq_transpose_DM_DMfreq)
        df_final=pd.concat(listdf)
        df_final_corr=pd.concat(listdf_corr)
        df_final_corr_p0_trans = pd.concat(listdf_corr_p0)


    columns_names_ordered = (sorted(df_final_corr_p0_trans.columns.values))
    columns_names_ordered.remove("DM")
    columns_names_ordered.remove("DMFreq")
    final_order = ["DM","DMFreq"]+columns_names_ordered
    final_order.append("Unassigned")

    df_final_corr_p0_trans_reindex = df_final_corr_p0_trans.reindex(columns = final_order)# new colum order is configured
    df_final_corr_p0_trans_reindex = df_final_corr_p0_trans_reindex.replace(np.nan, 0)
    df_final_corr.columns =position_list+["DMFreq","DM","aa"]
    df_final_corr_reindex = df_final_corr.reindex(columns =["DM","DMFreq","aa"]+position_list)# new colum order is configured
    df_final.columns =position_list+["DMFreq","DM","aa"]
    df_final_reindex = df_final.reindex(columns =["DM","DMFreq","aa"]+position_list)# new colum order is configured

    c = df_final_corr_p0_trans_reindex.shape[0]
    for DM in dic_NM:   
        c = c+1
        df_final_corr_p0_trans_reindex.loc[c]=[DM, dic_NM[DM]]+ ([0]*(len(df_final_corr_p0_trans_reindex.columns.values)-3))+[dic_NM[DM]]

    return df_final_corr_reindex,df_final_reindex,df_final_corr_p0_trans_reindex

def main(file,infile1):

    """
    Reading configuration file
    """
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(file)
 
    logging.info("Reading SiteListMaker configuration file") 
    
    seq_column_name = config["SiteListMaker_Parameters"].get("Sequence_column_name") # Sequence column name
    x = config["SiteListMaker_Parameters"].getint("x") # Parameter that indicates the extension (left and right) of the aminoacids wanted to be analyzed 
    Dm_column_name =  config["SiteListMaker_Parameters"].get("Calibrated_Delta_MH_column_name")#DM column name
    PeakAssignation_column_name =  config["SiteListMaker_Parameters"].get("PeakAssignation_column_name") #Name of column that contains peak assignation
    PeakNaming =  config["SiteListMaker_Parameters"].get("PeakNaming") #Parameter that indicates how peaks are named
    Clean_Frequency_Table =  config["SiteListMaker_Parameters"].get("Clean_Frequency_Table") #Name of the Clean Frequency Table ouput file
    Clean_P0_Frequency_Table =config["SiteListMaker_Parameters"].get("Clean_P0_Frequency_Table") #Name of the lean P0 Frequency Table ouput file
    Frequency_Table =  config["SiteListMaker_Parameters"].get("Frequency_Table") #Name of the Frequency Table ouput file
    
    logging.info("Processing input file")
    
    df = readInfile(infile1,Dm_column_name)
 
    dic_SL = {} # Dictionary for scanFreq is created
    dic_NM = {}
    seqlist = []
    cont = 0
    MaxPeptide = {}
    MaxScan={}

    for row in df.iterrows(): # Input file is processed row per row and dic_SL and dic_SL_Peptide are completed
        MassMod = df.loc[cont,Dm_column_name]
        if df.loc[cont,seq_column_name].find("_") !=-1 and df.loc[cont,PeakAssignation_column_name].upper()==PeakNaming.upper():
            aa = "Unassigned"
            if df.loc[cont,Dm_column_name] not in dic_NM.keys():
                dic_NM[df.loc[cont,Dm_column_name]]={}
                dic_NM[df.loc[cont,Dm_column_name]]=1
            else: 
                dic_NM[df.loc[cont,Dm_column_name]]=dic_NM[df.loc[cont,Dm_column_name]] +1  

        elif df.loc[cont,seq_column_name].find("_") ==-1 and df.loc[cont,PeakAssignation_column_name].upper()==PeakNaming.upper():
            if df.loc[cont,Dm_column_name] not in dic_SL.keys():
                dic_SL[df.loc[cont,Dm_column_name]]={}  

            dic_SL = SiteListMaker((df.loc[cont,seq_column_name]),dic_SL,x,df.loc[cont,Dm_column_name])
        cont = cont+1

    df_final_corr_reindex,df_final_reindex,df_final_corr_p0_trans_reindex=tablesMaker(dic_SL,dic_NM,x)
    logging.info("Writting output file")
    df_final_reindex.to_csv(infile1[:-4]+"_"+Frequency_Table+".txt", index=False, sep='\t', encoding='utf-8')
    df_final_corr_reindex.to_csv(infile1[:-4]+"_"+Clean_Frequency_Table+".txt", index=False, sep='\t', encoding='utf-8')
    df_final_corr_p0_trans_reindex.to_csv(infile1[:-4]+"_"+Clean_P0_Frequency_Table+".txt", index=False, sep='\t', encoding='utf-8')
    logging.info("end of the script")




if __name__ == '__main__':

    try:
        remove('Solver.ini')
    except:
        None

    # parse arguments
    parser = argparse.ArgumentParser(
        description='SiteListMaker',
        epilog='''
        Example:
            python SiteListMaker.py
        ''')
      
    # default DM0Solver configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
   
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)      
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + 'SL_log.txt'
    log_file_debug = outfile = args.infile[:-4] + 'SL_log_debug.txt'
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
    

        
    main(args.config, infile1)
