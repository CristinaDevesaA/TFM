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
import pandas as pd
from Bio import SeqIO
from os import remove
import math
import numpy as np
from optparse import OptionParser
import configparser
import argparse
import os
import logging
from pathlib import Path
import sys
import operator


###################
# Local functions #
###################
def readInfile(infile,cal_Dm_mh_colum_name):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile, sep="\t",                               
    float_precision='high',low_memory=False,dtype={cal_Dm_mh_colum_name:str})
    df[cal_Dm_mh_colum_name].astype("float64").dtypes
    return df


def Obtain_n (MasterProtein_column, dicc_fasta,clean_seq,m) :
    MasterProtein = MasterProtein_column.replace(" ","_") 
    MasterProtein=MasterProtein.strip("\n").split("_")
    MasterProtein= MasterProtein[0]+"_"+MasterProtein[1] # The id is extracted from the Master Protein name 
    final_q_pos = ""
    # The fasta sequence corresponding to this identifier is saved 
    for iden in dicc_fasta:
        if MasterProtein==iden:
            result=str(dicc_fasta[iden].seq.upper()).replace("X","L")
            break
    
    pattern=re.compile(clean_seq.replace("L","l").replace("I","[IL]").replace("l","[IL]")) # Problems that may exist with leucine and isoleucine are solved
    
    dicc_seqs={}
    pos = 0
    
    # The corresponding fasta sequence is rigorously scrutinized so that no chance is missed  
    while True:
        match = pattern.search(result, pos)
        if not match:
            break
        s = match.start()
        e = match.end()
        s = s+1
        if s-1 == 0:
            pos1 = 0
        else:
            pos1 = s-2
            
        try:
            p2 = result[e+1]
            pos2 = e+1
        except:
            pos2 = e
        s=s-1
        
        final_pos = s+1
        if final_q_pos == "":
            final_q_pos = str(final_pos+m-1)
            initial_q_pos = final_pos
        else: 
            final_q_pos = final_q_pos+","+str(final_pos+m-1)
            initial_q_pos = final_pos
            
        
        pos = e #  Move forward in text for the next search
    return final_q_pos,inital_q_pos

def ListMaker(df,seq,counts,dicc_fasta,MasterProtein):
    
    """
    ListMaker returns the frequency of the scan, the amino acid in the first position, cleaned sequence and m and l positions.
    """
    if seq.find("#")!=-1: 
        clean_seq = seq[seq.find("#")+1:] # Clean sequence is obtained.
        aa = "U" # Aminoacid first position
        freq=int(counts.loc[seq]) # Frecuency of the scan
        m = 0 # Position starting from N-term
        l = -1 # Position starting from C-term
        n = 0 # Protein position
        pd = clean_seq+":"+seq[seq.find("[")+1:seq.find("]")] # Amino acid positions and DM
        dqna = seq[seq.find("[")+1:seq.find("]")]
        x = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m)[1]
        o = x+len(clean_seq)
    elif seq.find(":")!=-1:
        clean_seq = seq[:seq.find(":")] # Clean sequence is obtained.
        aa = "U" # Aminoacid first position
        freq=int(counts.loc[seq]) # Frecuency of the scan
        m = -1 # Position starting from N-term
        l = 0 # Position starting from C-term
        pd = clean_seq+":"+seq[seq.find(":")+1]  # peptide and DM
        n = len(clean_seq)+1 # Protein position
        dqna = seq[seq.find("[")+1:seq.find("]")]
        x = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m)[1]
        o = x+len(clean_seq)
    elif seq.find("_") != -1 :
        clean_seq = seq[:seq.find("_")]
        aa = "U" # Aminoacid first position
        freq=int(counts.loc[seq]) # Frecuency of the scan
        m = np.nan # Position starting from N-term
        l = np.nan # Position starting from C-term
        pd = (seq[:seq.find("_")])+":"+(seq[seq.find("_")+1:]) # Peptide and DM
        dqna = (seq[seq.find("_")+1:])
        n = np.nan # Protein position
        x = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m)[1]
        o = x+len(clean_seq)
                                 
    else:
        clean_seq = seq[:seq.find("[")]+seq[seq.find("]")+1:] # Clean sequence is obtained.
        aa = seq[seq.find("[")-1] # Aminoacid first position
        freq=int(counts.loc[seq]) # Frecuency of the scan
        m = int(seq.find("[")) # Position starting from N-term
        l = int(len(clean_seq)-m+1) # Position starting from C-term
        pd = clean_seq+":"+seq[seq.find("[")+1:seq.find("]")] # peptide and DM
        dqna = seq[seq.find("[")+1:seq.find("]")] 
        n = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m)[0] # Protein positions
        x = n
        o = n
     
    return clean_seq,seq,aa,freq,m,l,pd,n,dqna,x,o


##################
# Main functions #
##################
def main(file,infile1,fastafile):
    
    """
    Reading configuration file
    """
    import pandas as pd
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(file)
    logging.info("Reading PDMTableMaker configuration file")  
    
    seq_column_name = config["PDMTableMaker_Parameters"].get("sequence_column_name") # Sequence with DM column name
    DM_column_name = config["PDMTableMaker_Parameters"].get("DM_column_name") # DM column name
    Theo_mh_column_name = config["PDMTableMaker_Parameters"].get("Theo_mh_column_name") # Theoretical column name
    MasterProtein_column_name =  config["PDMTableMaker_Parameters"].get("MasterProtein_column_name") # Master protien name column name 
    Outfile_suffix =  config["PDMTableMaker_Parameters"].get("Outfile_suffix") # Chosen suffix for output file
    nconditions =  config["PDMTableMaker_Conditions"].getint("number_of_conditions") # Number of conditions

    logging.info("Processing input file")
    
    dicc_fasta = SeqIO.index(fastafile, "fasta")

    
    df = readInfile(infile1,DM_column_name) # A dataframe is created based on input file
    counts = df[seq_column_name].value_counts() # The number of times that the species appear is saved in the variable count
    df2 = pd.DataFrame(columns=["p","q","pdm","pd","d","Theo_mh","ScanFreq","a","m","n","l","x","o","qdna","A","M","L","N","X","O","qdNA"],dtype=float) # Dataframe 2 is created with the aim of 

    cont = 0
    seqlist = [] # In this list it will be saved the sequences already analyzed
    dic_pd_M = {}
    for index, row in df.iterrows():
        meet = 0 
        if nconditions > 0: # if there are any conditions 
            for i in range(nconditions):
                i = i+1
                if row[config["PDMTableMaker_Conditions"].get("Condition"+str(i))] == config["PDMTableMaker_Conditions"].get("Value"+str(i)): # Each condition is checked
                    meet = meet+1
                                                                    
        if nconditions == meet: #If all conditions are met
            
            if  row[seq_column_name] not in seqlist:
                seqlist.append(row[seq_column_name])
                p,seq,aa,ScanFreq,m,l,pd,n,dqna,x,o= ListMaker(df,row[seq_column_name],counts,dicc_fasta,row[MasterProtein_column_name])              
                d = float(row[DM_column_name])
                qId = row[MasterProtein_column_name].split("|")[1]
                df2.loc[cont,"p"] = p
                df2.loc[cont,"q"] = qId
                df2.loc[cont,"pdm"] = row[seq_column_name]
                df2.loc[cont,"pd"] = pd
                df2.loc[cont,"d"] = d
                df2.loc[cont,"a"] = aa
                df2.loc[cont,"m"] = m
                df2.loc[cont,"l"] = l
                df2.loc[cont,"n"] = n
                df2.loc[cont,"x"] = x
                df2.loc[cont,"o"] = o
               
                df2.loc[cont,"ScanFreq"] =  int(ScanFreq)
                df2.loc[cont,"Theo_mh"] =  row[Theo_mh_column_name]
                
                if row[seq_column_name].find("_")==-1: 
                    qdna = qId+":"+str(dqna)+":"+str(n)+aa
                    df2.loc[cont,"qdna"] = qdna
                else:
                    qdna = qId+":"+str(dqna)+"::"+"U"
                    df2.loc[cont,"qdna"] = qdna
                    
                cont = cont+1
                # The M and L positions with maximun ScanFreq are saved ina dictionary
                if row[seq_column_name].find("_")==-1 and row[seq_column_name].find("#")==-1 and row[seq_column_name].find(":")==-1:
                    if p not in dic_pd_M.keys():
                        dic_pd_M[p] = {}
                    if d not in dic_pd_M[p].keys():
                        dic_pd_M[p][d] = {}
                        dic_pd_M[p][d] = ScanFreq,m,l,n,qdna,x,o
                    if ScanFreq > dic_pd_M[p][d][0]:
                        dic_pd_M[p][d] = ScanFreq,m,l,n,qdna,x,o
                    

    
    cont2 = 0
    for index, row in df2.iterrows(): # A M and  L columns are added 
        d = row["d"] 
        try:

            if row["pdm"].find("_") == -1:
                M = dic_pd_M[row["p"]][row["d"]][1]
                L = dic_pd_M[row["p"]][row["d"]][2]
                N = dic_pd_M[row["p"]][row["d"]][3]
                X = dic_pd_M[row["p"]][row["d"]][5]
                O = dic_pd_M[row["p"]][row["d"]][6]
                qdNA = dic_pd_M[row["p"]][row["d"]][4]
                df2.loc[cont2,"M"] = M
                df2.loc[cont2,"A"]=row["p"][M-1]
                df2.loc[cont2,"L"]= L
                df2.loc[cont2,"N"]= N
                df2.loc[cont2,"X"]= X
                df2.loc[cont2,"O"]= O
                df2.loc[cont2,"qdNA"]= qdNA

            else:
                df2.loc[cont2,"M"] = np.nan
                df2.loc[cont2,"A"] = "U"
                df2.loc[cont2,"L"] = np.nan
                df2.loc[cont2,"N"]= np.nan
                df2.loc[cont2,"qdNA"]= row["q"]+":"+str(round(row["d"],6))+"::"+"U"
            cont2 = cont2+1
        except: 
            pass


    logging.info("Writing output file")
    name = infile1[:-4]+Outfile_suffix+".txt"
    df2.to_csv(name, index=False, sep='\t', encoding='utf-8')
    logging.info('end script')
    return df2



if __name__ == '__main__':
    
    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='PDMTableMaker',
        epilog='''
        Example:
            python PDMTbleMaker.py
        ''')
      
    # default DM0Solver configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-f', '--fastafile', required=True, help='Path to input file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    
    args = parser.parse_args()
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + 'PDMTable_log.txt'
    log_file_debug = outfile = args.infile[:-4] + 'PDMTable_log_debug.txt'
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
    
    # configuration files are read      
    try:
        open('Solver.ini',"r")
        SiteSolverini='Solver.ini'
        logging.info("Modified SiteSolverini configuration file is going to be use")
        
    except:
        open("config/Solver.ini","r")
        SiteSolverini="config/Solver.ini"
        
    main(SiteSolverini, args.infile, args.fastafile)
