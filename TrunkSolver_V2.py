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
import math
import time
from pathlib import Path
from optparse import OptionParser
import configparser
from Bio import SeqIO
import numpy as np 
import pandas as pd 
from pandas import ExcelWriter
import argparse
import os
import logging
from pathlib import Path
import tkinter as tk





###################
# Local functions #
###################
def readInfile(infile):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile, sep="\t", float_precision='high', low_memory=False)
    return df

def get_fasta_report(file):
    '''
    Create the FASTA report
    '''
    def _create_key_id(rec):
        if (rec.startswith("sp") or rec.startswith("tr")) and "|" in rec:
            return rec.split("|")[1]
        else:
            return rec
    indb = SeqIO.index(file, "fasta", key_function=_create_key_id)
    return indb




def theoretical_mh_by_hand(subseq,label_mass,dic_mod,dic_aa,selectedaa,Mproton,Hydrogen,O2,decnum):
    
    """
    Theoretical mass is calculated taking into account fix modifications, label and subseq. This funtions returns theoretical
    mass fix modifications positions and the subsequence adding fix modifications.
    """
    decnum = int(decnum.replace("f","").replace(".",""))
    H20 = 2*Hydrogen+O2
    pattern = "[a-z]"
    
    sumall = 0
    newsequence = []
    c = 0
    mods_position = []
    
    # For each amino acid, depending on which is its condition (fix modifications, label, N-termina), mass is added  
    for aa in subseq:
        c = c+1
        
        if re.search(pattern,aa) != None: # If aa in N-terminal position  
            aa = aa.upper()
            
            if aa in dic_mod.keys(): # If it has a fix modification
                sumall = sumall+float(dic_mod[aa][1])
                sumall = sumall+label_mass
                newsequence.append(dic_mod[aa][0]+"TMT")
                mods_position.append(str(c)+"_"+aa+"_"+str(label_mass)+"_N")
                mods_position.append(str(c)+"_S_"+str(dic_mod[aa][2]))  
                    
            else: # If it has not a fix modification
          
                sumall = sumall+float(dic_aa[aa]) 
                sumall = sumall+label_mass
                cosa = float(dic_aa[aa])+float(label_mass)
                newsequence.append(aa+"-TMT")
                mods_position.append(str(c)+"_"+aa+"_"+str(label_mass)+"_N")
                
        else: # If it has not N-terminal position
         
            if aa in dic_mod.keys(): # If it has a fix modification
    
                sumall = sumall+float(dic_mod[aa][1])
                newsequence.append(dic_mod[aa][0])
                mods_position.append(str(c)+"_S_"+str(dic_mod[aa][2]))
   
                
            else: # If it has not a fix modification
                
                sumall = sumall+dic_aa[aa]
                newsequence.append(aa)
              
                
    mods_position = ",".join(mods_position)
    sumall = sumall+H20+Mproton 
    sumall = round(sumall,decnum)

    return sumall,newsequence,mods_position





def tag(seq,subseq):
    
    """
    Tag function returns a variable indicating if the cut is tryptic or not taking into account the 
    sequence plus one extra amino acid at both sides
    """

    Truncation = []
    number = subseq.find(seq)
    number2 = number+len(seq)
    loss = subseq.replace(seq,"")
    left = ""
    rigth = ""
    
    # Take into account if there is amino acid at both ends or just at one 
    if number2-number <= len(seq) and number != 0:
        left = subseq[0:number]
        rigth = subseq[number2:] 
    elif number != 0:
        left = loss    
    elif number2-number <= len(seq):
        rigth = loss
          

    # If the cut is tryptic or not is determined at both sides
    if left != "":

        if left[-1] == "K" or left[-1] == "R" and seq[0] != "P":   
            Truncation.append("YeS")
        else:
             Truncation.append("No") 
            
    if rigth != "":

        if seq[-1] == "K" or seq[-1] == "R" and rigth[0] != "P":
             Truncation.append("YeS")
        else:
             Truncation.append("No")
           
          
    # if one of the ends have a non tryptic cut means that there is a truncation
    if "No" in Truncation:
        Truncation1 = "No"
    else:
        Truncation1 = "YeS"

    return Truncation1





def Obtain_values(seq,MasterProtein_column,dic_fasta):
    
    """
    Taking in to account sequence, master protein and fasta dictionary of the input file this function returns 
    the entire sequence corresponding to the  master protein and all the initial and final positions that match with the 
    selected sequence (seq).
    """
    
    clean_seq = seq[:seq.find("[")]+seq[seq.find("]")+1:].upper() #The clean sequence is obtained.

    MasterProtein = MasterProtein_column.strip(" ")





    
    # The fasta sequence corresponding to this identifier is saved 
    for iden in dic_fasta:
        if MasterProtein == iden:
            result = str(dic_fasta[iden].seq.upper()).replace("X","L")
            break
   
    pattern = re.compile(clean_seq.replace("L","l").replace("I","[IL]").replace("l","[IL]")) # Problems that may exist with leucine and isoleucine are solved
    
    dic_seqs = {}
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
        
        cut = tag(result[s:e],result[pos1:pos2])
        dic_seqs[str(s)+":"+str(e)]=cut
        pos = e #  Move forward in text for the next search
    
    # If the same sequence is found in the same protein sequence those that have a tryptic cut will be selected     
    if "YeS" in dic_seqs.values():
        for key in list(dic_seqs):
            if dic_seqs[key] == "No":
                dic_seqs.pop(key)  
           
    return dic_seqs,result








def best_combination(subseq,Exp_Mh,cont,j,Error,label_mass,dic_mod,selectedaa,dic_CombList,dic_aa,decnum,Mproton,Hydrogen,O2, distanceDMsub,distanceDM,x):
    """
    Best_Combinations function returns the combinations of "Combination list" that give rise to an error less than or equal to the allowed and 
    variables that indicate  whether TrunkSolver should stop extending the length of the sequence being analyzed
    """
    # New subsequence mass is calulated by theoretical_mh_by_hand function

    ther,newsequence,mods_position = theoretical_mh_by_hand(subseq,label_mass,dic_mod,dic_aa,selectedaa,Mproton,Hydrogen,O2,decnum) 
  
    # initial parameters are set 
    ngreater = 0
    minimun_DiffMasa = 2*10**100
    minimun_DiffPPM = 2*10**100
    TrunkSequence = ""
    TrunkDM = ""
    TrunkLabel = ""
    Trunk_Label_ppm = ""
    DiffMasa2 = ""
    New_DM = ""
    New_Theo_MH = ""
    
    
    DiffMasa2 = Exp_Mh-ther # Mass difference between experimental an theoretical mass 
   
    # All posible options of "Combination List" are examined
    number_option = 0  

    for iden in dic_CombList: 
        number_option = number_option+1     
        mass = dic_CombList[iden]      
        total_value2 = ther+mass

        total_value2 = float(format (total_value2,decnum))

        DiffMasa = (total_value2-Exp_Mh) #Mass difference between experimental and theoretical plus one option of "Combination List"     

        # Mass difference is calculated in ppms  
        DiffPPM = abs(((total_value2-Exp_Mh)*1000000)/total_value2)

        # The lowest DiffPPM value is saved 
        if DiffPPM < minimun_DiffPPM:
            minimun_DiffMasa = DiffMasa2
            minimun_DiffPPM = DiffPPM
            DiffMasa_sequence = format( minimun_DiffMasa,decnum)
            TrunkSequence = (("".join(subseq)).upper()+"_"+str(DiffMasa_sequence))
            TrunkDM = minimun_DiffMasa
            TrunkLabel = iden
            Trunk_Label_ppm = minimun_DiffPPM
            New_DM = (Exp_Mh-total_value2)
            New_Theo_MH = total_value2 


      
    
    # if the distance between the original  DM site and the actual amino acid position is greater than the parameter x stablished by the user 
    # j variable indicates TrunkSolver function this subseq is not a possible option
    if abs(distanceDMsub-distanceDM) > int(x):
    
        j = "subseq1_stop"
   

    return minimun_DiffPPM,TrunkSequence,TrunkDM,TrunkLabel,mods_position,Trunk_Label_ppm,j,New_DM, New_Theo_MH





def TrunkSolver(seq,dic_seqs,Exp_mh,calibrated_delta_MH,result,Error,dic_aa,dic_CombList,dic_mod,NT_label,selectedaa,decnum,Mproton,Hydrogen,O2,Theo_mh,x):
    
   
    """
    This function creates all possible combinations of labels and truncations and returns the best option that
    matches all the parameters set by the user.e
    """
     
    # initial parameters
    minimun = 2*10**10
    final_Trunk_Label_ppm = 2*10**10 
    dic_result = {} # dictionary in which all results will be saved 
    match_number = 0
    x = int(x)
    addition = "" # Variable created with the aim of controlling if there are more than one possible subsequence  
   
    # All possible sequence found in protein sequence will be analyzed
    for key in dic_seqs:
        
        
        fields = key.split(":")
        position = int(fields[0]) # Initial position 
        final_position = int(fields[1]) # Final position 
        
        DMseqposition = seq.find("[")-1 
        DMresultposition = position + DMseqposition+1 # Position of oririginal DM site
        
        result=result.upper()
        result1 = list(result) 
        
        deltapeptide = seq[:seq.find("[")]+seq[seq.find("]")+1:] # Clean seqeunce
        
        # Total sequence lenght is obtained. 
        #Only those positions that are not further from the original site of the DM than the one indicated by the user will be taken into account.
        
        if len(seq[:seq.find("[")])<=x:
            firstpart = x
        else:
            
            firstpart = len(seq[:seq.find("[")])
        
        if len(seq[seq.find("]")+1:]) <= x:
            secondpart = x
        else:
            secondpart = len(seq[seq.find("]")+1:])
        
        totallength = secondpart+firstpart
        
        i = 1 #Position counter
        k = 1 #Position counter
        subseq1 = " " #Sequence starting form i position 
        subseq2 = " " #Sequence starting form k position 
        listseqs=1
     

        # While the extension of the subsequence can be elongated one amino acid will be added 
        while (len(subseq1) <= totallength or len(subseq2) <= totallength) and (len(subseq1) != 0 or len(subseq2) != 0) and listseqs != []:
            
            listseqs = [] # All possible subsequences will be saved in this list
            
            # C-term extension
            if len(subseq1) < len(seq)*3:
                distanceDMsub1 = position+1+i #Distance from the new position to the original DM site
                if position+i <= len(result1[position:])+position:
                    subseq1 = result1[position:position+i+1]
                    subseq1[0] = subseq1[0].lower()
                    listseqs.append(subseq1)
                   
                    i = i+1
                   

            # N-term extension
            if len(subseq2) < len(seq)*3:
                distanceDMsub2 = final_position-k
                if final_position-k >= 0:
                    subseq2 = result1[final_position-k:final_position]
                    subseq2[0] = subseq2[0].lower()
                    listseqs.append(subseq2)
                    
                    k = k+1
                   
          
 
            if len(listseqs) > 0:
                cont = 0
                for subseq in listseqs: # Best combination of this subsequence is saved

                    cont = cont+1  
                    j = ""
                    
                    # Name of distance variable dependidng of C-term o N-term extension
                    if cont == 1 : 
                        distanceDMsub = distanceDMsub1
                    elif cont == 2 : 
                        distanceDMsub = distanceDMsub2
                    
                    
                    minimun_DiffPPM,TrunkSequence,TrunkDM,TrunkLabel,mods_position,Trunk_Label_ppm,j,New_DM, New_Theo_MH = best_combination(subseq,Exp_mh,int(cont),j,Error,float(NT_label),dic_mod,selectedaa,dic_CombList,dic_aa,decnum,Mproton,Hydrogen,O2,distanceDMsub,DMresultposition,x)

                    if TrunkSequence != "" and j == "": # if best_combination  function finds a possible solution 
                       
                        # If the PPM difference is lower than the minimun the option will be saved in dic_result dictionary 
                        if abs(minimun_DiffPPM) <= minimun:
                            minimun = abs(minimun_DiffPPM)
                            final_TrunkSequence1 = TrunkSequence
                            final_TrunkSequence2 = final_TrunkSequence1[:final_TrunkSequence1.find("_")]
                            ptag = result[result.find(final_TrunkSequence2)-1:result.find(final_TrunkSequence2)+len(final_TrunkSequence2)+1]
                    
                            if deltapeptide==final_TrunkSequence2: 
                                cut=" "
                            else:
                                cut = tag(final_TrunkSequence2,ptag)
                            dic_result[minimun_DiffPPM] = TrunkDM,TrunkSequence,TrunkLabel,mods_position,Trunk_Label_ppm,cut,New_DM, New_Theo_MH



    # if there is more than one possibility those that have a trytic digestion will have preference

    if dic_result :
        for key in dic_result:
      
            if dic_result[key][4]<final_Trunk_Label_ppm: 

                    
                minimun = abs(key)
                final_TrunkDM = dic_result[key][0]
                final_TrunkSequence = dic_result[key][1]
                if dic_result[key][5] == "No":
                    final_TrunkLabel = "Truncation;"+dic_result[key][2]
                if dic_result[key][5] == " ":
                    final_TrunkLabel = dic_result[key][2]
                if dic_result[key][5].upper() == "YES":

                    final_TrunkLabel = "TrypticCut;"+dic_result[key][2]
                final_mods_position = dic_result[key][3]
                final_Trunk_Label_ppm = dic_result[key][4]
                final_New_DM = dic_result[key][6]
                final_New_Theo_MH = dic_result[key][7]


                    
                    
    # If no option is chosen, these variables will acquire their initial value
    else: 
        final_TrunkDM = calibrated_delta_MH
        final_TrunkSequence = seq
        final_TrunkLabel = ""
        final_mods_position = ""
        final_Trunk_Label_ppm = 0
        final_New_Theo_MH = Theo_mh 
        final_New_DM = calibrated_delta_MH
 
        
    if final_Trunk_Label_ppm > Error:
        final_TrunkDM = calibrated_delta_MH
        final_TrunkSequence = seq
        final_TrunkLabel = ""
        final_mods_position = ""
        final_New_Theo_MH = Theo_mh 
        final_New_DM = calibrated_delta_MH
        
       
    return final_TrunkSequence,final_TrunkDM,final_TrunkLabel,final_mods_position,minimun,final_Trunk_Label_ppm,addition,final_New_DM,final_New_Theo_MH







def applyTSSolver(row,MasterProtein_column_name,dic_fasta,Seq_column_name,Exp_mh_column_name,Delta_MH_cal_column_name,Error,dic_aa,dic_CombList,dic_mod,NT_label,selectedaa,decnum,Mproton,Hydrogen,O2,Theo_mh_column_name,x,TrunkSequence_output_column_name,TrunkDM_output_column_name,TrunkLabel_output_column_name,TrunkLabel_ppm_output_column_name,New_Theo_mh_output_column_name,New_Deltamass_output_column_name,Static_modifications_position_output_column_name,fix_mod_column_name,TrunkCleanPeptide_output_column_name,Missing_cleavages_output_column_name,Truncation_output_column_name):
    

        
    if row[Seq_column_name].find("_") != -1: # If another program has already corrected it, TrunkSolver does not annote anything new
        final_TrunkSequence = row[Seq_column_name]
        final_TrunkDM = row[Delta_MH_cal_column_name]
        final_TrunkLabel = " "
        final_TrunkCleanPeptide_output_column_name = final_TrunkSequence[:final_TrunkSequence.find("_")]
        final_mods_position = row[fix_mod_column_name]
        final_Trunk_Label_ppm = " "
        match_number = 0
        final_New_Theo_MH = row[Theo_mh_column_name]
        final_New_DM = row[Delta_MH_cal_column_name]
        addition = ""
                

    else:

        dic_seqs,result=Obtain_values(row[Seq_column_name],row[MasterProtein_column_name],dic_fasta)
        final_TrunkSequence,final_TrunkDM,final_TrunkLabel,final_mods_position,minimun,final_Trunk_Label_ppm,addition,final_New_DM,final_New_Theo_MH = TrunkSolver(row[Seq_column_name],dic_seqs,row[Exp_mh_column_name],row[Delta_MH_cal_column_name],result,Error,dic_aa,dic_CombList,dic_mod,NT_label,selectedaa,decnum,Mproton,Hydrogen,O2,row[Theo_mh_column_name],x)
        if final_TrunkSequence.find("_")!=-1: 
            final_TrunkCleanPeptide_output_column_name = final_TrunkSequence[:final_TrunkSequence.find("_")]
        else: 
            final_TrunkCleanPeptide_output_column_name = final_TrunkSequence[:final_TrunkSequence.find("[")]+final_TrunkSequence[final_TrunkSequence.find("]")+1:]


    number = final_TrunkCleanPeptide_output_column_name.count("K") + final_TrunkCleanPeptide_output_column_name.count("R")

    if final_TrunkCleanPeptide_output_column_name[-1]=="K" or final_TrunkCleanPeptide_output_column_name[-1]=="R": 
        number = number-1 
    number = number- final_TrunkCleanPeptide_output_column_name.count("KP")-final_TrunkCleanPeptide_output_column_name.count("RP")
    if number <0 : 
        number = 0 
    missingcleavage = number  
    if final_mods_position == " " or final_mods_position == "" :
        final_mods_position = row[fix_mod_column_name]
    if final_TrunkLabel.find("Trunc")!=-1:
        Truncated = 1
    else:
        Truncated = 0 
            
    row[TrunkSequence_output_column_name] = final_TrunkSequence
    row[TrunkDM_output_column_name] = final_TrunkDM
    row[TrunkLabel_output_column_name] = final_TrunkLabel
    row[TrunkCleanPeptide_output_column_name] =  final_TrunkCleanPeptide_output_column_name
    row[Missing_cleavages_output_column_name] =  missingcleavage
    row[TrunkLabel_ppm_output_column_name] = final_Trunk_Label_ppm
    row[New_Theo_mh_output_column_name]= final_New_Theo_MH
    row[New_Deltamass_output_column_name] = final_New_DM
    row[Static_modifications_position_output_column_name] = final_mods_position
    row[Truncation_output_column_name]=Truncated

    
    return row
    



##################
# Main functions #
##################

def main(file,file1,infile1, infilefasta):
    """
    Reading configuration file
    """
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(file)

    logging.info("Reading MassMod configuration file")
        
    Mproton = config["Masses"].getfloat("m_proton") #Proton mass
    Hydrogen = config["Masses"].getfloat("m_hydrogen") # Hydrogen mass
    O2 = config["Masses"].getfloat("m_oxygen") # Oxigen mass
        
    
    # Dictionary of aa at their masses is created
    dic_aa = {}       
    for option, value in config["Aminoacids"].items():
        dic_aa[option.strip("\n").upper()] = float(value.strip("\n").replace(",","."))
        
    
    # A  dictionary which fix modifications indicating their "new name", and the mass is created  
    dic_mod = {}   
    for option, value in config["Fix_modifications"].items():
        option = option.upper()
        new_name = value.split("\t")[1][1:]
        value = float(value.split("\t")[0].strip("\n").replace(",","."))
        if option.upper() == "NT":
            NT_label = value
        
        else:        
            if value == NT_label:
                selectedaa = option.upper() # The aa that has a fix label (ej ; K-TMT) is saved
        
            value2 = value+dic_aa[option]
            dic_mod[option[0]] = new_name,value2,value
            
    
    
    
    
    config.read(file1) # Read Trunksolver configuration file
    
    logging.info("Reading TrunkSolver configuration file")
              
    Error = config["TrunkSolver_Parameters"].getfloat("Relative_Error") # Relative error is saved     
    Exp_mh_column_name = config["TrunkSolver_Parameters"].get("exp_mh_column_name") # Experimental mh column name
    Theo_mh_column_name = config["TrunkSolver_Parameters"].get("theo_mh_column_name") # Theoretical mh name
    fix_mod_column_name = config["TrunkSolver_Parameters"].get("static_modifications_column_name") # Fix modifications column name  
    MasterProtein_column_name = config["TrunkSolver_Parameters"].get("MasterProtein_column_name") # Master protein column nam e
    Seq_column_name = config["TrunkSolver_Parameters"].get("Sequence_column_name") # Sequence colum name
    Delta_MH_cal_column_name = config["TrunkSolver_Parameters"].get("Calibrated_delta_mh_column_name") # DM column name
    New_Deltamass_output_column_name = config["TrunkSolver_Parameters"].get("New_Deltamass_output_column_name") # New deltamass  column name
    New_Theo_mh_output_column_name = config["TrunkSolver_Parameters"].get("New_Theo_mh_output_column_name") # New theoretical mh column name
    x =  config["TrunkSolver_Parameters"].get("x") # Number of positions to the right and left,  that the TrunkSolver is allowed to extend from the original DM site
    decnum =  config["TrunkSolver_Parameters"].get("decnum") # Number of decimales set by the user that will appear in the TrunkSequence column
    decnum = "."+str(decnum)+"f"
    
    TrunkSequence_output_column_name =  config["TrunkSolver_Parameters"].get("TrunkSequence_output_column_name") # Column name of the output where the chosen sequence is annotated
    TrunkCleanPeptide_output_column_name =  config["TrunkSolver_Parameters"].get("TrunkPlainPeptide_output_column_name") # Column name of the output where the chosen sequence is annotated
    Missing_cleavages_output_column_name =  config["TrunkSolver_Parameters"].get("Missing_cleavages_output_column_name") # Column name of the output where missing cleavages will be annotated
    Truncation_output_column_name =  config["TrunkSolver_Parameters"].get("Truncation_output_column_name") # Column name of the output where truncation are annotated
    TrunkDM_output_column_name =  config["TrunkSolver_Parameters"].get("TrunkDM_output_column_name") # Column name of the output where the recalculated DM is annotated, taking in to account the label 
    TrunkLabel_output_column_name =  config["TrunkSolver_Parameters"].get("TrunkLabel_output_column_name") # Column name of the output where the chosen label is annotated
    TrunkLabel_ppm_output_column_name =  config["TrunkSolver_Parameters"].get("TrunkLabel_ppm_output_column_name") # Column name of the output where the calculated error in ppm is annotated
    Static_modifications_position_output_column_name = config["TrunkSolver_Parameters"].get("Static_modifications_position_output_column_name") # Column name of the output where the  new fix modifications positions are annotated
    output_file_suffix = config["TrunkSolver_Parameters"].get("output_file_suffix") # Chosen suffix for output file 



    # All labels that want to be checked are save in dic_CombList dictionary
    dic_CombList = {}   
    for option, value in config["TrunkSolver_CombList"].items():
        dic_CombList[option.strip("\n").upper()] = float(value.strip("\n").replace(",","."))
        
    
    

    dic_fasta = get_fasta_report(infilefasta)
    df = readInfile(infile1)
    
    # Output columns are overwritten
    try: 
        df.drop([TrunkSequence_output_column_name,TrunkDM_output_column_name,TrunkLabel_ppm_output_column_name,TrunkLabel_ppm_output_column_name,New_Theo_mh_output_column_name,New_Deltamass_output_column_name,Static_modifications_position_output_column_name,Matchnumber_output_column_name], axis=1)
    except:
        pass
    
    df[TrunkSequence_output_column_name]=""
    df[TrunkDM_output_column_name]=np.nan
    df[TrunkLabel_output_column_name]=""
    df[TrunkLabel_ppm_output_column_name]=np.nan
    df[TrunkCleanPeptide_output_column_name] = ""
    df[TrunkCleanPeptide_output_column_name] = ""
    df[Missing_cleavages_output_column_name]=""
    df[New_Theo_mh_output_column_name]=np.nan
    df[New_Deltamass_output_column_name]=np.nan  
    df[Static_modifications_position_output_column_name]=""
    df[Truncation_output_column_name]=""

    
    
    logging.info("Processing input file")
    
    
    df = df.apply(lambda y: applyTSSolver(y,MasterProtein_column_name,dic_fasta,Seq_column_name,Exp_mh_column_name,Delta_MH_cal_column_name,Error,dic_aa,dic_CombList,dic_mod,NT_label,selectedaa,decnum,Mproton,Hydrogen,O2,Theo_mh_column_name,x,TrunkSequence_output_column_name,TrunkDM_output_column_name,TrunkLabel_output_column_name,TrunkLabel_ppm_output_column_name,New_Theo_mh_output_column_name,New_Deltamass_output_column_name,Static_modifications_position_output_column_name,fix_mod_column_name,TrunkCleanPeptide_output_column_name,Missing_cleavages_output_column_name,Truncation_output_column_name), axis = 1)
       
    # write outputfile
    logging.info("Writing output file")

    outfilename = infile1[:-4]
    outfile= outfilename+output_file_suffix+".txt"
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')


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
    
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-f', '--fastafile', required=True, help='Path to input fastafile')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    parser.add_argument('-Mm', '--Massconfig', default=defaultconfigM, help='Path to custom config.ini file')
    
    # these will overwrite the config if specified
    parser.add_argument('-r', '--relerror', help='Maximum allowable relative error (ppm)')
    parser.add_argument('-a', '--abserror', help='Maximum allowable absolute error')

    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
   
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.abserror is not None:
        config.set('TrunkSolver_Parameters', 'Absolute_Error', str(args.abserror))
        config.set('Logging', 'create_ini', '1')
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

    infile1=args.infile
    infilefasta=args.fastafile
    
    #start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    


    main(args.Massconfig,args.config, infile1,infilefasta)
  

