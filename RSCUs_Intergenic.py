# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 10:45:22 2021

@author: Catarina
"""
import os, csv, mysql.connector, sshtunnel, json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
import pandas as pd

####################################################################################
#                                 User Information                                 #
####################################################################################

# SQL Connection Information
Database = 'CAISlocalGC'
User = 'catherineweibel'
Host = '127.0.0.1'
Password ='savemefromsandhu'

# The tables where the coding sequences are stored
# ** These are stored as a list so that we can iterate through them **
# CodingTable = ['gene_clump_sequence_table']
# Verbose prints out the progress of the script for the user when set to True
Verbose = False

####################################################################################
#                             Program Executes Below                               #
####################################################################################

# A connection to the SQL server is established
cnx = mysql.connector.connect(user= User, password= Password,host= Host , database= Database, connect_timeout=1000)
cursor = cnx.cursor(buffered = True)
cursor2 = cnx.cursor()

# change to go through each species ID on full table
query = ("SELECT * FROM local_gc_table")
cursor.execute(query)
#iterate though each window in the species
for window in cursor:
    if Verbose == True:
        print('Extracting UID: %s'%window[0])
        sys.stdout.flush()
    
    RawCount = {'F':{'TTT':0,'TTC':0},
                'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
                'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
                'Y':{'TAT':0,'TAC':0},
                '*':{'TAA':0,'TAG':0,'TGA':0},
                'C':{'TGT':0,'TGC':0},
                'W':{'TGG':0},
                'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
                'H':{'CAT':0,'CAC':0},
                'Q':{'CAA':0,'CAG':0},
                'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
                'I':{'ATT':0,'ATC':0,'ATA':0},
                'M':{'ATG':0},
                'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
                'N':{'AAT':0,'AAC':0},
                'K':{'AAA':0,'AAG':0},
                'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
                'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
                'D':{'GAT':0,'GAC':0},
                'E':{'GAA':0,'GAG':0},
                'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}
    
    Sum = {'F':0,
          'L':0,
          'S':0,
          'Y':0,
          '*':0,
          'C':0,
          'W':0,
          'P':0,
          'H':0,
          'Q':0,
          'R':0,
          'I':0,
          'M':0,
          'T':0,
          'N':0,
          'K':0,
          'V':0,
          'A':0,
          'D':0,
          'E':0,
          'G':0} 
    
    ProbTable = {'F':{'TTT':0,'TTC':0},
                 'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
                 'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
                 'Y':{'TAT':0,'TAC':0},
                 '*':{'TAA':0,'TAG':0,'TGA':0},
                 'C':{'TGT':0,'TGC':0},
                 'W':{'TGG':0},
                 'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
                 'H':{'CAT':0,'CAC':0},
                 'Q':{'CAA':0,'CAG':0},
                 'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
                 'I':{'ATT':0,'ATC':0,'ATA':0},
                 'M':{'ATG':0},
                 'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
                 'N':{'AAT':0,'AAC':0},
                 'K':{'AAA':0,'AAG':0},
                 'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
                 'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
                 'D':{'GAT':0,'GAC':0},
                 'E':{'GAA':0,'GAG':0},
                 'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}
    
    Sum_ProbTable= {'F':0,
                   'L':0,
                   'S':0,
                   'Y':0,
                   '*':0,
                   'C':0,
                   'W':0,
                   'P':0,
                   'H':0,
                   'Q':0,
                   'R':0,
                   'I':0,
                   'M':0,
                   'T':0,
                   'N':0,
                   'K':0,
                   'V':0,
                   'A':0,
                   'D':0,
                   'E':0,
                   'G':0}

    Codon_freqTable_expected = {'F':{'TTT':0,'TTC':0},
                 'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
                 'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
                 'Y':{'TAT':0,'TAC':0},
                 '*':{'TAA':0,'TAG':0,'TGA':0},
                 'C':{'TGT':0,'TGC':0},
                 'W':{'TGG':0},
                 'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
                 'H':{'CAT':0,'CAC':0},
                 'Q':{'CAA':0,'CAG':0},
                 'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
                 'I':{'ATT':0,'ATC':0,'ATA':0},
                 'M':{'ATG':0},
                 'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
                 'N':{'AAT':0,'AAC':0},
                 'K':{'AAA':0,'AAG':0},
                 'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
                 'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
                 'D':{'GAT':0,'GAC':0},
                 'E':{'GAA':0,'GAG':0},
                 'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}

    Codon_freqTable_raw = {'F':{'TTT':0,'TTC':0},
                           'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
                           'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
                           'Y':{'TAT':0,'TAC':0},
                           '*':{'TAA':0,'TAG':0,'TGA':0},
                           'C':{'TGT':0,'TGC':0},
                           'W':{'TGG':0},
                           'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
                           'H':{'CAT':0,'CAC':0},
                           'Q':{'CAA':0,'CAG':0},
                           'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
                           'I':{'ATT':0,'ATC':0,'ATA':0},
                           'M':{'ATG':0},
                           'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
                           'N':{'AAT':0,'AAC':0},
                           'K':{'AAA':0,'AAG':0},
                           'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
                           'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
                           'D':{'GAT':0,'GAC':0},
                           'E':{'GAA':0,'GAG':0},
                           'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}
    
    local_RSCUS_table = {'F':{'TTT':0,'TTC':0},
                        'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
                        'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
                        'Y':{'TAT':0,'TAC':0},
                        '*':{'TAA':0,'TAG':0,'TGA':0},
                        'C':{'TGT':0,'TGC':0},
                        'W':{'TGG':0},
                        'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
                        'H':{'CAT':0,'CAC':0},
                        'Q':{'CAA':0,'CAG':0},
                        'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
                        'I':{'ATT':0,'ATC':0,'ATA':0},
                        'M':{'ATG':0},
                        'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
                        'N':{'AAT':0,'AAC':0},
                        'K':{'AAA':0,'AAG':0},
                        'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
                        'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
                        'D':{'GAT':0,'GAC':0},
                        'E':{'GAA':0,'GAG':0},
                        'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}
    
######################################PUTTING STUFF INTO THE DICTIONARIES

    totalCodonCount = 0
    totalnucleotidecount = 0
    totalGCcount = 0

    window_UID = window[0]
    Species_UID = int(window[1])
    gene_clump_nucleotide_sequence = window[3].upper()
    sequenceLength = len(gene_clump_nucleotide_sequence)
    GC_local_prob = float(window[2])/100
    notGC_local_prob = 1- GC_local_prob

    if Verbose == True:
        print('local GC extracted: %s'%GC_local_prob)
        sys.stdout.flush()

    # We partition each coding sequence into a list codons and then compare them to our dictionary
    codonList = [gene_clump_nucleotide_sequence[n:n+3] for n in range(0,sequenceLength,3)]

    #Count up all of the codons and update dictionaries
    if codonList != []:
        for AA in RawCount:
            for Codon in RawCount[AA]:
                CodonCount = codonList.count(Codon)
                RawCount[AA][Codon] += CodonCount
                totalCodonCount += CodonCount
                totalnucleotidecount += CodonCount*3
                
                if Codon in ['TTA','TAT','ATT','AAT','ATA','TAA','AAA','TTT']:
                    totalGCcount +=0
                elif Codon in ['GGG','CCC','GCG','CGC','GGC','CCG','CGG','GCC']:
                    totalGCcount += CodonCount*3
                elif Codon in ['GAT','CAT','AGT','ACT','ATG','ATC','GTA','CTA','TGA','TCA','TAG','TAC','TTC','AAC']:
                    totalGCcount += CodonCount
                elif Codon in ['TTG','AAG','ACA','AGA','TCT','TGT','GTT','CTT','GAA','CAA']:
                    totalGCcount += CodonCount
                elif Codon in ['GGA','GGT','CCA','CCT','GAG','GTG','CAC','CTC','AGC','TGC','ACG','GCT']:
                    totalGCcount +=CodonCount*2
                elif Codon in ['AGG','TGG','ACC','GAC','CAG','GTC','CTG','TCC','TCG','CGT','CGA','GCA']:
                    totalGCcount +=CodonCount*2
            Sum[AA] = sum(RawCount[AA].values())
        GC_genic_prob = totalGCcount/totalnucleotidecount
        notGC_genic_prob = 1 - GC_genic_prob
    
    #calculate expected codon probabilities given the local GC context
    for AA in RawCount:
        prob_sum = 0
        for Codon in RawCount[AA]:
            if Codon in ['TTA','TAT','ATT','AAT','ATA','TAA','AAA','TTT']:
                Prob = (notGC_local_prob*notGC_local_prob*notGC_local_prob)*1/8
            elif Codon in ['GGG','CCC','GCG','CGC','GGC','CCG','CGG','GCC']:
                Prob = (GC_local_prob*GC_local_prob*GC_local_prob)*1/8
            elif Codon in ['GAT','CAT','AGT','ACT','ATG','ATC','GTA','CTA','TGA','TCA','TAG','TAC','TTC','AAC']:
                Prob = (GC_local_prob*notGC_local_prob*notGC_local_prob)*1/8
            elif Codon in ['TTG','AAG','ACA','AGA','TCT','TGT','GTT','CTT','GAA','CAA']:
                Prob = (GC_local_prob*notGC_local_prob*notGC_local_prob)*1/8
            elif Codon in ['GGA','GGT','CCA','CCT','GAG','GTG','CAC','CTC','AGC','TGC','ACG','GCT']:
                Prob = (GC_local_prob*GC_local_prob*notGC_local_prob)*1/8
            elif Codon in ['AGG','TGG','ACC','GAC','CAG','GTC','CTG','TCC','TCG','CGT','CGA','GCA']:
                Prob = (GC_local_prob*GC_local_prob*notGC_local_prob)*1/8
            else:
                print("storm the castle")
            ProbTable[AA][Codon] = Prob
            prob_sum += Prob
        Sum_ProbTable[AA] = prob_sum

    if Verbose == True: 
        print("-------SumTable-----------------")    
        print(Sum)      
        print("-------ProbTable-----------------")    
        print(ProbTable)   
        print("-------SUMProbTable-----------------")    
        print(Sum_ProbTable)
        sys.stdout.flush()

    #calculate expected and observed frequency of each codon
    for AA in RawCount:
        for Codon in RawCount[AA]:
            Codon_freq_exp = ProbTable[AA][Codon]/Sum_ProbTable[AA]
            Codon_freqTable_expected[AA][Codon] = Codon_freq_exp
            if Sum[AA] != 0:
              Codon_freq_raw = RawCount[AA][Codon]/Sum[AA]
              Codon_freqTable_raw[AA][Codon] = Codon_freq_raw
            else: 
              Codon_freqTable_raw[AA][Codon] = Codon_freqTable_expected[AA][Codon]
    
    #Calculate RSCUS values for each codon (if an amino acid does not appear in this sequence, its RSCUS will be 1, so it will not impact the geometric mean when calculating CAIS)
    for AA in RawCount:
        for Codon in RawCount[AA]:
            local_RSCUS = Codon_freqTable_raw[AA][Codon]/Codon_freqTable_expected[AA][Codon]
            local_RSCUS_table[AA][Codon]=local_RSCUS 

    #upload RSCUS values to SQL database
    sqlUploadCodingDataStatement = "INSERT INTO RSCUS_by_gene_clump_table (SpeciesUID,gene_clump_uid,sequenceLength,gc,RSCUS_Table,RawCount_Table) VALUES (%s,%s,%s,%s,%s,%s)"
    uploadCodingData = (Species_UID,window_UID,sequenceLength,GC_local_prob,json.dumps(local_RSCUS_table),json.dumps(RawCount))
    cursor2.execute(sqlUploadCodingDataStatement,uploadCodingData)

cnx.commit()
cnx.close()

    # for AA in local_RSCUS_table:
    #     for Codon in local_RSCUS_table[AA]:
    #         if Codon == "TTT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TTT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TTC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TTC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TTA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TTA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TTG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TTG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CTT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CTT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CTC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CTC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CTA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CTA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CTG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CTG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TCT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TCT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TCC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TCC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TCA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TCA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TCG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TCG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "AGT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_AGT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "AGC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_AGC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TAT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TAT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TAC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TAC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TAA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TAA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TAG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TAG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TGA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TGA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TGT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TGT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TGC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TGC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "TGG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_TGG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CCT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CCT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CCC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CCC=%s WHERE gene_clump_uid=%s" 
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CCA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CCA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CCG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CCG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CAT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CAT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CAC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CAC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CAA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CAA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CAG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CAG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CGT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CGT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CGC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CGC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CGA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CGA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "CGG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_CGG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "AGA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_AGA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "AGG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_AGG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "ATT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_ATT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "ATC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_ATC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "ATA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_ATA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "ATG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_ATG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "ACT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_ACT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "ACC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_ACC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "ACA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_ACA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "ACG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_ACG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "AAT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_AAT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "AAC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_AAC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "AAA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_AAA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "AAG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_AAG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GTT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GTT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GTC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GTC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GTA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GTA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GTG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GTG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GCT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GCT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GCC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GCC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GCA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GCA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GCG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GCG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GAT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GAT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GAC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GAC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GAA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GAA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GAG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GAG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GGT":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GGT=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GGC":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GGC=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GGA":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GGA=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
    #         if Codon == "GGG":
    #           sqlUploadCodingDataStatement = "UPDATE RSCUS_by_gene_clump_table SET RSCUS_GGG=%s WHERE gene_clump_uid=%s"
    #           cursor2.execute(sqlUploadCodingDataStatement,(local_RSCUS_table[AA][Codon],window_UID))
