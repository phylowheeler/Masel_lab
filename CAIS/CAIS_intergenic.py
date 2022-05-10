# -*- coding: utf-8 -*-
"""
Created on Fri Aug 3 14:14:26 2021

@author: Catarina

for each species, CAIS controlled for local GC content will be a geometric mean of RSCU values across gene clumps 

for ease of computation, calculate as log of (sum of RSCUs across clumps across a species) then raise e to that power.... to get a tiny number to calculate CAIS
"""
import os, csv, mysql.connector, sshtunnel, json, math
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
CodingTable = ['']
# Verbose prints out the progress of the script for the user when set to True
Verbose = False


####################################################################################
#                             Program Executes Below                               #
####################################################################################

# A connection to the SQL server is established
speciesUID_list = [14]
CAIS_dict = {}

for species in speciesUID_list:
    cnx = mysql.connector.connect(user= User, password= Password,host= Host , database= Database, connect_timeout=1000)
    cursor = cnx.cursor(buffered = True)
    query = ("SELECT RSCUS_Table, RawCount_Table FROM RSCUS_by_gene_clump_table WHERE SpeciesUID = %s")
    cursor.execute(query,(species,))

    total_codons = 0
    log_RSCUS_sum = 0
    for window in cursor:
        local_RSCUS = json.loads(window[0])
        RawCount = json.loads(window[1])
        for AA in RawCount:
            for codon in RawCount[AA]:
                if RawCount[AA][codon] != 0:
                    total_codons += RawCount[AA][codon]
                    log_RSCUS_sum += math.log(local_RSCUS[AA][codon])*RawCount[AA][codon]
    CAIS = math.exp(log_RSCUS_sum/total_codons)
    CAIS_dict[species] = CAIS
    print(CAIS)
    cnx.close()

# for i in range(1,10):
#     if Verbose == True:
#         print('Extracting SpeciesUID: %s'%i)
#         sys.stdout.flush()
        
#     SelectionStatement = "SELECT SUM(*) as row_sum FROM %s WHERE SpeciesUID = %s"%(CodingTable,i)
#     # note "*" is WHOLE ROW OF ENTRIES
#     mycursor.execute(SelectionStatement)
#     speciesSpecificResults = mycursor.fetchall()
    
#     LogOfCAIS = 0
#         for k in range(0,len(row_sum)):
#             LogOfCAIS += math.log(RelativeAdaptednessTable[AA][Codon])      
#             # We divide by the total number of codons in all coding sequences
#             LogOfCAIS = (1/totalCodonCount)*LogOfCAIS
#             # and we invert the log to get the CAIS
#             CAIS = math.exp(LogOfCAIS)  
            
