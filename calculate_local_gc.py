import numpy as np
import pandas as pd  
from pathlib import Path
import csv, os, sys, mysql.connector, sshtunnel
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from Bio.Seq import Seq, MutableSeq
from BCBio import GFF

import operator
from functools import reduce


Database = 'CAISlocalGC'
User = 'catherineweibel'
Host = '127.0.0.1'
Password ='savemefromsandhu'

#Breaks the assembly into scaffolds and sorts all CDS features by which scaffold they are found in, and lists the start and stop coordinates of each feature
def parse_input(gff_file,assembly_file):
	scaffold_dict = {}

	for seq_record in SeqIO.parse(assembly_file, "fasta"):
		scaffold = seq_record.id
		seq = (seq_record.seq)
		scaffold_dict[scaffold] = [seq,[[],[],[]]]
		in_handle = open(gff_file)
		limit_info = dict(gff_id=[scaffold],gff_type=["CDS"])
		for rec in GFF.parse(in_handle,limit_info=limit_info,target_lines=1):
			for feature in rec.features:
				scaffold_dict[scaffold][1][0].append(feature.id)
				scaffold_dict[scaffold][1][1].append(feature.location)
				scaffold_dict[scaffold][1][2].append(int(feature.qualifiers["phase"][0]))
		in_handle.close()
	return scaffold_dict

#Masks the coding regions within a sequence
def mask_sequences(sequence,features):
	seq = sequence
	newseq = MutableSeq(str(seq))
	if not features[0]:
		pass
	else:
		for i in range(len(features[0])):
			newseq = GeneMasker(newseq,features[1][i].start,features[1][i].end)
	return newseq

#Masking subfunction
def GeneMasker(sequence,gene_start,gene_end):
	gene_index_list = range(int(gene_start), int(gene_end))
	for index in gene_index_list:
		sequence[index] = "-"
	return sequence

# Identifies the start and stop coordinates of local sequence 'windows' surrounding coding regions and calculate local GC content
def create_gene_dict(masked_sequence,sequence,species_UID,features):
	gene_dict = {}
	window_size = 1500
	index = 0
	UID = 0

	#divide the scaffold sequence into local windows around each gene clump
	while index < len(masked_sequence):
		letter = masked_sequence[index]
		if letter == "-":
			windowstart = index - window_size
			if windowstart < 0:
				windowstart = 0
			while letter == "-" and index < len(masked_sequence):
				letter = masked_sequence[index]
				index += 1
			count = 0
			while count < window_size and index < len(masked_sequence):
				letter = masked_sequence[index]
				if letter == "-":
					count = 0
				else:
					count += 1
				index += 1
			windowend = index

			#calculate intergenic GC content within the window
			window_intergenic_seq = str(masked_sequence[windowstart:windowend]).replace("-","")
			gc = GC(Seq(window_intergenic_seq))

			# identify genes that are within this window
			window_features = gene_assigner(windowstart,windowend,features)
			if not window_features[0]:
				print("no genes in window")

			#merge together features within the same window, trim off codons that are split between windows
			else:
				merging_dict = {}
				for i in range(len(window_features[0])):
					gene_id = window_features[0][i]
					start = window_features[1][i].start
					end = window_features[1][i].end
					strand = window_features[1][i].strand
					phase = window_features[2][i]

					if strand == -1:
						gene_seq = sequence[start:end].reverse_complement()
						if gene_id in merging_dict:
							merging_dict[gene_id] = [gene_seq + merging_dict[gene_id][0],phase]
						else:
							merging_dict[gene_id] = [gene_seq,phase]
					else:
						gene_seq = sequence[start:end]
						if gene_id in merging_dict:
							merging_dict[gene_id][0] = merging_dict[gene_id][0] + gene_seq
						else:
							merging_dict[gene_id] = [gene_seq,phase]

				window_coding_seq = Seq('')
				for gene in merging_dict:
					merging_dict[gene][0] = merging_dict[gene][0][merging_dict[gene][1]:]
					overhang = len(merging_dict[gene][0]) % 3
					if overhang != 0:
						merging_dict[gene][0] = merging_dict[gene][0][:len(merging_dict[gene][0])-overhang]
					window_coding_seq = window_coding_seq + merging_dict[gene][0]

				#create dictionary of the coding sequence within each window, and its local GC content
				gene_dict[UID] = [species_UID,gc,window_coding_seq]
				UID += 1
		index += 1
	return gene_dict

# Assigns features from a genomic feature table to the sequence windows they are found within
def gene_assigner(windowstart,windowend,features):
	window_features = [[],[],[]]
	if not features[1]:
		print("no genes in window")
	else:
		for i in range(len(features[1])):
			if features[1][i].start-1 > windowstart and features[1][i].start-1 < windowend:
				window_features[0].append(features[0][i])
				window_features[1].append(features[1][i])
				window_features[2].append(features[2][i])
	return window_features

def main(argv):
	scaffold_dict = parse_input(sys.argv[1],sys.argv[2])
	species_UID = sys.argv[3]
	# intergenic_gc_table = pd.DataFrame(columns = ['species_UID','GC','sequence'])

	for scaffold in scaffold_dict:
		sequence = scaffold_dict[scaffold][0]
		features = scaffold_dict[scaffold][1]
		masked_sequence = mask_sequences(sequence,features)
		gene_dict = create_gene_dict(masked_sequence,sequence,species_UID,features)

		# df = pd.DataFrame.from_dict(data=gene_dict, orient='index',columns = ['species_UID','GC','sequence'])	
		# intergenic_gc_table = pd.concat([intergenic_gc_table,df],ignore_index = True,sort=False)

		for i in range(len(gene_dict)):
			Species_UID = gene_dict[i][0]
			GC = gene_dict[i][1]
			Sequence = str(gene_dict[i][2])

			#connect to MySQL database
			cnx = mysql.connector.connect(user= User, password= Password,host= Host , database= Database, connect_timeout=1000)
			mycursor = cnx.cursor(buffered = True)

			#upload data
			sqlUploadCodingDataStatement = "INSERT INTO local_gc_table (Species_UID,GC,Sequence) VALUES (%s,%s,%s)"
			uploadCodingData = (Species_UID,GC,Sequence)
			mycursor.execute(sqlUploadCodingDataStatement,uploadCodingData)
			cnx.commit()
			cnx.close()

	# filename = Path(species_UID)
	# filename = str(filename.with_suffix('')) + "_local_gc.csv"
	# intergenic_gc_table.to_csv(filename)

if __name__ == "__main__":
	main(sys.argv[1:])

# returns a csv file with each unique gene, its coding sequence, and its local gc content
# Input must be in the form python3, calculate_local_gc.py, annotation file in gff3 format, assembly file in fasta format, species UID