import sys
import os 
import time
from YClon.util import *
from YClon.main import YClon
import pandas as pd
from rich.console import Console
from rich.markdown import Markdown


def main():
	version = "0.1 - Dec 4th 2023"
	clonotyped = False
	method = "AHAM"
	thr = 0.09
	metric = "hamming"
	cdr1 = "cdr1"
	cdr2 = "cdr2"
	sequence_column = "cdr3"
	all_cdrs = False
	vcolumn = "v_call"
	jcolumn = "j_call"
	seqID = "sequence_id"
	separator = "\t"
	ksize = 3
	short_output = False
	every_in_the_folder = False
	public = False
	analysis = False


	filename = ""
	filename_temp = ""
	out_filename = ""
	out_report_name = ""
	folder = ""
	if (any("--input" in i for i in sys.argv) == False) and (any("--folder" in i for i in sys.argv) == True):
		folder = sys.argv[sys.argv.index("--folder")+1]
	# if (any(".tsv" in i for i in sys.argv) == True) and (any("--input" in i for i in sys.argv) == False):
	# 	filename = sys.argv[1]
	# elif (any("--input" in i for i in sys.argv) == False) and (any("--folder" in i for i in sys.argv) == False):
	# 	print("Please provide a file")
	# 	exit()
	# elif (any("--input" in i for i in sys.argv) == False) and (any("--folder" in i for i in sys.argv) == True):
	# 	filename = sys.argv[sys.argv.index("--folder")+1]
	# else:
	# 	filename = sys.argv[sys.argv.index("--input")+1]
	# print(filename)

	# filename_temp = filename.split(".")

	for x in range(1,len(sys.argv)):
		if sys.argv[x].find("--input") != -1:
			filename = sys.argv[x+1]
		elif sys.argv[x].find("--method") != -1: 
			method = sys.argv[x+1]
		elif sys.argv[x].find("--metric") != -1: 
			metric = sys.argv[x+1]
		elif sys.argv[x].find("--thr") != -1:
			thr = float(sys.argv[x+1])
		elif sys.argv[x].find("--sequence") != -1:
			sequence_column = sys.argv[x+1]
		elif sys.argv[x].find("--v_gene") != -1:
			vcolumn = sys.argv[x+1]
		elif sys.argv[x].find("--j_gene") != -1:
			jcolumn = sys.argv[x+1]
		elif sys.argv[x].find("--seq_id") != -1:
			seqID = sys.argv[x+1]
		elif sys.argv[x].find("--sep") != -1:
			separator = sys.argv[x+1]
		elif sys.argv[x].find("--kmer_length") != -1:
			ksize = int(sys.argv[x+1])
		elif sys.argv[x].find("--dir_out") != -1:
			path = sys.argv[x+1]
			out_filename = os.path.join(path,os.path.basename(out_filename))
		elif sys.argv[x].find("--out") != -1:
			out_filename = os.path.join(path,sys.argv[x+1]+"_YPub_public_clones."+filename_temp[1])
		elif sys.argv[x].find("--short_output") != -1:
			short_output = True
		elif sys.argv[x].find("--analysis") != -1:
			analysis == True
		
		# elif sys.argv[x].find("--folder") != -1:
		# 	every_in_the_folder = True
		# 	filename = sys.argv[x+1]

	organise_repertoires_from_folder(folder, seqID, separator)
	filename = os.path.join(folder,"ypub_input.tsv")
	print(filename)
	filename_temp = filename.split(".")
	out_filename = filename_temp[0]+"_YPub_public_clones."+filename_temp[1]
	YClon(out_filename, filename, thr, sequence_column, vcolumn, jcolumn, seqID, separator, ksize, short_output, all_cdrs,metric)
	
	if analysis ==True:
		print("Analysing clones...")
		public_clones = pd.read_csv(out_filename,sep="\t")
		count_clones_in_more_than_one_rep = public_clones.origin_repertoire.groupby(public_clones['clone_id']).nunique()
		print("Counting public clones...")
		
		public_clones = public_clones[public_clones['clone_id'].map(count_clones_in_more_than_one_rep)>1]
		count_clones_in_more_than_one_rep=pd.DataFrame({"clone_id":count_clones_in_more_than_one_rep.index,
														"publicity":count_clones_in_more_than_one_rep.values})
		public_clones = pd.merge(public_clones,count_clones_in_more_than_one_rep)
		
		print("Saving file...")
		public_clones.to_csv(out_filename,quoting=False,index=False,sep="\t")
	os.remove(filename)
