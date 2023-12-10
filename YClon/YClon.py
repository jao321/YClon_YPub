import sys
import os 
import time
from YClon.util import *
from YClon.main import YClon
import pandas as pd
from rich.console import Console
from rich.markdown import Markdown


def main():
	version = "2.0.2 - Dec 4th 2023"
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

	if((("--version" in sys.argv) or ("-v" in sys.argv))):
		print(version)
	if((("--help" in sys.argv) or ("-h" in sys.argv))):
		console = Console()
		with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "docs","YClon_help.md"), "r+") as help_file:
			console.print(Markdown(help_file.read()))
		sys.exit(0)

	if("diversity" in sys.argv):
		clones = {}
		filename = sys.argv[2]
		print(filename)
		try:
			repertoire = pd.read_csv(filename,usecols = ['clone_id','clone_seq_count'], sep=separator)
		except:
			print("Please provide the path to a file with clone_id and clone_seq_count columns")
			exit()
		repertoire = repertoire.drop_duplicates(subset=["clone_id"])
		for x in repertoire["clone_id"]:
			clones[x]=repertoire[repertoire["clone_id"]==x]["clone_seq_count"].values[0]
			# exit()

		print("---------------------------------------------")
		print("DIVERSITY REPORT")
		print("\n")
		print("Simpson diversity:"+str(1-simpson_di(clones)))
		print("\n")
		print("Shannon diverity: "+str(shannon_di(clones)))
		print("\n")
		print("Shannon eveness: "+str(shannon_di(clones)/len(clones)))
		print("---------------------------------------------")
		print("\n")
		exit()



	filename = ""
	filename_temp = ""
	out_filename = ""
	out_report_name = ""
	if (any(".tsv" in i for i in sys.argv) == True) and (any("--input" in i for i in sys.argv) == False):
		filename = sys.argv[1]
	elif (any("--input" in i for i in sys.argv) == False) and (any("--folder" in i for i in sys.argv) == False):
		print("Please provide a file")
		exit()
	elif (any("--input" in i for i in sys.argv) == False) and (any("--folder" in i for i in sys.argv) == True):
		filename = sys.argv[sys.argv.index("--folder")+1]
	else:
		filename = sys.argv[sys.argv.index("--input")+1]

	print(filename)

	filename_temp = filename.split(".")

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
			out_filename = os.path.join(path,sys.argv[x+1]+"_YClon_clonotyped."+filename_temp[1])
		elif sys.argv[x].find("--short_output") != -1:
			short_output = True
		elif sys.argv[x].find("--folder") != -1:
			every_in_the_folder = True
			filename = sys.argv[x+1]


	
	path = directory_path(filename)

	if every_in_the_folder == False:
		filename_temp = filename.split(".")
		out_filename = filename_temp[0]+"_YClon_clonotyped."+filename_temp[1]
		YClon(out_filename, filename, thr, sequence_column, vcolumn, jcolumn, seqID, separator, ksize, short_output, all_cdrs,metric)
	# 	# clonotyping(filename, thr, sequence_column, vcolumn, jcolumn, seqID, separator, short_output,clonotyped)
	# else:
	# 	files = os.listdir(path)
	# 	for x in files:
	# 		if x.find(".tsv") != -1:
	# 			filename_temp = x.split(".")
	# 			out_filename = filename_temp[0]+"_YClon_clonotyped."+filename_temp[1]
	# # 			# clonotyping(path+"/"+x, thr, sequence_column, vcolumn, jcolumn, seqID, separator, short_output,clonotyped)
	
