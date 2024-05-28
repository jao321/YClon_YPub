'''
Functions and variables that will be used to auxiliate YClon processes
'''


import os



def directory_path(file_path):
	temp =	file_path.split(os.sep)
	file_path = file_path.replace(temp[len(temp)-1],"")
	return(file_path)

def simpson_di(data):
	#from https://gist.github.com/martinjc/f227b447791df8c90568


    ''' Given a hash { 'species': count } , returns the Simpson Diversity Index
    
    >>> simpson_di({'a': 10, 'b': 20, 'c': 30,})
    0.3888888888888889
    '''

    def p(n, N):
        """ Relative abundance """
        if n == 0:
            return 0
        else:
            return float(n)/N

    N = sum(data.values())
    
    return sum(p(n, N)**2 for n in data.values() if n != 0)


def shannon_di(data):
    #from https://gist.github.com/audy/783125

    ''' Given a hash { 'species': count } , returns the SDI
    
    >>> sdi({'a': 10, 'b': 20, 'c': 30,})
    1.0114042647073518 '''
    
    from math import log as ln
    
    def p(n, N):
        """ Relative abundance """
        if n == 0:
            return 0
        else:
            return (float(n)/N) * ln(float(n)/N)
            
    N = sum(data.values())
    
    return -sum(p(n, N) for n in data.values() if n != 0)




def get_columns_index(seqID,sequence_column,vcolumn,jcolumn,head,all_cdrs):
	try:
		seq_id_indx = head.index(seqID)
	except:
		print("\nWARNING\nThere is no column named "+seqID+"\n")
		exit()

	try:
		junc_indx = head.index(sequence_column)
	except:
		print("\nWARNING\nThere is no column named "+sequence_column+"\n")
		exit()

	try:
		vGene_indx = head.index(vcolumn)
	except:
		print("\nWARNING\nThere is no column named "+vcolumn+"\n")
		exit()

	try:
		jGene_indx = head.index(jcolumn)
	except:
		print("\nWARNING\nThere is no column named "+jcolumn+"\n")
		exit()
	
	if all_cdrs == True:
		try:
			cdr2_indx = head.index(cdr2)
		except:
			print("\nWARNING\nThere is no column named "+cdr2+"\n")
			exit()
		try:
			cdr1_indx = head.index(cdr1)
		except:
			print("\nWARNING\nThere is no column named "+cdr1+"\n")
			exit()
        
		return seq_id_indx, junc_indx, vGene_indx, jGene_indx, cdr2_indx, cdr1_indx
	else:
		return seq_id_indx, junc_indx, vGene_indx, jGene_indx


def parse_AIRR(filename, seqID, sequence_column, vcolumn, jcolumn, all_cdrs = False, separator = "\t"):
	f = open(filename, 'r')
	x = f.readline().strip()
	head = x.split(separator)
	number_of_columns = len(head)

	if all_cdrs == True:
		seq_id_indx, junc_indx, vGene_indx, jGene_indx, cdr2_indx, cdr1_indx = get_columns_index(seqID, sequence_column, vcolumn, jcolumn, head,all_cdrs)
	else:
		seq_id_indx, junc_indx, vGene_indx, jGene_indx = get_columns_index(seqID, sequence_column, vcolumn, jcolumn, head,all_cdrs)
	
	colunas = [head[seq_id_indx],head[junc_indx]]

	file_size = 0
	i=0
	fail = 0
	clonotypes = {}
	print("Processing file")
	for x in f:
		file_size +=1
		data = list(x.split(separator))
		if len(data)!= number_of_columns:
			print('The number of columns is different than the number of header names given the separator')
			fail+=1
			continue
		try:
			if len(data[junc_indx]) < 4:
				print('Less than 4 sequences that share the same V, J gene and CDR3 length')
				fail += 1
				continue
			cdr3len = str(len(data[junc_indx])).strip()
			vGene = data[vGene_indx].split(',')[0].split('*') #include all v gene alleles
			jGene = data[jGene_indx].split(',')[0].split('*') #include all j gene alleles
			if all_cdrs == True:
				cdr1len = str(len(data[cdr1_indx])).strip()
				cdr2len = str(len(data[cdr2_indx])).strip()
                
			if jGene != "" and vGene != "" and cdr3len != 0:
				if all_cdrs == True:
					if cdr1len != 0 and cdr2len != 0:
						key = vGene[0]+","+jGene[0]+","+cdr3len+","+cdr1len+","+cdr2len
					else:
						fail+=1
						continue 
				key = vGene[0]+","+jGene[0]+","+cdr3len
			else:
				print('Either V or J gene are not annotated for this row, or the cdr3 is empty')
				fail +=1
				continue
			clonotypes.setdefault(key, [])
			proclonotype = [data[seq_id_indx].strip(),data[junc_indx].strip().lower()]
			clonotypes[key].append(proclonotype)
		except:
			fail +=1
			continue
		
	f.close()
	return clonotypes, colunas, seq_id_indx, junc_indx, vGene_indx, jGene_indx, file_size, fail


def organise_repertoires_from_folder(folder, seqID, separator):
	rep_list = os.listdir(folder)
	ypub_input = open(os.path.join(folder,"ypub_input.tsv"), "w")
	header=False
	print("Organising input files")
	for repertoire in rep_list:
		print(repertoire)
		with open(os.path.join(folder,repertoire)) as f:
			for values in f:
				if(values.find(seqID) != -1) and (header==False):
					tmp = values.strip().split(separator)
					prod_indx = tmp.index("productive")
					for col_name in tmp:
						ypub_input.write(col_name+"\t")
					ypub_input.write("origin_repertoire\n")
					header = True
				elif(values.find(seqID) == -1) :
					tmp = values.strip().split(separator)
					if tmp[prod_indx].find("T")!= -1:
						for col_name in tmp:
							ypub_input.write(col_name+"\t")
						ypub_input.write(repertoire+"\n")

