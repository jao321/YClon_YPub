'''
Functions and variables that will be used to auxiliate YClon processes
'''


import os
import gzip


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
		print("\nWARNING\nThere is no column named "+seqID+"\n", flush=True)
		exit()

	try:
		junc_indx = head.index(sequence_column)
	except:
		print("\nWARNING\nThere is no column named "+sequence_column+"\n", flush=True)
		exit()

	try:
		vGene_indx = head.index(vcolumn)
	except:
		print("\nWARNING\nThere is no column named "+vcolumn+"\n", flush=True)
		exit()

	try:
		jGene_indx = head.index(jcolumn)
	except:
		print("\nWARNING\nThere is no column named "+jcolumn+"\n", flush=True)
		exit()
	
	if all_cdrs == True:
		try:
			cdr2_indx = head.index(cdr2)
		except:
			print("\nWARNING\nThere is no column named "+cdr2+"\n", flush=True)
			exit()
		try:
			cdr1_indx = head.index(cdr1)
		except:
			print("\nWARNING\nThere is no column named "+cdr1+"\n", flush=True)
			exit()
        
		return seq_id_indx, junc_indx, vGene_indx, jGene_indx, cdr2_indx, cdr1_indx
	else:
		return seq_id_indx, junc_indx, vGene_indx, jGene_indx


def parse_AIRR(f,head, seqID, sequence_column, vcolumn, jcolumn, all_cdrs = False, separator = "\t"):
	head = head.split(separator)
	number_of_columns = len(head)

	if all_cdrs == True:
		seq_id_indx, junc_indx, vGene_indx, jGene_indx, cdr2_indx, cdr1_indx = get_columns_index(seqID, sequence_column, vcolumn, jcolumn, head,all_cdrs)
	else:
		seq_id_indx, junc_indx, vGene_indx, jGene_indx = get_columns_index(seqID, sequence_column, vcolumn, jcolumn, head,all_cdrs)
	
	colunas = [head[seq_id_indx],head[junc_indx]]

	i=0
	fail = 0
	clonotypes = {}
	print("Processing file", flush=True)
	for x in f:
		data = list(x.split(separator))
		if len(data)!= number_of_columns:
			fail+=1
			continue
		try:
			if len(data[junc_indx]) < 4:
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
						key = vGene[0]+"."+jGene[0]+"."+cdr3len+"."+cdr1len+"."+cdr2len
					else:
						fail+=1
						continue 
				key = vGene[0]+"."+jGene[0]+"."+cdr3len
			else:
				fail +=1
				continue
			clonotypes.setdefault(key, [])
			proclonotype = [data[seq_id_indx].strip(),data[junc_indx].strip().lower()]
			clonotypes[key].append(proclonotype)
		except:
			fail +=1
			continue
	key_with_largest_list = max(clonotypes, key=lambda k: len(clonotypes[k]))
	length_of_largest_list = len(clonotypes[key_with_largest_list])

	print("There are "+str(len(clonotypes))+" unique combinations of VJ+cdr lengths", flush=True)
	print(f"The largest is {key_with_largest_list}  with {length_of_largest_list} sequences", flush=True)
	if type(f)!=list:
		f.close()
	return clonotypes, colunas, seq_id_indx, junc_indx, vGene_indx, jGene_indx, fail


def organise_repertoires_from_folder(folder, sequence_column, vcolumn, jcolumn, seqID, separator, all_cdrs, format):
	rep_list = os.listdir(folder)
	ypub_input = open(os.path.join(folder,"ypub_input.tsv"), "w")
	header=False
	print("Organising input files", flush=True)
	for repertoire in rep_list:
		print("sorting "+repertoire, flush=True)
		if format=='OAS':
			f = gzip.open(os.path.join(folder,repertoire),'rt')
			f.readline()
			seq_counter_OAS=1
		else:
			f = open(os.path.join(folder,repertoire))
		for values in f:
			if(values.find(sequence_column) != -1) and (header==False):
				head = values.strip().split(separator)
				# prod_indx = tmp.index("productive")
				# for col_name in tmp:
				if format=='OAS':
					head.append('sequence_id')
					seq_id_indx, junc_indx, vGene_indx, jGene_indx = get_columns_index(seqID, sequence_column, vcolumn, jcolumn, head,all_cdrs)
					ypub_input.write('sequence_id'+separator+head[junc_indx]+separator+head[vGene_indx]+separator+head[jGene_indx]+separator)
				else:
					seq_id_indx, junc_indx, vGene_indx, jGene_indx = get_columns_index(seqID, sequence_column, vcolumn, jcolumn, head,all_cdrs)
					ypub_input.write(head[seq_id_indx]+separator+head[junc_indx]+separator+head[vGene_indx]+separator+head[jGene_indx]+separator)
				ypub_input.write("origin_repertoire\n")
				header = True
			elif(values.find(sequence_column) == -1) :
				tmp = values.strip().split(separator)
				# if tmp[prod_indx].find("T")!= -1:
				if format=='OAS':
					ypub_input.write(repertoire+'_'+str(seq_counter_OAS)+separator+tmp[junc_indx]+separator+tmp[vGene_indx]+separator+tmp[jGene_indx]+separator)
					seq_counter_OAS+=1
				else:
					ypub_input.write(tmp[seq_id_indx]+separator+tmp[junc_indx]+separator+tmp[vGene_indx]+separator+tmp[jGene_indx]+separator)
				# for col_name in tmp:
				# 	ypub_input.write(col_name+separator)
				ypub_input.write(repertoire+"\n")
	ypub_input.close()


