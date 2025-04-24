import os
import time
import numpy as np
import pandas as pd
from alive_progress import alive_bar
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity, pairwise_distances
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial import distance
from scipy.cluster.hierarchy import linkage
# from multiprocessing import Pool
from YClon.util import parse_AIRR, directory_path
import multiprocessing
from functools import partial
import itertools
import tarfile


def calculate_hamming_distance(i, j):
    return distance.hamming(i, j)

def build_kmers_tf_idf(sequence, ksize=3): 
    ''' 
        Gets a sequence and returns a list with
        n-grrams of ksize available in that sequence (sliding window = 1)
    '''
    ngrams = zip(*[sequence[i:] for i in range(ksize)])
    return [''.join(ngram) for ngram in ngrams]

def colapse_unique(clonotypes_value, colunas,sequence_column):
    if len(clonotypes_value) > 1:
        pre_clone = pd.DataFrame(clonotypes_value)
        pre_clone.columns = colunas
        return pre_clone.drop_duplicates(sequence_column)
    else:
        return clonotypes_value

def most_common(lst):
    return max(set(lst), key=lst.count)

def clonotype(pre_clone, clonotypes_vjcdr, key, seqID, sequence_column,  
              thr=0.09, ksize=3, metric="hamming",seq_clone_id=[]):
        starting_time = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        print(f"[{starting_time}] clustering {key} {len(pre_clone)} unique junction sequences", flush=True)
        junc_seq = pre_clone[sequence_column]
        if metric == "kmer":
            clusterer = AgglomerativeClustering(distance_threshold = thr, n_clusters= None, metric="precomputed",linkage='single')
            vectorizer = CountVectorizer(min_df=1, analyzer=lambda x: build_kmers_tf_idf(x, ksize))
            tf_idf_matrix = vectorizer.fit_transform(junc_seq)	
            dist = 1 - cosine_similarity(tf_idf_matrix)
        elif metric == "hamming":
            input_pair = [list(map(int,list(x.lower().replace("a","0").replace("c","1").replace("t","2").replace("g","3").replace("n","4")))) for x in junc_seq]
            clusterer = AgglomerativeClustering(distance_threshold = thr, n_clusters= None,metric="precomputed",linkage='single')
            dist = pairwise_distances(input_pair,metric="hamming", n_jobs=-1)
            
        cluster_pre_clone = clusterer.fit(dist)
        clone = cluster_pre_clone.labels_
        junc_seq = junc_seq.reset_index(drop=True)
        junc_seq =  junc_seq.to_frame()
        junc_seq['clone_id'] = clone
        clonotypes_df = pd.DataFrame(clonotypes_vjcdr)
        clonotypes_df.columns = [seqID,sequence_column]

        seq_clone_id = clonotypes_df.merge(junc_seq)
        seq_clone_id = seq_clone_id[['sequence_id','clone_id']].to_dict('split')['data']
        # for i in range(0, len(clone)):
        #     if clone[i] != -1:
        #         clone_id = key+'_'+str(int(clone[i])+1)
        #     else:
        #         clone_id = key+'_'+str(int(clone[i]))

        #     all_again = list(tmp_df[tmp_df[sequence_column] == junc_seq[i]][seqID].unique())
        #     for k in all_again:
        #         seq_cluster = k+","+str(clone_id)
        #         if seq_cluster not in seq_clone_id:
        #             seq_clone_id.append(seq_cluster)
        #     # print(len(seq_clone_id))	
        return seq_clone_id



def write_output(in_airr, seqID, out_filename, 
                 clonotipo, separator, seq_id_indx, vGene_indx, 
                 jGene_indx, junc_indx, short_output=False):
    most_common_cdr3 = {}
    most_common_seq_id = {}
    if short_output ==  True:
        out_small_name = out_filename.replace("_YClon_clonotyped.","_YClon_clonotyped_only_essential_columns.")
        out = open(out_small_name, 'w+')
    else:
        out = open(out_filename, 'w+')
    for x in in_airr:
        if x.find(seqID) == -1:
            data = x.strip().split(separator)
            if data[seq_id_indx] in clonotipo:
                if short_output == True:
                    if data[seq_id_indx] in clonotipo:
                        out.write(data[seq_id_indx]+separator+data[vGene_indx]+separator+data[jGene_indx]+separator+data[junc_indx]+separator)
                        out.write(clonotipo[data[seq_id_indx]])
                else:
                    for i in range(0, len(data)):
                        out.write(data[i].strip()+separator)
                    out.write(clonotipo[data[seq_id_indx]]+"\n")
                if clonotipo[data[seq_id_indx]] not in most_common_cdr3:
                    most_common_cdr3[clonotipo[data[seq_id_indx]]] = []
                    most_common_seq_id[clonotipo[data[seq_id_indx]]] = []
                    most_common_cdr3[clonotipo[data[seq_id_indx]]].append(data[junc_indx].strip())
                    most_common_seq_id[clonotipo[data[seq_id_indx]]].append(data[seq_id_indx])
                else:
                    most_common_cdr3[clonotipo[data[seq_id_indx]]].append(data[junc_indx].strip())
                    most_common_seq_id[clonotipo[data[seq_id_indx]]].append(data[seq_id_indx])
        else:
            data = x.strip().split(separator)
            seq_id_indx = data.index(seqID)
            if short_output == True:
                out.write(seqID+separator+data.index(vGene_indx)+separator+data.index(jGene_indx)+separator+data.index(junc_indx)+separator+"clone_id\n")
            else:
                if x.find("clone_id") == -1:
                    out.write(x.strip()+separator+"clone_id\n")
                    clonotyped = False
                else:
                    out.write(x.strip()+separator+"clone_id_YClon\n")
                    clonotyped = True
    return most_common_seq_id, most_common_cdr3, clonotyped



def add_seq_count(path,seqID,clonotyped,maior,separator,out_filename):
    temp_filename = path+"YClon_temp.txt"
    temp = open(temp_filename, 'w')
    out = open(out_filename, 'r')
    for x in out:
        if x.find(seqID) != -1:
            temp.write(x.strip()+separator+"clone_seq_count\n")
            if clonotyped == True:
                i = x.strip().split(separator).index("clone_id_YClon")
            else:
                i = x.strip().split(separator).index("clone_id")
        else:
            temp.write(x.strip()+separator+str(len(maior[x.strip().split(separator)[i]]))+"\n")
    temp.close()
    out.close()
    os.remove(out_filename)
    os.rename(temp_filename,out_filename)
    

def write_report(most_common_cdr3,most_common_seq_id,maior,out_filename):
    if out_filename.find("_YClon_clonotyped")!= -1:
        out_report_name=out_filename.replace("_YClon_clonotyped.","_YClon_clonotyped_report.")
    else:
        out_report_name=out_filename.replace(".tsv","report.tsv")
    out_report = open(out_report_name, 'w+')
    out_report.write("sequence_id\tseq_count\tmost_common_cdr3\tclone_id\n")
    for i in most_common_cdr3:
        cdr3 = most_common(most_common_cdr3[i])
        out_report.write(most_common_seq_id[i][most_common_cdr3[i].index(cdr3)]+"\t"+str(len(maior[i]))+"\t"+cdr3+"\t"+i+"\n")

def process_key(key, clonotypes, colunas, sequence_column, seqID, thr, ksize, metric):
    seq_clone_id = []
    pre_clone = colapse_unique(clonotypes[key], colunas, sequence_column)
    
    if len(pre_clone) > 1:
        seq_clone_id = clonotype(pre_clone, clonotypes[key], key, seqID, sequence_column, thr, ksize, metric)
    else:
        pre_clone = pd.DataFrame(clonotypes[key])
        pre_clone.columns = colunas
        seq_id = pre_clone[seqID]
        for i in range(0, len(clonotypes[key])):
            seq_clone_id=[str(seq_id[i])+","+key+"_1"]
    return seq_clone_id

def YClon(out_filename,filename, thr, sequence_column, vcolumn, jcolumn, seqID, separator, ksize, short_output, all_cdrs,metric): 
    start_time = time.time()
    print("Opening and reading "+filename, flush=True)
    if tarfile.is_tarfile(filename):
        with tarfile.open(filename) as tar:
            binary = tar.extractfile(filename.replace('.tar.gz',''))
            f = binary.read().decode('utf-8').split('\n')
            head = f[0]
            f = f[1:]
    else: 
        f = open(filename, 'r')
        head = f.readline().strip()
    in_airr = open(filename, 'r')
    clonotypes, colunas, seq_id_indx, junc_indx, vGene_indx, jGene_indx, fail = parse_AIRR(f,head, seqID, sequence_column, vcolumn, jcolumn, all_cdrs, separator)
    path = directory_path(filename)
    temp_filename = path+"YClon_temp.txt"

    
    # worker = partial(process_key, 
    #             clonotypes=clonotypes,
    #             colunas=colunas,
    #             sequence_column=sequence_column,
    #             seqID=seqID,
    #             thr=thr,
    #             ksize=ksize,
    #             metric=metric)
    
    # num_processes = multiprocessing.cpu_count()
    # with multiprocessing.Pool(processes=num_processes) as pool:
    #     results = pool.map(worker, clonotypes.keys())
    
    # # Combine results
    # results = (itertools.chain.from_iterable(results))
    # results = []
    print('creating temporary file...', flush=True)
    temp = open(temp_filename, 'w')
    count = 0
    for key in clonotypes: #each key is the combination of V gene, J gene and the length of cdr3, the values are the sequence ID and cdr3 sequence
    # bar()
        pre_clone = colapse_unique(clonotypes[key],colunas, sequence_column)
        if len(pre_clone) > 1:
            results=clonotype(pre_clone, clonotypes[key], key, seqID, sequence_column, thr, ksize, metric)
            # print(len(results))
            for x in results:
                temp.write(x[0]+','+key+'_'+str(x[1])+'\n')
        else:
            ct=1
            pre_clone = pd.DataFrame(clonotypes[key])
            pre_clone.columns = colunas
            seq_id = pre_clone[seqID]
            for i in range(0, len(clonotypes[key])):
                temp.write(str(seq_id[i])+","+key+"_"+str(ct))
        # unico_pq_VJLen +=1
        # count += 1
        # total_clust += 1
        # pre_clone = pd.DataFrame(clonotypes[key])
        # pre_clone.columns = colunas
        # seq_id = pre_clone[seqID]
        # for i in range(0, len(clonotypes[key])):
        #     temp.write(str(seq_id[i])+","+str(count)+"\n")
    # results = (itertools.chain.from_iterable(results))
    

    print('Assigning clonotypes...', flush=True)
    # for x in results:
    #     temp.write(x+'\n')
    

    temp.close()
    in_temp = open(temp_filename, 'r')
    out = open(out_filename, 'w+')

    clonotipo = {}
    maior ={}
    maximo = 0


    for x in in_temp:
        data = x.strip().split(",")
        clonotipo[','.join(data[0:-1])] = data[-1].strip()
        if data[-1] not in maior:
            maior[data[-1]] = []
            maior[data[-1]].append(','.join(data[0:-1]))
        else:
            maior[data[-1]].append(','.join(data[0:-1]))

    seq_list = []
    if tarfile.is_tarfile(filename):
        with tarfile.open(filename) as tar:
            binary = tar.extractfile(filename.replace('.tar.gz',''))
            in_airr = binary.read().decode('utf-8').split('\n')
    else: 
        in_airr = open(filename, 'r')
    most_common_seq_id, most_common_cdr3, clonotyped = write_output(in_airr, seqID, out_filename, clonotipo, separator, seq_id_indx, vGene_indx, jGene_indx, junc_indx, short_output)

    out.close()
    in_temp.close()
    # os.remove(temp_filename)

    add_seq_count(path, seqID, clonotyped, maior,separator,out_filename)

    write_report(most_common_cdr3,most_common_seq_id,maior,out_filename)

    current_time = time.time()
    elapsed_time = current_time - start_time

    print("The work was completed in: " + "%.3f" % int(elapsed_time) + " seconds", flush=True)
    print(str(fail)+ " sequences could not be assigned to any clones", flush=True)
