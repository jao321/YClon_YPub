import os
import time
import numpy as np
import pandas as pd
from alive_progress import alive_bar
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity, pairwise_distances
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial import distance
from multiprocessing import Pool
from YClon.util import parse_AIRR, directory_path


def calculate_hamming_distance(i, j):
    return distance.hamming(i, j)

def distance_matrix_hamming(input_pair, n_jobs=None):
    pool = Pool(processes=n_jobs)
    results = pool.starmap(calculate_hamming_distance, [(input_pair[i], input_pair[j]) for i in range(len(input_pair)) for j in range(len(input_pair))])
    distance_matrix = np.zeros((len(input_pair), len(input_pair)))
    ct=0
    for i in range(len(input_pair)):
        for j in range(len(input_pair)):
            distance_matrix[i, j] = results[ct]
            ct+=1
        return distance_matrix
    
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

def clonotype(pre_clone, clonotypes, seqID, sequence_column, unico_pq_VJLen, total_clust, count, temp, thr=0.09, ksize=3, metric="hamming"):
        seq_id = pre_clone[seqID]

        junc_seq = pre_clone[sequence_column]
        if metric == "kmer":
            clusterer = AgglomerativeClustering(distance_threshold = thr, n_clusters= None, metric="precomputed",linkage='complete')
            vectorizer = CountVectorizer(min_df=1, analyzer=lambda x: build_kmers_tf_idf(x, ksize))
            tf_idf_matrix = vectorizer.fit_transform(junc_seq)	
            dist = 1 - cosine_similarity(tf_idf_matrix)
        elif metric == "hamming":
            # input_pair = [list(x.lower()) for x in junc_seq]
            input_pair = [list(map(int,list(x.lower().replace("a","0").replace("c","1").replace("t","2").replace("g","3").replace("n","4")))) for x in junc_seq]
            clusterer = AgglomerativeClustering(distance_threshold = thr, n_clusters= None,metric="precomputed",linkage='complete')
            # dist = distance_matrix_hamming(input_pair, n_jobs=None)
            dist = pairwise_distances(input_pair,metric="hamming",n_jobs=-1)
            
        cluster_pre_clone = clusterer.fit(dist)
        clone = cluster_pre_clone.labels_
        maior = count + cluster_pre_clone.labels_.max() + 1
        junc_seq = junc_seq.reset_index(drop=True)
        key = list(pre_clone["v_call"].unique())[0].split("*")[0]+","+list(pre_clone["j_call"].unique())[0].split("*")[0]+","+len(list(pre_clone[sequence_column].unique())[0])
        clonotypes_df = pd.DataFrame(clonotypes[key])
        clonotypes_df.columns = [seqID,sequence_column]
        tmp_df = pre_clone
        tmp_df.columns = ["tmp_seq",sequence_column]
        pre_clone = pd.merge(tmp_df,clonotypes_df,on=[sequence_column])[[seqID,sequence_column]]
        del tmp_df
        del clonotypes_df
        for i in range(0, len(clone)):
            if clone[i] != -1:
                clone_id = count + int(clone[i]) +1
            else:
                clone_id = maior + 1
                maior += 1
                total_clust += 1

            all_again = pre_clone.loc[pre_clone[sequence_column] == junc_seq[i]][seqID]
            for k in all_again:
                temp.write(k+","+str(clone_id)+"\n")	
        count = maior
        total_clust += cluster_pre_clone.labels_.max() + 1
        return unico_pq_VJLen, total_clust, maior


def write_output(in_airr, seqID, out_filename, clonotipo, separator, seq_id_indx, vGene_indx, jGene_indx, junc_indx, short_output=False):
    most_common_cdr3 = {}
    most_common_seq_id = {}
    if short_output ==  True:
        out_small_name = out_filename.replace("_YClon_clonotyped.","_YClon_clonotyped_only_essential_columns.")
        out = open(out_small_name, 'w+')
    else:
        out = open(out_filename, 'w+')
    for x in in_airr:
        if x.find(seqID) == -1:
            data = x.split(separator)
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
            data = x.split(separator)
            seq_id_indx = data.index(seqID)
            if short_output == True:
                out.write(seqID+separator+vcolumn+separator+jcolumn+separator+sequence_column+separator+"clone_id\n")
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
            # print(x.strip().split(separator)[i])
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
	  # print(i)
        cdr3 = most_common(most_common_cdr3[i])
	  # print(most_common_seq_id[i])
        out_report.write(most_common_seq_id[i][most_common_cdr3[i].index(cdr3)]+"\t"+str(len(maior[i]))+"\t"+cdr3+"\t"+i+"\n")


def YClon(out_filename,filename, thr, sequence_column, vcolumn, jcolumn, seqID, separator, ksize, short_output, all_cdrs,metric): 
    start_time = time.time()
    clonotypes, colunas, seq_id_indx, junc_indx, vGene_indx, jGene_indx, file_size, fail = parse_AIRR(filename, seqID, sequence_column, vcolumn, jcolumn, all_cdrs, separator)
    path = directory_path(filename)
    temp_filename = path+"YClon_temp.txt"
    temp = open(temp_filename, 'w')
    count = 0
    total_clust = 0
    unico_pq_VJLen = 0
    maior = 0
    a = 0

    with alive_bar(len(clonotypes), title="Clonotyping") as bar: 
        for key in clonotypes: #each key is the combination of V gene, J gene and the length of cdr3, the values are the sequence ID and cdr3 sequence
            bar()
            
            pre_clone = colapse_unique(clonotypes[key],colunas, sequence_column)
            if len(pre_clone) > 1:
                unico_pq_VJLen, total_clust, maior = clonotype(pre_clone, clonotypes[key], seqID, sequence_column, unico_pq_VJLen, total_clust, count, temp, thr,ksize,metric)
            else:
                unico_pq_VJLen +=1
                count += 1
                total_clust += 1
                pre_clone = pd.DataFrame(clonotypes[key])
                pre_clone.columns = colunas
                seq_id = pre_clone[seqID]
                for i in range(0, len(clonotypes[key])):
                    temp.write(str(seq_id[i])+","+str(count)+"\n")

    pre_clone = []

    temp.close()
    in_temp = open(temp_filename, 'r')
    in_airr = open(filename, 'r')
    filename_temp = filename.split(".")
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

    most_common_seq_id, most_common_cdr3, clonotyped = write_output(in_airr, seqID, out_filename, clonotipo, separator, seq_id_indx, vGene_indx, jGene_indx, junc_indx, short_output)

    out.close()
    in_temp.close()
    os.remove(temp_filename)

    add_seq_count(path, seqID, clonotyped, maior,separator,out_filename)

    write_report(most_common_cdr3,most_common_seq_id,maior,out_filename)

    current_time = time.time()
    elapsed_time = current_time - start_time

    print("The work was completed in: " + "%.3f" % int(elapsed_time) + " seconds")
    print(str(fail)+ " sequences could not be assigned to any clones")
