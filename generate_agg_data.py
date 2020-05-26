import numpy as np
import pdb
import pandas as pd
import commands
import multiprocessing

from scipy.stats import pearsonr,spearmanr,kendalltau

def cal_sim(z_features):
    pdb.set_trace()
    z_feature = pd.read_csv(z_features,header=0,index_col=0,sep='\t')
    z_feature = z_feature.as_matrix()
    pdb.set_trace()
    sim_all = np.zeros((z_feature.shape[0],z_feature.shape[0]),dtype=float)
    for ix in range(z_feature.shape[0]):
        if ix%1000 == 0:
            print('now sim %d/2034'%(ix))
        for jx in range(z_feature.shape[0]):
            sim_all[ix][jx] = cos_sim(z_feature[ix],z_feature[jx])
    return sim_all

def sort_find_neighbors(sim_all,cell_name_list,k_neighbor=50):
    sort_matrix = np.argsort(-sim_all)
    pdb.set_trace()
    for index,item in enumerate(cell_name_list):
        #sel_name_list = ' '.join(["bigwig_all/"+cell_name_list[tmp]+".bw" for tmp in list(sort_matrix[index][:k_neighbor])])
        sel_name_list = '/data2/fly/PBMCs/raw_data/call_peaks/'+cell_name_list[index]+'/'+cell_name_list[index]+'.bw '+' '.join(["/data2/fly/PBMCs/raw_data/call_peaks/"+cell_name_list[tmp]+'/'+cell_name_list[tmp]+".bw" for tmp in list(sort_matrix[index][:k_neighbor])])
        aa = commands.getoutput('mkdir bedgraph_PBMC')
        aa = commands.getoutput('bigWigMerge '+sel_name_list+' bedgraph_PBMC/'+cell_name_list[index]+'.bedGraph')
        cc = commands.getoutput('sort -k1,1 -k2,2n '+'bedgraph_PBMC/'+cell_name_list[index]+'.bedGraph > bedgraph_PBMC/'+cell_name_list[index]+'.bedgraph')
        dd = commands.getoutput('rm -f bedgraph_PBMC/'+cell_name_list[index]+'.bedGraph')
        dd = commands.getoutput('mkdir new_folder_cisTopics_PBMC_agg')
        bb = commands.getoutput('bedGraphToBigWig bedgraph_PBMC/'+cell_name_list[index]+'.bedgraph'+' /data2/fly/fly_scFAN_bak_20190603/resources/hg19.autoX.chrom.sizes new_folder_cisTopics_PBMC_agg/'+cell_name_list[index]+'.bw')
        if index%100 == 0:
            print('now num: %d/2034'%(index))
        if bb != '':
            print item

def pearsonrSim(x,y):
    return pearsonr(x,y)[0]


def cos_sim(vector_a, vector_b):
    """
    :return: sim
    """
    vector_a = np.mat(vector_a)
    vector_b = np.mat(vector_b)
    num = float(vector_a * vector_b.T)
    denom = np.linalg.norm(vector_a) * np.linalg.norm(vector_b)
    cos = num / denom
    sim = 0.5 + 0.5 * cos
    return sim


if __name__=='__main__':
    all_sim = cal_sim('/data2/fly/PBMCs/raw_data/compared_models/feature_cisTopic20_PBMC_scale.txt')
    #all_sim = cal_sim('cis_matrix_topic')
    #all_sim = cal_sim('/data2/fly/PBMCs/raw_data/compared_models/Tsne_matrix_cisTopic20')
    cell_names = open('/data2/fly/PBMCs/raw_data/all_cell_name_no_dir','r').readlines()
    #cell_names = open('/data2/fly/PBMCs/raw_data/all_cell_test2','r').readlines()
    cell_name_list = [item.strip() for item in cell_names]
    sort_find_neighbors(all_sim,cell_name_list,k_neighbor=100)

