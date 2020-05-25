
# scFAN

Predicting transcription factor binding in single cells through deep learning

# Introduction

Characterizing genome-wide binding profiles of transcription factor (TF) is essential for understanding many biological processes. Although techniques have been developed to assess binding profiles within a population of cells, determining binding profiles at a single cell level remains elusive. Here we report scFAN (Single Cell Factor Analysis), a deep learning model that predicts genome-wide TF binding profiles in individual cells. scFAN is pre-trained on genome-wide bulk ATAC-seq, DNA sequence and ChIP-seq data, and utilizes single-cell ATAC-seq to predict TF binding in individual cells. We demonstrate the efficacy of scFAN by studying sequence motifs enriched within predicted binding peaks and investigating the effectiveness of predicted TF peaks for discovering cell types. We develop a new metric called "TF activity score" to characterize the state of each cell, and show that the activity scores can reliably capture cell identities. The method allows us to discover and study cellular identities and heterogeneity based on chromatin accessibility profiles.

# Prerequisites

-Python (2.7). Python 2.7.13 is recommended. we have also updated our codes for Python3.7, will upload in a separate subfolder, stay tuned. 

-Numpy

-Keras(2.0.2). For Python 3.7 version, it's Keras(2.3.1).

-Scipy

-[collections](https://docs.python.org/2.7/library/collections.html#)

-[pybedtools](https://daler.github.io/pybedtools/main.html)

-[pyfasta](https://pypi.org/project/pyfasta/)

# Data

Bulk ATAC-seq data was collected from Greenleaf et at(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47753), Schimidt et al(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70482) and Liu et al(https://www.ahajournals.org/doi/full/10.1161/CIRCRESAHA.116.310456?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub%3Dpubmed)

Single cell ATAC-seq data was retrieved from buenrostro et al(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65360) and Corces et al(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74310) which is same as chromVAR.

As for the three cell types related Chip-seq data, we gained GM12878 based Chip-seq data from ENCODE-DREAM Challenge dataset which used in FactorNet(https://www.sciencedirect.com/science/article/pii/S1046202318303293). We also retrieved the K562 and H1ESC Chip-seq data from Li's work(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1642-2#Decs).We also got other feature data such as mapability data from FactorNet.

Before training, you also need the copy of hg19 genome fasta data and put it into the resources folder. Because of the space limitation, the data is not included in the current folder. You can download it to your drive and unpress it in the following commands:  
 <pre><code>$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz 
$ tar zxvf chromFa.tar.gz   
$ cat chr*.fa > hg19.fa 
</code></pre>
# Usage

## Training:
<pre><code>$ python scFAN_train.py -i Datafolder -e 5 -oc outputdir
</code></pre>
* Parameters:  
\-- `-i`: input train data folder. e.g. Data/GM12787  
\-- `-e`: epoch times. eg. 5  
\-- `-oc`: model save path. e.g. model_out  
## Prediction on single cells:
<pre><code>$ python scFAN_predict.py -i Datafolder -scindir scATAC-seq_data_folder -moname motifname -pb True -oc modeldir
</code></pre>
* Parameters:  
\-- `-i`: input single cell data folder. e.g. /data2/fly/PBMCs/raw_data/gz_files/new_folder/LMPP  
\-- `-scindir`: scATAC-seq data folder. e.g. /data2/fly/scFAN_data/new_folder_PBMC_agg
\-- `-moname`: motif name. e.g. LMPP
\-- `-pb`: whether process batch effect. default: True
\-- `-oc`: model saved path. e.g. model_out
* Example:
\-- `python scFAN_predict.py -i /data2/fly/PBMCs/raw_data/gz_files/new_folder/LMPP -scindir /data2/fly/scFAN_data/new_folder_cisTopics_PBMC_agg -moname LMPP -pb True -oc multiTask_H1hESC_add_ATAC_moreTFs multiTask_GM12878_add_ATAC_moreTFs multiTask_K562_ATAC_more_chipdata`

## Acknowledgement
We referenced some of the codes from our lab's previous work [FactorNet](https://www.sciencedirect.com/science/article/pii/S1046202318303293) (here is the link to its [codes](https://github.com/uci-cbcl/FactorNet)) in the data preprocessing part and data generating structure, so we would like to thank for Daniel Quang's work and contribution. 
