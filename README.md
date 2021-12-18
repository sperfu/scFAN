
# scFAN: Predicting transcription factor binding in single cells through deep learning

# Introduction

Characterizing genome-wide binding profiles of transcription factor (TF) is essential for understanding many biological processes. Although techniques have been developed to assess binding profiles within a population of cells, determining binding profiles at a single cell level remains elusive. Here we report scFAN (Single Cell Factor Analysis), a deep learning model that predicts genome-wide TF binding profiles in individual cells. scFAN is pre-trained on genome-wide bulk ATAC-seq, DNA sequence and ChIP-seq data, and utilizes single-cell ATAC-seq to predict TF binding in individual cells. We demonstrate the efficacy of scFAN by studying sequence motifs enriched within predicted binding peaks and investigating the effectiveness of predicted TF peaks for discovering cell types. We develop a new metric called "TF activity score" to characterize the state of each cell, and show that the activity scores can reliably capture cell identities. The method allows us to discover and study cellular identities and heterogeneity based on chromatin accessibility profiles.

# Prerequisites

-Python (2.7). Python 2.7.13 is recommended. we have also updated our codes for Python3.6, will upload in a separate subfolder. 

-Numpy

-Keras(2.0.2). For Python 3.6 version, it's Keras(2.0.6).

-theano(1.0.5)

-Scipy

-[collections](https://docs.python.org/2.7/library/collections.html#)

-[pybedtools](https://daler.github.io/pybedtools/main.html) (0.7.10)

-[pyfasta](https://pypi.org/project/pyfasta/)

-[ucsc-bigwigmerge](https://anaconda.org/bioconda/ucsc-bigwigmerge)

# Conda installation(python3 version)  

For users' convenience, we have also create conda env for users to clone and use.  
## Setup

Clone the repository.

```
git clone https://github.com/sperfu/scFAN.git
```

Navigate to the root of this repo and setup the conda environment.

```
conda env create -f scFAN.yml
```

Activate conda environment.

```
conda activate scFAN
```


# Data  
## Bulk ATAC-seq data  
Bulk ATAC-seq data was collected from [Greenleaf et at](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47753), [Schimidt et al](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70482) and [Liu et al](https://www.ahajournals.org/doi/full/10.1161/CIRCRESAHA.116.310456?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub%3Dpubmed)
## scATAC-seq data  
As for Corces (Buenrostro 2015 and Corces) dataset of 2210 cells, single cell ATAC-seq data was retrieved from [Buenrostro 2015 et al](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65360) and [Corces et al](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74310). 
As for PBMCs dataset from [Buenrostro_2018](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96772), it contains 2805 cells, the filtered 2034 cells can be downloaded from [here](https://www.dropbox.com/sh/8o8f0xu6cvr46sm/AAB6FMIDvHqnG6h7athgcm5-a/Buenrostro_2018.tar.gz?dl=0).
## Chip-seq data  
As for the three cell types related Chip-seq data, we retrieved GM12878, K562 and H1ESC Chip-seq data from [Li's work](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1642-2#Decs).We also got other feature data such as mapability data from FactorNet.
## other resources data
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
### Training on newly added TF:
Our model is capable of adopting other TF related Chip-seq data to train a new model on new TFs, details are in [Data folder](https://github.com/sperfu/scFAN/blob/master/Data).  

## Prediction on bulk cells:
Our pretrained model is capable of predicting TF binding using bulk data, the output will give you the probability of different TF binding using different bulk ATAC-seq data.
<pre><code>$ python predict_bulk_TF.py -i data/GM12878 -oc mul
</code></pre>
* Parameters:  
\-- `-i`: input train data folder. e.g. Data/GM12787  
\-- `-oc`: pretrained model save path. e.g. multiTask_K562_ATAC_more_chipdata

## Aggregating single cells to smooth the scATAC-seq data:
After retrieving scATAC-seq data from [here](https://drive.google.com/drive/folders/1R9V53HgpdrjYdFJ04nF_BxjaUfVI7LI1), users can perform aggregating procedure to smooth the scATAC-seq data by calculating the similarity between cells, and new bigwig files are stored in a new folder, see details in [generate_agg_data.py](https://github.com/sperfu/scFAN/blob/master/generate_agg_data.py).  (You need to change the directory in the script to your own before running the code).  
## Prediction on single cells:
<pre><code>$ python scFAN_predict.py -i Datafolder -scindir scATAC-seq_data_folder -moname motifname -pb True -oc modeldir
</code></pre>
* Parameters:  
\-- `-i`: input single cell data folder, which contains peaks and other resources. e.g. /data2/fly/PBMCs/raw_data/gz_files/new_folder/LMPP  (These can be retrieved from our [google drive](https://drive.google.com/drive/folders/1R9V53HgpdrjYdFJ04nF_BxjaUfVI7LI1),under single cell data subfolder, named as celltype.tar.gz, e.g. LMPP.tar.gz etc.)  
\-- `-scindir`: scATAC-seq data folder, which contains all the processed bigwig files. e.g. /data2/fly/scFAN_data/new_folder_PBMC_agg (if you want to use scATAC-seq without aggregation, just change this parameter to the raw scATAC-seq folder) (These can be retrieved from our [google drive](https://drive.google.com/drive/folders/1R9V53HgpdrjYdFJ04nF_BxjaUfVI7LI1),under single cell data subfolder, named as PBMC_scATAC_bigwig.tar.gz & Corces_scATAC_bigwig.tar.gz, note that these two files are quite large.)   
\-- `-moname`: single cell type name. e.g. LMPP  
\-- `-pb`: whether process batch effect. default: True  
\-- `-oc`: model saved path. e.g. model_out  
* Example:  
 `python scFAN_predict.py -i /data2/fly/PBMCs/raw_data/gz_files/new_folder/LMPP -scindir /data2/fly/scFAN_data/new_folder_cisTopics_PBMC_agg -moname LMPP -pb True -oc multiTask_H1hESC_add_ATAC_moreTFs multiTask_GM12878_add_ATAC_moreTFs multiTask_K562_ATAC_more_chipdata`

## Some path information
<strong><em>Please</em></strong> be aware that some of the path in the codes may need to altered accorrding to user's own path. I listed all the paths that need to alter in the following:

* change this path in [utils.py](https://github.com/sperfu/scFAN/blob/4efa63381702676e44e9582eed5903fda428ec20/utils.py#L14) (as well as py3 version [utils.py](https://github.com/sperfu/scFAN/blob/4efa63381702676e44e9582eed5903fda428ec20/python3_codes/utils.py#L14)) to your own temp folder to allow pybedtool to create some temp files.
* In generate_agg_data.py, [This line](https://github.com/sperfu/scFAN/blob/4efa63381702676e44e9582eed5903fda428ec20/generate_agg_data.py#L56) (as well as py3 version [line](https://github.com/sperfu/scFAN/blob/4efa63381702676e44e9582eed5903fda428ec20/python3_codes/generate_agg_data_py3.py#L56)) contains the precalculated latent embeddings from Cistopic(using latent dim number of 20). We use this embedding to precalculate cell-cell similarities as mentioned in our paper. Please refer to [Cistopic](https://github.com/aertslab/cisTopic/tree/76ba23dedb60042bf9610537a0727100c2d4c486) git respository to calculate the embedding matrix.
* In generate_agg_data.py, [This line](https://github.com/sperfu/scFAN/blob/4efa63381702676e44e9582eed5903fda428ec20/generate_agg_data.py#L59) (as well as py3 version [line](https://github.com/sperfu/scFAN/blob/4efa63381702676e44e9582eed5903fda428ec20/python3_codes/generate_agg_data_py3.py#L59)) contains all the single cell names, you can change this path to your own, all the cell names (named: all_cell_name_2210.txt & all_cell_name_no_dir) are stored in our [google drive](https://drive.google.com/drive/folders/1R9V53HgpdrjYdFJ04nF_BxjaUfVI7LI1), please download them from two single cell data subfolders(Corces dataset & PBMCs dataset). As for the [line 26 & 32](https://github.com/sperfu/scFAN/blob/4efa63381702676e44e9582eed5903fda428ec20/generate_agg_data.py#L26) (py3 version [line 26 & 32](https://github.com/sperfu/scFAN/blob/4efa63381702676e44e9582eed5903fda428ec20/python3_codes/generate_agg_data_py3.py#L26)), you also need to change to your own path which contains all the single cell bigwig(bw) file(These can also be downloaded through google drive from single cell data subfolders mentioned above, you might need to unzip them into your own folder after downloading).

## Citation
If you use our tool, please cite our work: Fu L, Zhang L, Dollinger E, Peng Q, Nie Q, Xie X. [Predicting transcription factor binding in single cells through deep learning](https://advances.sciencemag.org/content/6/51/eaba9031). Sci Adv. 2020 Dec 18;6(51):eaba9031. doi: 10.1126/sciadv.aba9031. PMID: 33355120.

## Acknowledgement
We referenced some of the codes from our lab's previous work [FactorNet](https://www.sciencedirect.com/science/article/pii/S1046202318303293) (here is the link to its [codes](https://github.com/uci-cbcl/FactorNet)) in the data preprocessing part and data generating structure, so we would like to thank for Daniel Quang's work and contribution.  

## Other git reference  
scFAN is also available at UCI Xie Lab, see website at [scFAN](https://github.com/uci-cbcl/scFAN).
 

