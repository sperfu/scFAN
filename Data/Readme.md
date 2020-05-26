# Data
The Data folder consists of most of the training bulk data and all the the single cell data.
## Training data
As described in the main text, in the training phase, we trained three models based on different TFs combination, which is generally from three cell lines: GM12878, K562 and H1ESC.   
There are these files listed below in the training data:  
  `-- bigwig.txt` consisting of all the epigenomic file you use in the training procedure.  
  `-- chip.txt` consisting of all the TF-related chip-seq file names.  
  `--  *.gz` files consisting of each TF-related chip-seq bed file.  
  `**  `to save the disk space, we recommend you to create a symbolic link to all the bigwig files, as shown below:  
GM12878:  
  <pre><code>$ ln -s ../SRR891269.forward.1x2.bw SRR891269.forward.1x2.bw
$ ln -s ../wgEncodeDukeMapabilityUniqueness35bp.bigWig wgEncodeDukeMapabilityUniqueness35bp.bigWig
</code></pre>  
K562:  
  <pre><code>$ ln -s ../K562_rep2.bw K562_rep2.bw
$ ln -s ../wgEncodeDukeMapabilityUniqueness35bp.bigWig wgEncodeDukeMapabilityUniqueness35bp.bigWig
</code></pre>  
H1ESC:  
  <pre><code>$ ln -s ../SRR4011946.forward.bw SRR4011946.forward.bw
$ ln -s ../wgEncodeDukeMapabilityUniqueness35bp.bigWig wgEncodeDukeMapabilityUniqueness35bp.bigWig
</code></pre>  

For each cell line, we trained a separate model, right now we collected three cell lines data due to the available open datasets: GM12878, K562 and H1ESC. So there are three subfolders in the Data folder.   

Correspondingly, we also have three pre-trained models in the model folder. Due to the space limit, we provide a link for you to download if you would like to use a pre-trained model. Specifically, wgEncodeDukeMapabilityUniqueness35bp.bigWig can be downloaded from [here](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/) or from the data folder we provide in our [cloud drive](https://drive.google.com/open?id=1R9V53HgpdrjYdFJ04nF_BxjaUfVI7LI1).  
### Train on new dataset 
Our pre-trained model is fully capable of training new TF binding related data. Users can download TF related Chip-seq data and run the bash script add_new_TF.sh, then a new model that can be applied to new TFs is available to use. One example is shown in [add_new_TF.sh](https://github.com/sperfu/scFAN/blob/master/Data/add_new_TF.sh), please refer to that.


## Prediction data
In the prediction phase, we utilized over 2000 cells retrieved from public datasets, which is shown in the main Readme.md and in the paper. All the single cell data we retrieved covered about 21 subtypes, so we collect each subtype data in one folder and predict them once for all.  
The single cell dataset is predicted based on the pre-trained model, so you need to put the single cell data into corresponding folder, we provided the all the single cell peak files in the link, you can download them and unzip them into the corresponding folder(GM12878):  
<pre><code>$ cd GM12878
$ tar zxvf cell_name.tar.gz
</code></pre>   
Also, to save the disk space, we also recommend you to create a symbolic link to all the bigwig files, as shown below:
<pre><code>$ ln -s ../../wgEncodeDukeMapabilityUniqueness35bp.bigWig wgEncodeDukeMapabilityUniqueness35bp.bigWig
</code></pre>  
For the record, due to the size of the disk, for each single cell data, I only uploaded all the needed bed files into GM12878 subfolder, you can simply copy those tar.gz files into either K562 or H1ESC subfolders, for example(K562):  
<pre><code>$ cp *tar.gz ../K562/ && tar zxvf *tar.gz
$ cd gz_files_BJ
</code></pre>   
because all the bed files and bigwig files are the same for each single cell. The only places you need to change are the link. All the single cell bigwig files provided here were aggregated version, I have uploaded all the files in the google drive, you can download from there.   
 

Considering the file size in the repository. We uploaded all the files in google drive, you can get the data from [here](https://drive.google.com/open?id=1R9V53HgpdrjYdFJ04nF_BxjaUfVI7LI1).  

If you have any questions, please feel free to contact me.
Email: fulaiyi@gmail.com or laiyif@uci.edu
