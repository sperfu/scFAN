# Data
The Data folder consists of most of the training bulk data and all the the single cell data.
## Training data
As described in the main text, in the training phase, we trained three models based on different TFs combination, which is generally from three cell lines: GM12878, K562 and H1ESC.   
There are these files listed below in the training data:  
  `-- bigwig.txt` consisting of all the epigenomic file you use in the training procedure.  
  `-- chip.txt` consisting of all the TF-related chip-seq file names.  
  `--  *.gz` files consisting of each TF-related chip-seq bed file.  
  `**  `to save the disk space, we recommend you to create a symbolic link to all the bigwig files, as shown below:
  <pre><code>$ ln -s ../SRR891269.forward.1x2.bw SRR891269.forward.1x2.bw
$ ln -s ../wgEncodeDukeMapabilityUniqueness35bp.bigWig wgEncodeDukeMapabilityUniqueness35bp.bigWig
</code></pre>  

For each cell line, we trained a separate model, right now we collected three cell lines data due to the available open datasets: GM12878, K562 and H1ESC. So there are three subfolders in the Data folder.  
Correspondingly, we also have three pre-trained models in the model folder. Due to the space limit, we also provide a link for you to download if you would like to use a pre-trained model.  

## Prediction data
In the prediction phase, we utilized over 2000 cells retrieved from public datasets, which is shown in the main Readme.md and in the paper. All the single cell data we retrieved covered about 21 subtypes, so we collect each subtype data in one folder and predict them once for all.  
The single cell dataset is predicted based on the pre-trained model, so you need to put the single cell data into corresponding folder, we provided the all the single cell peak files in the link, you can download them and unzip them into the corresponding folder:  
<pre><code>$ cd GM12878
$ mkdir gz_files_GM12878
$ cd gz_files_GM12878
$ wget link
$ tar zxvf cell_name.tar.gz
</code></pre>   
Also, to save the disk space, we also recommend you to create a symbolic link to all the bigwig files, as shown below:
<pre><code>$ ln -s ../../SRR891269.forward.1x2.bw SRR891269.forward.1x2.bw
$ ln -s ../../wgEncodeDukeMapabilityUniqueness35bp.bigWig wgEncodeDukeMapabilityUniqueness35bp.bigWig
</code></pre>  
Considering the file size in the repository. We uploaded all the files in google drive, please contact us for the data link.

Email: laiyif@uci.edu
