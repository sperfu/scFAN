# Usage
This document shows an example of using calculated TF activity score to cluster cells, it contains two subfolders,
which is the results of Corces dataset and PBMCs dataset respectively.

When you enter each subfolder, just run the command: 
<pre><code>$ Rscript draw_tsne.R 
</code></pre>  
Then you will get Tsne plot of clustering these cell. We loaded the precalculated Tsne results to let you see the result in
our paper, you could also perform a new calculation by uncommenting some codes in the draw_tsne.R script. But due to the initialization
difference in the tsne function, the final plot may present a little different. 
