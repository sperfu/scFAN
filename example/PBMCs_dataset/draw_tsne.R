#library(RcppCNPy)
library(tsne)
library(ggplot2)
setEPS()
## load files
print("loading data...")
load_tsne_result = TRUE
if (load_tsne_result == TRUE) {
    tsne_result = read.csv('Tsne_matrix_PBMCs_save',header=TRUE)[,2:3]
    all_label = read.table('label_PBMCs_dataset.txt',stringsAsFactors=FALSE)
    tsne_result = cbind(tsne_result,all_label)
    names(tsne_result) = c("tSNE1","tSNE2","cell_type")
    print("draw ori plot...")
    p = ggplot(tsne_result,aes(x=tSNE1, y=tSNE2, colour=cell_type)) + theme_classic() + geom_point()
    postscript("tsne_PBMCs_dataset.eps")
    print(p)
    } else {
       load("TF_activity_score_PBMCs_pretrained_using_GM.Rdata")
       load("TF_activity_score_PBMCs_pretrained_using_K562.Rdata")
       load("TF_activity_score_PBMCs_pretrained_using_H1.Rdata")
       
       all_label = read.table('label_PBMCs_dataset.txt',stringsAsFactors=FALSE)

       print("process data...")
       merge_three = cbind(TF_activity_score_PBMCs_pretrained_using_GM,TF_activity_score_PBMCs_pretrained_using_K562,TF_activity_score_PBMCs_pretrained_using_H1)

       scaledData = merge_three
       pcs = prcomp(scaledData, scale. = F, center = F)
       pcs$tSNE = tsne(pcs$x[, 1:20])
       pcs$tSNEProj = data.frame(pcs$tSNE)
       pcs$tSNEProj2 = cbind(pcs$tSNEProj,all_label)
       names(pcs$tSNEProj2) = c("tSNE1", "tSNE2","cell_type")
       print("draw plot...")
       setEPS()
       p=ggplot(pcs$tSNEProj2,aes(x=tSNE1, y=tSNE2, colour=cell_type)) + theme_classic() + geom_point()
       postscript("tsne_PBMCs_dataset.eps")
       print(p)
       write.csv(pcs$tSNEProj2[,1:2],'Tsne_matrix_PBMCs')
    }

