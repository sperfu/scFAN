library(RcppCNPy)
library(tsne)
library(ggplot2)
setEPS()
## load files
print("loading data...")
load_tsne_result = TRUE
if (load_tsne_result == TRUE) {
    tsne_result = read.csv('Tsne_matrix_Corces_save',header=TRUE)[,2:3]
    all_label = read.table('label_Corces_dataset.txt',stringsAsFactors=FALSE)
    tsne_result = cbind(tsne_result,all_label)
    names(tsne_result) = c("tSNE1","tSNE2","cell_type")
    print("draw ori plot...")
    filter_index = read.table('all_filtered_cell_index.txt') ## filter cells
    p = ggplot(tsne_result[filter_index[,1]+1,],aes(x=tSNE1, y=tSNE2, colour=cell_type)) + theme_classic() + geom_point()
    postscript("tsne_Corces_dataset_ori.eps")
    print(p)
    } else {
            all_data_GM = npyLoad("TF_activity_score_Corces_pretrained_using_GM.npy")
            all_data_K562 = npyLoad("TF_activity_score_Corces_pretrained_using_K562.npy")
            all_data_H1 = npyLoad("TF_activity_score_Corces_pretrained_using_H1.npy")
            
            all_label = read.table('label_Corces_dataset.txt',stringsAsFactors=FALSE)
            
            filter_index = read.table('all_filtered_cell_index.txt') ## filter cells
            
            print("process data...")
            
            merge_three = cbind(all_data_GM,all_data_K562,all_data_H1)
            
            scaledData = merge_three
            pcs = prcomp(scaledData, scale. = F, center = F)
            pcs$tSNE = tsne(pcs$x[, 1:30])
            pcs$tSNEProj = data.frame(pcs$tSNE)
            pcs$tSNEProj2 = cbind(pcs$tSNEProj,all_label)
            names(pcs$tSNEProj2) = c("tSNE1", "tSNE2","cell_type")
            print("draw plot...")
            setEPS()
            p=ggplot(pcs$tSNEProj2[filter_index[,1]+1,],aes(x=tSNE1, y=tSNE2, colour=cell_type)) + theme_classic() + geom_point()
            postscript("tsne_Corces_dataset.eps")
            print(p)
            setEPS()
            write.csv(pcs$tSNEProj2[,1:2],'Tsne_matrix_Corces')
    }
