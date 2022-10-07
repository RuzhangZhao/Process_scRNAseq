## sucovid_cite
## download from:
## https://doi.org/10.1016/j.cell.2020.10.037
## Author: Ruzhang Zhao<rzhao@jhu.edu>, Zhicheng Ji <zhicheng.ji@duke.edu>
# Process scRNAseq 
library(Matrix)
library(hdf5r)

file <- hdf5r::h5file(filename = 'for_sharing_RAW_gex_data_w_annos.h5ad', mode = 'r')
data <- file[['X']][['data']][]
indices <- file[['X']][['indices']][]
indptr <- file[['X']][['indptr']][]
obs <- file[['obs']][['_index']][]
var <- file[['var']][['_index']][]

x <- sparseMatrix(i=indices + 1, p=indptr, x=data)

colnames(x) <- obs
rownames(x) <- var[-length(var)]

saveRDS(x,file='count.rds')


### Process ADT(surface protein)

library(data.table)
af <- list.files('pro/')
m <- sapply(af,function(f) {
    d <- fread(paste0('pro/',f),data.table = F)  
    rownames(d) <- d[,1]
    d <- as.matrix(d[,-1])
    rownames(d) <- paste0(gsub('heathlab_dc_9_17_pbmc_pro_|.txt.gz','',f),':',sub(':.*','',rownames(d)))
    d
})
m <- do.call(rbind,m)
m <- t(m)
saveRDS(m,file='pro.rds')



expr <- readRDS("count.rds")
pro<-readRDS("pro.rds")
intersect_name<-intersect(colnames(expr),colnames(pro))
expr<-expr[,intersect_name]
pro<-pro[,intersect_name]

# do cell filteration step 1
index1<-which(colSums(expr > 0) >= 1000)
expr <- expr[,index1]
pro <- pro[,index1]  

# step 2
if (length(grep('^MT-',rownames(expr))) > 0 ){
    index2<-colSums(expr[grep('^MT-',rownames(expr)),])/colSums(expr) <= 0.1
    expr <- expr[,index2]
    pro <- pro[,index2]
}
# step 3
expr <- expr[rowSums(expr > 0)/ncol(expr) >= 0.001,]
saveRDS(expr,"count.rds")
saveRDS(pro,"pro.rds")

