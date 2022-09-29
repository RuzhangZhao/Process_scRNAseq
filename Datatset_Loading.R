## CITEseq Dataset1 
## bmcite
## download link:
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128639
library(SeuratData)
dataset_name = "bmcite"
bm <- LoadData(ds = "bmcite")
expr<-bm@assays$RNA@counts
adt<-bm@assays$ADT@data

## CITEseq Dataset2
## cbmc8k_cite
## download link: 
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866
expr<-read.table(file =file.path("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz") , header = T, row.names=1,sep=",", as.is=T)
adt<-read.table(file =file.path("GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz") , header = T, row.names=1,sep=",", as.is=T)



## CITEseq Dataset3
## pbmc_cite
## download link: 
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866
expr<-read.table(file =file.path("GSE100866_PBMC_vs_flow_10X-RNA_umi.csv.gz") , header = T, row.names=1, sep=",", as.is=T)
adt<-read.table(file =file.path("GSE100866_PBMC_vs_flow_10X-ADT_umi.csv.gz") , header = T, row.names=1, sep=",", as.is=T)

## CITEseq Dataset4
## FLiver_cite
## download link: 
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166895
RNA_file<-"GSE166895_postQC_mRNAraw_FL-FBM-CB.csv.gz"
expr<-t(read.table(file =RNA_file, header = T, row.names=1,sep=",", as.is=T))
colnames(expr)<-sub("\\-.*", "",colnames(expr))
colnames(expr)<-make.unique(colnames(expr),"DELETE")
expr<-expr[,-which(grepl("DELETE",colnames(expr)))]

ADT_file<-"GSE166895_overlap_ADTonmRNAraw_FL-FBM-CB.csv.gz"
adt<-t(read.table(file =ADT_file, header = T, row.names=1,sep=",", as.is=T))
colnames(adt)<-sub("\\-.*", "",colnames(adt))
colnames(adt)<-make.unique(colnames(adt),"DELETE")
adt<-adt[,-which(grepl("DELETE",colnames(adt)))]

expr<-expr[,which(colnames(expr)%in%colnames(adt))]
adt<-adt[,colnames(expr)]


## CITEseq Dataset5
## FBM_cite
## download link: 
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166895
RNA_file<-"/dcs04/hongkai/data/rzhao/pbmc/cellatlas/GSE166895_postQC_mRNAraw_FBM-MNCs.csv.gz"
expr<-t(read.table(file =RNA_file, header = T, row.names=1,sep=",", as.is=T))
ADT_file<-"/dcs04/hongkai/data/rzhao/pbmc/cellatlas/GSE166895_postQC_ADT_raw_FBM-MNCs.csv.gz"
adt<-t(read.table(file =ADT_file, header = T, row.names=1,sep=",", as.is=T))
colnames(adt)<-paste0(colnames(adt),'-1')
expr<-expr[,colnames(expr)%in%colnames(adt)]
adt<-adt[,colnames(expr)]

## CITEseq Dataset6
## seurat_cite
## download link: 
## https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
library(SeuratDisk)
reference<-LoadH5Seurat("pbmc_multimodal.h5seurat")
expr<-reference@assays$SCT@counts
adt<-reference@assays$ADT@counts

## CITEseq Dataset7
## sucovid_cite
## download from:
## https://doi.org/10.1016/j.cell.2020.10.037
## The count and adt file should be processed from the 
expr <- readRDS("count.rds")
pro<-readRDS("adt.rds")

## cell sorting Dataset1
## duo4_pbmc
library(DuoClustering2018)
sce_names<-c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
sn<-sce_names[1]
sce<-do.call(paste0("sce_full_",sn),list())
duocount<-sce@assays$data@listData$counts
duolabel<-sce$phenoid
saveRDS(duocount,"duo4_expr.rds")
saveRDS(duolabel,"duo4_label.rds")

## cell sorting Dataset2
## duo4un_pbmc
library(DuoClustering2018)
sce_names<-c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
sn<-sce_names[2]
sce<-do.call(paste0("sce_full_",sn),list())
duocount<-sce@assays$data@listData$counts
duolabel<-sce$phenoid
saveRDS(duocount,"duo4un_expr.rds")
saveRDS(duolabel,"duo4un_label.rds")

## cell sorting Dataset3
## duo8_pbmc
library(DuoClustering2018)
sce_names<-c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
sn<-sce_names[3]
sce<-do.call(paste0("sce_full_",sn),list())
duocount<-sce@assays$data@listData$counts
duolabel<-sce$phenoid
saveRDS(duocount,"duo8_expr.rds")
saveRDS(duolabel,"duo8_label.rds")

## cell sorting Dataset4
## GBM_sd
## download link:
## https://doi.org/10.1016/j.celrep.2017.10.030
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465
expr<-read.table("GBM_raw_gene_counts.csv",header = T)
cell_label<-read.table("GBM_metadata.csv",header = T)$Location

## cell sorting Dataset5
## homo_tissue
## download link:
## https://www.cell.com/cell/fulltext/S0092-8674(18)30116-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418301168%3Fshowall%3Dtrue
## https://github.com/RuzhangZhao/Process_scRNAseq/blob/main/DataPreprocessing_homo_tissue.R

## cell sorting Dataset6
## mus_tissue
## download link:
## https://www.cell.com/cell/fulltext/S0092-8674(18)30116-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418301168%3Fshowall%3Dtrue
## https://github.com/RuzhangZhao/Process_scRNAseq/blob/main/DataPreprocessing_mus_tissue.R

## cell sorting Dataset7
## mus_tissue_large
## download link:
## https://www.cell.com/cell/fulltext/S0092-8674(18)30116-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418301168%3Fshowall%3Dtrue
## https://github.com/RuzhangZhao/Process_scRNAseq/blob/main/DataPreprocessing_mus_tissue.R
## Change the cutoff '1000' in 44 & 61 to be '500'

## cell sorting Dataset8
## zheng_pbmc
## download link:
## https://www.10xgenomics.com/resources/datasets?query=&page=1&configure%5Bfacets%5D%5B0%5D=chemistryVersionAndThroughput&configure%5Bfacets%5D%5B1%5D=pipeline.version&configure%5BhitsPerPage%5D=500&menu%5Bproducts.name%5D=Single%20Cell%20Gene%20Expression
## Refer to Single Cell 3' Paper: Zheng et al. 2017 (v1 Chemistry)
## Also refer to https://github.com/RuzhangZhao/Process_scRNAseq/blob/main/DataPreprocessing_zheng_pbmc.R

## multiomeATAC Dataset1
## snareseq
## download link:
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126074

library(Signac)
library(Seurat)
# load processed data matrices for each assay
rna <- Read10X("SNAREseq/GSE126074_AdBrainCortex_rna/", gene.column = 1)
atac <- Read10X("SNAREseq/GSE126074_AdBrainCortex_atac/", gene.column = 1)

# create a Seurat object and add the assays
chrom_assay <- CreateChromatinAssay(
    counts = atac,
    sep = c(":", "-"),
    min.cells = 10,
    min.features = 1500
)
snare<-CreateSeuratObject(chrom_assay)
snare <- FindTopFeatures(snare, min.cutoff = 'q0')
snare <- RunTFIDF(snare)
snare <- RunSVD(snare)
snare_lsi<-snare@reductions$lsi@cell.embeddings[,2:30]
snare_lsi<-t(snare_lsi)
saveRDS(snare_lsi,"SNAREseq/snare_lsi.rds")


## multiomeATAC Dataset2
## pbmc3k_multi
## download link:
## https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0

## multiomeATAC Dataset3
## pbmc10k_multi
## download link:
## https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0

## multiomeATAC Dataset4
## homo_brain3k
## download link:
## https://www.10xgenomics.com/resources/datasets/frozen-human-healthy-brain-tissue-3-k-1-standard-2-0-0

## multiomeATAC Dataset5
## mus_brain5k
## download link:
## https://www.10xgenomics.com/resources/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0

## multiomeATAC Dataset6
## lymphoma
## download link:
## https://www.10xgenomics.com/resources/datasets/fresh-frozen-lymph-node-with-b-cell-lymphoma-14-k-sorted-nuclei-1-standard-2-0-0

## For Dataset 2,3,4,5,6, change the dataname below, where they share the same template.  
dataname<-"pbmc3k"
library(Seurat)
library(Matrix)
library(SeuratObject)

folder_path = paste0(dataname,"/")
mtx <- paste0(folder_path,"filtered_feature_bc_matrix/matrix.mtx.gz")
cells <- paste0(folder_path,"filtered_feature_bc_matrix/barcodes.tsv.gz")
features <- paste0(folder_path,"filtered_feature_bc_matrix/features.tsv.gz")
MultiomeATAC<-Matrix::readMM(mtx)
MultiomeATAC_barcode<-fread(cells,header=FALSE)
MultiomeATAC_features<-fread(features,header=FALSE)
MultiomeATAC_rna<-MultiomeATAC[which(MultiomeATAC_features$V3 == "Gene Expression"),]
MultiomeATAC_rna_meta<-MultiomeATAC_features[which(MultiomeATAC_features$V3 == "Gene Expression"),]
rownames(MultiomeATAC_rna)<-MultiomeATAC_rna_meta$V2
colnames(MultiomeATAC_rna)<-MultiomeATAC_barcode$V1
MultiomeATAC_atac<-MultiomeATAC[which(MultiomeATAC_features$V3 == "Peaks"),]
MultiomeATAC_atac_meta<-MultiomeATAC_features[which(MultiomeATAC_features$V3 == "Peaks"),]
rownames(MultiomeATAC_atac)<-MultiomeATAC_atac_meta$V2
colnames(MultiomeATAC_atac)<-MultiomeATAC_barcode$V1
library(Signac)
library(Signac)
library(Seurat)
chrom_assay <- CreateChromatinAssay(
    counts = MultiomeATAC_atac,
    sep = c(":", "-"),
    min.cells = 10,
    min.features = 1500
)
MultiomeATAC<-CreateSeuratObject(chrom_assay)
MultiomeATAC <- RunTFIDF(MultiomeATAC)
MultiomeATAC <- FindTopFeatures(MultiomeATAC, min.cutoff = 'q0')
MultiomeATAC <- RunSVD(MultiomeATAC)
MultiomeATAC_lsi<-MultiomeATAC@reductions$lsi@cell.embeddings
MultiomeATAC_lsi<-MultiomeATAC_lsi[,2:30]
MultiomeATAC_lsi<-t(MultiomeATAC_lsi)
pro<-MultiomeATAC_lsi
saveRDS(MultiomeATAC_rna,paste0(folder_path,dataname,"_rna.mat.rds"))
saveRDS(pro,paste0(folder_path,dataname,"_lsi.rds"))
