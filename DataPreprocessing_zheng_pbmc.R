# Download the sorted pbmc dataset from the link below
# https://support.10xgenomics.com/single-cell-gene-expression/datasets
# Use the dataset under Single Cell 3' Paper: Zheng et al. 2017 
# Use the following cell types
#(1) CD19+ B Cells
#(2) CD14+ Monocytes
#(3) CD34+ Cells
#(4) CD56+ Natural Killer Cells
#(5) CD8+ Cytotoxic T cells
#(6) CD4+/CD45RO+ Memory T Cells
#(7) CD8+/CD45RA+ Naive Cytotoxic T Cells
#(8) CD4+/CD45RA+/CD25- Naive T cells
#(9) CD4+/CD25+ Regulatory T Cells
# Save all the unzip folders under the same upper folder
# named as "Filtered_Folders"
library(Matrix)
# Load all folder names
file_names<-list.files("Filtered_Folders/")
# Check if there are non-overlap genes: 0: all are the same
#for (i in 1:8){
#  gene <- readLines(paste0("Filtered_Folders/",file_names[i], "/hg19/genes.tsv"))
#  geneplus <- readLines(paste0("Filtered_Folders/",file_names[i+1], "/hg19/genes.tsv"))
#   print(sum(gene!=geneplus))
#}
# All the genes are identical

# Define the cell type names
unique_cell_name<-c("B cells","Monocytes","CD34+","NK cells","Cytotoxic","CD4+ Memory","Naive Cytotoxic","CD4+ Naive T","Regulatory")
# Read the first file 
gene <- readLines(paste0("Filtered_Folders/",file_names[1], "/hg19/genes.tsv"))
expr <- readMM(paste0("Filtered_Folders/",file_names[1], "/hg19/matrix.mtx"))
barcode<- readLines(paste0("Filtered_Folders/",file_names[1], "/hg19/barcodes.tsv"))
rownames(expr)<-gene
colnames(expr)<-barcode
cell_label<-rep(unique_cell_name[1],ncol(expr))
# Loop to read all the other files
for ( i in 2:9){
  barcode<- readLines(paste0("Filtered_Folders/",file_names[i], "/hg19/barcodes.tsv"))
  expr_tmp <- readMM(paste0("Filtered_Folders/",file_names[i], "/hg19/matrix.mtx"))
  rownames(expr_tmp)<-gene
  colnames(expr_tmp)<-barcode
  cell_label<-c(cell_label,rep(unique_cell_name[i],ncol(expr_tmp)))
  expr<-cbind(expr,expr_tmp)
  print(ncol(expr) == length(cell_label))
}

# Mitochondria reads check 
# Notice that there is no mitochondria read
# No need to deal with this case in preprocessing
length(grep('^MT-',gene)) > 0

# Retain cells with positive expression in at least 500 genes
index1<-which(colSums(expr > 0) >= 500)
expr <- expr[,index1]
cell_label<-cell_label[index1]

# Retain genes that are expressed in at least 0.1% of cells.
index2<-rowSums(expr > 0)/ncol(expr) >= 0.001
expr <- expr[index2,]

saveRDS(expr,"zheng_cells.rds")
saveRDS(cell_label,"zheng_label.rds")

