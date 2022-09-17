# Download HCL DGE Data from the link below
# https://figshare.com/articles/dataset/HCL_DGE_Data/7235471
# Use dge_rmbatch_data.tar
# Rename the downloaded folder name to be "hcl_rmbatch_data"
# This is the downloaded path
foldpath<-"dge_rmbatch_data/"
# Get file names from the folder
fnames<-list.files(foldpath)
# Remove filenames associated with fetal
fnames<-fnames[-which(substr(fnames,1,5)%in%c("Fetal","Neona","Place"))]
# There should be 62 files, where most start with Adult
length(fnames)

# First get the overlap of gene names
lst_genename = vector("list", length(fnames))
# Loop to read every files and save its gene names
.<-sapply(seq_along(fnames), function(i){
  print(i)
  expr <- read.table(paste0(foldpath,fnames[i]), stringsAsFactors = F)
  lst_genename[[i]]<<-rownames(expr)
  1
})

# Get the union genes stored in final_genename (suggested by paper authors)
# We got the following reply from Dr. Guo's Lab 
# Dear Ruzhang,
# The gene A wasn't detected in BM or expressed in less than three cells 
# in BM, so we filtered these genes like gene A. Besides, we recommend to
# use union genes in all samples for those tissue-specific genes.
# Best,
# Lijiang.

final_genename<-lst_genename[[1]]
for ( i in c(1:length(lst_genename)) ){
  final_genename<-union(final_genename,lst_genename[[i]])
  print(length(final_genename))
}

# Mitochondria reads check 
# Notice that there are mitochondria reads
# Need to deal with this case in preprocessing
length(grep('^MT-',final_genename)) > 0

# Read the first file as start
i<-1
expr_new<- read.table(paste0(foldpath,fnames[i]), stringsAsFactors = F)    
# Retain cells with positive expression in at least 500 genes
index1<-colSums(expr_new > 0) >= 500
expr_new<-expr_new[,index1]
# Use the union of genes
# Some genes only appear in final_genename
# but not in current tissue 
# have expression 0 
expr_add<-matrix(0,final_gene_length-nrow(expr_new),ncol(expr_new))
rownames(expr_add)<-final_genename[which(!final_genename%in%rownames(expr_new))]
colnames(expr_add)<-colnames(expr_new)
expr_combine<-rbind(expr_new,expr_add)
expr_combine<-expr_combine[final_genename,]

# record the label
label<-c(rep(fnames[i],ncol(expr_new)))
# record the final expression matrix
expr<-expr_combine
# The same process is conducted for all files
.<-lapply(2:length(fnames),function(i){
  message(i)
  expr_new<- read.table(paste0(foldpath,fnames[i]), stringsAsFactors = F)    
  index1<-colSums(expr_new > 0) >= 500
  expr_new<-expr_new[,index1]
  expr_add<-matrix(0,final_gene_length-nrow(expr_new),ncol(expr_new))
  rownames(expr_add)<-final_genename[which(!final_genename%in%rownames(expr_new))]
  colnames(expr_add)<-colnames(expr_new)
  expr_combine<-rbind(expr_new,expr_add)
  expr_combine<-expr_combine[final_genename,]
  
  label<<-c(label,rep(fnames[i],ncol(expr_new)))
  expr<<-cbind(expr,expr_combine)
  print(ncol(expr))
  1
})

## filter genes
# Retain genes that are expressed in at least 0.1% of cells
NCOLEXPR<-ncol(expr)
index3<-(Rfast::rowsums(expr > 0)/NCOLEXPR>= 0.001)

expr <- expr[index3,]

# Refine the label 
fnames<-label
# Remove the ".rmbatchdge.txt.gz" at the end of file names
unique(substr(fnames,nchar(fnames)-17, nchar(fnames) ))
fnames<-substr(fnames,1,nchar(fnames)-18)
# Check current labels
# unique(fnames)

# Remove "Adult" at the head of some file names
index_adult<-which((substr(fnames,1, 5 ))== "Adult")
# unique(substr(fnames[index_adult],1, 5 ))
fnames[index_adult]<-substr(fnames[index_adult],6, nchar(fnames[index_adult]) )

# Remove the count at end of file names
new_fnames<-substr(fnames,1, nchar(fnames)-1 )
unique(new_fnames)

# Use abbreviation to save the label 
# Give some long names abbreviation
new_fnames[which(new_fnames == "AdrenalGland")]<-"A.Gland"
new_fnames[which(new_fnames == "AscendingColon")]<-"A.Colon"
new_fnames[which(new_fnames == "BoneMarrow")]<-"BM"
new_fnames[which(new_fnames == "Cerebellum")]<-"Cer"
new_fnames[which(new_fnames == "Duodenum")]<-"Duo"
new_fnames[which(new_fnames == "Epityphlon")]<-"Epit"
new_fnames[which(new_fnames == "Esophagus")]<-"Eso"
new_fnames[which(new_fnames == "Fallopiantube")]<-"Fall"
new_fnames[which(new_fnames == "Gallbladder")]<-"Gallbla"
new_fnames[which(new_fnames == "Pancreas")]<-"Pan"
new_fnames[which(new_fnames == "PeripheralBlood")]<-"Per.Blood"
new_fnames[which(new_fnames == "SigmoidColon")]<-"S.Colon"
new_fnames[which(new_fnames == "TemporalLobe")]<-"Temp"
new_fnames[which(new_fnames == "TransverseColon")]<-"Trans.Colon"
new_fnames[which(new_fnames == "ChorionicVillus")]<-"Cho.Vil"
new_fnames[which(new_fnames == "CordBloodCD34P")]<-"C.BloodCD34P"
new_fnames[which(new_fnames == "CordBlood")]<-"C.Blood"

saveRDS(expr,"human_tissue_expr.rds")
saveRDS(new_fnames,"human_tissue_label.rds")
