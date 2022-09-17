# Download MCA DGE Data from the link below
# https://figshare.com/articles/dataset/MCA_DGE_Data/5435866
# Use MCA_BatchRemove_dge.zip
# Rename the downloaded folder name to be "mca_rmbatch_data"
# This is the downloaded path
foldpath<-"dge_rmbatch_data/"
# Get file names from the folder
fnames<-list.files(foldpath)

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
# Record the length of final union genes
final_gene_length<-length(final_genename)

# Mitochondria reads check 
# Notice that there is no mitochondria read
# No need to deal with this case in preprocessing
length(grep('^MT-',final_genename)) > 0

# Read the first file as start
expr_new<- read.table(paste0(foldpath,fnames[1]), stringsAsFactors = F)    
# Retain cells with positive expression in at least 1000 genes
index1<-colSums(expr_new > 0) >= 1000
expr_new<-expr_new[,index1]
# Use the union genes fill the gene not appear in this sample
# with 0 (those genes are not detected or have been filtered)
expr_add<-matrix(0,final_gene_length-nrow(expr_new),ncol(expr_new))
rownames(expr_add)<-final_genename[which(!final_genename%in%rownames(expr_new))]
colnames(expr_add)<-colnames(expr_new)
expr_combine<-rbind(expr_new,expr_add)
# Follow the order if final genes
expr<-expr_combine[final_genename,]
# Save the label
label<-rep(fnames[1],ncol(expr_new))
print(ncol(expr))
# The same process is conducted for all files
.<-lapply(2:length(fnames),function(i){
  print(i)
  expr_new<- read.table(paste0(foldpath,fnames[i]), stringsAsFactors = F)    
  index1<-colSums(expr_new > 0) >= 1000 # This cutoff can be changed to be 500
  if(sum(index1)==1){
    expr_combine<-  rep(0,final_gene_length)
    names(expr_combine)<-final_genename
    expr_combine[rownames(expr_new)]<-expr_new[,index1]
    label<<-c(label,rep(fnames[i],1))
    expr<<-cbind(expr,expr_combine)
  }else{
    expr_new<-expr_new[,index1]
    expr_add<-matrix(0,final_gene_length-nrow(expr_new),ncol(expr_new))
    rownames(expr_add)<-final_genename[which(!final_genename%in%rownames(expr_new))]
    colnames(expr_add)<-colnames(expr_new)
    expr_combine<-rbind(expr_new,expr_add)
    expr_combine<-expr_combine[final_genename,]
    label<<-c(label,rep(fnames[i],ncol(expr_new)))
    expr<<-cbind(expr,expr_combine)
  }
  print(ncol(expr))
  1
})


## filter genes
# Retain genes that are expressed in at least 0.1% of cells
index3<-rowSums(expr > 0)/ncol(expr)>= 0.001
expr <- expr[index3,]

# Refine the label
fnames<-label
# Remove the ".txt.gz" at the end of file names
fnames<-substr(fnames,1, nchar(fnames)-7 )
# unique(substr(fnames,nchar(fnames)-9, nchar(fnames) ))
# Remove the ".batch_dge" at the end of file names
fnames<-substr(fnames,1, nchar(fnames)-10 )
# length(unique(fnames))
# unique(substr(fnames,nchar(fnames)-2, nchar(fnames) ))
# Remove the "_rm" at the end of file names
fnames<-substr(fnames,1, nchar(fnames)-3 )

# Remove the count at end of file names
new_fnames<-fnames
uni_fnames<-unique(fnames)
# unique(substr(uni_fnames,nchar(uni_fnames), nchar(uni_fnames)))
.<-sapply(1:length(uni_fnames), function(i){
  index_i<-which(fnames == uni_fnames[i])  
  if(substr(uni_fnames[i],nchar(uni_fnames[i]), nchar(uni_fnames[i]))%in% c("1","2","3","4","5","6")){
    new_fnames[index_i]<<- substr(uni_fnames[i],1, nchar(uni_fnames[i])-1)   
  }
  1
})

# Use abbreviation to save the label 
# Give some long names abbreviation

new_fnames[which(new_fnames == "BoneMarrowcKit")]<-"BMcKit"
new_fnames[which(new_fnames == "BoneMarrow")]<-"BM"
new_fnames[which(new_fnames == "CJ7.EB14.Ezh2.")]<-"CJ7.Ezh2.1"
new_fnames[which(new_fnames == "CJ7.EB14.WT.")]<-"CJ7.WT"
new_fnames[which(new_fnames == "EmbryonicMesenchymeE14.")]<-"E.Mes"
new_fnames[which(new_fnames == "EmbryonicStemCells" )]<-"E.SC"
new_fnames[which(new_fnames == "FetalBrain")]<-"F.Brain"
new_fnames[which(new_fnames == "FetalFemaleGonad")]<-"F.Fgonad"
new_fnames[which(new_fnames == "FetalIntestine")]<-"F.Int"
new_fnames[which(new_fnames == "FetalKidney" )]<-"F.Kid"
new_fnames[which(new_fnames == "FetalLiverE14.")]<-"F.Liv"
new_fnames[which(new_fnames == "FetalLung")]<-"F.Lun"
new_fnames[which(new_fnames == "FetalMaleGonad")]<-"F.Mgonad"
new_fnames[which(new_fnames == "FetalPancreas" )]<-"F.Pan"
new_fnames[which(new_fnames == "FetalStomach" )]<-"F.Sto"
new_fnames[which(new_fnames ==
    "MammaryGland.Involution.CD45." )]<-"M.I.CD45"
new_fnames[which(new_fnames == "MammaryGland.Involution" )]<-"M.I"
new_fnames[which(new_fnames == "MammaryGland.Pregnancy" )]<-"M.P"
new_fnames[which(new_fnames == 
    "MammaryGland.Virgin.CD45."  )]<-"M.V.CD45"
new_fnames[which(new_fnames == "MammaryGland.Virgin" )]<-"M.V"
new_fnames[which(new_fnames == "MesenchymalStemCells" )]<-"Mes.SC"
new_fnames[which(new_fnames == 
    "MesenchymalStemCellsPrimary" )]<-"Mes.SCP"
new_fnames[which(new_fnames == "NeonatalCalvaria" )]<-"N.Cal"
new_fnames[which(new_fnames == "NeonatalHeart" )]<-"N.Heart"
new_fnames[which(new_fnames == "NeonatalMuscle" )]<-"N.Mus"
new_fnames[which(new_fnames == "NeonatalPancreas" )]<-"N.Pan"
new_fnames[which(new_fnames == "NeonatalRib" )]<-"N.Rib"
new_fnames[which(new_fnames == "NeonatalSkin" )]<-"N.Skin"
new_fnames[which(new_fnames == "NeontalBrain" )]<-"N.Brain"
new_fnames[which(new_fnames == "PeripheralBlood" )]<-"Per.Blood"
new_fnames[which(new_fnames == "Pancreas" )]<-"Pan"
new_fnames[which(new_fnames == "PlacentaE14."  )]<-"P.E14"
new_fnames[which(new_fnames == "SmallIntestine.CD4" )]<-"SI.CD4"
new_fnames[which(new_fnames == "SmallIntestine" )]<-"SI"
new_fnames[which(new_fnames == "TrophoblastStemCells" )]<-"T.SC"

saveRDS(expr,"mouse_tissue_expr.rds")
saveRDS(new_fnames,"mouse_tissue_label.rds")