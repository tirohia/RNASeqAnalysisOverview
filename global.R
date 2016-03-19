# load gene annotations, data and experimental design information from the working directory. 

# load annotation data for GSEA from the file annotations.txt
filename<-"annotations.txt"
annDat <-data.frame(read.delim(filename, skip=0, as.is=TRUE, sep="\t")) #read in data

# load raw read counts from Data.csv
filename<-"data.csv"
exp.l.Dat<-data.frame(read.csv(filename,row.names=1)) #read in data

# load the experimental design from Design.matrix.txt, a csv file with a row for each condition, i.e. group/donor/cell type.  
# then construct the design matrix - maodel.matrix(as.formula(paste0(etc etc))) doesn't work. It should though, why it doesn't, I don't know
filename<-"Design.matrix.txt"
ct<-read.csv(filename,row.names = 1)
design<-data.frame(samples = 1:length(ct[1,]))
dn<-NULL
for (i in 1:length(rownames(ct))){
  if (i == 1 ){
    design<-cbind(design, model.matrix(~0 + factor(as.numeric(ct[i,]))))
    groups<-unique(as.matrix(ct)[i,])
    dn <- c(dn, paste0(rownames(ct)[i],(groups)))
  }
  else {
    interimMatrix<-model.matrix(~ factor(as.numeric(ct[i,])))
    interimMatrix <- interimMatrix[,-1]
    design<-cbind(design, interimMatrix)
    groups<-unique(as.matrix(ct)[i,])
    groups<-groups[-1]
    dn <- c(dn, paste0(rownames(ct)[i],(groups)))
  }
}
design<-design[,-1]
colnames(design) <- dn


## define the list of possible comparisons. Though why it's done here I have no idea. 
#ct<-data.frame(t(ct))
#dn<-unlist(lapply(rownames(ct),function(x){paste0(x,unique(as.numeric(ct[x,])))}))
eg<-expand.grid(dn,dn); eg2<-eg[eg[,1]!=eg[,2],]
eg3<-paste(eg2[,1],"-",eg2[,2])
eg4<-c("all possible", eg3, "bespoke")
