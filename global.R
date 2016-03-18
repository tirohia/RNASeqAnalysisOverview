# load gene annotations, data and experimental design information from the working directory. 

# load annotation data for GSEA from the file annotations.txt
filename<-"annotations.txt"
annDat <-data.frame(read.delim(filename, skip=0, as.is=TRUE, sep="\t")) #read in data

# load raw read counts from Data.csv
filename<-"data.csv"
exp.l.Dat<-data.frame(read.csv(filename,row.names=1)) #read in data

# load the experimental design from Design.matrix.txt, a csv file with a row for each condition, i.e. group/donor/cell type.  
filename<-"Design.matrix.txt"
ct<-read.csv(filename,row.names = 1)
formula<-as.formula(paste(c("~ 0",rownames(ct)),collapse=" + "))
ct<-data.frame(t(ct))
design<-model.matrix(formula, ct)


## define the list of possible comparisons. Though why it's done here I have no idea. 
ct<-data.frame(t(ct))
dn<-unlist(lapply(rownames(ct),function(x){paste0(x,unique(as.numeric(ct[x,])))}))
eg<-expand.grid(dn,dn); eg2<-eg[eg[,1]!=eg[,2],]
eg3<-paste(eg2[,1],"-",eg2[,2])
eg4<-c("all possible", eg3, "bespoke")
