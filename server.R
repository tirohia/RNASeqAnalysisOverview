library(shiny)
library(edgeR)
library(matrixStats)
library(gplots)
library(ggplot2)
library(colorRamps)
library(rms)
library(fastICA)
library(limma)
library(Biobase)
source('heatmap-mik.R')

options(contrasts=c("contr.treatment","contr.treatment")) # to allow Design to work properly with ordered factors


# Define server logic
shinyServer(function(input, output) {
  genes <- reactive({
    naList<-!geneLists[input$list,]==""
    g<-as.vector(unlist(geneLists[input$list,!geneLists[input$list,]==""]))
    
    if(input$list2=="Differential_expression"){
      g<-gg()
    }

    if(input$list2=="user"){
      g<-strsplit(gsub(":"," ",gsub(","," ",gsub(";"," ",input$genes,fixed=T),
                                   fixed=T),fixed=T)," ")[[1]]
      g<-g[g!=""]
      g<-toupper(g)
    }
   return(g)
  })
	

eDat <- reactive({
    g<-genes() # if not looking at correlated genes just get the gene list

    if(input$list=="gcor"){ # if looking at correlated genes, just get the list of these genes
      pCor<-cbind(expDat()[match(unlist(annDat[match(g,annDat[,3]),1]),rownames(expDat())),])
      if(nrow(pCor)==1) pCor<-t(pCor)
      corDat<-cor(pCor,t(expDat()))
      p<-rownames(expDat())[order(abs(corDat),decreasing=T)[sort(abs(corDat),decreasing=T)>=input$corThres]]
      g<-unlist(annDat[match(p,annDat[,1]),3])
    }

    # then, whether looking at correlated genes or not ...
    mt<-as.vector(na.omit(match(g,annDat[,3])))
    g<-g[!is.na(match(g,annDat[,3]))]
    p<-unlist(annDat[mt,1])
    x<-expDat()[match(p,rownames(expDat())),]
    xn<-names(x)

    if(input$list=="gsingle"){  x<-matrix(x,nrow = 1,ncol = length(x));rownames(x)<-g; colnames(x)<-xn}

    if(input$list!="gsingle"){rownames(x)<-g; z<-x[apply(is.na(x),1,sum)==0,]}
    if(input$list=="gsingle"){x[is.na(x)]<-0; z<-x}

    return(z) # the function eDat returns a matrix of data for all tumours but only the genes you are focussig on
  })
  

  
output$diag.expl <- renderText(
if(input$expl=="No"){
"All plots reflect any trimming or normalisation you do in this tab. Operations are done in order they appear in the sidebar, from top to bottom."}
else
if(input$expl=="Yes"){
"All plots reflect any trimming or normalisation you do in this tab. Operations are done in order they appear in the sidebar, from top to bottom. Trimming and Normalisation affect both limma and EdgeR data processing but log-transformation afects only data that is used for limma, not EdgeR. The normalisation methods available are: (1) QUantile normalisation uses the normalizeQuantiles function from limmma, mapping all columns to the same distribution. The remaining methods ceom from EdgeR. (2) The method TMM is the weighted trimmed mean of M-values (to the reference) proposed by Robinson and Oshlack (2010), where the weights are from the delta method on Binomial data. If refColumn is unspecified, the library whose upper quartile is closest to the mean upper quartile is used. (3) The method RLE is the scaling factor method proposed by Anders and Huber (2010). A median library is calculated from the geometric mean of all columns and the median ratio of each sample to the median library is taken as the scale factor. (4) The  method upperquartile is the upper-quartile normalization method of Bullard et al (2010), in which the scale factors are calculated from the 75% quantile of the counts for each library, after removing genes which are zero in all libraries. (5) For method none the normalization factors are set to 1. Normalization factors are adjusted to multiply to 1 so the effective library size is the original library size multiplied by the scaling factor. Rows that have zero counts for all columns are trimmed before normalization factors are computed. Note that some of these methods log2 transform the data while others do not."
}
)


output$dif.expl <- renderText(
"Tab-delimited 'Data.txt' file and 'Design.matrix.txt' files will have been read in from this dashboard directory. It assumes your expression data (RNAseq or microarray) has already been normalised in an appropriate way, however visulaistaion and some top and bottom trimming are possible in the 'Diagnostic Plots and Trimming' tab before you perform differential expression analysis")

##################################################################################
## individual reactive functions that will be used to build up the DE gene list ##
##################################################################################

# trimming - reads in originally loaded data in golbally defined object e<-exp.l.Dat
trim.e <- reactive({
  e<-exp.l.Dat[rowMeans(exp.l.Dat)>=input$trimmer.min,]
  e<-e[rowMeans(e)<=input$trimmer.max,]
  return(e)
})

#calculate normalising factors for EdgeR methods 
norm.factor <- reactive({
e<-trim.e()
if(input$norm=="t"){y<-calcNormFactors(e, method="TMM")}
if(input$norm=="r"){y<-calcNormFactors(e, method="RLE")}
if(input$norm=="uq"){y<-calcNormFactors(e, method="upperquartile")}
return(y)
})

#normalise trimmed data - reads in trim.e()
pre.expDat <- reactive({
if(input$norm=="q"){f<-normalizeQuantiles(trim.e(), ties=TRUE)} # quantie normalisation
if(input$norm=="t" | input$norm=="r" | input$norm=="uq"){f<-trim.e()*norm.factor()} # one of the EdgeR normalisation methods
if(input$norm=="none"){f<-trim.e()} # no normalisation
return(f)
})

# +/- log2 Transform  - reads in pre.expDat()
######### This is the surce of expDat() {trimmed and normalised expression data} that is used extensively later #########
expDat <- reactive({
	if(input$l2=="Yes"){g<-log2(pre.expDat()+1)} else {g<-pre.expDat()}
	  return(g)
})

# Define contrast matrices (or not), reads in glabally defined object design
cm <- reactive({
  if(input$wc=="bespoke"){
    contrast.matrix <- makeContrasts(contrasts=input$bsc, levels=design)	
  }
  if(input$wc!="all possible" && input$wc!="bespoke") {
    print(design)
    print(input$wc)
    contrast.matrix <- makeContrasts(contrasts=input$wc, levels=design)	
  }
  return(contrast.matrix)
})

###############################################
gg <- reactive({
###############################################

######### LIMMA workflow #########
	
if(input$Anal.Type=="Limma"){
eset <-new("ExpressionSet",exprs= as.matrix(expDat()))
fit <- lmFit(eset, design)

if(input$wc=="all possible"){
fit2<-eBayes(fit)
g.out <- rownames(topTable(fit2, number=input$max.number, genelist=fit2$genes, adjust.method=input$mtc, sort.by="F", resort.by=NULL, p.value=input$pSlider, lfc=input$lfcSlider))
}
if(input$wc!="all possible") {
fit <- contrasts.fit(fit, cm())
fit2<-eBayes(fit)
g.out <- rownames(topTable(fit2, number=input$max.number, genelist=fit2$genes, adjust.method=input$mtc, sort.by="p", resort.by=NULL, p.value=input$pSlider, lfc=input$lfcSlider))
		}
}

######### LIMMA voom workflow #########
if(input$Anal.Type=="LV"){
  dge <- DGEList(counts=trim.e()) # don't use pre-normalised data in the Limma voom pipeline of course

  dge <-calcNormFactors(dge, method="TMM") # set the default, this is also what is done if input$norm ="t"
  if(input$norm=="r"){dge <-calcNormFactors(dge, method="RLE")}
  if(input$norm=="uq"){dge <-calcNormFactors(dge, method="upperquartile")}

  v <- voom(dge, design, plot=F) 
  fit <- lmFit(v, design)

  if(input$wc=="all possible"){
    fit2<-eBayes(fit)
    DEgenes<-rownames(topTable(fit2, number=input$max.number, genelist=fit2$genes, adjust.method=input$mtc, sort.by="F", resort.by=NULL, p.value=input$pSlider, lfc=input$lfcSlider))
    if (length(DEgenes > 0)){
      g.out <- DEgenes
    }
    else {
      g.out<-"No genes found to be differentially expressed. Altering required p-values might be neccesary"
    }
  }
  if(input$wc!="all possible") {
    fit <- contrasts.fit(fit, cm())
    fit2<-eBayes(fit)
    DEgenes <- rownames(topTable(fit2, number=input$max.number, genelist=fit2$genes, adjust.method=input$mtc, sort.by="p", resort.by=NULL, p.value=input$pSlider, lfc=input$lfcSlider))
    if (length(DEgenes > 0)){
      g.out <- DEgenes
    }
    else {
      g.out<-"No genes found to be differentially expressed. Altering required p-values might be neccesary"
    }
  }
}
	
######### EdgeR workflow #########
  
if(input$Anal.Type=="EdgeR"){
  counts<-trim.e()
  counts<-counts[rowSums(counts)>0,] #rows of all zeros cause problems for  glmTreat
  y <- DGEList(counts=counts) # don't use pre-normalised data for EdgeR of course
  y<-calcNormFactors(y, method="TMM") # set the default, this is also what is done if input$norm ="t"
  if(input$norm=="r"){y<-calcNormFactors(y, method="RLE")}
  if(input$norm=="uq"){y<-calcNormFactors(y, method="upperquartile")}

  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)

  if(input$wc=="all possible"){
    tr <- glmTreat(fit, contrast=NULL , lfc=input$lfcSlider)
    DEgenes <- rownames(topTags(tr, n=input$max.number, adjust.method="BH", sort.by="PValue", p.value=input$pSlider))
    if (length(DEgenes > 0)){
      g.out <- DEgenes
    }
    else {
      g.out<-"No genes found to be differentially expressed. Altering required p-values might be neccesary"
    }
    
  }
  
  if(input$wc!="all possible") {
    DEgenes <- rownames(topTags(tr, n=input$max.number, adjust.method="BH", sort.by="PValue", p.value=input$pSlider))
    tr <- glmTreat(fit, contrast=cm() , lfc=input$lfcSlider)   	
    if (length(DEgenes > 0)){
      g.out <- DEgenes
    }
    else {
      g.out<-"No genes found to be differentially expressed. Altering required p-values might be neccesary"
    }
  }	 
}
  
return(g.out)	
})


#exectue save file dialogue
observeEvent(input$sv, {	
# generate date-named directory into which all snp output files will be written	
dn<-paste("Saved_",date(),sep="")
		dir.create(file.path(getwd(), dn))
		setwd(file.path(getwd(), dn))
		
write.table(eDat(), paste0("Gene.names_and_Data_",date(),input$fn,".txt"), 
			row.names = T, col.names = T, append = F, 
			quote = F, sep = "\t")
#restore wd
setwd(file.path(getwd()))			
	})


 output$s.c <- renderPlot({ 
   output$s.c <- renderPlot({ 
     if(input$ctype=="Spearman"){
       c.c<-cor(expDat(), use="na.or.complete", method="spearman")}
     else if(input$ctype=="Pearson"){
       c.c<-cor(expDat(), use="na.or.complete", method="pearson")}
     heatmap.2(c.c,trace="none",keysize = 1, key.title = "Correlation", mar=c(15,15), cexRow=0.7, cexCol=0.7, Rowv=T, Colv=T, main="Correlation structure of samples")
     #counts.pca <- prcomp(t(log(e+1)),retx=TRUE)
     #colList<-c("red","red","red","red","blue","blue","blue","blue","green","green","green","green")
     #colList<-c("red","blue","green","orange","red","blue","green","orange","red","blue","green","orange")
     #plot(counts.pca$x[,1],counts.pca$x[,3],pch=19,col =colList, xlab="PC 1",ylab="PC2")
     
   })
    #c.c<-cor(expDat(), use="na.or.complete") 
    #heatmap.2(c.c,trace="none",keysize = 1, key.title = "Correlation", mar=c(15,15), cexRow=0.7, cexCol=0.7, Rowv=T, Colv=T, main="Correlation structure of samples")
})

 output$m.sd <- renderPlot({   
#xlim=range(0:input$xslid2/100*input$xslid/100*max(rowMeans(expDat())))
xlim<-range(0:max(rowMeans(expDat())))
plot(rowMeans(expDat()),(rowSds(as.matrix(expDat()))/rowMeans(expDat())), xlim=xlim, pch=".", col="blue", main="Coefficient of variation relationship to mean", xlab="Gene Mean", ylab="Coefficient of Variation")
})

 output$m.hist <- renderPlot({   
   
 	#PexpDat<-expDat()[rowMeans(expDat())<=input$xslid2/100*input$xslid/100*max(rowMeans(expDat())),]
   PexpDat<-expDat()[rowMeans(expDat())<=max(rowMeans(expDat())),]
plot(hist(rowMeans(PexpDat), breaks=100), main="Histogram of gene means", xlab="Gene Mean", col="blue")
})

  output$gene.cor <- renderPlot({    
nr<-nrow(expDat())
s<-sample(nr, nr/100,replace = FALSE)
c.g<-cor(t(expDat()[s,]))
plot(density(c.g, na.rm =T), col="blue", xlab="Pearson's correlation between genes", main="Correlation structure of a randon 1% of genes"); abline(v=0,lty=2,col="grey")
})

  output$distrs <- renderPlot({ 
nr<-nrow(expDat())
f<-rowSums(expDat()); g<-length(f[f==0])

ff<-unlist(expDat())
s<-sample(1:length(ff),100100,replace = FALSE)
ff<-sort(ff[s])
ff<-ff[1:100000]#trim off top 100 values
par(mfrow=c(3,1))
plot(1:length(ff),sort(ff), pch=".", xlab="linear uniform distribution", col="blue", ylab="sorted data", main=paste("Compare the data to data ranks  - ",round(length(ff[ff==0])/length(ff)*100,2),"% of data = 0  - ",round(g/nr*100,2),"% of genes = 0 in all samples"), cex.main=1.5)
plot(sort(rnorm(100000, mean = 0, sd = 1)),sort(ff), pch=".", xlab="normal distribution", col="blue", ylab="sorted data", main="Compare the data to a normal distribution", cex.main=1.5)
plot(sort(rnbinom(100000,1,0.1)),sort(ff), pch=".", xlab="negative binomial distribution", col="blue", ylab="sorted data", main="Compare the data to a negative binomial distribution", cex.main=1.5)
})


  output$box.pl <- renderPlot({    
  	par(mar=c(20, 5, 5, 5) + 0.1)
boxplot(expDat(), las=3, main="Boxplot of read counts", col="lightblue", ylab="read counts")
})

  output$box.pll <- renderPlot({    
  	par(mar=c(20, 5, 5, 5) + 0.1)
boxplot(expDat() +0.1, log="y", las=3, main="Log10 Boxplot of read counts", col="lightblue", ylab="read counts")
})

  output$trend.linear <- renderPlot({    
  	y<-unlist()
  	yl<-length(y)
plot(1:yl,sort(y))
})

gp <- reactive({
t<-paste(genes(),collapse = "+")   
v_GO<-paste("http://gather.genome.duke.edu/?cmd=report&gene_box=",t,"&tax_id=9606&annot_type=gene_ontology&network=0&homologs=0",sep = "")
glg<-data.frame(read.delim(v_GO))
glg
})

tp <- reactive({
  t<-paste(genes(),collapse = "+") 
  print(paste("http://gather.genome.duke.edu/?cmd=report&gene_box=",t,"&tax_id=9606&annot_type=transfac&network=0&homologs=0",sep = ""))
  v_TF<-paste("http://gather.genome.duke.edu/?cmd=report&gene_box=",t,"&tax_id=9606&annot_type=transfac&network=0&homologs=0",sep = "")
  glt<-data.frame(read.delim(v_TF))
  glt[,2]<-gsub("^..", "", glt[,2],perl=TRUE)
  glt
})

output$t.go <- renderDataTable({
	data.frame(gp())
  })
  
output$t.tf <- renderDataTable({
	data.frame(tp())
  })

?renderPlot
output$expMap <- renderPlot({    
    # Metagene based on genes in heatmap
    ## Below Cris added lines to calculate both the centroid and second-fourth PCs of the gene set shown in the heatmap
    # Note that what we are using is Z-transformed data from the outset
    gdat<-t(apply(eDat(),1,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
    gdat[is.nan(gdat)] = 0 # replace NaNs with zeros
    ss<-rank(svd(gdat)$v[,1])/ncol(gdat)
    ss_centroid<-colMeans(gdat)
    # correct directions
    if(cor(apply(gdat,2,mean),ss)<0) ss<-abs(1-ss)

   #generate object of correlations to centroid and PCs
    CP<-data.frame(cbind(cor(ss_centroid,t(gdat)),cor(ss,t(gdat))))
   
    # Set bounds for heat map values
    gdat[gdat< -3]<- -3
    gdat[gdat>3]<-3
    
    
    if(input$colors=="greenred"){
      cc<-rbind(
        greenred(length(ss))[rank(ss)],
        greenred(length(ss_centroid))[rank(ss_centroid)])
        ccx<-greenred(50)
    } 
    if(input$colors=="bluered"){
      cc<-rbind(
      bluered(length(ss))[rank(ss)],
                 bluered(length(ss_centroid))[rank(ss_centroid)])
      ccx<-bluered(50)
    }
    if(input$colors=="blueyellow"){
      cc<-rbind(
      blue2yellow(length(ss))[rank(ss)],
                blue2yellow(length(ss_centroid))[rank(ss_centroid)])
      ccx<-blue2yellow(50)
    }
    if(input$colors=="matlab"){
      cc<-rbind(
      matlab.like(length(ss))[rank(ss)],
                matlab.like(length(ss_centroid))[rank(ss_centroid)])
      ccx<-matlab.like(50)
    }
    

    oo<-order(ss)
                
# order heatmap by ...

  if(input$order=="No") oo<-1:length(ss)
  if(input$order=="HM1") oo<-order(ss)
  if(input$order=="HMC") oo<-order(ss_centroid)
  if(input$order=="group") oo<-order(ss_centroid)
  if(input$order=="donor") oo<-order(ss_centroid)

  #up.panel<-rbind(colorpanel(ncol(gdat),"green","black","red")[rank(ss_centroid)][oo], colorpanel(ncol(gdat),"green","black","red")[rank(ss)][oo],as.character(ct+4)[oo])
  up.panel<-rbind(colorpanel(ncol(gdat),"green","black","red")[rank(ss_centroid)][oo], colorpanel(ncol(gdat),"green","black","red")[rank(ss)][oo],as.character(ct[1,]+3)[oo],as.character(ct[2,]+6)[oo])
  rownames(up.panel)<-c("Centroid","1st PC","Group","Donor")
  ct
  #print(as.character(ct[1,]*4)[oo])
  #print("asldnfasndf")
  print(up.panel)

  rowsep=1:nrow(gdat); colsep=1:ncol(gdat)

  if(input$order!="hc"){
    suppressWarnings(
    heatmap.mik(gdat[,oo], trace='none',col=ccx, keysize=.5, key=T, ColSideColors=up.panel,mar=c(20,15),Colv=F,cexCol=1.2,scale=input$radioScale, rowsep=rowsep, colsep=colsep, sepwidth=c(0.001, 0.001), sepcol="black",distfun=function (y) dist(y,method = input$radioDist), hclustfun=function (y) hclust(y, method = input$radioCLM))        )
  }
  if(input$order=="hc"){
      suppressWarnings(
      heatmap.mik(gdat, trace='none', col=ccx, keysize=.5, key=T, ColSideColors=up.panel,mar=c(20,15),cexCol=1.2,scale=input$radioScale,distfun=function (y) dist(y,method = input$radioDist),hclustfun=function (y) hclust(y, method = input$radioCLM), rowsep= rowsep, colsep= colsep))
        
  }

},height=1200,width=1000)


  
})

