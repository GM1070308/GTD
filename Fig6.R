library(BayesPrism)
library(gplots)
### Quality control of input scRNA data

norm.to.one <- function(ref, 
                        pseudo.min){
  
  G <- ncol(ref)
  
  phi <- ref/rowSums(ref) * (1-pseudo.min*G) + pseudo.min
  
  #if the minimum value is greater than zero. simply normalize by total depth
  min.value <- apply(ref,1,min)
  which.row <- min.value>0
  if(any(which.row)){
    #cat("One or more cell types have all genes with non-zero expression. pseudo.min is not applied to these cell types. \n")
    phi[which.row,] <- ref[which.row,,drop=F]/rowSums(ref[which.row,,drop=F])
  }
  
  
  return(phi)
  
}

collapse <- function(ref, labels){
  
  stopifnot(nrow(ref) == length(labels))
  
  #remove NA in labels
  non.na.idx <- !is.na(labels)
  if(sum(!non.na.idx)>0) print("Warning: NA found in the cell type/state labels. These cells will be excluded!")
  labels <- labels[non.na.idx]
  ref <- ref[non.na.idx,]
  
  labels.uniq <- unique(labels)
  
  ref.collapsed <- do.call(rbind,
                           lapply(labels.uniq,
                                  function(label.i) 
                                    colSums(ref[labels==label.i,,drop=F])
                           )
  )
  
  rownames(ref.collapsed) <- labels.uniq
  
  return(ref.collapsed)
}
plot.cor.phi <- function(input, 
                         input.labels,
                         min.exp=3,
                         pseudo.min=1E-8, 
                         my_palette= colorRampPalette(c("#3C5488FF", "white", "#DC0000FF")),
                         title ='',
                         symkey=F,
                         symbreaks=F,
                         pdf.prefix=NULL,
                         ...){
  
  input <- input[,colSums(input) >= min.exp]
  
  ref.ct <- collapse(ref = input, labels = input.labels)
  ref.ct <- norm.to.one (ref = ref.ct, pseudo.min = pseudo.min)
  
  ref.ct <- scale(log2(ref.ct),center=T,scale=F)
  ref.ct <- t(ref.ct)
  cor.mat <- cor(ref.ct)
  
  hc.cols <- hclust(as.dist(1-cor.mat),method="ward.D2")
  hc.rows <- hc.cols
  
  if(!is.null(pdf.prefix)) pdf(paste(pdf.prefix,"_cor_phi.pdf",sep=""),pointsize=8,useDingbats=FALSE)
  heatmap.2(cor.mat, 
            col=my_palette, 
            density.info="none", 
            trace="none",
            dendrogram='both', 
            symm=TRUE, 
            symkey= symkey, 
            symbreaks= symbreaks,
            scale="none" , 
            Rowv=as.dendrogram(hc.rows), 
            Colv=as.dendrogram(hc.cols),
            keysize = 1,
            key.title = NULL,
            key.xlab = "Pearson's correlation",
            ...)
  
  if(!is.null(pdf.prefix)) dev.off()	
}
setwd('~/thesis/BayesPrism/')
plot.cor.phi (input=sc.dat,
              input.labels=cell.state.labels,
              title="cell type correlation",
              #specify pdf.prefix if need to output to pdf
              pdf.prefix="gbm.cor.state.cs", 
              cexRow=0.9, cexCol=0.9,
              margins=c(10,10))

compute.specificity <- function(input.matrix,
                                pseudo.min = 1E-8){
  
  ref.ct <- norm.to.one (ref = input.matrix, pseudo.min = pseudo.min)
  
  #compute gene expression specificity score
  exp.spec <- t(ref.ct) /  colSums(ref.ct)
  max.spec <- apply(exp.spec,1,max)
  
  return(max.spec)
}

plot.scRNA.outlier <- function(input,
                               cell.type.labels,
                               species,
                               pdf.prefix=NULL,
                               return.raw=FALSE){
  
  input <- collapse(ref = input, labels = cell.type.labels)
  ref.ct <- norm.to.one (ref = input, pseudo.min = 1E-8)
  
  #compute mean expression
  exp.mean <- colMeans(ref.ct)
  exp.mean.log <- log(exp.mean)
  
  #compute gene expression specificity score
  max.spec <- compute.specificity(input)
  
  category.matrix <- assign.category(input.genes = colnames(ref.ct), 
                                     species= species)
  
  #make plots	
  if(!is.null(pdf.prefix)) pdf(paste(pdf.prefix,"_scRNA_outlier.pdf",sep=""),pointsize=8,useDingbats=FALSE , width = 10, height = 6)
  
  cex=1
  plot(y= max.spec, x= exp.mean.log, cex=0.5, pch=16, col="#dbd9d9", 
       xlab="Log mean expression", ylab="Maximum expression specificity",
       xlim = c(-19,-3), ylim = c(0,1.05), xaxt='n',
       xaxs = 'i', yaxs = 'i', cex.lab = 1.5)
  axis(1, at=c(-19, -17,-15,-13,-11,-9,-7,-5,-3))
  points(y= max.spec[category.matrix[,"Rb"]], x= exp.mean.log[category.matrix[,"Rb"]],cex= cex,col="#E64B35FF")
  points(y= max.spec[category.matrix[,"Mrp"]], x= exp.mean.log[category.matrix[,"Mrp"]],cex= cex,col="#3C5488FF")
  points(y= max.spec[category.matrix[,"other_Rb"]], x= exp.mean.log[category.matrix[,"other_Rb"]],cex= 1,col="#F39B7FFF")
  
  points(y= max.spec[category.matrix[,"chrM"]], x= exp.mean.log[category.matrix[,"chrM"]],cex= cex,col="#4DBBD5FF")
  
  points(y= max.spec[category.matrix[,"act"]], x= exp.mean.log[category.matrix[,"act"]],cex= cex,col="#91D1C2FF")
  points(y= max.spec[category.matrix[,"hb"]], x= exp.mean.log[category.matrix[,"hb"]],cex= cex,col="#8491B4FF")
  
  if(species=="hs")
    points(y= max.spec[category.matrix[,"MALAT1"]], x= exp.mean.log[category.matrix[,"MALAT1"]],cex=1,col="#00A087FF")
  
  if(!is.null(pdf.prefix)) dev.off()
  
  if(return.raw) return(cbind.data.frame(exp.mean.log, max.spec, category.matrix))
}

assign.category <- function(input.genes,
                            species=c("mm","hs")){
  
  stopifnot(length(species)==1 & species %in% c("hs","mm"))
  
  #load gene list
  if(species=="hs") gene.list <- read.table(system.file("extdata", "genelist.hs.new.txt", package="BayesPrism"),sep="\t",header=F,stringsAsFactors=F)
  if(species=="mm") gene.list <- read.table(system.file("extdata", "genelist.mm.new.txt", package="BayesPrism"),sep="\t",header=F,stringsAsFactors=F)
  
  
  #detect if EMSEMBLE ID (starts with ENS) or gene symbol is used
  if( sum(substr(input.genes,1,3)=="ENS")> length(input.genes)*0.8  ){
    # use 80% of colnames to avoid sometimes there are manually added gene name such as GFP-XXX.
    #use EMSEMBLE ID
    #strip the "." from ENSXXX.X
    cat("EMSEMBLE IDs detected.\n")
    input.genes.short <- unlist(lapply(input.genes, function(gene.id) strsplit(gene.id,split="\\.")[[1]][1]))
    gene.df <- gene.list[,c(1,2)]
  }
  else{
    #use gene symbols
    cat("Gene symbols detected. Recommend to use EMSEMBLE IDs for more unique mapping.\n")
    input.genes.short <- input.genes
    gene.df <- gene.list[,c(1,3)]
  }
  
  gene.group.matrix <- do.call(cbind.data.frame, lapply(unique(gene.df[,1]), 
                                                        function(gene.group.i) input.genes.short %in% gene.df[gene.df[,1]== gene.group.i,2]))
  colnames(gene.group.matrix) <- unique(gene.df[,1])
  rownames(gene.group.matrix) <- input.genes
  return(gene.group.matrix)
}

plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=FALSE, #return the data used for plotting. 
  pdf.prefix="gbm.sc.stat.type"
)

plot.bulk.outlier <- function(bulk.input,
                              sc.input,
                              cell.type.labels,
                              species,
                              pdf.prefix=NULL,
                              return.raw=FALSE){
  
  bulk.norm <- norm.to.one (ref = bulk.input, pseudo.min = 1E-8)
  
  #compute mean expression
  exp.mean <- colMeans(bulk.norm)
  exp.mean.log <- log(exp.mean)
  
  #compute gene expression specificity score
  sc.input <- collapse(ref = sc.input, labels = cell.type.labels)
  max.spec <- compute.specificity(sc.input)
  
  max.spec <- max.spec[match(colnames(bulk.norm),names(max.spec))]
  
  category.matrix <- assign.category(input.genes = colnames(bulk.norm), 
                                     species= species)
  
  
  #make plots	
  if(!is.null(pdf.prefix)) pdf(paste(pdf.prefix,"_bulk_outlier.pdf",sep=""),pointsize=8,useDingbats=FALSE, width = 10, height = 6 )
  
  cex=1
  plot(y= max.spec, x= exp.mean.log, cex=0.5, pch=16, col="#dbd9d9", 
       xlab="Log mean expression", ylab="Maximum expression specificity",
       xlim = c(-19,-3), ylim = c(0,1.05), xaxt='n', xaxs = 'i', yaxs = 'i', cex.lab = 1.5)
  axis(1, at=c(-19, -17,-15,-13,-11,-9,-7,-5,-3))
  points(y= max.spec[category.matrix[,"Rb"]], x= exp.mean.log[category.matrix[,"Rb"]],cex= cex,col="#E64B35FF")
  points(y= max.spec[category.matrix[,"Mrp"]], x= exp.mean.log[category.matrix[,"Mrp"]],cex= cex,col="#3C5488FF")
  points(y= max.spec[category.matrix[,"other_Rb"]], x= exp.mean.log[category.matrix[,"other_Rb"]],cex= 0.3,col="#F39B7FFF")
  
  points(y= max.spec[category.matrix[,"chrM"]], x= exp.mean.log[category.matrix[,"chrM"]],cex= cex,col="#4DBBD5FF")
  
  points(y= max.spec[category.matrix[,"act"]], x= exp.mean.log[category.matrix[,"act"]],cex= cex,col="#91D1C2FF")
  points(y= max.spec[category.matrix[,"hb"]], x= exp.mean.log[category.matrix[,"hb"]],cex= cex,col="#8491B4FF")
  
  if(species=="hs")
    points(y= max.spec[category.matrix[,"MALAT1"]], x= exp.mean.log[category.matrix[,"MALAT1"]],cex=1,col="#00A087FF")
  
  if(!is.null(pdf.prefix)) dev.off()
  
  if(return.raw) return(cbind.data.frame(exp.mean.log, max.spec, category.matrix))
}
bk.dat <- t(read.csv('~/EPIC/bk2.txt', sep = ' '))

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=FALSE,
  pdf.prefix="gbm.bk.stat"
)

plot.bulk.vs.sc <- function(sc.input,
                            bulk.input,
                            pdf.prefix=NULL,
                            return.value=FALSE){
  
  gene.shared <- intersect(colnames(sc.input), colnames(bulk.input))
  if(length(gene.shared)==0)
    stop("Error: gene names of reference and mixture do not match!") 
  if(length(gene.shared) < 100)
    warning("Warning: very few gene from reference and mixture match! Please double check your gene names.")
  
  #align reference and mixture
  bulk.input <- bulk.input[, gene.shared]
  sc.input <- sc.input[, gene.shared]
  
  bulk.norm <- norm.to.one (ref = matrix(colSums(bulk.input),nrow=1),
                            pseudo.min = 1E-8)[1,]
  sc.norm <- norm.to.one (ref = matrix(colSums(sc.input),nrow=1),
                          pseudo.min = 1E-8)[1,]
  
  gene.tab.path <- system.file("extdata","gencode.v22.broad.category.txt", package="BayesPrism")	
  gene.list <- read.table(gene.tab.path, sep="\t",header=F,stringsAsFactors=F)
  
  #decide if emsembl ID or gene symbol is used
  if( sum(substr(gene.shared,1,3)=="ENS")> length(gene.shared)*0.8  ){
    # use 80% of colnames to avoid sometimes there are manually added gene name such as GFP-XXX.
    #use EMSEMBLE ID
    #strip the "." from ENSXXX.X
    cat("EMSEMBLE IDs detected.\n")
    gene.shared <- unlist(lapply(gene.shared, function(gene.id) strsplit(gene.id,split="\\.")[[1]][1]))
    gene.df <- gene.list[match(gene.shared, gene.list[,8]),c(8,9)]
  }
  else{
    #use gene symbols
    cat("Gene symbols detected. Recommend to use EMSEMBLE IDs for more unique mapping.\n")
    gene.df <- gene.list[match(gene.shared, gene.list[,5]),c(5,9)]
  }
  colnames(gene.df) <- c("gene.name", "category")
  
  #plotting three groups: lincRNA, protein_coding and pseudogene
  
  plot.df <- data.frame(log2.bulk = log2(bulk.norm),
                        log2.sc = log2(sc.norm),
                        gene.df)
  
  selected.gene.types <- c("lincRNA", "protein_coding", "pseudogene")
  plot.df <- plot.df[plot.df$category %in% selected.gene.types,]
  
  #make plots	
  if(!is.null(pdf.prefix)) pdf(paste(pdf.prefix,"_cor_by_geneType.pdf",sep=""),pointsize=8,useDingbats=FALSE )
  
  pc.idx <- plot.df$category=="protein_coding"
  lnc.idx <- plot.df$category=="lincRNA"
  psg.idx <- plot.df$category=="pseudogene"
  
  par(mfrow = c(1,3), pty = "s")
  
  xlim <- range(plot.df[,"log2.sc"])
  ylim <- range(plot.df[,"log2.bulk"])
  plot.each(y.value = plot.df[pc.idx,"log2.bulk"], 
            x.value = plot.df[pc.idx,"log2.sc"], 
            xlab ="log2 total expression in scRNA-seq", 
            ylab ="log2 total expression in bulk",
            title = "protein coding", 
            x.lim = xlim, y.lim = ylim)
  
  plot.each(y.value = plot.df[lnc.idx,"log2.bulk"], 
            x.value = plot.df[lnc.idx,"log2.sc"], 
            xlab ="log2 total expression in scRNA-seq", 
            ylab ="log2 total expression in bulk",
            title = "lncRNA", 
            x.lim = xlim, y.lim = ylim)
  
  plot.each(y.value = plot.df[psg.idx,"log2.bulk"], 
            x.value = plot.df[psg.idx,"log2.sc"], 
            xlab ="log2 total expression in scRNA-seq", 
            ylab ="log2 total expression in bulk",
            title = "pseudogene", 
            x.lim = xlim, y.lim = ylim)
  
  if(!is.null(pdf.prefix)) dev.off()
  
  if(return.value) return(plot.df)
}

plot.each <- function(y.value, 
                      x.value,
                      cex=0.2,
                      pch=16,
                      col=adjustcolor("black",alpha.f=0.5),
                      title=NULL,
                      x.lim=NULL,
                      y.lim=NULL,
                      ...){
  
  if(is.null(x.lim)) x.lim <- range(x.value, na.rm=T)
  if(is.null(y.lim)) y.lim <- range(y.value, na.rm=T)
  
  plot(x= x.value, y= y.value, cex= cex,col= col,pch=16, 
       ylim=y.lim,xlim=x.lim, main= title, ...)
  abline(a=0,b=1,lty=2,col="red")
  
  cor.sp <- cor(x.value, y.value, method="spearman", use="complete")
  cor.ps <- cor(x.value, y.value, use="complete" )
  mse <-  mean((x.value- y.value)^2, na.rm=T)
  
  text(x=0.75*(x.lim[2]-x.lim[1])+x.lim[1], y=0.40*(y.lim[2]-y.lim[1])+y.lim[1], labels=paste("R=", round(cor.ps,3) ,sep=""))
  text(x=0.75*(x.lim[2]-x.lim[1])+x.lim[1], y=0.30*(y.lim[2]-y.lim[1])+y.lim[1], labels=paste("rho=", round(cor.sp,3),sep=""))
  text(x=0.75*(x.lim[2]-x.lim[1])+x.lim[1], y=0.20*(y.lim[2]-y.lim[1])+y.lim[1], labels=paste("MSE=", signif(mse,3),sep=""))
  
  NULL		   	
}

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat,
                 pdf.prefix="gbm.bk.vs.sc"
)
