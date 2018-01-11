library(monocle)
library(reshape2)
library(ggplot2)
library(stringr)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(shiny)


#### Functions needed ####

lookupGeneId<-function(eset,gene_names){
  res <- rownames(fData(eset))[fData(eset)$gene_short_name %in% gene_names]
  res <- c(res,rownames(fData(eset))[rownames(fData(eset)) %in% gene_names])
  res <- unique(res)
  res
}

lookupGeneName<-function(eset,gene_id){
  res <- fData(eset[gene_id,])$gene_short_name
  #res <- unique(res)
  res
}


meltCDS<-function(cds,geneset,logMode=F){
  sub<-cds[lookupGeneId(cds,geneset),]
  sub.expr<-as.matrix(exprs(sub))
  if(logMode){
    sub.expr<-log10(sub.expr+1)
  }
  sub.expr.melt<-melt(sub.expr)
  colnames(sub.expr.melt)<-c("gene_id","cell_id","value")
  res<-merge(sub.expr.melt,pData(sub),by.x="cell_id",by.y="cell_id")
  res<-merge(res,fData(sub),by.x="gene_id",by.y="gene_id")
  res
}

myBarMap<-function(cds,geneset,facet_by="kmeans_tSNE_cluster",color_by="factor(kmeans_tSNE_cluster)",cluster="both",showSummary=T,...){
  sub.melt<-meltCDS(cds,geneset,...)
  facet_by_melt<-strsplit(facet_by,"\\+")[[1]]
  sub.melt.summary<-sub.melt %>%
    dplyr::group_by_(.dots=c("gene_short_name",facet_by_melt)) %>%
    dplyr::summarise(mean=mean(value),median=median(value),sd=sd(value),upper_bound=mean+sd,lower_bound=max(mean-sd,0))
  
  if(cluster %in% c("row","both",T)){
    sub.sum.mat<-sub.melt.summary %>%
      recast(as.formula(paste("gene_short_name ~",facet_by)),measure.var="mean",fun.aggregate=mean)
    sub.sum.hclust<-hclust2(dist(sub.sum.mat[,-1]))
    gene.order.idx<-order.dendrogram(as.dendrogram(sub.sum.hclust))
    gene.order<-sub.sum.mat$gene_short_name[gene.order.idx]
    sub.melt$gene_short_name<-factor(sub.melt$gene_short_name, levels=gene.order)
  }
  
  if(cluster %in% c("column","both",T)){
    sub.mat<-sub.melt %>%
      recast(as.formula("gene_short_name ~ cell_id"),measure.var="value",fun.aggregate=mean)
    sub.hclust<-hclust2(dist(t(sub.mat[,-1])))
    cell.order.idx<-order.dendrogram(as.dendrogram(sub.hclust))
    cell.order<-colnames(sub.mat[,-1])[cell.order.idx]
    #print(cell.order)
    sub.melt$cell_id<-factor(sub.melt$cell_id,levels=cell.order)
  }
  
  p<-ggplot(sub.melt)
  p<-p + geom_bar(aes_string(x="cell_id",y="value",fill=color_by,color=color_by),stat="identity")
  
  if(showSummary){
    p<-p + geom_hline(aes(yintercept=mean),data=sub.melt.summary,size=1.0)
    p<-p + geom_hline(aes(yintercept=upper_bound),data=sub.melt.summary,linetype="dashed")
    p<-p + geom_hline(aes(yintercept=lower_bound),data=sub.melt.summary,linetype="dashed")
  }
  p<-p +
    facet_grid(as.formula(paste("gene_short_name ~", facet_by)),scale="free",space="free_x",labeller=labeller(.default=label_both,gene_short_name=label_value)) +
    theme_bw() + guides(color=FALSE) + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          strip.text.y = element_text(angle=0,hjust=0),
          strip.background = element_blank(),
          panel.background = element_rect(fill="white"),
          panel.margin = unit(0, "lines"),
          panel.grid = element_blank()
    )
  p
}

myBoxplot.subset <- function(cds,markers=NULL,logMode=T,color_by="color", scaled = FALSE){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-as.matrix(exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)]))
    if(logMode){
      genes<-log2(genes+1)
    }
    geneMeans<-rowMax(genes)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    genes<-melt(genes)
    colnames(genes)<-c("cell_id","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="gene_id",all.x=TRUE,sort=FALSE)
    #print(head(genes))
    tmp<-merge(tmp,genes,by.x=0,by.y="cell_id")
    tmp$gene_short_name <- factor(tmp$gene_short_name, levels = markers)
    #print(head(tmp))
    
    p<-ggplot(tmp,aes(factor(subset.cluster),value, fill = subset.cluster))
    p + geom_boxplot() + facet_wrap('gene_short_name', scales = "free_y") + theme_bw() + ylab("log2(Transcripts+1)")
  }
}

myBoxplot.cluster <- function(cds,markers=NULL,logMode=T,color_by="color", scaled = FALSE){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-as.matrix(exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)]))
    if(logMode){
      genes<-log2(genes+1)
    }
    geneMeans<-rowMax(genes)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    genes<-melt(genes)
    colnames(genes)<-c("cell_id","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="gene_id",all.x=TRUE,sort=FALSE)
    #print(head(genes))
    tmp<-merge(tmp,genes,by.x=0,by.y="cell_id")
    tmp$gene_short_name <- factor(tmp$gene_short_name, levels = markers)
    #print(head(tmp))
    
    p<-ggplot(tmp,aes(factor(kmeans_tSNE_cluster),value, fill = kmeans_tSNE_cluster))
    p + geom_boxplot() + facet_wrap('gene_short_name', scales = "free_y") + theme_bw() + ylab("log2(Transcripts+1)")
  }
}

cbindPad <- function(...){
  args <- list(...)
  n <- sapply(args,nrow)
  mx <- max(n)
  pad <- function(x, mx){
    if (nrow(x) < mx){
      nms <- colnames(x)
      padTemp <- matrix("", mx - nrow(x), ncol(x))
      colnames(padTemp) <- nms
      if (ncol(x)==0) {
        return(padTemp)
      } else {
        return(rbind(x,padTemp))
      }
    }
    else{
      return(x)
    }
  }
  rs <- lapply(args,pad,mx)
  return(do.call(cbind,rs))
}

#Adding needed function
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
