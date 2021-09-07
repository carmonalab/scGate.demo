## UTILS
my.summary <- function(actual,pred){
  ct <- table(actual,pred)
  
  ### basic statistics
  measure <- ROCit::measureit(score = pred, class = actual,
                              measure = c("FNR","FPR","SPEC","SENS","PPV"))
  st <- as.data.frame((print(measure)%>%as.matrix()))[2,]
  res.sum<- st[,c("SENS","SPEC","FPR","FNR","PPV")]
  
  ## add mcc coefficient  
  ch2 <- chisq.test(x=actual,y=pred)
  mcc <- sqrt(ch2$statistic/length(actual))
  
  res.sum$MCC <-  mcc
  res.sum$p.value <- ch2$p.value
  
  return(list('conting' = ct,
              'summary' = res.sum ))
}


plot_summary_list <- function(summary.list, title = '', variable.name = 'impurity'){
  impurities <- names(summary.list)
  toshow <-do.call("cbind",summary.list%>%lapply(function(x){x$summary%>%as.matrix()%>%t()}))
  colnames(toshow) <- names(summary.list)
  toshow[,as.character(impurities)]
  
  ts <- (toshow[c('SENS','SPEC','PPV','MCC'),as.character(impurities)])%>%t()%>%as.data.frame()
  ts[,variable.name] <- names(summary.list)
  tsm <-ts%>%melt(id.var = variable.name)
  
  plot <- ggplot(tsm,aes(x=get(variable.name),y=value,group = variable,colour = variable)) +
    geom_point() + geom_line(aes(lty=variable)) + ggtitle(title) +  theme_bw() + labs(x = variable.name)
  return(plot)
}


plot_positives <- function(obj,binary.prediction,feature,label = 1,title = '',min.thr = 10){
  pos <- binary.prediction==label #
  positives <- obj@meta.data[pos,feature]%>%table%>%as.data.frame()
  colnames(positives) <- c(feature,'Freq') 
  plot <- ggplot(positives%>%subset(Freq>min.thr),aes(x=get(feature),y=Freq)) + geom_bar(stat="identity") + theme_bw() + ggtitle(title) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + labs(x = feature)
  return(plot)
}


process_seurat <- function(obj,nfeat = 2000, 
                           selection.method = "vst",pca.dims=5,umap.dims = 30,
                           blacklist =NULL,
                           seed = 123,
                           scaling =T
){
  obj <- NormalizeData(obj, verbose = FALSE);
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeat, verbose = FALSE);
  
  if(!isFALSE(blacklist)){ 
    if(is.null(blacklist)){ 
      blacklist <- scGate::genes.blacklist.Hs
    }
    genes.remove <- unlist(blacklist) #non interesting and/or technical artifact-associated genes
    removing <- intersect(obj@assays$RNA@var.features, genes.remove)
    cat(paste('removing: ',round(100*removing%>%length()/nfeat,2),'% genes\n' ))#, file = con, append =T)
    obj@assays$RNA@var.features <- setdiff(obj@assays$RNA@var.features, genes.remove)
  }
  
  set.seed(seed)
  if(scaling){
    obj <- ScaleData(obj);
  }
  obj <- RunPCA(obj, features = obj@assays$RNA@var.features, ndims.print = 1:pca.dims, nfeatures.print = pca.dims);
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:umap.dims, seed.use=seed); 
  return(obj)
}


pr_ggplot <- function(pr.object,score_title = "Threshold", title = ""){
  plt <- ggplot(data.frame(pr.object$curve),aes(x=X1,y=X2,color=X3)) +
    geom_line() +
    labs(x="Recall", y="Precision",subtitle=paste("AUC-PR: ",format(pr.object$auc.integral,digits=3))
         ,colour=score_title) +
    #scale_colour_gradient2(low="red", mid="orange",high="yellow")
    scale_color_gradientn(colours = rainbow(7)) +
    theme_bw()+
    ggtitle(title)
  return(plt)
}


# thr {NULL: , "extreme", vector}
wrap_pr_curve<-function(zscore_obj,gt,score_method = "sd.in",thr = NULL,drop.na =F){
  require(PRROC)
  if(drop.na){
    isna <- zscore_obj$summary_zscores%>%apply(1,function(x){any(is.na(x))})%>%as.vector()
    zscore_obj$summary_zscores <- aver<- zscore_obj$summary_zscores[!isna,] 
    gt <- gt[!isna]
  }
  
  if(score_method=="sd.in"){
    scoring <- zscore_obj$summary_zscores$max_positive_zscore
    
    if(!length(thr)>1){
      if(is.null(thr)){
        thr <- 7
      }else if(thr == 'extreme'){
        thr <- max(zscore_obj$summary_zscores$max_negative_zscore)          
      }
      ## impute minimum score to cells with contamination zscore above threshold (default sd.out = 7)
      mmin <- min(scoring);
      scoring[zscore_obj$summary_zscores$max_negative_zscore > thr] <- mmin 
      prcurve <- pr.curve(scores.class0 = scoring[gt], scores.class1 = scoring[!gt], curve = TRUE)
      return(prcurve)
    }else if(length(thr)>1){
      mmin <- min(scoring);
      
      aucprs <- lapply(thr,function(qscoring,score = scoring){
        score[zscore_obj$summary_zscores$max_negative_zscore > qscoring] <- mmin 
        auc <- pr.curve(scores.class0 = score[gt], scores.class1 = score[!gt], curve = F,dg.compute = F)$auc.integral
        return(auc)            
      })
      return(data.frame("sd.out"=thr,"aucpr"=unlist(aucprs)))
    }else{
      print("Warning, thr must be one of NULL, extreme or a vector")
      rturn(NULL)
    }        
  }
  
  if(score_method=="sd.out"){
    scoring <- -zscore_obj$summary_zscores$max_negative_zscore  # invert sign (lower zscores in contamination signatures is beter)
    
    if(!length(thr)>1){
      if(is.null(thr)){
        thr <- -4
      }else if(thr == 'extreme'){
        thr <- min(zscore_obj$summary_zscores$max_positive_zscore)          
      }
      ## impute minimum score to cells with a signature scoring below the threshold (default sd.in = -4)
      mmin <- min(scoring);
      scoring[zscore_obj$summary_zscores$max_positive_zscore < thr] <-mmin
      prcurve <- pr.curve(scores.class0 = scoring[gt], scores.class1 = scoring[!gt], curve = TRUE)
      return(prcurve)
      
    }else if(length(thr)>1){
      mmin <- min(scoring);
      
      aucprs <- lapply(thr,function(qscoring,score = scoring){
        score[zscore_obj$summary_zscores$max_positive_zscore < qscoring] <- mmin 
        auc <- pr.curve(scores.class0 = score[gt], scores.class1 = score[!gt], curve = F,dg.compute = F)$auc.integral
        return(auc)            
      })
      return(data.frame("sd.in"=thr,"aucpr"=unlist(aucprs)))
    }else{
      print("Warning, thr must be one of NULL, extreme or a vector")
      rturn(NULL)
    }        
  }
  
}

#markerlist: list of cell signatures for model training
#@positive_celltypes must be found in markerlist names


#@max.impurity : result must not be sensible to 
get_zscores <- function(object,gating.model, max.impurity=0.5, vervose =T, return_seurat =F,gated_object =F,ncores = 1){  

  if(!gated_object){  
  gated.object <- scGate(object, gating.model = gating.model ,
                         max.impurity = max.impurity, 
                         return_signature_scores = T,verbose = vervose, ncores = ncores)
  }else{
    gated.object <- object
  }
  zscore_obj <- list();  # output object
  
  if(class(gating.model)!='list'){   # one layer models
    zscore_cols <- paste0(names(gating.model@markers),"_Zscore")
    zpos_cols <- paste0(gating.model@positive_celltypes,"_Zscore")
    scoring_methods <- gating.model@scoring_method
  }else{                            
    # Multiple layered models
    zcols_list <- lapply(gating.model,function(model){
      zcols <- paste0(names(model@markers),"_Zscore")
      return(zcols)
    })

    zpos_list <- lapply(gating.model,function(model){
      zpos_cols <- paste0(model@positive_celltypes,"_Zscore")
      return(zpos_cols)
    })

    zscore_cols <- unlist(zcols_list)
    zpos_cols <- unlist(zpos_list)

    
    scoring_methods <- lapply(gating.model,function(model){
      scoring_meth <- model@scoring_method
      return(scoring_meth)
    })%>%unlist()%>%unique()

    }
  
  
  zscores_Ucell <- gated.object@meta.data[,zscore_cols]   

  zscore_obj$z <- zscores_Ucell;
  zscore_obj$positive_celltypes_z <- zpos_cols;
  zscore_obj$scoring <- scoring_methods
  
  
  max_neg <- zscore_obj$z%>%select(-all_of(zpos_cols))%>%apply(1,max)
  max_pos <- zscore_obj$z%>%select(zpos_cols)%>%apply(1,max)
  summary_zscores <- data.frame("max_positive_zscore" = max_pos, "max_negative_zscore" = max_neg)
  
  zscore_obj$summary_zscores <- summary_zscores
  
  if(return_seurat == T){  
    gated.object <- AddMetaData(gated.object,metadata = max_pos, col.name = '"max_positive_zscore"')
    gated.object <- AddMetaData(gated.object,metadata = max_neg, col.name = '"max_negative_zscore"')
    return(gated.object)
  }else{
    return(zscore_obj)
  }
}

basic_qc_seurat <- function(data,nFeature_RNA.colname = "nFeature_RNA",Total.counts.colname = "nCount_RNA",
                            min.features = 500, max.features = 6000, min.counts = 600,max.counts = 20000,
                            max.ribo = 50, max.mito = 15){
  percent.ribo.dv <- PercentageFeatureSet(data, pattern = "^RP[LS]")
  percent.mito.dv <- PercentageFeatureSet(data, pattern = "^MT-")
  
  data <- AddMetaData(data, metadata = percent.ribo.dv, col.name = "percent.ribo")
  data <- AddMetaData(data, metadata = percent.mito.dv, col.name = "percent.mito")
  
  mtd <- data@meta.data
  ind <- (mtd[,nFeature_RNA.colname]>min.features & mtd[,nFeature_RNA.colname]<max.features &
            mtd[,Total.counts.colname]>min.counts & mtd[,Total.counts.colname]< max.counts &
            mtd$percent.ribo < max.ribo  &
            mtd$percent.mito < max.mito)
  
  cellnames = (mtd%>%rownames)[ind]
  data <- data%>%subset(cells = cellnames)# & 
  
  return(data)  
  
}
