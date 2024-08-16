#' This function performs enrichment scoring.
#'
#' @param rbp_interested The name of the RBP interested.
#' @param Events_DS The difference matrix of the clipping events, obtained by the MRAS function DS_matrix().
#' @param Event_DS_sig Events with significant differences clipped by threshold filtering in Events_DS.
#' @param DS_pvalue Significance level of differential splicing analysis.
#' @param DS_dPSI Threshold for differential splicing.
#' @param RBP_use RBP expression matrix, the last column is the log2FC value expressed under two conditions
#' @param result_type The types of return values are "top10", "tab_simple", and "tab_all".
#'                    "top10" returns the top 10 MRAS results;
#'                    "tab_simple" returns a simplified version of the MRAS scoring matrix (without the events corresponding to the RBP);
#'                    "tab_all" returns the full version of the MRAS scoring matrix. matrix (including the events corresponding to the RBP).
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 2.
#' @param rbp_event_deal_all_total Regulatory matrix of RBP-Event.
#' @param rbp_event_deal_all Significant regulatory matrix of RBP-Event.
#' @param path_use The path to the file used to store the output.
#'
#' @return Results of MRAS analysis.
#' @export
#' @import foreach
#' @import doParallel
#' @import utils
#' @import data.table
#' @import tidyr
#' @import fgsea
#' @import BiocParallel
#' @import parallel
#' @import fastmatch
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise left_join n
#' @importFrom stats p.adjust na.omit fisher.test binomial ave
#'
#'

score_matrix<-function(rbp_interested,Events_DS,Event_DS_sig,RBP_use,DS_pvalue,DS_dPSI,
                       result_type = c("Top10","tab_simple","tab_all"),threads = 2,
                       rbp_event_deal_all_total,rbp_event_deal_all,path_use){
  #RBP_use<-rbp_expr
  Events=OR=RBP=dpsi_score=i=logFC=modeFraction=nes1_es=nes1_nes=nes1_p=nes1_size=
    nes2_es=nes2_nes=nes2_p=nes2_size=pval=score3=total_size=NULL
  result_type<-match.arg(result_type)
  rbp_expr_use<-cbind(name=rownames(RBP_use),RBP_use)

  uni_use<-merge(as.data.frame(rbp_event_deal_all),rbp_expr_use[,c(1,ncol(rbp_expr_use))],by.x="RBP",by.y="name")
  uni_use$score<-as.numeric(uni_use$score)
  uni_use$logFC<-as.numeric(uni_use$logFC)
  uni_use$aver<-uni_use$score*uni_use$logFC
  uni_use_tab<-as.matrix(table(uni_use$RBP))
  uni_use_tab<-cbind(rownames(uni_use_tab),uni_use_tab)
  colnames(uni_use_tab)<-c("RBP","count")
  #unique
  uni_DS<-merge(Event_DS_sig,uni_use,by.x="Events",by.y="ID")
  uni_DS$aver<-as.numeric(uni_DS$aver)
  uni_DS$abs_dPSI<-as.numeric(uni_DS$abs_dPSI)
  uni_DS$dpsi_score<-(uni_DS$abs_dPSI)*(uni_DS$aver)
  uni_DS %>% dplyr::group_by(RBP) %>%
    dplyr::summarise(logFC=unique(logFC),Events=tj(Events),score1=sum(dpsi_score),m1=n()) -> rbp_uni_mat
  rbp_uni_mat$score1<-abs(as.numeric(rbp_uni_mat$score1))
  rbp_uni_mat$score1_nor<-(rbp_uni_mat$score1-min(rbp_uni_mat$score1))/(max(rbp_uni_mat$score1)-min(rbp_uni_mat$score1))

  #weight matrix 1
  w1<-as.data.frame(table(rbp_event_deal_all[,1]))
  colnames(w1)<-c("ID","weight1")
  w1$weight1<-1+sqrt((max(w1$weight1)-w1$weight1)/(max(w1$weight1)-min(w1$weight1)))
  w1$weight1<-round(w1$weight1,3)


  #weight matrix 2
  if (file.exists(paste0(path_use,"w2.txt"))) file.remove(paste0(path_use,"w2.txt"))

  for (rr in unique(rbp_event_deal_all_total$RBP)) {
    w2_part<-rbp_event_deal_all_total[which(rbp_event_deal_all_total$RBP==rr),]
    w2_part$weight<-1+sqrt((w2_part$dPSI-min(w2_part$dPSI))/(max(w2_part$dPSI)-min(w2_part$dPSI)))
    data.table::fwrite(as.data.frame(w2_part),file = paste0(path_use,"w2.txt"),col.names = F,row.names = F,
           quote = F,sep = "\t",append = T)
    # cat(rr,"\t")
  }
  w2<-data.table::fread(paste0(path_use,"w2.txt"))
  colnames(w2)<-c("ID","score","RBP","label","weight2")

  # cat("weight all","\t")
  w3<-left_join(w2,w1,by = "ID")
  w3[is.na(w3)]<-1
  w3$weight_all<-w3$weight1*w3$weight2*w3$score
  w3$if_DS<-ifelse(w3$ID %in% Event_DS_sig[,1],1,0)

  RBP_candidate<-unique(w3$RBP)
  # cat("xunhuan","\t")
  # t1=proc.time()
  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  gsea<-function(pathways , stats ,maxSize = length(stats)-1,absEps = NULL,sampleSize  = 101,minSize = 1,
                 eps = 1e-50,scoreType   = c("std", "pos", "neg"),BPPARAM = NULL,nproc  = 0,nPermSimple = 1000){
    preparePathwaysAndStats <- function(pathways, stats, minSize, maxSize, gseaParam, scoreType){
      # Error if pathways is not a list
      if (!is.list(pathways)) {
        stop("pathways should be a list with each element containing names of the stats argument")
      }

      # Error if stats is not named
      if (is.null(names(stats))) {
        stop("stats should be named")
      }

      # Error if stats are non-finite
      if (any(!is.finite(stats))){
        stop("Not all stats values are finite numbers")
      }

      # Warning message for ties in stats
      ties <- sum(duplicated(stats[stats != 0]))
      if (ties != 0) {
        warning("There are ties in the preranked stats (",
                paste(round(ties * 100 / length(stats), digits = 2)),
                "% of the list).\n",
                "The order of those tied genes will be arbitrary, which may produce unexpected results.")
      }

      # Warning message for duplicate gene names
      if (any(duplicated(names(stats)))) {
        warning("There are duplicate gene names, fgsea may produce unexpected results.")
      }

      if (all(stats > 0) & scoreType == "std"){
        warning("All values in the stats vector are greater than zero and scoreType is \"std\", ",
                "maybe you should switch to scoreType = \"pos\".")
      }

      stats <- sort(stats, decreasing=TRUE)
      stats <- abs(stats) ^ gseaParam

      minSize <- max(minSize, 1)
      maxSize <- min(maxSize, length(stats)-1)

      pathwaysFiltered <- lapply(pathways, function(p) { unique(na.omit(fastmatch::fmatch(p, names(stats)))) })
      pathwaysSizes <- sapply(pathwaysFiltered, length)

      toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)

      pathwaysFiltered <- pathwaysFiltered[toKeep]
      pathwaysSizes <- pathwaysSizes[toKeep]

      list(filtered=pathwaysFiltered,
           sizes=pathwaysSizes,
           stats=stats)
    }
    pp <- preparePathwaysAndStats(pathways , stats , minSize = 1,
                                  maxSize = maxSize, gseaParam = 1, scoreType)
    pathwaysFiltered <- pp$filtered
    pathwaysSizes <- pp$sizes
    stats <- pp$stats
    m <- length(pathwaysFiltered)
    if (m == 0) {
      return(data.table(pathway=character(),
                        pval=numeric(),
                        padj=numeric(),
                        log2err=numeric(),
                        ES=numeric(),
                        NES=numeric(),
                        size=integer(),
                        leadingEdge=list()))
    }

    # Warning message for deprecated absEps parameter
    if (!is.null(absEps)){
      warning("You are using deprecated argument `absEps`. ",
              "Use `eps` argument instead. ",
              "`absEps` was assigned to `eps`.")
      eps <-  absEps
    }

    # Warning message for to small value for sampleSize
    if (sampleSize < 3){
      warning("sampleSize is too small, so sampleSize = 3 is set.")
      sampleSize <- max(3, sampleSize)
    }

    #To avoid warnings during the check
    log2err=nMoreExtreme=pathway=pval=padj=NULL
    nLeZero=nGeZero=leZeroMean=geZeroMean=nLeEs=nGeEs=isCpGeHalf=NULL
    ES=NES=size=leadingEdge=NULL
    .="damn notes"

    minSize <- max(minSize, 1)
    eps <- max(0, min(1, eps))

    if (sampleSize %% 2 == 0){
      sampleSize <-  sampleSize + 1
    }

    gseaStatRes <- do.call(rbind,
                           lapply(pathwaysFiltered,
                                  fgsea::calcGseaStat,
                                  stats             = stats,
                                  returnLeadingEdge = TRUE,
                                  scoreType         = scoreType))
    leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
    pathwayScores <- unlist(gseaStatRes[, "res"])

    setUpBPPARAM <- function(nproc=0, BPPARAM=NULL){
      if (is.null(BPPARAM)) {
        if (nproc != 0) {
          if (nproc == 1) {
            result <- BiocParallel::SerialParam(progressbar = TRUE)
          }
          if (.Platform$OS.type == "windows") {
            # windows doesn't support multicore, using snow instead
            result <- BiocParallel::SnowParam(workers = nproc, progressbar = TRUE)
          } else {
            result <- BiocParallel::MulticoreParam(workers = nproc, progress = TRUE)
          }
        } else {
          result <- BiocParallel::bpparam()
        }
        return(result)
      }
      else {
        return(BPPARAM)
      }
    }

    seeds <- sample.int(10^9, 1)
    BPPARAM <- setUpBPPARAM(nproc=nproc, BPPARAM=BPPARAM)
    fgseaSimpleImpl <- function(pathwayScores, pathwaysSizes,
                                pathwaysFiltered, leadingEdges,
                                permPerProc, seeds,toKeepLength,
                                stats, BPPARAM, scoreType){
      K <- max(pathwaysSizes)
      universe <- seq_along(stats)

      counts <- BiocParallel::bplapply(seq_along(permPerProc), function(i) {
        nperm1 <- permPerProc[i]
        leEs <- rep(0, toKeepLength)
        geEs <- rep(0, toKeepLength)
        leZero <- rep(0, toKeepLength)
        geZero <- rep(0, toKeepLength)
        leZeroSum <- rep(0, toKeepLength)
        geZeroSum <- rep(0, toKeepLength)
        if (toKeepLength == 1) {
          set.seed(seeds[i])
          for (j in seq_len(nperm1)) {
            randSample <- sample.int(length(universe), K)
            randEsP <- fgsea::calcGseaStat(
              stats = stats,
              selectedStats = randSample,
              gseaParam = 1,
              scoreType = scoreType)
            leEs <- leEs + (randEsP <= pathwayScores)
            geEs <- geEs + (randEsP >= pathwayScores)
            leZero <- leZero + (randEsP <= 0)
            geZero <- geZero + (randEsP >= 0)
            leZeroSum <- leZeroSum + pmin(randEsP, 0)
            geZeroSum <- geZeroSum + pmax(randEsP, 0)
          }
        } else {
          aux <- calcGseaStatCumulativeBatch(
            stats = stats,
            gseaParam = 1,
            pathwayScores = pathwayScores,
            pathwaysSizes = pathwaysSizes,
            iterations = nperm1,
            seed = seeds[i],
            scoreType = scoreType)
          leEs = get("leEs", aux)
          geEs = get("geEs", aux)
          leZero = get("leZero", aux)
          geZero = get("geZero", aux)
          leZeroSum = get("leZeroSum", aux)
          geZeroSum = get("geZeroSum", aux)
        }
        data.table::data.table(pathway=seq_len(toKeepLength),
                               leEs=leEs, geEs=geEs,
                               leZero=leZero, geZero=geZero,
                               leZeroSum=leZeroSum, geZeroSum=geZeroSum
        )
      }, BPPARAM=BPPARAM)

      counts <- rbindlist(counts)

      # Getting rid of check NOTEs
      leEs=leZero=geEs=geZero=leZeroSum=geZeroSum=NULL
      pathway=padj=pval=ES=NES=geZeroMean=leZeroMean=NULL
      nMoreExtreme=nGeEs=nLeEs=nLeZero=nGeZero=size=NULL
      leadingEdge=NULL
      .="damn notes"

      pvals <- counts[, list(leZeroMean = sum(leZeroSum) / sum(leZero),
                             geZeroMean = sum(geZeroSum) / sum(geZero),
                             nLeZero = sum(leZero),
                             nGeZero = sum(geZero),
                             nLeEs = sum(leEs),
                             nGeEs = sum(geEs)),
                      by = .(pathway)]

      pvals[, ES := pathwayScores[pathway]]

      pvals[, NES := as.numeric(NA)]

      switch(scoreType,
             std = pvals[(ES > 0 & geZeroMean != 0) | (ES <= 0 & leZeroMean != 0),
                         NES := ES / ifelse(ES > 0, geZeroMean, abs(leZeroMean))],
             pos = pvals[(ES >= 0 & geZeroMean != 0), NES := ES / geZeroMean],
             neg = pvals[(ES <= 0 & leZeroMean != 0), NES := ES / abs(leZeroMean)])

      pvals[, pval := as.numeric(NA)]
      pvals[!is.na(NES), pval := pmin((1+nLeEs) / (1 + nLeZero),
                                      (1+nGeEs) / (1 + nGeZero))]


      pvals[, padj := as.numeric(NA)]
      pvals[!is.na(pval), padj := p.adjust(pval, method = "BH")]

      switch(scoreType,
             std = pvals[, nMoreExtreme :=  ifelse(ES > 0, nGeEs, nLeEs)],
             pos = pvals[, nMoreExtreme :=  nGeEs],
             neg = pvals[, nMoreExtreme :=  nLeEs])

      pvals[, size := pathwaysSizes[pathway]]
      pvals[, pathway := names(pathwaysFiltered)[pathway]]
      pvals[, leadingEdge := .(leadingEdges)]

      pvals
    }
    calcGseaStatCumulativeBatch <- function(stats, gseaParam, pathwayScores, pathwaysSizes, iterations, seed, scoreType) {
      .Call('_fgsea_calcGseaStatCumulativeBatch', PACKAGE = 'fgsea', stats, gseaParam, pathwayScores, pathwaysSizes, iterations, seed, scoreType)
    }

    simpleFgseaRes <- fgseaSimpleImpl(pathwayScores=pathwayScores, pathwaysSizes=pathwaysSizes,
                                      pathwaysFiltered=pathwaysFiltered, leadingEdges=leadingEdges,
                                      permPerProc=nPermSimple, seeds=seeds, toKeepLength=m,
                                      stats=stats, BPPARAM=BiocParallel::SerialParam(), scoreType=scoreType)

    switch(scoreType,
           std = simpleFgseaRes[, modeFraction := ifelse(ES >= 0, nGeZero, nLeZero)],
           pos = simpleFgseaRes[, modeFraction := nGeZero],
           neg = simpleFgseaRes[, modeFraction := nLeZero])

    simpleFgseaRes[, leZeroMean := NULL]
    simpleFgseaRes[, geZeroMean := NULL]
    simpleFgseaRes[, nLeEs := NULL]
    simpleFgseaRes[, nGeEs := NULL]
    simpleFgseaRes[, nLeZero := NULL]
    simpleFgseaRes[, nGeZero := NULL]

    simpleFgseaRes[modeFraction < 10, pval := as.numeric(NA)]
    simpleFgseaRes[modeFraction < 10, padj := as.numeric(NA)]
    simpleFgseaRes[modeFraction < 10, NES := as.numeric(NA)]

    if (any(simpleFgseaRes$modeFraction < 10)){
      warning("There were ",
              paste(sum(simpleFgseaRes$modeFraction < 10)),
              " pathways for which P-values were not calculated properly due to ",
              "unbalanced (positive and negative) gene-level statistic values. ",
              "For such pathways pval, padj, NES, log2err are set to NA. ",
              "You can try to increase the value of the argument nPermSimple (for example set it nPermSimple = ",
              paste0(format(nPermSimple * 10, scientific = FALSE), ")"))
    }

    # Storing NA fgseaSimple results in a separate data.table
    naSimpleRes <- simpleFgseaRes[is.na(pval)]
    naSimpleRes[, padj := as.numeric(NA)]
    naSimpleRes[, log2err := as.numeric(NA)]
    naSimpleRes[, modeFraction := NULL]

    simpleFgseaRes <- simpleFgseaRes[!is.na(pval)]

    res<-do.call(cbind,simpleFgseaRes)

    return(res)
  }

  # t1=proc.time()
  ES_mat_RBP<-foreach::foreach(rr = RBP_candidate, .packages = c('data.table','dplyr','stats'),.errorhandling = 'remove',.combine = 'rbind') %dopar% {
    # cat(which(target==rr),"\n")
    # library(data.table)
    w3_part<-w3[which(w3$RBP == rr),]
    w3_part<-w3_part[which(w3_part$ID %in% Events_DS[,1]),]

    Events_DS_RBP<-as.data.frame(Events_DS[w3_part$ID,])
    Events_DS_RBP$w1<-w3_part$weight1
    Events_DS_RBP$w2<-w3_part$weight2
    Events_DS_RBP$rank<-as.numeric(Events_DS_RBP$abs_dPSI)*as.numeric(Events_DS_RBP$w1)*as.numeric(Events_DS_RBP$w2)
    Events_DS_RBP<-Events_DS_RBP[order(Events_DS_RBP$rank,decreasing = T),]
    w3_part<-w3_part[order(as.numeric(w3_part$weight_all),decreasing = T),]

    #NES1
    w3_part1<-w3_part[which(as.numeric(w3_part$score)>0.5),]
    gmt_nes1<-list(unlist(w3_part1$ID))
    names(gmt_nes1)[1]<-paste0(rr,"_NES1")
    ranks_nes1<-Events_DS_RBP$rank
    names(ranks_nes1)<-Events_DS_RBP$Events

    #NES2
    Events_DS_RBP_2<-Events_DS_RBP[which((as.numeric(Events_DS_RBP$pvalue)<DS_pvalue)&(as.numeric(Events_DS_RBP$abs_dPSI)>DS_dPSI)),]
    gmt_nes2<-list(unlist(Events_DS_RBP_2$Events))
    names(gmt_nes2)[1]<-paste0(rr,"_NES2")
    ranks_nes2<-w3_part$weight_all
    names(ranks_nes2)<-w3_part$ID

    if (length(gmt_nes1[[1]])==0){
        stop("At this threshold, no ",rr,"-regulated events were found.")
    }else{
      if (length(gmt_nes2[[1]])==0){
        stop("At this threshold, no Differential splicing events were found.")
      }else{
        nes1<-gsea(pathways =gmt_nes1, stats =ranks_nes1,scoreType = "pos")
        nes2<-gsea(pathways =gmt_nes2, stats =ranks_nes2,scoreType = "pos")

        nes<-data.table::data.table(NA)
        nes[,RBP:=strsplit(as.character(nes1[1,1]),"_")[[1]][1]]
        nes[,nes1_size:=nes1[1,7]]
        nes[,nes1_es:=nes1[1,2]]
        nes[,nes1_nes:=nes1[1,3]]
        nes[,nes1_p:=nes1[1,4]]
        nes[,nes2_size:=nes2[1,7]]
        nes[,nes2_es:=nes2[1,2]]
        nes[,nes2_nes:=nes2[1,3]]
        nes[,nes2_p:=nes2[1,4]]
        overlap<-intersect(unlist(gmt_nes1),unlist(gmt_nes2))
        nes[,overlap:=length(overlap<-intersect(unlist(gmt_nes1),unlist(gmt_nes2)))]
        nes[,total_size:=nrow(w3_part)]
        nes<-nes[,-1]
        n1<-nes[,overlap]
        n2<-nes[,nes1_size]-nes[,overlap]
        n3<-nes[,nes2_size]-nes[,overlap]
        n4<-nes[,total_size]+nes[,overlap]-nes[,nes1_size]-nes[,nes2_size]
        tableR<-matrix(c(n1,n2,n3,n4),nrow = 2,ncol = 2)
        fish<-fisher.test(tableR,alternative = "greater")
        if ((is.na(fish$estimate))|(is.infinite(as.numeric(fish$estimate)))){
          if (is.na(fish$estimate)){
            nes[,OR:=0]
            nes[,pval:=0]
          }else{
            nes[,OR:=200]
            nes[,pval:=0]
          }

        }else{
          if (fish$estimate<200){
            nes[,OR:=fish$estimate]
            nes[,pval:=fish$p.value]
          }else{
            nes[,OR:=200]
            nes[,pval:=fish$p.value]
          }
        }
        nes[,score3:=as.numeric(nes[,nes1_nes])*as.numeric(nes[,nes2_nes])*sqrt(as.numeric(nes[,OR]))]
        nes<-as.matrix(nes)
        # nes<-do.call(cbind,nes)
        return(nes)
      }
    }
  }

  parallel::stopCluster(cll)

  rbp_uni_mat<-do.call(cbind,rbp_uni_mat)

  rbp_uni_mat_2<-merge(as.matrix(rbp_uni_mat),as.matrix(ES_mat_RBP))

  rbp_uni_mat_2$score3<-as.numeric(rbp_uni_mat_2$score1_nor)*as.numeric(rbp_uni_mat_2$score3)

  rbp_uni_mat_2<-rbp_uni_mat_2[order(as.numeric(rbp_uni_mat_2[,ncol(rbp_uni_mat_2)]),decreasing = T),]

  rbp_uni_mat_3<-(rbp_uni_mat_2[,-3])

  if (file.exists(paste0(path_use,"w2.txt"))) file.remove(paste0(path_use,"w2.txt"))
  score_matrix_all<-matrix(NA,nrow = 1,ncol = 22)
  if (is.null(rbp_interested)){
    score_matrix_all[1,1]<-NA
    score_matrix_all[1,2]<-NA
  }else{
    score_matrix_all[1,1]<-rbp_interested
    score_matrix_all[1,2]<-ifelse(rbp_interested %in% rbp_uni_mat_3[,1],which(rbp_uni_mat_3[,1]==rbp_interested),NA)
     }
  for (q in 1:10) {
    score_matrix_all[1,(2*q+1)]<-rbp_uni_mat_3[q,1]
    score_matrix_all[1,(2*q+2)]<-rbp_uni_mat_3[q,ncol(rbp_uni_mat_3)]
  }
  colnames(score_matrix_all)<-c("rbp_interested","rank",paste0(c("RBP","rank"),rep(1:10,each=2)))
  data.table::fwrite(as.data.frame(score_matrix_all),file = paste0(path_use,"result_top10.txt"),
                     row.names = F,col.names = T,quote = F,sep = "\t")
  data.table::fwrite(as.data.frame(rbp_uni_mat_3),file = paste0(path_use,"result_tab_simple.txt"),
         row.names = F,col.names = T,quote = F,sep = "\t")
  save(rbp_uni_mat_2,file = paste0(path_use,"result_tab_all.RData"))
  if (result_type == "Top10"){
    return(score_matrix_all)
  }else {
    if (result_type == "tab_simple"){
      tmp<-rbp_uni_mat_3[,c("RBP","score1_nor","nes1_nes","nes2_nes","OR","score3")]
      colnames(tmp)<-c("RBP","D","NES1","NES2","odds","MRAS_Score")
      return(tmp)
    }else{
      tmp<-rbp_uni_mat_2[,c("RBP","Events","score1_nor","nes1_nes","nes2_nes","OR","score3")]
      colnames(tmp)<-c("RBP","Targets","D","NES1","NES2","odds","MRAS_Score")
      return(tmp)
    }
  }
}

