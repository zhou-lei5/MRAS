
#' MRAS_net_single
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param num1,num2 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'                  Input requirements: if the input is less than 1, the input is considered as the sample proportion, and greater than 1, the input is considered as the number of samples.
#'                  When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param cor_cutoff The cutoff of co-regulation of RBP to event j, with a default value of 0.3.
#' @param cor_p_cutoff The cutoff of p-value of co-regulation of RBP to event j, with a default value of 0.05.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param BS if input the binding data.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param path_use The path to the file used to store the output.
#' @param method The method to calculate correlation.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network.
#' @export
#'
#'
MRAS_net_single<-function(expr,psi,num1 = 0.1,num2 = 0.1,method="spearman",
                          cor_cutoff=0.3,cor_p_cutoff=0.05,
                          dpsi_network_threshold = 0.1,BS = NULL,
                          threads = 2,path_use = path_use){
  if ((num1>0.5)&(num1<1)) stop("Wrong input : num!")
  if ((num2>0.5)&(num2<1)) stop("Wrong input : num!")
  if (num1<1) num1<-ceiling(ncol(expr)*num1)
  if (num2<1) num2<-ceiling(ncol(expr)*num2)
  use_mat1<-t(apply(expr, 1, function(q){
    use<-c(order(q,decreasing=TRUE)[1:num1],
           order(q,decreasing=FALSE)[1:num1])
    return(use)
  }))
  use_mat2<-t(apply(expr, 1, function(q){
    use<-c(order(q,decreasing=TRUE)[1:num2],
           order(q,decreasing=FALSE)[1:num2])
    return(use)
  }))
  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  tj<-function(gene_ID_use,sep=","){
    if (is.na(gene_ID_use[1])){
      return(NA)
    }else{
      gene_ID_use<-gene_ID_use[complete.cases(gene_ID_use)]
      name<-paste(gene_ID_use,collapse = sep)
      return(name)
    }
  }
  i=0
  re<-foreach (i = 1:nrow(expr), .combine = "c") %dopar% {
    t1<-matrix(as.numeric(unlist(expr[i,use_mat1[i,]])),nrow = 1)
    expr_p<-t.test(t1[1:(ncol(t1)/2)],t1[((ncol(t1)/2)+1):(ncol(t1))])
    t2<-psi[,use_mat1[i,]]
    rr<-cor(t(t1),t(t2),method = method)
    df<-ncol(t1)-2
    rownames(rr)<-rownames(expr)[i]
    pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
    pp[which(pp>1)]<-2-pp[which(pp>1)]
    t3<-psi[,use_mat2[i,]]
    psi1<-apply(t3[,1:(ncol(t3)/2)], 1, mean)
    psi2<-apply(t3[,((ncol(t3)/2)+1):(ncol(t3))], 1, mean)
    dpsi<-psi1-psi2
    abs_dpsi<-abs(as.numeric(dpsi))
    abs_dpsi[which(abs_dpsi==1)]<-0
    logp<-(-log10(pp))
    if_cor<-ifelse((abs(rr)>cor_cutoff)&(pp<cor_p_cutoff),1,0)
    P<-apply(t2,1,function(x){
      x1<-as.numeric(x[1:(ncol(t2)/2)])
      x2<-as.numeric(x[((ncol(t2)/2)+1):(ncol(t2))])
      va1<-var(x1)
      va2<-var(x2)
      if (va1==0) x1[1]<-x1[1]+0.001
      if (va2==0) x2[1]<-x2[1]+0.001
      p<-t.test(x1,x2)
      return(p$p.value)
    })
    if_P<-(P<0.05)
    if_real<-ifelse((abs_dpsi>dpsi_network_threshold)&(abs_dpsi<1)&(if_P),1,0)
    use<-as.numeric(abs(rr)*if_cor*sqrt(abs_dpsi))
    corr_mat<-data.frame(rbp = rownames(expr)[i],
                         event = colnames(rr)[col(rr)],
                         cor = c(unlist(rr)),
                         p = c(unlist(pp)),
                         abs_dpsi = c(abs_dpsi),
                         logp = c(logp),
                         if_cor = c(if_cor),
                         if_real = c(if_real),
                         use = c(use))
    if (!is.null(BS)) {
      BS_part<-BS[which(BS[,1]==rownames(expr)[i]),]
      colnames(BS_part)<-c("rbp","event","BS_score")
      corr_mat<-dplyr::left_join(corr_mat,BS_part)
      # a<-stats::glm(if_real ~ scale(use) + BS_score ,family = binomial(link = "logit") ,data =corr_mat)
      a<-stats::glm(if_real ~ use + BS_score ,family = binomial(link = "logit") ,data =corr_mat)
      rbp_net_mat<-data.frame(rbp = rownames(expr)[i],
                              expr_p = expr_p$p.value,
                              BS = "T",
                              logit1 = a$coefficients[1],
                              logit2 = a$coefficients[2],
                              logit3 = a$coefficients[3])
    }else{
      # a<-stats::glm(if_real ~ scale(use) ,family = binomial(link = "logit") ,data =corr_mat)
      a<-stats::glm(if_real ~ use ,family = binomial(link = "logit") ,data =corr_mat)
      rbp_net_mat<-data.frame(rbp = rownames(expr)[i],
                              expr_p = expr_p$p.value,
                              BS = "F",
                              logit1 = a$coefficients[1],
                              logit2 = a$coefficients[2],
                              logit3 = 0)

    }
    return(list(rbp_net_mat=rbp_net_mat,corr_mat=corr_mat))
  }
  rbp_net_mat <- do.call(rbind, lapply(re, function(x){
    if (dim(x)[2]==6) return(x)
  }))
  corr_mat <- do.call(rbind, lapply(re, function(x){
    if (dim(x)[2]==9) return(x)
  }))
  parallel::stopCluster(cll)
  return(list(rbp_net_mat=rbp_net_mat,corr_mat=corr_mat))
}
#' MRAS_net_group
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param num1,num2 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'                  Input requirements: if the input is less than 1, the input is considered as the sample proportion, and greater than 1, the input is considered as the number of samples.
#'                  When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param cor_cutoff The cutoff of co-regulation of RBP to event j, with a default value of 0.3.
#' @param cor_p_cutoff The cutoff of p-value of co-regulation of RBP to event j, with a default value of 0.05.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param BS if input the binding data.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param path_use The path to the file used to store the output.
#' @param MRAS_net_single_re The output of function `MRAS_net_single()`.
#' @param string_net The relationship between RBPs.
#' @param method The method to calculate correlation.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network (After considering the relationship between RBPs).
#' @export
#'
MRAS_net_group<-function(expr = expr,psi = psi,num1 = 0.1,num2 = 0.1,method="spearman",BS=NULL,
                         cor_cutoff=0.3,cor_p_cutoff=0.05,
                         MRAS_net_single_re = MRAS_net_single_re,
                         dpsi_network_threshold = 0.1,
                         string_net = string_net,
                         threads = 2,path_use = path_use){
  rbp_net_mat<-MRAS_net_single_re$rbp_net_mat
  rbp_corr<-MRAS_net_single_re$corr_mat

  rbp_corr_deal<-rbp_corr[which((abs(rbp_corr$cor)>cor_cutoff)&(rbp_corr$p<cor_p_cutoff)),]
  rbp_corr_deal$cor<-abs(rbp_corr_deal$cor)
  rbp_corr_deal_tab<-table(rbp_corr_deal$col)
  rbp_corr_deal_remove<-names(rbp_corr_deal_tab)[which(rbp_corr_deal_tab == 1)]
  rbp_corr_deal<-rbp_corr_deal[!(rbp_corr_deal$event %in% rbp_corr_deal_remove),]

  rbp_corr_mat<-reshape2::dcast(rbp_corr_deal,rbp~event,value.var = "cor",fill = 0)
  rownames(rbp_corr_mat)<-rbp_corr_mat[,1]
  rbp_corr_mat<-rbp_corr_mat[,-1]
  rbp_corr_mat<-rbp_corr_mat[intersect(rownames(string_net),rownames(rbp_corr_mat)),]
  rbp_corr_mat<-as.matrix(rbp_corr_mat)
  expr_rbp<-as.matrix(expr[intersect(rownames(string_net),rownames(rbp_corr_mat)),])

  string_net_use<-string_net[intersect(rownames(string_net),rownames(rbp_corr_mat)),]
  string_net_use<-string_net_use[,intersect(rownames(string_net),rownames(rbp_corr_mat))]
  #+1
  expr_rbp_log<-as.matrix(apply(expr_rbp,2,function(x){
    return(log2(as.numeric(x)+1))
  }))
  rownames(expr_rbp_log)<-rownames(expr_rbp)
  if ((num1>0.5)&(num1<1)) stop("Wrong input : num!")
  if ((num2>0.5)&(num2<1)) stop("Wrong input : num!")
  if (num1<1) num1<-ceiling(ncol(expr)*num1)
  if (num2<1) num2<-ceiling(ncol(expr)*num2)
  use_mat1<-t(apply(expr, 1, function(q){
    use<-c(order(q,decreasing=TRUE)[1:num1],
           order(q,decreasing=FALSE)[1:num1])
    return(use)
  }))
  use_mat2<-t(apply(expr, 1, function(q){
    use<-c(order(q,decreasing=TRUE)[1:num2],
           order(q,decreasing=FALSE)[1:num2])
    return(use)
  }))

  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  filter_expr<-function(expr_matrix,n){
    results<-matrix(TRUE,nrow = nrow(expr_matrix),ncol = 1)
    for (i in 1:nrow(expr_matrix)) {
      length0<-length(which(as.numeric(expr_matrix[i,])<1))
      results[i,1]<-ifelse((length0)>n,FALSE,TRUE)
      # cat(round(i/nrow(expr_matrix),2),"\t")
    }
    return(results[,1])
  }
  MRAS_net_group_pre1<-function(rbp,psi_deal,use_mat1 = use_mat1,use_mat2 = use_mat2,method=method,
                                cor_cutoff=cor_cutoff,cor_p_cutoff=cor_p_cutoff,
                                num1 = 0.5,num2 = 0.5,dpsi_network_threshold = 0.1,
                                rbp_bind_event){

    t1<-rbp_bind_event[,use_mat1[rbp,]]
    t2<-psi_deal[,use_mat1[rbp,]]
    t1_2<-cbind(t1,t2)
    rr<-apply(t1_2, 1, function(x){
      return(cor(x[1:ncol(t1)],x[(ncol(t1)+1):ncol(t1_2)],method = method))
    })
    df<-ncol(t1)-2
    pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
    pp[which(pp>1)]<-2-pp[which(pp>1)]
    t3<-psi_deal[,use_mat2[rbp,]]
    psi1<-apply(t3[,1:(ncol(t3)/2)], 1, mean)
    psi2<-apply(t3[,((ncol(t3)/2)+1):(ncol(t3))], 1, mean)
    dpsi<-psi1-psi2
    abs_dpsi<-abs(as.numeric(dpsi))
    abs_dpsi<-abs(as.numeric(dpsi))
    abs_dpsi[which(abs_dpsi==1)]<-0
    logp<-(-log10(pp))
    if_cor<-ifelse((abs(rr)>cor_cutoff)&(pp<cor_p_cutoff),1,0)
    P<-apply(t2,1,function(x){
      x1<-as.numeric(x[1:(ncol(t2)/2)])
      x2<-as.numeric(x[((ncol(t2)/2)+1):(ncol(t2))])
      va1<-var(x1)
      va2<-var(x2)
      if (va1==0) x1[1]<-x1[1]+0.001
      if (va2==0) x2[1]<-x2[1]+0.001
      p<-t.test(x1,x2)
      return(p$p.value)
    })
    if_P<-(P<0.05)
    if_real<-ifelse((abs_dpsi>dpsi_network_threshold)&(abs_dpsi<1)&(if_P),1,0)
    use<-as.numeric(abs(rr)*if_cor*sqrt(abs_dpsi))
    corr_mat<-data.frame(rbp = rbp,
                         event = rownames(rbp_bind_event),
                         cor_new = c(unlist(rr)),
                         p_new = c(unlist(pp)),
                         abs_dpsi_new = c(abs_dpsi),
                         logp_new = c(logp),
                         if_cor_new = c(if_cor),
                         if_real_new = c(if_real),
                         use_new = c(use))
    return(corr_mat)
  }

  MRAS_net_group_pre2<-function(expr,re_new_all_part,rbp_corr,BS=NULL){
    re_new_all_part<-as.data.frame(re_new_all_part)
    rbp_corr<-as.data.frame(rbp_corr)
    network_re<-dplyr::left_join(rbp_corr,re_new_all_part,by = c("rbp", "event"))
    network_re$use_new[is.na(network_re$use_new)]<-network_re$use[is.na(network_re$use_new)]
    tj<-function(gene_ID_use,sep=","){
      if (is.na(gene_ID_use[1])){
        return(NA)
      }else{
        gene_ID_use<-gene_ID_use[complete.cases(gene_ID_use)]
        name<-paste(gene_ID_use,collapse = sep)
        return(name)
      }
    }
    network_re$use_new<-as.numeric(network_re$use_new)
    i=0
    re_final<-foreach::foreach (i = 1:nrow(expr), .combine = "rbind") %dopar% {
      network_re_part<-network_re[which(network_re[,1]==rownames(expr)[i]),]
      if (!is.null(BS)){
        BS_part<-BS[which(BS[,1]==rownames(expr)[i]),]
        colnames(BS_part)<-c("rbp","event","BS_score")
        network_re_part<-dplyr::left_join(network_re_part,BS_part)
        network_re_part[is.na(network_re_part)]<-0
        # a<-stats::glm(if_real ~ scale(use) + scale(use_new) + BS_score ,family = binomial(link = "logit") ,data =network_re_part)
        a<-stats::glm(if_real ~ use + use_new + BS_score ,family = binomial(link = "logit") ,data =network_re_part)
        rbp_net_mat_group<-data.frame(rbp = rownames(expr)[i],
                                      BS = "T",
                                      logit1_group = a$coefficients[1],
                                      logit2_group = a$coefficients[2],
                                      logit3_group = a$coefficients[3],
                                      logit4_group = a$coefficients[4])
      }else{
        # a<-stats::glm(if_real ~ scale(use) + scale(use_new) ,family = binomial(link = "logit") ,data =network_re_part)
        a<-stats::glm(if_real ~ use + use_new ,family = binomial(link = "logit") ,data =network_re_part)
        rbp_net_mat_group<-data.frame(rbp = rownames(expr)[i],
                                      BS = "F",
                                      logit1_group = a$coefficients[1],
                                      logit2_group = a$coefficients[2],
                                      logit3_group = a$coefficients[3],
                                      logit4_group = 0)
      }
      return(rbp_net_mat_group)
    }
    return(list(network_re = network_re,re_final = re_final))
  }
  rrr=0
  re_new_all<-foreach (rrr = 1:ncol(string_net_use),.combine = "rbind",.errorhandling = "pass") %dopar% {
    # for(rr in 1:ncol(string_net)){
    rbp<-colnames(string_net_use)[rrr]
    string_net_part<-string_net_use[,rrr]
    rbp_corr_mat_part<-rbp_corr_mat*string_net_part
    string_net_value<-unlist(apply(rbp_corr_mat_part, 2, function(x){
      return(prod(x[which(x!=0)]))
    }))
    rbp_corr_mat_part[which(rbp_corr_mat_part!=0)]<-1
    rbp_bind_event<-t(rbp_corr_mat_part) %*% expr_rbp_log
    rbp_bind_event<-as.matrix(rbp_bind_event)
    new_matrix <- 2^rbp_bind_event
    new_matrix<-new_matrix -1
    rbp_bind_event<-new_matrix*string_net_value
    ff1<-filter_expr(rbp_bind_event,ncol(rbp_bind_event)*0.8)
    rbp_bind_event<-rbp_bind_event[ff1,]
    psi_deal<-psi[rownames(rbp_bind_event),]
    if (nrow(psi_deal)>0){
      re_new<-MRAS_net_group_pre1(rbp = rbp,psi_deal = psi_deal,use_mat1 = use_mat1,use_mat2 = use_mat2,method = method,
                                  cor_cutoff=cor_cutoff,cor_p_cutoff=cor_p_cutoff,
                                  num1 = num1,num2 = num2,dpsi_network_threshold = dpsi_network_threshold,
                                  rbp_bind_event = rbp_bind_event)
      re_new<-as.matrix(re_new)
      return(re_new)
    }else{
      err<-matrix(0,1,9)
      colnames(err)<-c("rbp","event","cor_new","p_new","abs_dpsi_new","logp_new","if_cor_new","if_real_new","use_new")
      return(err)
    }
    return(re_new)
  }
  re_new_all<-as.matrix(re_new_all)
  re_new_all<-re_new_all[!(is.na(re_new_all[,9])),]
  re_new_all<-re_new_all[which(as.numeric(re_new_all[,9]) != 0),]
  network_re_final<-MRAS_net_group_pre2(expr = expr,re_new_all_part = re_new_all[,c(1,2,9)],
                                        rbp_corr = rbp_corr,BS = BS)
  parallel::stopCluster(cll)
  rbp_net_mat_final<-dplyr::left_join(rbp_net_mat,network_re_final$re_final)
  return(list(rbp_net_mat_group = rbp_net_mat_final,rbp_corr_group = network_re_final$network_re))
}
#' MRAS_net_single_work
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param num1,num2 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'                  Input requirements: if the input is less than 1, the input is considered as the sample proportion, and greater than 1, the input is considered as the number of samples.
#'                  When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential splicing events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param cor_cutoff The cutoff of co-regulation of RBP to event j, with a default value of 0.3.
#' @param cor_p_cutoff The cutoff of p-value of co-regulation of RBP to event j, with a default value of 0.05.
#' @param BS if input the binding data.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param path_use The path to the file used to store the output.
#' @param method The method to calculate correlation.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network.
#' @export
#'

MRAS_net_single_work<-function(expr,psi,num1 = 0.5,num2 = 0.5,method="spearman",
                               cor_cutoff=0.3,cor_p_cutoff=0.05,
                               dpsi_network_threshold = 0.1,BS = NULL,
                               threads = 2,path_use = path_use){
  if ((num1>0.5)&(num1<1)) stop("Wrong input : num!")
  if ((num2>0.5)&(num2<1)) stop("Wrong input : num!")
  if (num1<1) num1<-ceiling(ncol(expr)*num1)
  if (num2<1) num2<-ceiling(ncol(expr)*num2)
  use_mat1<-t(apply(expr, 1, function(q){
    use<-c(order(q,decreasing=TRUE)[1:num1],
           order(q,decreasing=FALSE)[1:num1])
    return(use)
  }))
  use_mat2<-t(apply(expr, 1, function(q){
    use<-c(order(q,decreasing=TRUE)[1:num2],
           order(q,decreasing=FALSE)[1:num2])
    return(use)
  }))
  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  tj<-function(gene_ID_use,sep=","){
    if (is.na(gene_ID_use[1])){
      return(NA)
    }else{
      gene_ID_use<-gene_ID_use[complete.cases(gene_ID_use)]
      name<-paste(gene_ID_use,collapse = sep)
      return(name)
    }
  }
  i=0
  re<-foreach (i = 1:nrow(expr), .combine = "rbind") %dopar% {
    t1<-matrix(unlist(expr[i,use_mat1[i,]]),nrow = 1)
    expr_p<-t.test(t1[1:(ncol(t1)/2)],t1[((ncol(t1)/2)+1):(ncol(t1))])
    t2<-psi[,use_mat1[i,]]
    rr<-cor(t(t1),t(t2),method = method)
    df<-ncol(t1)-2
    rownames(rr)<-rownames(expr)[i]
    pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
    pp[which(pp>1)]<-2-pp[which(pp>1)]
    t3<-psi[,use_mat2[i,]]
    psi1<-apply(t3[,1:(ncol(t3)/2)], 1, mean)
    psi2<-apply(t3[,((ncol(t3)/2)+1):(ncol(t3))], 1, mean)
    dpsi<-psi1-psi2
    abs_dpsi<-abs(as.numeric(dpsi))
    abs_dpsi[which(abs_dpsi==1)]<-0
    logp<-(-log10(pp))
    if_cor<-ifelse((abs(rr)>cor_cutoff)&(pp<cor_p_cutoff),1,0)
    P<-apply(t2,1,function(x){
      x1<-as.numeric(x[1:(ncol(t2)/2)])
      x2<-as.numeric(x[((ncol(t2)/2)+1):(ncol(t2))])
      va1<-var(x1)
      va2<-var(x2)
      if (va1==0) x1[1]<-x1[1]+0.001
      if (va2==0) x2[1]<-x2[1]+0.001
      p<-t.test(x1,x2)
      return(p$p.value)
    })
    if_P<-(P<0.05)
    if_real<-ifelse((abs_dpsi>dpsi_network_threshold)&(abs_dpsi<1)&(if_P),1,0)
    use<-as.numeric(abs(rr)*if_cor*sqrt(abs_dpsi))
    corr_mat<-data.frame(rbp = rownames(expr)[i],
                         event = colnames(rr)[col(rr)],
                         cor = c(unlist(rr)),
                         p = c(unlist(pp)),
                         abs_dpsi = c(abs_dpsi),
                         logp = c(logp),
                         if_cor = c(if_cor),
                         if_real = c(if_real),
                         use = c(use))
    return(corr_mat)
  }
  parallel::stopCluster(cll)
  return(re)
}
#' MRAS_net_group_work
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param num1,num2 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'                  Input requirements: if the input is less than 1, the input is considered as the sample proportion, and greater than 1, the input is considered as the number of samples.
#'                  When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param cor_cutoff The cutoff of co-regulation of RBP to event j, with a default value of 0.3.
#' @param cor_p_cutoff The cutoff of p-value of co-regulation of RBP to event j, with a default value of 0.05.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param BS if input the binding data.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param path_use The path to the file used to store the output.
#' @param string_net The relationship between RBPs.
#' @param rbp_corr_work The output of function `MRAS_net_single_work()`.
#' @param method The method to calculate correlation.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network (After considering the relationship between RBPs).
#' @export
#'



MRAS_net_group_work<-function(expr = expr,psi = psi,num1 = 0.5,num2 = 0.5,
                              cor_cutoff=0.3,cor_p_cutoff=0.05,
                              method="spearman",BS=NULL,
                              rbp_corr_work = rbp_corr_work,
                              dpsi_network_threshold = 0.1,
                              string_net = string_net,
                              threads = 2,path_use = path_use){
  rbp_corr<-rbp_corr_work

  rbp_corr_deal<-rbp_corr[which((abs(rbp_corr$cor)>cor_cutoff)&(rbp_corr$p<cor_p_cutoff)),]
  rbp_corr_deal$cor<-abs(rbp_corr_deal$cor)
  rbp_corr_deal_tab<-table(rbp_corr_deal$col)
  rbp_corr_deal_remove<-names(rbp_corr_deal_tab)[which(rbp_corr_deal_tab == 1)]
  rbp_corr_deal<-rbp_corr_deal[!(rbp_corr_deal$event %in% rbp_corr_deal_remove),]

  rbp_corr_mat<-reshape2::dcast(rbp_corr_deal,rbp~event,value.var = "cor")
  rbp_corr_mat[is.na(rbp_corr_mat)]<-0
  rownames(rbp_corr_mat)<-rbp_corr_mat[,1]
  rbp_corr_mat<-rbp_corr_mat[,-1]
  rbp_corr_mat<-rbp_corr_mat[intersect(rownames(string_net),rownames(rbp_corr_mat)),]
  rbp_corr_mat<-as.matrix(rbp_corr_mat)
  expr_rbp<-as.matrix(expr[intersect(rownames(string_net),rownames(rbp_corr_mat)),])

  string_net_use<-string_net[intersect(rownames(string_net),rownames(rbp_corr_mat)),]
  string_net_use<-string_net_use[,intersect(rownames(string_net),rownames(rbp_corr_mat))]
  #+1
  expr_rbp_log<-as.matrix(apply(expr_rbp,2,function(x){
    return(log2(as.numeric(x)+1))
  }))
  rownames(expr_rbp_log)<-rownames(expr_rbp)
  if ((num1>0.5)&(num1<1)) stop("Wrong input : num!")
  if ((num2>0.5)&(num2<1)) stop("Wrong input : num!")
  if (num1<1) num1<-ceiling(ncol(expr)*num1)
  if (num2<1) num2<-ceiling(ncol(expr)*num2)
  use_mat1<-t(apply(expr, 1, function(q){
    use<-c(order(q,decreasing=TRUE)[1:num1],
           order(q,decreasing=FALSE)[1:num1])
    return(use)
  }))
  use_mat2<-t(apply(expr, 1, function(q){
    use<-c(order(q,decreasing=TRUE)[1:num2],
           order(q,decreasing=FALSE)[1:num2])
    return(use)
  }))

  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  filter_expr<-function(expr_matrix,n){
    results<-matrix(TRUE,nrow = nrow(expr_matrix),ncol = 1)
    i<-0
    for (i in 1:nrow(expr_matrix)) {
      length0<-length(which(as.numeric(expr_matrix[i,])<1))
      results[i,1]<-ifelse((length0)>n,FALSE,TRUE)
      # cat(round(i/nrow(expr_matrix),2),"\t")
    }
    return(results[,1])
  }
  MRAS_net_group_pre1<-function(rbp,psi_deal,use_mat1 = use_mat1,use_mat2 = use_mat2,method = method,
                                cor_cutoff=cor_cutoff,cor_p_cutoff=cor_p_cutoff,
                                num1 = 0.5,num2 = 0.5,dpsi_network_threshold = 0.1,
                                rbp_bind_event){
    t1<-rbp_bind_event[,use_mat1[rbp,]]
    t2<-psi_deal[,use_mat1[rbp,]]
    t1_2<-cbind(t1,t2)
    rr<-apply(t1_2, 1, function(x){
      return(cor(x[1:ncol(t1)],x[(ncol(t1)+1):ncol(t1_2)],method = method))
    })
    df<-ncol(t1)-2
    pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
    pp[which(pp>1)]<-2-pp[which(pp>1)]
    t3<-psi_deal[,use_mat2[rbp,]]
    psi1<-apply(t3[,1:(ncol(t3)/2)], 1, mean)
    psi2<-apply(t3[,((ncol(t3)/2)+1):(ncol(t3))], 1, mean)
    dpsi<-psi1-psi2
    abs_dpsi<-abs(as.numeric(dpsi))
    abs_dpsi[which(abs_dpsi==1)]<-0
    logp<-(-log10(pp))
    if_cor<-ifelse((abs(rr)>cor_cutoff)&(pp<cor_p_cutoff),1,0)
    P<-apply(t2,1,function(x){
      x1<-as.numeric(x[1:(ncol(t2)/2)])
      x2<-as.numeric(x[((ncol(t2)/2)+1):(ncol(t2))])
      va1<-var(x1)
      va2<-var(x2)
      if (va1==0) x1[1]<-x1[1]+0.001
      if (va2==0) x2[1]<-x2[1]+0.001
      p<-t.test(x1,x2)
      return(p$p.value)
    })
    if_P<-(P<0.05)
    if_real<-ifelse((abs_dpsi>dpsi_network_threshold)&(abs_dpsi<1)&(if_P),1,0)
    use<-as.numeric(abs(rr)*if_cor*sqrt(abs_dpsi))
    corr_mat<-data.frame(rbp = rbp,
                         event = rownames(rbp_bind_event),
                         cor_new = c(unlist(rr)),
                         p_new = c(unlist(pp)),
                         abs_dpsi_new = c(abs_dpsi),
                         logp_new = c(logp),
                         if_cor_new = c(if_cor),
                         if_real_new = c(if_real),
                         use_new = c(use))
    return(corr_mat)
  }
  rrr=0
  re_new_all<-foreach (rrr = 1:ncol(string_net_use),.combine = "rbind",.errorhandling = "pass") %dopar% {
    # for(rr in 1:ncol(string_net)){
    rbp<-colnames(string_net_use)[rrr]
    # cat(rbp,"\t")
    string_net_part<-string_net_use[,rrr]
    rbp_corr_mat_part<-rbp_corr_mat*string_net_part
    string_net_value<-unlist(apply(rbp_corr_mat_part, 2, function(x){
      return(prod(x[which(x!=0)]))
    }))
    rbp_corr_mat_part[which(rbp_corr_mat_part!=0)]<-1
    rbp_bind_event<-t(rbp_corr_mat_part) %*% expr_rbp_log
    rbp_bind_event<-as.matrix(rbp_bind_event)
    new_matrix <- 2^rbp_bind_event
    new_matrix<-new_matrix -1
    rbp_bind_event<-new_matrix*string_net_value
    ff1<-filter_expr(rbp_bind_event,ncol(rbp_bind_event)*0.8)
    rbp_bind_event<-rbp_bind_event[ff1,]
    psi_deal<-psi[rownames(rbp_bind_event),]
    if (nrow(psi_deal)>0){
      re_new<-MRAS_net_group_pre1(rbp = rbp,psi_deal = psi_deal,use_mat1 = use_mat1,use_mat2 = use_mat2,method = method,
                                  cor_cutoff=cor_cutoff,cor_p_cutoff=cor_p_cutoff,
                                  num1 = num1,num2 = num2,dpsi_network_threshold = dpsi_network_threshold,
                                  rbp_bind_event = rbp_bind_event)
      re_new<-as.matrix(re_new)
      return(re_new)
    }else{
      err<-matrix(0,1,9)
      colnames(err)<-c("rbp","event","cor_new","p_new","abs_dpsi_new","logp_new","if_cor_new","if_real_new","use_new")
      return(err)
    }
    return(re_new)
  }
  parallel::stopCluster(cll)
  re_new_all<-as.matrix(re_new_all)
  re_new_all<-re_new_all[which(as.numeric(re_new_all[,9]) != 0),]
  re_new_all_part<-as.data.frame(re_new_all[,c(1,2,9)])
  rbp_corr<-as.data.frame(rbp_corr)
  network_re<-dplyr::left_join(rbp_corr,re_new_all_part,by = c("rbp", "event"))
  network_re$use_new[is.na(network_re$use_new)]<-network_re$use[is.na(network_re$use_new)]
  return(network_re)
}

#' MRAS_net_single_sc
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param cor_cutoff The cutoff of co-regulation of RBP to event j, with a default value of 0.
#' @param cor_p_cutoff The cutoff of p-value of co-regulation of RBP to event j, with a default value of 0.05.
#' @param num1,num2 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'                  Input requirements: if the input is less than 1, the input is considered as the sample proportion, and greater than 1, the input is considered as the number of samples.
#'                  When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param BS if input the binding data.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param path_use The path to the file used to store the output.
#' @param method The method to calculate correlation.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network.
#' @export
#'
MRAS_net_single_sc<-function(expr,psi,num1 = 0.1,num2 = 0.1,method="spearman",
                             cor_cutoff=0,cor_p_cutoff=0.05,
                             dpsi_network_threshold = 0.1,BS = NULL,
                             threads = 2,path_use = path_use){
  # browser()
  if ((num1>0.5)&(num1<1)) stop("Wrong input : num!")
  if ((num2>0.5)&(num2<1)) stop("Wrong input : num!")
  # if (num1<1) num1<-ceiling(ncol(expr)*num1)
  # if (num2<1) num2<-ceiling(ncol(expr)*num2)
  use_mat1<-apply(expr, 1, function(q){
    len<-length(q[q!=0])
    num_use<-ceiling(len*num1)
    if (num_use < 5) stop("Please filter the RBP expr mat.")
    use<-c(order(q[q!=0],decreasing=TRUE)[1:num_use],
           order(q[q!=0],decreasing=FALSE)[1:num_use])
    return(which(q!=0)[use])
  })
  use_mat2<-apply(expr, 1, function(q){
    len<-length(q[q!=0])
    num_use<-ceiling(len*num2)
    if (num_use < 5) stop("Please filter the RBP expr mat.")
    use<-c(order(q[q!=0],decreasing=TRUE)[1:num_use],
           order(q[q!=0],decreasing=FALSE)[1:num_use])
    return(which(q!=0)[use])
  })
  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  tj<-function(gene_ID_use,sep=","){
    if (is.na(gene_ID_use[1])){
      return(NA)
    }else{
      gene_ID_use<-gene_ID_use[complete.cases(gene_ID_use)]
      name<-paste(gene_ID_use,collapse = sep)
      return(name)
    }
  }
  i=0
  re<-foreach (i = 1:nrow(expr), .combine = "c") %dopar% {
    t1<-matrix(as.numeric(unlist(expr[i,use_mat1[[i]]])),nrow = 1)
    expr_p<-t.test(t1[1:(ncol(t1)/2)],t1[((ncol(t1)/2)+1):(ncol(t1))])
    t2<-psi[,use_mat1[[i]]]
    rr<-cor(t(t1),t(t2),method = method)
    df<-ncol(t1)-2
    rownames(rr)<-rownames(expr)[i]
    pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
    pp[which(pp>1)]<-2-pp[which(pp>1)]
    t3<-psi[,use_mat2[[i]]]
    psi1<-apply(t3[,1:(ncol(t3)/2)], 1, mean)
    psi2<-apply(t3[,((ncol(t3)/2)+1):(ncol(t3))], 1, mean)
    dpsi<-psi1-psi2
    abs_dpsi<-abs(as.numeric(dpsi))
    abs_dpsi[which(abs_dpsi==1)]<-0
    logp<-(-log10(pp))
    # if_cor<-ifelse((abs(rr)>0.3)&(pp<0.05),1,0)
    if_cor<-ifelse((abs(rr)>cor_cutoff)&(pp<cor_p_cutoff),1,0)
    P<-apply(t2,1,function(x){
      x1<-as.numeric(x[1:(ncol(t2)/2)])
      x2<-as.numeric(x[((ncol(t2)/2)+1):(ncol(t2))])
      va1<-var(x1)
      va2<-var(x2)
      if (va1==0) x1[1]<-x1[1]+0.001
      if (va2==0) x2[1]<-x2[1]+0.001
      p<-t.test(x1,x2)
      return(p$p.value)
    })
    if_P<-(P<0.05)
    if_real<-ifelse((abs_dpsi>dpsi_network_threshold)&(abs_dpsi<1)&(if_P),1,0)
    use<-as.numeric(abs(rr)*if_cor*sqrt(abs_dpsi))
    corr_mat<-data.frame(rbp = rownames(expr)[i],
                         event = colnames(rr)[col(rr)],
                         cor = c(unlist(rr)),
                         p = c(unlist(pp)),
                         abs_dpsi = c(abs_dpsi),
                         logp = c(logp),
                         if_cor = c(if_cor),
                         if_real = c(if_real),
                         use = c(use))
    if (!is.null(BS)) {
      BS_part<-BS[which(BS[,1]==rownames(expr)[i]),]
      colnames(BS_part)<-c("rbp","event","BS_score")
      corr_mat<-dplyr::left_join(corr_mat,BS_part)
      # a<-stats::glm(if_real ~ scale(use) + BS_score ,family = binomial(link = "logit") ,data =corr_mat)
      a<-stats::glm(if_real ~ use + BS_score ,family = binomial(link = "logit") ,data =corr_mat)
      rbp_net_mat<-data.frame(rbp = rownames(expr)[i],
                              expr_p = expr_p$p.value,
                              BS = "T",
                              logit1 = a$coefficients[1],
                              logit2 = a$coefficients[2],
                              logit3 = a$coefficients[3])
    }else{
      # a<-stats::glm(if_real ~ scale(use) ,family = binomial(link = "logit") ,data =corr_mat)
      a<-stats::glm(if_real ~ use ,family = binomial(link = "logit") ,data =corr_mat)
      rbp_net_mat<-data.frame(rbp = rownames(expr)[i],
                              expr_p = expr_p$p.value,
                              BS = "F",
                              logit1 = a$coefficients[1],
                              logit2 = a$coefficients[2],
                              logit3 = 0)

    }
    return(list(rbp_net_mat=rbp_net_mat,corr_mat=corr_mat))
  }
  rbp_net_mat <- do.call(rbind, lapply(re, function(x){
    if (dim(x)[2]==6) return(x)
  }))
  corr_mat <- do.call(rbind, lapply(re, function(x){
    if (dim(x)[2]==9) return(x)
  }))
  parallel::stopCluster(cll)
  return(list(rbp_net_mat=rbp_net_mat,corr_mat=corr_mat))
}

#' MRAS_net_group_sc
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param num1,num2 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'                  Input requirements: if the input is less than 1, the input is considered as the sample proportion, and greater than 1, the input is considered as the number of samples.
#'                  When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param cor_cutoff The cutoff of co-regulation of RBP to event j, with a default value of 0
#' @param cor_p_cutoff The cutoff of p-value of co-regulation of RBP to event j, with a default value of 0.05.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param BS if input the binding data.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param path_use The path to the file used to store the output.
#' @param MRAS_net_single_re The output of function `MRAS_net_single_sc()`.
#' @param method The method to calculate correlation.
#' @param string_net The relationship between RBPs.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network.
#' @export
#'

MRAS_net_group_sc<-function(expr = expr,psi = psi,num1 = 0.1,num2 = 0.1,
                            cor_cutoff=0,cor_p_cutoff=0.05,
                            method="spearman",BS=NULL,
                            MRAS_net_single_re = MRAS_net_single_re,
                            dpsi_network_threshold = 0.1,
                            string_net = string_net,
                            threads = 2,path_use = path_use){
  rbp_net_mat<-MRAS_net_single_re$rbp_net_mat
  rbp_corr<-MRAS_net_single_re$corr_mat

  rbp_corr_deal<-rbp_corr[which((abs(rbp_corr$cor)>cor_cutoff)&(rbp_corr$p<cor_p_cutoff)),]
  rbp_corr_deal$cor<-abs(rbp_corr_deal$cor)
  rbp_corr_deal_tab<-table(rbp_corr_deal$event)
  rbp_corr_deal_remove<-names(rbp_corr_deal_tab)[which(rbp_corr_deal_tab == 1)]
  rbp_corr_deal<-rbp_corr_deal[!(rbp_corr_deal$event %in% rbp_corr_deal_remove),]

  rbp_corr_mat<-reshape2::dcast(rbp_corr_deal,rbp~event,value.var = "cor")
  rbp_corr_mat[is.na(rbp_corr_mat)]<-0
  rownames(rbp_corr_mat)<-rbp_corr_mat[,1]
  rbp_corr_mat<-rbp_corr_mat[,-1]
  rbp_corr_mat<-rbp_corr_mat[intersect(rownames(string_net),rownames(rbp_corr_mat)),]
  rbp_corr_mat<-as.matrix(rbp_corr_mat)
  expr_rbp<-as.matrix(expr[intersect(rownames(string_net),rownames(rbp_corr_mat)),])

  string_net_use<-string_net[intersect(rownames(string_net),rownames(rbp_corr_mat)),]
  string_net_use<-string_net_use[,intersect(rownames(string_net),rownames(rbp_corr_mat))]
  #+1
  expr_rbp_log<-as.matrix(apply(expr_rbp,2,function(x){
    return(log2(as.numeric(x)+1))
  }))
  rownames(expr_rbp_log)<-rownames(expr_rbp)
  if ((num1>0.5)&(num1<1)) stop("Wrong input : num!")
  if ((num2>0.5)&(num2<1)) stop("Wrong input : num!")
  use_mat1<-apply(expr, 1, function(q){
    len<-length(q[q!=0])
    num_use<-ceiling(len*num1)
    if (num_use < 5) stop("Please filter the RBP expr mat.")
    use<-c(order(q[q!=0],decreasing=TRUE)[1:num_use],
           order(q[q!=0],decreasing=FALSE)[1:num_use])
    return(which(q!=0)[use])
  })
  use_mat2<-apply(expr, 1, function(q){
    len<-length(q[q!=0])
    num_use<-ceiling(len*num2)
    if (num_use < 5) stop("Please filter the RBP expr mat.")
    use<-c(order(q[q!=0],decreasing=TRUE)[1:num_use],
           order(q[q!=0],decreasing=FALSE)[1:num_use])
    return(which(q!=0)[use])
  })

  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  filter_expr<-function(expr_matrix,n){
    results<-matrix(TRUE,nrow = nrow(expr_matrix),ncol = 1)
    i<-0
    for (i in 1:nrow(expr_matrix)) {
      length0<-length(which(as.numeric(expr_matrix[i,])<1))
      results[i,1]<-ifelse((length0)>n,FALSE,TRUE)
      # cat(round(i/nrow(expr_matrix),2),"\t")
    }
    return(results[,1])
  }
  MRAS_net_group_pre1<-function(rbp,psi_deal,use_mat1 = use_mat1,use_mat2 = use_mat2,method=method,
                                cor_cutoff=cor_cutoff,cor_p_cutoff=cor_p_cutoff,
                                num1 = 0.5,num2 = 0.5,dpsi_network_threshold = 0.1,
                                rbp_bind_event){

    t1<-rbp_bind_event[,use_mat1[[rbp]]]
    t2<-psi_deal[,use_mat1[[rbp]]]
    t1_2<-cbind(t1,t2)
    rr<-apply(t1_2, 1, function(x){
      return(cor(x[1:ncol(t1)],x[(ncol(t1)+1):ncol(t1_2)],method = method))
    })
    df<-ncol(t1)-2
    pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
    pp[which(pp>1)]<-2-pp[which(pp>1)]
    t3<-psi_deal[,use_mat2[[rbp]]]
    psi1<-apply(t3[,1:(ncol(t3)/2)], 1, mean)
    psi2<-apply(t3[,((ncol(t3)/2)+1):(ncol(t3))], 1, mean)
    dpsi<-psi1-psi2
    abs_dpsi<-abs(as.numeric(dpsi))
    abs_dpsi[which(abs_dpsi==1)]<-0
    logp<-(-log10(pp))
    if_cor<-ifelse((abs(rr)>cor_cutoff)&(pp<cor_p_cutoff),1,0)
    # if_cor<-ifelse((pp<0.05),1,0)
    P<-apply(t2,1,function(x){
      x1<-as.numeric(x[1:(ncol(t2)/2)])
      x2<-as.numeric(x[((ncol(t2)/2)+1):(ncol(t2))])
      va1<-var(x1)
      va2<-var(x2)
      if (va1==0) x1[1]<-x1[1]+0.001
      if (va2==0) x2[1]<-x2[1]+0.001
      p<-t.test(x1,x2)
      return(p$p.value)
    })
    if_P<-(P<0.05)
    if_real<-ifelse((abs_dpsi>dpsi_network_threshold)&(abs_dpsi<1)&(if_P),1,0)
    use<-as.numeric(abs(rr)*if_cor*sqrt(abs_dpsi))
    corr_mat<-data.frame(rbp = rbp,
                         event = rownames(rbp_bind_event),
                         cor_new = c(unlist(rr)),
                         p_new = c(unlist(pp)),
                         abs_dpsi_new = c(abs_dpsi),
                         logp_new = c(logp),
                         if_cor_new = c(if_cor),
                         if_real_new = c(if_real),
                         use_new = c(use))
    corr_mat<-as.matrix(corr_mat)
    return(corr_mat)
  }

  MRAS_net_group_pre2<-function(expr,re_new_all_part,rbp_corr,BS=NULL){
    re_new_all_part<-as.data.frame(re_new_all_part)
    rbp_corr<-as.data.frame(rbp_corr)
    network_re<-dplyr::left_join(rbp_corr,re_new_all_part,by = c("rbp", "event"))
    network_re$use_new[is.na(network_re$use_new)]<-network_re$use[is.na(network_re$use_new)]
    tj<-function(gene_ID_use,sep=","){
      if (is.na(gene_ID_use[1])){
        return(NA)
      }else{
        gene_ID_use<-gene_ID_use[complete.cases(gene_ID_use)]
        name<-paste(gene_ID_use,collapse = sep)
        return(name)
      }
    }
    network_re$use_new<-as.numeric(network_re$use_new)
    i=0
    re_final<-foreach (i = 1:nrow(expr), .combine = "rbind") %dopar% {
      network_re_part<-network_re[which(network_re[,1]==rownames(expr)[i]),]
      if (!is.null(BS)){
        BS_part<-BS[which(BS[,1]==rownames(expr)[i]),]
        colnames(BS_part)<-c("rbp","event","BS_score")
        network_re_part<-dplyr::left_join(network_re_part,BS_part)
        network_re_part[is.na(network_re_part)]<-0
        # a<-stats::glm(if_real ~ scale(use) + scale(use_new) + BS_score ,family = binomial(link = "logit") ,data =network_re_part)
        a<-stats::glm(if_real ~ use + use_new + BS_score ,family = binomial(link = "logit") ,data =network_re_part)
        rbp_net_mat_group<-data.frame(rbp = rownames(expr)[i],
                                      BS = "T",
                                      logit1_group = a$coefficients[1],
                                      logit2_group = a$coefficients[2],
                                      logit3_group = a$coefficients[3],
                                      logit4_group = a$coefficients[4])
      }else{
        # a<-stats::glm(if_real ~ scale(use) + scale(use_new) ,family = binomial(link = "logit") ,data =network_re_part)
        a<-stats::glm(if_real ~ use + use_new ,family = binomial(link = "logit") ,data =network_re_part)
        rbp_net_mat_group<-data.frame(rbp = rownames(expr)[i],
                                      BS = "F",
                                      logit1_group = a$coefficients[1],
                                      logit2_group = a$coefficients[2],
                                      logit3_group = a$coefficients[3],
                                      logit4_group = 0)
      }
      return(rbp_net_mat_group)
    }
    return(list(network_re = network_re,re_final = re_final))
  }
  rrr=0
  re_new_all<-foreach (rrr = 1:ncol(string_net_use),.combine = "rbind",.errorhandling = "pass") %dopar% {
    # for(rr in 1:ncol(string_net)){
    rbp<-colnames(string_net_use)[rrr]
    string_net_part<-string_net_use[,rrr]
    rbp_corr_mat_part<-rbp_corr_mat*string_net_part
    string_net_value<-unlist(apply(rbp_corr_mat_part, 2, function(x){
      return(prod(x[which(x!=0)]))
    }))
    rbp_corr_mat_part[which(rbp_corr_mat_part!=0)]<-1
    rbp_corr_mat_part_2<-as.data.frame(apply(rbp_corr_mat_part,2,as.numeric))
    rownames(rbp_corr_mat_part_2)<-rownames(rbp_corr_mat_part)
    rbp_bind_event<-t(rbp_corr_mat_part_2) %*% as.matrix(expr_rbp_log)
    rbp_bind_event<-as.matrix(rbp_bind_event)
    new_matrix <- 2^rbp_bind_event
    new_matrix<-new_matrix -1
    rbp_bind_event<-new_matrix*string_net_value
    ff1<-filter_expr(rbp_bind_event,ncol(rbp_bind_event)*0.8)
    rbp_bind_event<-rbp_bind_event[ff1,]
    psi_deal<-psi[rownames(rbp_bind_event),]
    if (nrow(psi_deal)>0){
      re_new<-MRAS_net_group_pre1(rbp = rbp,psi_deal = psi_deal,use_mat1 = use_mat1,use_mat2 = use_mat2,method=method,
                                  cor_cutoff=cor_cutoff,cor_p_cutoff=cor_p_cutoff,
                                  num1 = num1,num2 = num2,dpsi_network_threshold = dpsi_network_threshold,
                                  rbp_bind_event = rbp_bind_event)
      re_new<-as.matrix(re_new)
      return(re_new)
    }else{
      err<-matrix(0,1,9)
      colnames(err)<-c("rbp","event","cor_new","p_new","abs_dpsi_new","logp_new","if_cor_new","if_real_new","use_new")
      return(err)
      # re_new<-data.frame(rbp = 0,
      #                    event = 0,
      #                    cor_new = 0,
      #                    p_new = 0,
      #                    abs_dpsi_new = 0,
      #                    logp_new = 0,
      #                    if_cor_new = 0,
      #                    if_real_new = 0,
      #                    use_new = 0)
      # re_new<-as.matrix(re_new)
      # return(re_new)
    }
    # return(re_new)
  }

  re_new_all<-as.matrix(re_new_all)
  re_new_all<-re_new_all[!(is.na(re_new_all[,9])),]
  re_new_all<-re_new_all[which(as.numeric(re_new_all[,9]) != 0),]
  network_re_final<-MRAS_net_group_pre2(expr = expr,re_new_all_part = re_new_all[,c(1,2,9)],
                                        rbp_corr = rbp_corr,BS = BS)
  parallel::stopCluster(cll)
  rbp_net_mat_final<-dplyr::left_join(rbp_net_mat,network_re_final$re_final)
  return(list(rbp_net_mat_group = rbp_net_mat_final,rbp_corr_group = network_re_final$network_re))
}
#' MRAS_net_single_sc_work
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param num1,num2 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'                  Input requirements: if the input is less than 1, the input is considered as the sample proportion, and greater than 1, the input is considered as the number of samples.
#'                  When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param cor_cutoff The cutoff of co-regulation of RBP to event j, with a default value of 0.3.
#' @param cor_p_cutoff The cutoff of p-value of co-regulation of RBP to event j, with a default value of 0.05.
#' @param BS if input the binding data.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param path_use The path to the file used to store the output.
#' @param method The method to calculate correlation.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network.
#' @export
#'
MRAS_net_single_sc_work<-function(expr,psi,num1 = 0.5,num2 = 0.5,method="spearman",
                                  cor_cutoff=0,cor_p_cutoff=0.05,
                                  dpsi_network_threshold = 0.1,BS = NULL,
                                  threads = 2,path_use = path_use){
  if ((num1>0.5)&(num1<1)) stop("Wrong input : num!")
  if ((num2>0.5)&(num2<1)) stop("Wrong input : num!")
  use_mat1<-apply(expr, 1, function(q){
    len<-length(q[q!=0])
    num_use<-ceiling(len*num1)
    if (num_use < 5) stop("Please filter the RBP expr mat.")
    use<-c(order(q[q!=0],decreasing=TRUE)[1:num_use],
           order(q[q!=0],decreasing=FALSE)[1:num_use])
    return(which(q!=0)[use])
  })
  use_mat2<-apply(expr, 1, function(q){
    len<-length(q[q!=0])
    num_use<-ceiling(len*num2)
    if (num_use < 5) stop("Please filter the RBP expr mat.")
    use<-c(order(q[q!=0],decreasing=TRUE)[1:num_use],
           order(q[q!=0],decreasing=FALSE)[1:num_use])
    return(which(q!=0)[use])
  })
  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  tj<-function(gene_ID_use,sep=","){
    if (is.na(gene_ID_use[1])){
      return(NA)
    }else{
      gene_ID_use<-gene_ID_use[complete.cases(gene_ID_use)]
      name<-paste(gene_ID_use,collapse = sep)
      return(name)
    }
  }
  i=0
  re<-foreach (i = 1:nrow(expr), .combine = "rbind") %dopar% {
    t1<-matrix(unlist(expr[i,use_mat1[[i]]]),nrow = 1)
    expr_p<-t.test(t1[1:(ncol(t1)/2)],t1[((ncol(t1)/2)+1):(ncol(t1))])
    t2<-psi[,use_mat1[[i]]]
    rr<-cor(t(t1),t(t2),method = method)
    df<-ncol(t1)-2
    rownames(rr)<-rownames(expr)[i]
    pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
    pp[which(pp>1)]<-2-pp[which(pp>1)]
    t3<-psi[,use_mat2[[i]]]
    psi1<-apply(t3[,1:(ncol(t3)/2)], 1, mean)
    psi2<-apply(t3[,((ncol(t3)/2)+1):(ncol(t3))], 1, mean)
    dpsi<-psi1-psi2
    abs_dpsi<-abs(as.numeric(dpsi))
    abs_dpsi[which(abs_dpsi==1)]<-0
    logp<-(-log10(pp))
    if_cor<-ifelse((abs(rr)>cor_cutoff)&(pp<cor_p_cutoff),1,0)
    # if_cor<-ifelse((pp<0.05),1,0)
    P<-apply(t2,1,function(x){
      x1<-as.numeric(x[1:(ncol(t2)/2)])
      x2<-as.numeric(x[((ncol(t2)/2)+1):(ncol(t2))])
      va1<-var(x1)
      va2<-var(x2)
      if (va1==0) x1[1]<-x1[1]+0.001
      if (va2==0) x2[1]<-x2[1]+0.001
      p<-t.test(x1,x2)
      return(p$p.value)
    })
    if_P<-(P<0.05)
    if_real<-ifelse((abs_dpsi>dpsi_network_threshold)&(abs_dpsi<1)&(if_P),1,0)
    use<-as.numeric(abs(rr)*if_cor*sqrt(abs_dpsi))
    corr_mat<-data.frame(rbp = rownames(expr)[i],
                         event = colnames(rr)[col(rr)],
                         cor = c(unlist(rr)),
                         p = c(unlist(pp)),
                         abs_dpsi = c(abs_dpsi),
                         logp = c(logp),
                         if_cor = c(if_cor),
                         if_real = c(if_real),
                         use = c(use))
    return(corr_mat)
  }
  parallel::stopCluster(cll)
  return(re)
}
#' MRAS_net_group_sc_work
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param num1,num2 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'                  Input requirements: if the input is less than 1, the input is considered as the sample proportion, and greater than 1, the input is considered as the number of samples.
#'                  When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param cor_cutoff The cutoff of co-regulation of RBP to event j, with a default value of 0.3.
#' @param cor_p_cutoff The cutoff of p-value of co-regulation of RBP to event j, with a default value of 0.05.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param BS if input the binding data.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param path_use The path to the file used to store the output.
#' @param method The method to calculate correlation.
#' @param rbp_corr_work The output of function `MRAS_net_single_sc_work()`.
#' @param string_net The relationship between RBPs.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network.
#' @export
#'
MRAS_net_group_sc_work<-function(expr = expr,psi = psi,num1 = 0.5,num2 = 0.5,
                                 cor_cutoff=0,cor_p_cutoff=0.05,
                                 method="spearman",BS=NULL,
                                 rbp_corr_work = rbp_corr_work,
                                 dpsi_network_threshold = 0.1,
                                 string_net = string_net,
                                 threads = 2,path_use = path_use){
  rbp_corr<-rbp_corr_work

  rbp_corr_deal<-rbp_corr[which((abs(rbp_corr$cor)>cor_cutoff)&(rbp_corr$p<cor_p_cutoff)),]
  rbp_corr_deal$cor<-abs(rbp_corr_deal$cor)
  rbp_corr_deal_tab<-table(rbp_corr_deal$col)
  rbp_corr_deal_remove<-names(rbp_corr_deal_tab)[which(rbp_corr_deal_tab == 1)]
  rbp_corr_deal<-rbp_corr_deal[!(rbp_corr_deal$event %in% rbp_corr_deal_remove),]

  rbp_corr_mat<-reshape2::dcast(rbp_corr_deal,rbp~event,value.var = "cor")
  rbp_corr_mat[is.na(rbp_corr_mat)]<-0
  rownames(rbp_corr_mat)<-rbp_corr_mat[,1]
  rbp_corr_mat<-rbp_corr_mat[,-1]
  rbp_corr_mat<-rbp_corr_mat[intersect(rownames(string_net),rownames(rbp_corr_mat)),]
  rbp_corr_mat<-as.matrix(rbp_corr_mat)
  expr_rbp<-as.matrix(expr[intersect(rownames(string_net),rownames(rbp_corr_mat)),])

  string_net_use<-string_net[intersect(rownames(string_net),rownames(rbp_corr_mat)),]
  string_net_use<-string_net_use[,intersect(rownames(string_net),rownames(rbp_corr_mat))]
  #+1
  expr_rbp_log<-as.matrix(apply(expr_rbp,2,function(x){
    return(log2(as.numeric(x)+1))
  }))
  rownames(expr_rbp_log)<-rownames(expr_rbp)

  if ((num1>0.5)&(num1<1)) stop("Wrong input : num!")
  if ((num2>0.5)&(num2<1)) stop("Wrong input : num!")
  use_mat1<-apply(expr, 1, function(q){
    len<-length(q[q!=0])
    num_use<-ceiling(len*num1)
    if (num_use < 5) stop("Please filter the RBP expr mat.")
    use<-c(order(q[q!=0],decreasing=TRUE)[1:num_use],
           order(q[q!=0],decreasing=FALSE)[1:num_use])
    return(which(q!=0)[use])
  })
  use_mat2<-apply(expr, 1, function(q){
    len<-length(q[q!=0])
    num_use<-ceiling(len*num2)
    if (num_use < 5) stop("Please filter the RBP expr mat.")
    use<-c(order(q[q!=0],decreasing=TRUE)[1:num_use],
           order(q[q!=0],decreasing=FALSE)[1:num_use])
    return(which(q!=0)[use])
  })

  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  filter_expr<-function(expr_matrix,n){
    results<-matrix(TRUE,nrow = nrow(expr_matrix),ncol = 1)
    i<-0
    for (i in 1:nrow(expr_matrix)) {
      length0<-length(which(as.numeric(expr_matrix[i,])<1))
      results[i,1]<-ifelse((length0)>n,FALSE,TRUE)
      # cat(round(i/nrow(expr_matrix),2),"\t")
    }
    return(results[,1])
  }
  MRAS_net_group_pre1<-function(rbp,psi_deal,use_mat1 = use_mat1,use_mat2 = use_mat2,method=method,
                                cor_cutoff=cor_cutoff,cor_p_cutoff=cor_p_cutoff,
                                num1 = 0.5,num2 = 0.5,dpsi_network_threshold = 0.1,
                                rbp_bind_event){

    t1<-rbp_bind_event[,use_mat1[[rbp]]]
    t2<-psi_deal[,use_mat1[[rbp]]]
    t1_2<-cbind(t1,t2)
    rr<-apply(t1_2, 1, function(x){
      return(cor(x[1:ncol(t1)],x[(ncol(t1)+1):ncol(t1_2)],method = method))
    })
    df<-ncol(t1)-2
    pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
    pp[which(pp>1)]<-2-pp[which(pp>1)]
    t3<-psi_deal[,use_mat2[[rbp]]]
    psi1<-apply(t3[,1:(ncol(t3)/2)], 1, mean)
    psi2<-apply(t3[,((ncol(t3)/2)+1):(ncol(t3))], 1, mean)
    dpsi<-psi1-psi2
    abs_dpsi<-abs(as.numeric(dpsi))
    abs_dpsi[which(abs_dpsi==1)]<-0
    logp<-(-log10(pp))
    if_cor<-ifelse((abs(rr)>cor_cutoff)&(pp<cor_p_cutoff),1,0)
    # if_cor<-ifelse((pp<0.05),1,0)
    P<-apply(t2,1,function(x){
      x1<-as.numeric(x[1:(ncol(t2)/2)])
      x2<-as.numeric(x[((ncol(t2)/2)+1):(ncol(t2))])
      va1<-var(x1)
      va2<-var(x2)
      if (va1==0) x1[1]<-x1[1]+0.001
      if (va2==0) x2[1]<-x2[1]+0.001
      p<-t.test(x1,x2)
      return(p$p.value)
    })
    if_P<-(P<0.05)
    if_real<-ifelse((abs_dpsi>dpsi_network_threshold)&(abs_dpsi<1)&(if_P),1,0)
    use<-as.numeric(abs(rr)*if_cor*sqrt(abs_dpsi))
    corr_mat<-data.frame(rbp = rbp,
                         event = rownames(rbp_bind_event),
                         cor_new = c(unlist(rr)),
                         p_new = c(unlist(pp)),
                         abs_dpsi_new = c(abs_dpsi),
                         logp_new = c(logp),
                         if_cor_new = c(if_cor),
                         if_real_new = c(if_real),
                         use_new = c(use))
    return(corr_mat)
  }
  rrr=0
  re_new_all<-foreach (rrr = 1:ncol(string_net_use),.combine = "rbind",.errorhandling = "pass") %dopar% {
    # for(rr in 1:ncol(string_net)){
    rbp<-colnames(string_net_use)[rrr]
    cat(rbp,"\t")
    string_net_part<-string_net_use[,rrr]
    rbp_corr_mat_part<-rbp_corr_mat*string_net_part
    string_net_value<-unlist(apply(rbp_corr_mat_part, 2, function(x){
      return(prod(x[which(x!=0)]))
    }))
    rbp_corr_mat_part[which(rbp_corr_mat_part!=0)]<-1
    rbp_bind_event<-t(rbp_corr_mat_part) %*% expr_rbp_log
    rbp_bind_event<-as.matrix(rbp_bind_event)
    new_matrix <- 2^rbp_bind_event
    new_matrix<-new_matrix -1
    rbp_bind_event<-new_matrix*string_net_value
    ff1<-filter_expr(rbp_bind_event,ncol(rbp_bind_event)*0.8)
    rbp_bind_event<-rbp_bind_event[ff1,]
    psi_deal<-psi[rownames(rbp_bind_event),]
    if (nrow(psi_deal)>0){
      re_new<-MRAS_net_group_pre1(rbp = rbp,psi_deal = psi_deal,use_mat1 = use_mat1,use_mat2 = use_mat2,method=method,
                                  cor_cutoff=cor_cutoff,cor_p_cutoff=cor_p_cutoff,
                                  num1 = num1,num2 = num2,dpsi_network_threshold = dpsi_network_threshold,
                                  rbp_bind_event = rbp_bind_event)
      re_new<-as.matrix(re_new)
      return(re_new)
    }else{
      err<-matrix(0,1,9)
      colnames(err)<-c("rbp","event","cor_new","p_new","abs_dpsi_new","logp_new","if_cor_new","if_real_new","use_new")
      return(err)
    }
    return(re_new)
  }
  parallel::stopCluster(cll)
  re_new_all<-as.matrix(re_new_all)
  re_new_all<-re_new_all[which(as.numeric(re_new_all[,9]) != 0),]
  re_new_all_part<-as.data.frame(re_new_all[,c(1,2,9)])
  rbp_corr<-as.data.frame(rbp_corr)
  network_re<-dplyr::left_join(rbp_corr,re_new_all_part,by = c("rbp", "event"))
  network_re$use_new[is.na(network_re$use_new)]<-network_re$use[is.na(network_re$use_new)]
  return(network_re)
}

#' MRAS_net_single_fc
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param num1 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'             Input requirements: if the input is less than 1, the input is considered as the sample proportion.
#'             When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param cor_cutoff The cutoff of co-regulation of RBP to event j, with a default value of 0.3.
#' @param cor_p_cutoff The cutoff of p-value of co-regulation of RBP to event j, with a default value of 0.05.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param BS if input the binding data.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param path_use The path to the file used to store the output.
#' @param fc_mat fold change matrix.
#' @param dpsi_mat dpsi matrix.
#' @param method The method to calculate correlation.
#' @param real real dpsi.
#' @param sample batch matrix.

#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network.
#'
MRAS_net_single_fc<-function(expr,psi,fc_mat,dpsi_mat,num1 = 0.1,method="spearman",real=rbp_event_0.1,
                             cor_cutoff=0.3,cor_p_cutoff=0.05,
                             sample=batch_matrix,dpsi_network_threshold = 0.1,BS = NULL,
                             threads = 2,path_use = path_use){
  if ((num1>0.5)&(num1<1)) stop("Wrong input : num!")
  if (num1<1) num1<-ceiling(ncol(fc_mat)*num1)
  use_mat1<-t(apply(fc_mat, 1, function(q){
    use<-c(order(q,decreasing=TRUE)[1:num1],
           order(q,decreasing=FALSE)[1:num1])
    return(use)
  }))
  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  tj<-function(gene_ID_use,sep=","){
    if (is.na(gene_ID_use[1])){
      return(NA)
    }else{
      gene_ID_use<-gene_ID_use[complete.cases(gene_ID_use)]
      name<-paste(gene_ID_use,collapse = sep)
      return(name)
    }
  }
  i=0
  re<-foreach (i = 1:nrow(fc_mat), .combine = "rbind") %dopar% {
    t1<-matrix(as.numeric(unlist(fc_mat[i,use_mat1[i,]])),nrow = 1)
    expr_p<-t.test(t1[1:(ncol(t1)/2)],t1[((ncol(t1)/2)+1):(ncol(t1))])
    t2<-dpsi_mat[,use_mat1[i,]]
    rr<-cor(t(t1),t(t2),method = method)
    df<-ncol(t1)-2
    rownames(rr)<-rownames(expr)[i]
    pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
    pp[which(pp>1)]<-2-pp[which(pp>1)]

    # name<-colnames(fc_mat)[use_mat1[i,]]
    t3<-dpsi_mat[,use_mat1[i,]]
    psi1<-apply(t3[,1:(ncol(t3)/2)], 1, mean)
    psi2<-apply(t3[,((ncol(t3)/2)+1):(ncol(t3))], 1, mean)
    dpsi<-psi1-psi2
    abs_dpsi<-abs(as.numeric(dpsi))
    abs_dpsi[which(abs_dpsi==1)]<-0
    logp<-(-log10(pp))
    if_cor<-ifelse((abs(rr)>cor_cutoff)&(pp<cor_p_cutoff),1,0)
    P<-apply(t2,1,function(x){
      x1<-as.numeric(x[1:(ncol(t2)/2)])
      x2<-as.numeric(x[((ncol(t2)/2)+1):(ncol(t2))])
      va1<-var(x1)
      va2<-var(x2)
      if (va1==0) x1[1]<-x1[1]+0.001
      if (va2==0) x2[1]<-x2[1]+0.001
      p<-t.test(x1,x2)
      return(p$p.value)
    })
    if_P<-(P<0.05)
    real_part<-real[which(real[,1] == rownames(fc_mat)[i]),]
    rownames(real_part)<-real_part[,2]
    real_part<-real_part[rownames(psi),]

    # real_part<-real_part[which(as.numeric(real_part[,3])>0.1),]
    if_real<-ifelse(as.numeric(real_part[,3])>0.1,1,0)

    # if_real<-ifelse(rownames(psi) %in% real_part[,2],1,0)
    # if_real<-ifelse((abs_dpsi>dpsi_network_threshold)&(abs_dpsi<1)&(if_P),1,0)
    use<-as.numeric(abs(rr)*if_cor*sqrt(abs_dpsi))
    corr_mat<-data.frame(rbp = rownames(expr)[i],
                         event = colnames(rr)[col(rr)],
                         cor = c(unlist(rr)),
                         p = c(unlist(pp)),
                         abs_dpsi = c(abs_dpsi),
                         logp = c(logp),
                         if_cor = c(if_cor),
                         if_real = c(if_real),
                         use = c(use))
    if (!is.null(BS)) {
      BS_part<-BS[which(BS[,1]==rownames(expr)[i]),]
      colnames(BS_part)<-c("rbp","event","BS_score")
      corr_mat<-dplyr::left_join(corr_mat,BS_part)
      # a<-stats::glm(if_real ~ scale(use) + BS_score ,family = binomial(link = "logit") ,data =corr_mat)
      a<-stats::glm(if_real ~ use + BS_score ,family = binomial(link = "logit") ,data =corr_mat)
      rbp_net_mat<-data.frame(rbp = rownames(expr)[i],
                              expr_p = expr_p$p.value,
                              BS = "T",
                              logit1 = a$coefficients[1],
                              logit2 = a$coefficients[2],
                              logit3 = a$coefficients[3])
    }else{
      # a<-stats::glm(if_real ~ scale(use) ,family = binomial(link = "logit") ,data =corr_mat)
      a<-stats::glm(if_real ~ use ,family = binomial(link = "logit") ,data =corr_mat)
      rbp_net_mat<-data.frame(rbp = rownames(expr)[i],
                              expr_p = expr_p$p.value,
                              BS = "F",
                              logit1 = a$coefficients[1],
                              logit2 = a$coefficients[2],
                              logit3 = 0)
    }
    return(list(rbp_net_mat=rbp_net_mat,corr_mat=corr_mat))
  }
  rbp_net_mat <- do.call(rbind, lapply(re, function(x){
    if (dim(x)[2]==6) return(x)
  }))
  corr_mat <- do.call(rbind, lapply(re, function(x){
    if (dim(x)[2]==9) return(x)
  }))
  parallel::stopCluster(cll)
  return(list(rbp_net_mat=rbp_net_mat,corr_mat=corr_mat))
}


#' MRAS_net_single_fc2
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param num1 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'             Input requirements: if the input is less than 1, the input is considered as the sample proportion.
#'             When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param cor_cutoff The cutoff of co-regulation of RBP to event j, with a default value of 0.3.
#' @param cor_p_cutoff The cutoff of p-value of co-regulation of RBP to event j, with a default value of 0.05.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param BS if input the binding data.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param path_use The path to the file used to store the output.
#' @param fc_mat fold change matrix.
#' @param dpsi_mat dpsi matrix.
#' @param method The method to calculate correlation.
#' @param real real dpsi.
#' @param sample batch matrix.

#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network.
#'
MRAS_net_single_fc2<-function(expr,psi,fc_mat,dpsi_mat,num1 = 0.1,method=c("pearson","spearman"),real=rbp_event_0.1,
                              cor_cutoff = 0.3,cor_p_cutoff = 0.05,
                              sample=batch_matrix,dpsi_network_threshold = 0.1,BS = NULL,
                              threads = 2,path_use = path_use){
  if ((num1>0.5)&(num1<1)) stop("Wrong input : num!")
  if (num1<1) num1<-ceiling(ncol(fc_mat)*num1)
  use_mat1<-t(apply(fc_mat, 1, function(q){
    use<-c(order(q,decreasing=TRUE)[1:num1],
           order(q,decreasing=FALSE)[1:num1])
    return(use)
  }))
  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  tj<-function(gene_ID_use,sep=","){
    if (is.na(gene_ID_use[1])){
      return(NA)
    }else{
      gene_ID_use<-gene_ID_use[complete.cases(gene_ID_use)]
      name<-paste(gene_ID_use,collapse = sep)
      return(name)
    }
  }
  i=0
  re<-foreach (i = 1:nrow(fc_mat), .combine = "rbind",.errorhandling = "stop") %dopar% {
    t1<-matrix(as.numeric(unlist(fc_mat[i,use_mat1[i,]])),nrow = 1)
    expr_p<-t.test(t1[1:(ncol(t1)/2)],t1[((ncol(t1)/2)+1):(ncol(t1))])
    t2<-dpsi_mat[,use_mat1[i,]]
    rr<-cor(t(t1),t(t2),method = method)
    df<-ncol(t1)-2
    rownames(rr)<-rownames(expr)[i]
    pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
    pp[which(pp>1)]<-2-pp[which(pp>1)]

    # name<-colnames(fc_mat)[use_mat1[i,]]
    # name<-rownames(fc_mat)[i]
    # rbp_sample<-sample$sample_name[which(sample$File_target %in% name)]
    # ctrl_sample<-sample$sample_name[which((sample$type=="ctrl")&(sample$LEVEL %in% unique(sample$LEVEL[which(sample$File_target %in% name)])))]
    #
    # t3<-dpsi_mat[,use_mat1[i,]]
    # psi1<-apply(t3[,1:(ncol(t3)/2)], 1, mean)
    # psi2<-apply(t3[,((ncol(t3)/2)+1):(ncol(t3))], 1, mean)
    # dpsi<-psi1-psi2
    # abs_dpsi<-abs(as.numeric(dpsi))
    # abs_dpsi[which(abs_dpsi==1)]<-0
    logp<-(-log10(pp))
    if_cor<-ifelse((abs(rr)>cor_cutoff)&(pp<cor_p_cutoff),1,0)
    P<-apply(t2,1,function(x){
      x1<-as.numeric(x[1:(ncol(t2)/2)])
      x2<-as.numeric(x[((ncol(t2)/2)+1):(ncol(t2))])
      va1<-var(x1)
      va2<-var(x2)
      if (va1==0) x1[1]<-x1[1]+0.001
      if (va2==0) x2[1]<-x2[1]+0.001
      p<-t.test(x1,x2)
      return(p$p.value)
    })
    if_P<-(P<0.05)
    real_part<-real[which(real[,1] == rownames(fc_mat)[i]),]
    rownames(real_part)<-real_part[,2]
    real_part<-real_part[rownames(psi),]
    abs_dpsi<-abs(as.numeric(real_part[,3]))
    abs_dpsi[which(abs_dpsi==1)]<-0
    # real_part<-real_part[which(as.numeric(real_part[,3])>0.1),]

    if_real<-ifelse(as.numeric(real_part[,3])>0.1,1,0)
    # if_real<-ifelse((abs_dpsi>dpsi_network_threshold)&(abs_dpsi<1)&(if_P),1,0)
    use<-as.numeric(abs(rr)*if_cor*sqrt(abs_dpsi))
    corr_mat<-data.frame(rbp = rownames(expr)[i],
                         event = colnames(rr)[col(rr)],
                         cor = c(unlist(rr)),
                         p = c(unlist(pp)),
                         abs_dpsi = c(abs_dpsi),
                         logp = c(logp),
                         if_cor = c(if_cor),
                         if_real = c(if_real),
                         use = c(use))
    if (!is.null(BS)) {
      BS_part<-BS[which(BS[,1]==rownames(expr)[i]),]
      colnames(BS_part)<-c("rbp","event","BS_score")
      rownames(BS_part)<-BS_part[,2]
      BS_part<-BS_part[rownames(psi),]
      corr_mat<-cbind(corr_mat,BS_score=BS_part[,3])
      # a<-stats::glm(if_real ~ scale(use) + BS_score ,family = binomial(link = "logit") ,data =corr_mat)
      a<-stats::glm(if_real ~ use + BS_score ,family = binomial(link = "logit") ,data =corr_mat)
      rbp_net_mat<-data.frame(rbp = rownames(expr)[i],
                              expr_p = expr_p$p.value,
                              BS = "T",
                              logit1 = a$coefficients[1],
                              logit2 = a$coefficients[2],
                              logit3 = a$coefficients[3])
    }else{
      # a<-stats::glm(if_real ~ scale(use) ,family = binomial(link = "logit") ,data =corr_mat)
      a<-stats::glm(if_real ~ use ,family = binomial(link = "logit") ,data =corr_mat)
      rbp_net_mat<-data.frame(rbp = rownames(expr)[i],
                              expr_p = expr_p$p.value,
                              BS = "F",
                              logit1 = a$coefficients[1],
                              logit2 = a$coefficients[2],
                              logit3 = 0)
    }
    return(list(rbp_net_mat=rbp_net_mat,corr_mat=corr_mat))
  }
  rbp_net_mat <- do.call(rbind, lapply(re, function(x){
    if (dim(x)[2]==6) return(x)
  }))
  bs_num<-ifelse(is.null(BS),9,10)
  corr_mat <- do.call(rbind, lapply(re, function(x){
    if (dim(x)[2]==bs_num) return(x)
  }))
  parallel::stopCluster(cll)
  return(list(rbp_net_mat=rbp_net_mat,corr_mat=corr_mat))
}
#' MRAS_net_single_fc_work
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param num1 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'             Input requirements: if the input is less than 1, the input is considered as the sample proportion.
#'             When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param cor_cutoff The cutoff of co-regulation of RBP to event j, with a default value of 0.3.
#' @param cor_p_cutoff The cutoff of p-value of co-regulation of RBP to event j, with a default value of 0.05.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param BS if input the binding data.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param path_use The path to the file used to store the output.
#' @param fc_mat fold change matrix.
#' @param dpsi_mat dpsi matrix.
#' @param method The method to calculate correlation.
#' @param sample batch matrix.

#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network.
#'
MRAS_net_single_fc_work<-function(expr,psi,fc_mat,dpsi_mat,num1 = 0.1,method=c("pearson","spearman"),
                                  cor_cutoff = 0.3,cor_p_cutoff = 0.05,
                                  sample=batch_matrix,dpsi_network_threshold = 0.1,BS = NULL,
                                  threads = 2,path_use = path_use){
  if ((num1>0.5)&(num1<1)) stop("Wrong input : num!")
  if (num1<1) num1<-ceiling(ncol(fc_mat)*num1)
  use_mat1<-t(apply(fc_mat, 1, function(q){
    use<-c(order(q,decreasing=TRUE)[1:num1],
           order(q,decreasing=FALSE)[1:num1])
    return(use)
  }))

  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  tj<-function(gene_ID_use,sep=","){
    if (is.na(gene_ID_use[1])){
      return(NA)
    }else{
      gene_ID_use<-gene_ID_use[complete.cases(gene_ID_use)]
      name<-paste(gene_ID_use,collapse = sep)
      return(name)
    }
  }
  i=0
  re<-foreach (i = 1:nrow(expr), .combine = "rbind") %dopar% {
    t1<-matrix(unlist(fc_mat[i,use_mat1[i,]]),nrow = 1)
    expr_p<-t.test(t1[1:(ncol(t1)/2)],t1[((ncol(t1)/2)+1):(ncol(t1))])
    t2<-psi[,use_mat1[i,]]
    rr<-cor(t(t1),t(t2),method = method)
    df<-ncol(t1)-2
    rownames(rr)<-rownames(expr)[i]
    pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
    pp[which(pp>1)]<-2-pp[which(pp>1)]

    name<-colnames(fc_mat)[use_mat1[i,]]
    rbp_sample<-sample$sample_name[which(sample$File_target %in% name)]
    ctrl_sample<-sample$sample_name[which((sample$type=="ctrl")&(sample$LEVEL %in% unique(sample$LEVEL[which(sample$File_target %in% name)])))]

    t3<-psi[,c(rbp_sample,ctrl_sample)]
    psi1<-apply(t3[,1:length(rbp_sample)], 1, mean)
    psi2<-apply(t3[,(length(rbp_sample)+1):(length(rbp_sample)+length(ctrl_sample))], 1, mean)
    dpsi<-psi1-psi2
    abs_dpsi<-abs(as.numeric(dpsi))
    abs_dpsi[which(abs_dpsi==1)]<-0
    logp<-(-log10(pp))
    if_cor<-ifelse((abs(rr)>cor_cutoff)&(pp<cor_p_cutoff),1,0)
    P<-apply(t2,1,function(x){
      x1<-as.numeric(x[1:(ncol(t2)/2)])
      x2<-as.numeric(x[((ncol(t2)/2)+1):(ncol(t2))])
      va1<-var(x1)
      va2<-var(x2)
      if (va1==0) x1[1]<-x1[1]+0.001
      if (va2==0) x2[1]<-x2[1]+0.001
      p<-t.test(x1,x2)
      return(p$p.value)
    })
    if_P<-(P<0.05)
    if_real<-ifelse((abs_dpsi>dpsi_network_threshold)&(abs_dpsi<1)&(if_P),1,0)
    use<-as.numeric(abs(rr)*if_cor*sqrt(abs_dpsi))
    corr_mat<-data.frame(rbp = rownames(expr)[i],
                         event = colnames(rr)[col(rr)],
                         cor = c(unlist(rr)),
                         p = c(unlist(pp)),
                         abs_dpsi = c(abs_dpsi),
                         logp = c(logp),
                         if_cor = c(if_cor),
                         if_real = c(if_real),
                         use = c(use))
    return(corr_mat)
  }
  parallel::stopCluster(cll)
  return(re)
}

#' get_MRAS_net
#'
#' @param rbp_net_mat_group first net which can get from "MRAS".
#' @param rbp_corr_group_work second net.
#' @param Regulate_threshold P cutoff.
#' @param BS RBP binding matrix.
#' @param threads threads.
#' @param path_use The path to the file used to store the output.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network
#' @export
#'

get_MRAS_net<-function(rbp_net_mat_group,rbp_corr_group_work,
                       Regulate_threshold = 0.5,BS = NULL,
                       threads = 2,path_use){
  rbp_net_mat_group[is.na(rbp_net_mat_group)]<-0
  rbp_corr_group_work$use_new<-as.numeric(rbp_corr_group_work$use_new)
  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  i<-0
  network_work<-foreach (i = 1:nrow(rbp_net_mat_group),.combine = "rbind") %dopar% {
    rbp<-rbp_net_mat_group$rbp[i]
    network_work_part<-rbp_corr_group_work[which(rbp_corr_group_work$rbp == rbp),]
    logit1<-rbp_net_mat_group$logit1[i]
    logit2<-rbp_net_mat_group$logit2[i]
    logit1_group<-rbp_net_mat_group$logit1_group[i]
    logit2_group<-rbp_net_mat_group$logit2_group[i]
    logit3_group<-rbp_net_mat_group$logit3_group[i]
    if (is.null(BS)){
      # network_work_part$score_single<-exp(logit1+logit2*scale(as.numeric(network_work_part$use)))/(1+exp(logit1+logit2*scale(as.numeric(network_work_part$use))))
      # network_work_part$score_group<-exp(logit1_group+logit2_group*scale(as.numeric(network_work_part$use))+logit3_group*scale(as.numeric(network_work_part$use_new)))/(1+exp(logit1_group+logit2_group*scale(as.numeric(network_work_part$use))+logit3_group*scale(as.numeric(network_work_part$use_new))))
      network_work_part$score_single<-exp(logit1+logit2*as.numeric(network_work_part$use))/(1+exp(logit1+logit2*as.numeric(network_work_part$use)))
      network_work_part$score_group<-exp(logit1_group+logit2_group*as.numeric(network_work_part$use)+logit3_group*as.numeric(network_work_part$use_new))/(1+exp(logit1_group+logit2_group*as.numeric(network_work_part$use)+logit3_group*as.numeric(network_work_part$use_new)))

    }else{
      logit3<-rbp_net_mat_group$logit3[i]
      logit4_group<-rbp_net_mat_group$logit4_group[i]
      # network_work_part$score_single<-exp(logit1+logit2*scale(as.numeric(network_work_part$use))+logit3*network_work_part$BS_score)/(1+exp(logit1+logit2*scale(as.numeric(network_work_part$use))+logit3*network_work_part$BS_score))
      # network_work_part$score_group<-exp(logit1_group+logit2_group*scale(as.numeric(network_work_part$use))+logit3_group*scale(as.numeric(network_work_part$use_new))+logit4*network_work_part$BS_score)/(1+exp(logit1_group+logit2_group*scale(as.numeric(network_work_part$use))+logit3_group*scale(as.numeric(network_work_part$use_new))+logit4*network_work_part$BS_score))
      #
      network_work_part$score_single<-exp(logit1+logit2*as.numeric(network_work_part$use)+logit3*network_work_part$BS_score)/(1+exp(logit1+logit2*as.numeric(network_work_part$use)+logit3*network_work_part$BS_score))
      network_work_part$score_group<-exp(logit1_group+logit2_group*as.numeric(network_work_part$use)+logit3_group*as.numeric(network_work_part$use_new)+logit4_group*network_work_part$BS_score)/(1+exp(logit1_group+logit2_group*as.numeric(network_work_part$use)+logit3_group*as.numeric(network_work_part$use_new)+logit4_group*network_work_part$BS_score))
    }
    return(network_work_part)
  }
  parallel::stopCluster(cll)
  rbp_event_single<-cbind(network_work$event,network_work$score_single,network_work$rbp,network_work$abs_dpsi)
  colnames(rbp_event_single)<-c("ID","score","RBP","dpsi")
  rbp_event_single<-as.matrix(rbp_event_single)
  if (!dir.exists(paste0(path_use,"/deal"))) dir.create(paste0(path_use,"/deal/"))
  if (!dir.exists(paste0(path_use,"/deal/single"))) dir.create(paste0(path_use,"/deal/single/"))
  data.table::fwrite(as.data.frame(rbp_event_single),file = paste0(path_use,"deal/single/rbp_event_deal_all_total.txt"),
                     append = F,col.names = F,row.names = F,sep = "\t",quote = F)
  rbp_event_single_deal<-rbp_event_single[which(as.numeric(rbp_event_single[,2])>Regulate_threshold),1:3]
  rbp_event_single_deal<-matrix(rbp_event_single_deal,ncol = 3)
  data.table::fwrite(as.data.frame(rbp_event_single_deal),file = paste0(path_use,"deal/single/rbp_event_deal_all.txt"),
                     append = F,col.names = F,row.names = F,sep = "\t",quote = F)
  rbp_event_group<-cbind(network_work$event,network_work$score_group,network_work$rbp,network_work$abs_dpsi)
  colnames(rbp_event_group)<-c("ID","score","RBP","dpsi")
  rbp_event_group<-as.matrix(rbp_event_group)
  if (!dir.exists(paste0(path_use,"/deal/group"))) dir.create(paste0(path_use,"/deal/group/"))
  data.table::fwrite(as.data.frame(rbp_event_group),file = paste0(path_use,"deal/group/rbp_event_deal_all_total.txt"),
                     append = F,col.names = F,row.names = F,sep = "\t",quote = F)
  rbp_event_group_deal<-rbp_event_group[which(as.numeric(rbp_event_group[,2])>Regulate_threshold),1:3]
  rbp_event_group_deal<-matrix(rbp_event_group_deal,ncol = 3)
  data.table::fwrite(as.data.frame(rbp_event_group_deal),file = paste0(path_use,"deal/group/rbp_event_deal_all.txt"),
                     append = F,col.names = F,row.names = F,sep = "\t",quote = F)
  return(network_work)
}

#' get_MRAS_net_smooth
#'
#' @param MRAS_net_group_re first net.
#' @param Regulate_threshold P cutoff.
#' @param BS RBP binding matrix.
#' @param threads threads.
#' @param path_use The path to the file used to store the output.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network
#' @export
#'

get_MRAS_net_smooth<-function(MRAS_net_group_re,
                       Regulate_threshold = 0.5,BS = NULL,
                       threads = 2,path_use){
  rbp_net_mat_group<-MRAS_net_group_re$rbp_net_mat_group
  rbp_net_mat_group[is.na(rbp_net_mat_group)]<-0
  rbp_corr_group_work<-MRAS_net_group_re$rbp_corr_group
  rbp_corr_group_work$use_new<-as.numeric(rbp_corr_group_work$use_new)
  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  i<-0
  network_work<-foreach (i = 1:nrow(rbp_net_mat_group),.combine = "rbind") %dopar% {
    rbp<-rbp_net_mat_group$rbp[i]
    network_work_part<-rbp_corr_group_work[which(rbp_corr_group_work$rbp == rbp),]
    logit1<-rbp_net_mat_group$logit1[i]
    logit2<-rbp_net_mat_group$logit2[i]
    logit1_group<-rbp_net_mat_group$logit1_group[i]
    logit2_group<-rbp_net_mat_group$logit2_group[i]
    logit3_group<-rbp_net_mat_group$logit3_group[i]
    if (is.null(BS)){
      # network_work_part$score_single<-exp(logit1+logit2*scale(as.numeric(network_work_part$use)))/(1+exp(logit1+logit2*scale(as.numeric(network_work_part$use))))
      # network_work_part$score_group<-exp(logit1_group+logit2_group*scale(as.numeric(network_work_part$use))+logit3_group*scale(as.numeric(network_work_part$use_new)))/(1+exp(logit1_group+logit2_group*scale(as.numeric(network_work_part$use))+logit3_group*scale(as.numeric(network_work_part$use_new))))
      network_work_part$score_single<-exp(logit1+logit2*as.numeric(network_work_part$use))/(1+exp(logit1+logit2*as.numeric(network_work_part$use)))
      network_work_part$score_group<-exp(logit1_group+logit2_group*as.numeric(network_work_part$use)+logit3_group*as.numeric(network_work_part$use_new))/(1+exp(logit1_group+logit2_group*as.numeric(network_work_part$use)+logit3_group*as.numeric(network_work_part$use_new)))

    }else{
      logit3<-rbp_net_mat_group$logit3[i]
      logit4_group<-rbp_net_mat_group$logit4_group[i]
      # network_work_part$score_single<-exp(logit1+logit2*scale(as.numeric(network_work_part$use))+logit3*network_work_part$BS_score)/(1+exp(logit1+logit2*scale(as.numeric(network_work_part$use))+logit3*network_work_part$BS_score))
      # network_work_part$score_group<-exp(logit1_group+logit2_group*scale(as.numeric(network_work_part$use))+logit3_group*scale(as.numeric(network_work_part$use_new))+logit4*network_work_part$BS_score)/(1+exp(logit1_group+logit2_group*scale(as.numeric(network_work_part$use))+logit3_group*scale(as.numeric(network_work_part$use_new))+logit4*network_work_part$BS_score))
      #
      network_work_part$score_single<-exp(logit1+logit2*as.numeric(network_work_part$use)+logit3*network_work_part$BS_score)/(1+exp(logit1+logit2*as.numeric(network_work_part$use)+logit3*network_work_part$BS_score))
      network_work_part$score_group<-exp(logit1_group+logit2_group*as.numeric(network_work_part$use)+logit3_group*as.numeric(network_work_part$use_new)+logit4_group*network_work_part$BS_score)/(1+exp(logit1_group+logit2_group*as.numeric(network_work_part$use)+logit3_group*as.numeric(network_work_part$use_new)+logit4_group*network_work_part$BS_score))
    }
    return(network_work_part)
  }
  parallel::stopCluster(cll)
  rbp_event_single<-cbind(network_work$event,network_work$score_single,network_work$rbp,network_work$abs_dpsi)
  colnames(rbp_event_single)<-c("ID","score","RBP","dpsi")
  rbp_event_single<-as.matrix(rbp_event_single)
  if (!dir.exists(paste0(path_use,"/deal"))) dir.create(paste0(path_use,"/deal/"))
  if (!dir.exists(paste0(path_use,"/deal/single"))) dir.create(paste0(path_use,"/deal/single/"))
  data.table::fwrite(as.data.frame(rbp_event_single),file = paste0(path_use,"deal/single/rbp_event_deal_all_total.txt"),
                     append = F,col.names = F,row.names = F,sep = "\t",quote = F)
  rbp_event_single_deal<-rbp_event_single[which(as.numeric(rbp_event_single[,2])>Regulate_threshold),1:3]
  rbp_event_single_deal<-matrix(rbp_event_single_deal,ncol = 3)
  data.table::fwrite(as.data.frame(rbp_event_single_deal),file = paste0(path_use,"deal/single/rbp_event_deal_all.txt"),
                     append = F,col.names = F,row.names = F,sep = "\t",quote = F)
  rbp_event_group<-cbind(network_work$event,network_work$score_group,network_work$rbp,network_work$abs_dpsi)
  colnames(rbp_event_group)<-c("ID","score","RBP","dpsi")
  rbp_event_group<-as.matrix(rbp_event_group)
  if (!dir.exists(paste0(path_use,"/deal/group"))) dir.create(paste0(path_use,"/deal/group/"))
  data.table::fwrite(as.data.frame(rbp_event_group),file = paste0(path_use,"deal/group/rbp_event_deal_all_total.txt"),
                     append = F,col.names = F,row.names = F,sep = "\t",quote = F)
  rbp_event_group_deal<-rbp_event_group[which(as.numeric(rbp_event_group[,2])>Regulate_threshold),1:3]
  rbp_event_group_deal<-matrix(rbp_event_group_deal,ncol = 3)
  data.table::fwrite(as.data.frame(rbp_event_group_deal),file = paste0(path_use,"deal/group/rbp_event_deal_all.txt"),
                     append = F,col.names = F,row.names = F,sep = "\t",quote = F)
  return(network_work)
}

#' get_MRAS_net_fc
#'
#' @param rbp_net_mat first net.
#' @param rbp_corr_work second net.
#' @param Regulate_threshold P cutoff.
#' @param if_BS if use BS.
#' @param threads threads.
#' @param path_use The path to the file used to store the output.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @return network
#' @export
#'

get_MRAS_net_fc<-function(rbp_net_mat,rbp_corr_work,
                          Regulate_threshold = 0.5,if_BS = FALSE,
                          threads = 2,path_use){
  # rbp_net_mat<-MRAS_net_single_re$rbp_net_mat
  rbp_net_mat[is.na(rbp_net_mat)]<-0
  rbp_corr_work$use<-as.numeric(rbp_corr_work$use)
  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  i<-0
  network_work<-foreach (i = 1:nrow(rbp_net_mat),.combine = "rbind") %dopar% {
    rbp<-rbp_net_mat$rbp[i]
    network_work_part<-rbp_corr_work[which(rbp_corr_work$rbp == rbp),]
    logit1<-rbp_net_mat$logit1[i]
    logit2<-rbp_net_mat$logit2[i]
    if (!if_BS){
      # network_work_part$score_single<-exp(logit1+logit2*scale(as.numeric(network_work_part$use)))/(1+exp(logit1+logit2*scale(as.numeric(network_work_part$use))))
      # network_work_part$score_group<-exp(logit1_group+logit2_group*scale(as.numeric(network_work_part$use))+logit3_group*scale(as.numeric(network_work_part$use_new)))/(1+exp(logit1_group+logit2_group*scale(as.numeric(network_work_part$use))+logit3_group*scale(as.numeric(network_work_part$use_new))))
      network_work_part$score_single<-exp(logit1+logit2*as.numeric(network_work_part$use))/(1+exp(logit1+logit2*as.numeric(network_work_part$use)))
    }else{
      logit3<-rbp_net_mat$logit3[i]
      # network_work_part$score_single<-exp(logit1+logit2*scale(as.numeric(network_work_part$use))+logit3*network_work_part$BS_score)/(1+exp(logit1+logit2*scale(as.numeric(network_work_part$use))+logit3*network_work_part$BS_score))
      # network_work_part$score_group<-exp(logit1_group+logit2_group*scale(as.numeric(network_work_part$use))+logit3_group*scale(as.numeric(network_work_part$use_new))+logit4*network_work_part$BS_score)/(1+exp(logit1_group+logit2_group*scale(as.numeric(network_work_part$use))+logit3_group*scale(as.numeric(network_work_part$use_new))+logit4*network_work_part$BS_score))
      #
      network_work_part$score_single<-exp(logit1+logit2*as.numeric(network_work_part$use)+logit3*network_work_part$BS_score)/(1+exp(logit1+logit2*as.numeric(network_work_part$use)+logit3*network_work_part$BS_score))
    }
    return(network_work_part)
  }
  parallel::stopCluster(cll)
  rbp_event_single<-cbind(network_work$event,network_work$score_single,network_work$rbp,network_work$abs_dpsi)
  colnames(rbp_event_single)<-c("ID","score","RBP","dpsi")
  rbp_event_single<-as.matrix(rbp_event_single)
  if (!dir.exists(paste0(path_use,"/deal"))) dir.create(paste0(path_use,"/deal/"))
  if (!dir.exists(paste0(path_use,"/deal/single"))) dir.create(paste0(path_use,"/deal/single/"))
  data.table::fwrite(as.data.frame(rbp_event_single),file = paste0(path_use,"deal/single/rbp_event_deal_all_total.txt"),
                     append = F,col.names = F,row.names = F,sep = "\t",quote = F)
  rbp_event_single_deal<-rbp_event_single[which(as.numeric(rbp_event_single[,2])>Regulate_threshold),1:3]
  rbp_event_single_deal<-matrix(rbp_event_single_deal,ncol = 3)
  data.table::fwrite(as.data.frame(rbp_event_single_deal),file = paste0(path_use,"deal/single/rbp_event_deal_all.txt"),
                     append = F,col.names = F,row.names = F,sep = "\t",quote = F)
  return(network_work)
}

