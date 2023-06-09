#' This function constructs the RBP-EVents regulatory relationship.
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param num2 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'                  Input requirements: if the input is less than 1, the input is considered as the sample proportion, and greater than 1, the input is considered as the number of samples.
#'                  When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param rbp_corr The correlation matrix of RBP and events has 4 columns,
#'                 the first column is the row name of input x,
#'                 the second column is the row name of input y,
#'                 the third column is the magnitude of the correlation,
#'                 the fourth column is the significance level P-value.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 2.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param Regulate_threshold A threshold value to measure whether the RBP regulates a splicing event, greater than which the RBP is considered to regulate the splicing event, with a default threshold value of 0.5.
#' @param path_use The path to the file used to store the output.
#'
#' @return No return value, the result is stored under the path provided.
#' @export
#' @import doParallel
#' @import parallel
#'

MRAS_network<-function(expr,psi,num2 = 0.5,rbp_corr = rbp_corr,
                       dpsi_network_threshold = 0.1,
                       Regulate_threshold = 0.5,threads = 2,path_use = path_use ){
  expres_f2<-as.matrix(expr)
  # if (num1<1) num1<-ceiling(ncol(expr)*num1)
  if (num2<1) num2<-ceiling(ncol(expr)*num2)

  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  rbp_event_0.1<-foreach::foreach(i = rownames(expres_f2), .packages = "dplyr", .errorhandling = "pass",.combine = "rbind") %dopar% {
    sample_use<-colnames(expres_f2)[c(order(expres_f2[i,],decreasing=TRUE)[1:num2],
                                      order(expres_f2[i,],decreasing=FALSE)[1:num2])]
    rbp_expr<-expres_f2[,sample_use]
    rbp_psi<-psi[,sample_use]
    rbp_psi<-as.data.frame(rbp_psi)
    dpsi<-unlist(apply(rbp_psi, 1, function(x){
      d<-x[1:(ncol(rbp_psi)/2)]-x[((ncol(rbp_psi)/2)+1):(ncol(rbp_psi))]
      return(mean(d))
    }))
    rbp_psi$dd<-dpsi
    matt<-cbind(i,rownames(rbp_psi),rbp_psi$dd)
    colnames(matt)<-c("RBP","Event","label")
    return(matt)
  }
  parallel::stopCluster(cll)
  if (!dir.exists(paste0(path_use,"/deal"))) dir.create(paste0(path_use,"/deal/"))
  if (file.exists(paste0(path_use,"/deal/bb.txt"))) file.remove(paste0(path_use,"/deal/bb.txt"))
  if (file.exists(paste0(path_use,"/deal/rbp_event_deal_all.txt"))) file.remove(paste0(path_use,"/deal/rbp_event_deal_all.txt"))
  if (file.exists(paste0(path_use,"/deal/rbp_event_deal_all.RData")))  file.remove(paste0(path_use,"/deal/rbp_event_deal_all.RData"))
  if (file.exists(paste0(path_use,"/deal/rbp_event_deal_all_total.txt")))  file.remove(paste0(path_use,"/deal/rbp_event_deal_all_total.txt"))
  if (file.exists(paste0(path_use,"/deal/rbp_event_deal_all_total.RData")))  file.remove(paste0(path_use,"/deal/rbp_event_deal_all_total.RData"))

  for (i in rownames(expres_f2)){
    sample_use<-colnames(expres_f2)[c(order(expres_f2[i,],decreasing=TRUE)[1:num2],
                                      order(expres_f2[i,],decreasing=FALSE)[1:num2])]
    rbp_expr<-expres_f2[,sample_use]
    rbp_psi<-psi[,sample_use]
    rbp_psi<-as.data.frame(rbp_psi)
    rbp_event_real<-rbp_event_0.1[which(rbp_event_0.1[,1]==i),]
    rbp_FC_part<-rbp_corr[which(rbp_corr[,1]==i),]
    rbp_event_real<-rbp_event_real[which(rbp_event_real[,2] %in% unique(rbp_FC_part[,2])),]
    rbp_use<-merge(rbp_event_real,rbp_FC_part,by.x=c("RBP","Event"),by.y=c("row","col"))
    rbp_use$label<-abs(as.numeric(rbp_use$label))
    rbp_use$label<-ifelse(rbp_use$label==1,0,rbp_use$label)
    rbp_use$cor<-abs(rbp_use$cor)
    rbp_use<-rbp_use[which(rbp_use$p>0),]
    rbp_use$logp<-(-log10(rbp_use$p))
    rbp_use$if_cor<-ifelse((rbp_use$cor>0.3)&(rbp_use$p<0.05),1,0)
    rbp_use$if_real<-ifelse((rbp_use$label>dpsi_network_threshold)&(rbp_use$label<1),1,0)

    rbp_use$use<-rbp_use$cor*rbp_use$if_cor*sqrt(rbp_use$label)
    rbp_use$use<-as.numeric(rbp_use$use)
    a<-stats::glm(if_real ~ use ,family = binomial(link = "logit") ,data =rbp_use)
    bb<-tj(c(i,a$coefficients[1]))
    write.table(as.data.frame(bb),file = paste0(path_use,"deal/bb.txt"),append = T,col.names = F,row.names = F,sep = "\t",quote = F)
    # rbp_use$score<-exp(a$coefficients[1]+a$coefficients[2]*rbp_use$label+a$coefficients[3]*rbp_use$value+a$coefficients[4]*rbp_use$logp)/(1+exp(a$coefficients[1]+a$coefficients[2]*rbp_use$label+a$coefficients[3]*rbp_use$value+a$coefficients[4]*rbp_use$logp))
    rbp_use$score<-exp(a$coefficients[1]+a$coefficients[2]*rbp_use$use)/(1+exp(a$coefficients[1]+a$coefficients[2]*rbp_use$use))

    rbp_event<-cbind(rbp_use$Event,rbp_use$score,i,rbp_use$label)
    colnames(rbp_event)<-c("ID","score","RBP","dpsi")
    rbp_event<-as.matrix(rbp_event)
    data.table::fwrite(as.data.frame(rbp_event),file = paste0(path_use,"deal/rbp_event_deal_all_total.txt"),
           append = T,         col.names = F,row.names = F,sep = "\t",quote = F)
    rbp_event_deal<-rbp_event[which(as.numeric(rbp_event[,2])>Regulate_threshold),1:3]
    rbp_event_deal<-matrix(rbp_event_deal,ncol = 3)
    data.table::fwrite(as.data.frame(rbp_event_deal),file = paste0(path_use,"deal/rbp_event_deal_all.txt"),
           append = T,  col.names = F,row.names = F,sep = "\t",quote = F)
  }
}

