#' Master Regulator analysis of Alternative Splicing
#'
#' @param input_type input type, different options for user's different needs.
#'                   "1": input differential splicing events and use pre-constructed regulatory network;
#'                   "2": input data and inferred relationships by MRAS.
#'                   "3": input data and construct a network using the user's own data
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param rbp_interested The name of the RBP interested.
#' @param m Number of the first condition.
#' @param n Number of the second condition.
#' @param DS_pvalue Significance level of differential splicing analysis.
#' @param DS_dPSI Threshold for differential splicing.
#' @param rbp_event_deal_all_total,rbp_event_deal_all the network built by MRAS.
#' @param method Method of finding correlation.
#' @param num1,num2 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'                  Input requirements: if the input is less than 1, the input is considered as the sample proportion, and greater than 1, the input is considered as the number of samples.
#'                  When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param BS RBP binding matrix.
#' @param smooth network building type.User can choose "TRUE" or "FALSE" in the third input type.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param Regulate_threshold A threshold value to measure whether the RBP regulates a splicing event, greater than which the RBP is considered to regulate the splicing event, with a default threshold value of 0.5.
#' @param rbp_net_mat_group first net which can get from "MRAS".
#' @param path_use The path to the file used to store the output.
#' @param group TRUE or FALSE.If consider the synergy between RBPs.
#' @param sc TRUE or FALSE. If single cell data.
#' @param result_type The types of return values are "top10", "tab_simple", and "tab_all".
#'                    "top10" returns the top 10 MRAS results;
#'                    "tab_simple" returns a simplified version of the MRAS scoring matrix (without the events corresponding to the RBP);
#'                    "tab_all" returns the full version of the MRAS scoring matrix. matrix (including the events corresponding to the RBP).
#'
#' @return MRAS result. If you have any questions,you can get help here (https://github.com/zhou-lei5/MRAS).
#' @export
#'
#'

MRAS<-function(input_type,
               expr,psi,rbp_interested = NULL,m = 0,n = 0,DS_pvalue = 0.05,DS_dPSI = 0.1,
               rbp_event_deal_all_total = NULL,rbp_event_deal_all = NULL,
               method = "spearman",BS = NULL,group = T,
               dpsi_network_threshold = 0.1,Regulate_threshold = 0.5,rbp_net_mat_group = NULL,
               num1 = 0.5,num2 = 0.5,sc = F,
               threads = 1,path_use ,smooth = FALSE,
               result_type = c("Top10","tab_simple","tab_all")){
  if (input_type == "1"){
    cat("Step1:Performing differential splicing analysis...\n")
    Events_DS<-DS_matrix(psi,m = m,n = n)
    Event_DS_sig<-get_Event_DS_sig(Events_DS,DS_pvalue = DS_pvalue,DS_dPSI = DS_dPSI)
    data.table::fwrite(as.data.frame(Event_DS_sig),file = paste0(path_use,"DS_mat.txt"),
                       row.names = F,col.names = T,quote = F,sep = "\t")
    cat("Step2:Preparing data...\n")
    RBP_use<-get_RBP_use(expr,m = m,n = n)

    cat("Go directly to step4!\n")
    cat("Step4:Performing enrichment analysis...\n")
    result<-score_matrix(rbp_interested = rbp_interested,
                         Events_DS = Events_DS,Event_DS_sig = Event_DS_sig,
                         RBP_use = RBP_use,result_type = result_type,threads = threads,
                         rbp_event_deal_all_total = rbp_event_deal_all_total,
                         rbp_event_deal_all = rbp_event_deal_all,path_use = path_use)
    cat("Finish!\n")
    return(result)
  }
  if (input_type == "2"){
    cat("Step1:Performing differential splicing analysis...\n")
    Events_DS<-DS_matrix(psi,m = m,n = n)
    Event_DS_sig<-get_Event_DS_sig(Events_DS,DS_pvalue = DS_pvalue,DS_dPSI = DS_dPSI)
    data.table::fwrite(as.data.frame(Event_DS_sig),file = paste0(path_use,"DS_mat.txt"),
                       row.names = F,col.names = T,quote = F,sep = "\t")

    cat("Step2:Preparing data...\n")
    RBP_use<-get_RBP_use(expr,m = m,n = n)


    cat("Step3:Constructing RBP-Event regulatory relationship network...\n")
    rbp_corr_work<-MRAS_net_single_work(expr = expr,psi = psi,num1 = 0.5,num2 = 0.5,
                                        method=method,BS = BS,
                                        dpsi_network_threshold = dpsi_network_threshold,
                                        threads = threads,path_use = path_use)

    rbp_corr_group_work<-MRAS_net_group_work(expr = expr,psi = psi,num1 = 0.5,num2 = 0.5,
                                             method=method,BS=BS,
                                             rbp_corr_work = rbp_corr_work,
                                             dpsi_network_threshold = dpsi_network_threshold,
                                             string_net = string_net,
                                             threads = threads,path_use = path_use)
    network_work<-get_MRAS_net(rbp_net_mat_group = rbp_net_mat_group,
                               rbp_corr_group_work = rbp_corr_group_work,
                               Regulate_threshold = Regulate_threshold,BS = BS,
                               threads = threads,path_use = path_use)
    if (group){
      rbp_event_deal_all_total<-get_rbp_event_deal_all_total(path_useful = paste0(path_use,"deal/group"))
      rbp_event_deal_all<-get_rbp_event_deal_all(path_useful = paste0(path_use,"deal/group"))
    }else{
      rbp_event_deal_all_total<-get_rbp_event_deal_all_total(path_useful = paste0(path_use,"deal/single"))
      rbp_event_deal_all<-get_rbp_event_deal_all(path_useful = paste0(path_use,"deal/single"))
    }
    cat("Step4:Performing enrichment analysis...\n")
    result<-score_matrix(rbp_interested = rbp_interested,
                         Events_DS = Events_DS,Event_DS_sig = Event_DS_sig,
                         RBP_use = RBP_use,result_type = result_type,threads = threads,
                         rbp_event_deal_all_total = rbp_event_deal_all_total,
                         rbp_event_deal_all = rbp_event_deal_all,path_use = path_use)
    cat("Finish!\n")
    return(result)
  }
  if (input_type == "3"){
    cat("Step1:Performing differential splicing analysis...\n")
    Events_DS<-DS_matrix(psi,m = m,n = n)
    Event_DS_sig<-get_Event_DS_sig(Events_DS,DS_pvalue = DS_pvalue,DS_dPSI = DS_dPSI)
    data.table::fwrite(as.data.frame(Event_DS_sig),file = paste0(path_use,"DS_mat.txt"),
                       row.names = F,col.names = T,quote = F,sep = "\t")

    cat("Step2:Preparing data...\n")
    RBP_use<-get_RBP_use(expr,m = m,n = n)


    cat("Step3:Constructing RBP-Event regulatory relationship network...\n")
    if (!sc){
      MRAS_net_single_re<-MRAS_net_single(expr = expr,psi = psi,num1 = num1,num2 = num2,
                                          method=method,BS = BS,
                                          dpsi_network_threshold = dpsi_network_threshold,
                                          threads = threads,path_use = path_use)

      MRAS_net_group_re<-MRAS_net_group(expr = expr,psi = psi,num1 = num1,num2 = num2,
                                        method=method,BS=BS,
                                        MRAS_net_single_re = MRAS_net_single_re,
                                        dpsi_network_threshold = dpsi_network_threshold,
                                        string_net = string_net,
                                        threads = threads,path_use = path_use)
      if (!smooth){
        rbp_corr_work<-MRAS_net_single_work(expr = expr,psi = psi,num1 = 0.5,num2 = 0.5,
                                            method=method,BS = BS,
                                            dpsi_network_threshold = dpsi_network_threshold,
                                            threads = threads,path_use = path_use)

        rbp_corr_group_work<-MRAS_net_group_work(expr = expr,psi = psi,num1 = 0.5,num2 = 0.5,
                                                 method=method,BS=BS,
                                                 rbp_corr_work = rbp_corr_work,
                                                 dpsi_network_threshold = dpsi_network_threshold,
                                                 string_net = string_net,
                                                 threads = threads,path_use = path_use)
        network_work<-get_MRAS_net(rbp_net_mat_group = MRAS_net_group_re$rbp_net_mat_group,
                                   rbp_corr_group_work = rbp_corr_group_work,
                                   Regulate_threshold = Regulate_threshold,BS = BS,
                                   threads = threads,path_use = path_use)
      }else{
        network_work<-get_MRAS_net_smooth(MRAS_net_group_re,
                                          Regulate_threshold = 0.5,BS = NULL,
                                          threads = 6,path_use = path_use)
      }
    }else{
      MRAS_net_single_re<-MRAS_net_single_sc(expr = expr,psi = psi,num1 = num1,num2 = num2,
                                             dpsi_network_threshold = dpsi_network_threshold,BS = BS,
                                             threads = threads,path_use = path_use)

      MRAS_net_group_re<-MRAS_net_group_sc(expr = expr,psi = psi,num1 = num1,num2 = num2,BS=BS,
                                           MRAS_net_single_re = MRAS_net_single_re,
                                           dpsi_network_threshold = dpsi_network_threshold,
                                           string_net = string_net,
                                           threads = threads,path_use = path_use)
      if (!smooth){
        rbp_corr_work<-MRAS_net_single_sc_work(expr = expr,psi = psi,num1 = 0.5,num2 = 0.5,
                                               dpsi_network_threshold = dpsi_network_threshold,BS = BS,
                                               threads = threads,path_use = path_use)

        rbp_corr_group_work<-MRAS_net_group_sc_work(expr = expr,psi = psi,num1 = 0.5,num2 = 0.5,BS=BS,
                                                    rbp_corr_work = rbp_corr_work,
                                                    dpsi_network_threshold = dpsi_network_threshold,
                                                    string_net = string_net,
                                                    threads = threads,path_use = path_use)
        network_work<-get_MRAS_net(rbp_net_mat_group = MRAS_net_group_re$rbp_net_mat_group,
                                   rbp_corr_group_work = rbp_corr_group_work,
                                   Regulate_threshold = Regulate_threshold,BS = BS,
                                   threads = threads,path_use = path_use)
      }else{
        network_work<-get_MRAS_net_smooth(MRAS_net_group_re,
                                          Regulate_threshold = 0.5,BS = NULL,
                                          threads = 6,path_use = path_use)
      }
    }
    if (group){
      rbp_event_deal_all_total<-get_rbp_event_deal_all_total(path_useful = paste0(path_use,"deal/group"))
      rbp_event_deal_all<-get_rbp_event_deal_all(path_useful = paste0(path_use,"deal/group"))
    }else{
      rbp_event_deal_all_total<-get_rbp_event_deal_all_total(path_useful = paste0(path_use,"deal/single"))
      rbp_event_deal_all<-get_rbp_event_deal_all(path_useful = paste0(path_use,"deal/single"))
    }
    cat("Step4:Performing enrichment analysis...\n")
    result<-score_matrix(rbp_interested = rbp_interested,
                         Events_DS = Events_DS,Event_DS_sig = Event_DS_sig,
                         RBP_use = RBP_use,result_type = result_type,threads = threads,
                         rbp_event_deal_all_total = rbp_event_deal_all_total,
                         rbp_event_deal_all = rbp_event_deal_all,path_use = path_use)
    cat("Finish!\n")
    return(result)
  }
}
