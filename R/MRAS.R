#' Master Regulator analysis of Alternative Splicing
#'
#' @param expr RBP expression matrix.
#' @param psi Events splicing matrix.
#' @param rbp_interested The name of the RBP interested.
#' @param m Number of the first condition.
#' @param n Number of the second condition.
#' @param DS_pvalue Significance level of differential splicing analysis.
#' @param DS_dPSI Threshold for differential splicing.
#' @param method Method of finding correlation.
#' @param num1,num2 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'                  Input requirements: if the input is less than 1, the input is considered as the sample proportion, and greater than 1, the input is considered as the number of samples.
#'                  When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 1.
#' @param dpsi_network_threshold A threshold value that measures the absolute value of the degree of differential spliicng events in two different condition, above which variable clipping is considered to have occurred, with a default threshold of 0.1
#' @param Regulate_threshold A threshold value to measure whether the RBP regulates a splicing event, greater than which the RBP is considered to regulate the splicing event, with a default threshold value of 0.5.
#' @param path_use The path to the file used to store the output.
#' @param result_type The types of return values are "top10", "tab_simple", and "tab_all".
#'                    "top10" returns the top 10 MRAS results;
#'                    "tab_simple" returns a simplified version of the MRAS scoring matrix (without the events corresponding to the RBP);
#'                    "tab_all" returns the full version of the MRAS scoring matrix. matrix (including the events corresponding to the RBP).
#'
#' @return MRAS result.
#' @export
#'
#'

MRAS<-function(expr,psi,rbp_interested = NULL,m = 0,n = 0,DS_pvalue = 0.05,DS_dPSI = 0.1,
               method = c("pearson","spearman"),num1 = 0.5,num2 = 0.5,
               dpsi_network_threshold = 0.1,Regulate_threshold = 0.5,threads = 1,path_use,
               result_type = c("Top10","tab_simple","tab_all")){
  cat("Step1:Performing differential splicing analysis...\n")
  Events_DS<-DS_matrix(psi,m = m,n = n)
  Event_DS_sig<-get_Event_DS_sig(Events_DS,DS_pvalue = DS_pvalue,DS_dPSI = DS_dPSI)

  cat("Step2:Constructing RBP-Event regulatory relationship network...\n")
  rbp_corr<-cor_spearman_percent(expr,psi,method = method,num1 = num1,threads = threads)

  MRAS_network(expr,psi,num1 = num1,num2 = num2,rbp_corr = rbp_corr,
               dpsi_network_threshold = dpsi_network_threshold,
               Regulate_threshold = Regulate_threshold,threads = threads,path_use = path_use)
  rbp_event_deal_all<-get_rbp_event_deal_all(path_use = path_use)
  rbp_event_deal_all_total<-get_rbp_event_deal_all_total(path_use = path_use)

  cat("Step3:Preparing for enrichment analysis...\n")
  RBP_use<-get_RBP_use(expr,m = m,n = n)

  cat("Step4:Performing enrichment analysis...\n")
  result<-score_matrix(rbp_interested = rbp_interested,
                       Events_DS = Events_DS,Event_DS_sig = Event_DS_sig,
                       RBP_use = RBP_use,result_type = result_type,threads = threads,
                       rbp_event_deal_all_total = rbp_event_deal_all_total,
                       rbp_event_deal_all = rbp_event_deal_all,path_use = path_use)
  return(result)
}
