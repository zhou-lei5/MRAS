# run_mras.R
library(optparse)

# Define the command-line options
option_list <- list(
  make_option(c("-i", "--input_type"), type = "character", default = NULL, help = "Input type"),
  make_option(c("-e", "--expr"), type = "character", default = NULL, help = "Path to the expression data"),
  make_option(c("-p", "--psi"), type = "character", default = NULL, help = "Path to the PSI data"),
  make_option(c("--rbp_interested"), type = "character", default = NULL, help = "RBP of interest (optional)"),
  make_option(c("-m", "--m"), type = "numeric", default = 0, help = "Parameter m"),
  make_option(c("-n", "--n"), type = "numeric", default = 0, help = "Parameter n"),
  make_option(c("--RBP_cutoff"), type = "numeric", default = 0.05, help = "RBP cutoff"),
  make_option(c("--DS_pvalue"), type = "numeric", default = 0.05, help = "DS p-value"),
  make_option(c("--DS_dPSI"), type = "numeric", default = 0.1, help = "DS delta PSI"),
  make_option(c("--rbp_event_deal_all_total"), type = "character", default = NULL, help = "Total RBP event deal"),
  make_option(c("--rbp_event_deal_all"), type = "character", default = NULL, help = "RBP event deal"),
  make_option(c("-M", "--method"), type = "character", default = "spearman", help = "Method for analysis (default is 'spearman')"),
  make_option(c("--BS"), type = "character", default = NULL, help = "Path to the Binding Strength data"),
  make_option(c("--group"), type = "logical", default = TRUE, help = "Group parameter (default TRUE)"),
  make_option(c("--dpsi_network_threshold"), type = "numeric", default = 0.1, help = "Delta PSI network threshold"),
  make_option(c("--Regulate_threshold"), type = "numeric", default = 0.5, help = "Regulate threshold"),
  make_option(c("--rbp_net_mat_group"), type = "character", default = NULL, help = "RBP network matrix group"),
  make_option(c("--num1"), type = "numeric", default = 0.5, help = "Parameter num1"),
  make_option(c("--num2"), type = "numeric", default = 0.5, help = "Parameter num2"),
  make_option(c("--sc"), type = "logical", default = FALSE, help = "Single-cell parameter (default FALSE)"),
  make_option(c("-t", "--threads"), type = "numeric", default = 1, help = "Number of threads"),
  make_option(c("-o", "--path_use"), type = "character", default = NULL, help = "Path to use for the analysis"),
  make_option(c("--smooth"), type = "logical", default = FALSE, help = "Smooth data (default FALSE)"),
  make_option(c("-r", "--result_type"), type = "character", default = "Top10", help = "Result type: Top10, tab_simple, or tab_all")
)

# Parse the command-line arguments
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Load data files (example for expr and psi)
expr_data <- read.table(args$expr, header = TRUE, row.names = 1)
psi_data <- read.table(args$psi, header = TRUE, row.names = 1)


library(MRAS)
# Call the MRAS function with the command-line arguments
MRAS(input_type = args$input_type,
     expr = expr_data,
     psi = psi_data,
     rbp_interested = args$rbp_interested,
     m = args$m,
     n = args$n,
     RBP_cutoff = args$RBP_cutoff,
     DS_pvalue = args$DS_pvalue,
     DS_dPSI = args$DS_dPSI,
     rbp_event_deal_all_total = args$rbp_event_deal_all_total,
     rbp_event_deal_all = args$rbp_event_deal_all,
     method = args$method,
     BS = NULL,  # Adjust accordingly
     group = args$group,
     dpsi_network_threshold = args$dpsi_network_threshold,
     Regulate_threshold = args$Regulate_threshold,
     rbp_net_mat_group = args$rbp_net_mat_group,
     num1 = args$num1,
     num2 = args$num2,
     sc = args$sc,
     threads = args$threads,    
     path_use = args$path_use,
     smooth = args$smooth,
     result_type = args$result_type)