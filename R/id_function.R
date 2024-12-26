#' id_find
#'
#' @param file_path the input files path.
#' @param software Commonly used software for detecting differential splicing events (rMTAS, JUM, SUPPA).
#' @param subtype if `software` = "rMATS",subtype should be "JCEC";if `software` = "JUM",subtype should be "simplified" or "detailed";
#'
#' @return Returns a PSI matrix with canonical splicing event names
#' @export
#'
#'
id_find<-function(file_path,software,subtype = NULL){

  ######rMATS#########
  if (software == "rMATS") {
    path_rmats<-file_path
    file_list<-list.files(path_rmats)
    if (subtype == "JCEC"){
      if (file.exists(paste0(path_rmats,"/PSI_JCEC.txt"))) file.remove(paste0(path_rmats,"/PSI_JCEC.txt"))
      if ("SE.MATS.JCEC.txt" %in% file_list){
        SE<-read.delim(paste0(path_rmats,"/SE.MATS.JCEC.txt"),header = T)
        SE<-SE[,c(1:11,17,18,13,14,21)]
        SE<-as.data.frame(t(apply(SE, 1, function(x){
          x[6]<-as.numeric(x[6])+1
          x[8]<-as.numeric(x[8])+1
          x[10]<-as.numeric(x[10])+1
          return(x)
        })))

        SE[,1]<-paste0(SE[,3],"_","ES","_",SE[,4],"_",SE[,5],"_",SE[,6],"_",SE[,7],"_",SE[,8],"_",SE[,9],"_",SE[,10],"_",SE[,11])
        write.table(SE[,c(1,16)],file = paste0(path_rmats,"/PSI_JCEC.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      if ("RI.MATS.JCEC.txt" %in% file_list){
        RI<-read.delim(paste0(path_rmats,"/RI.MATS.JCEC.txt"),header = T)
        RI<-RI[,c(1:11,17,18,13,14,21)]
        RI<-as.data.frame(t(apply(RI, 1, function(x){
          x[6]<-as.numeric(x[6])+1
          x[8]<-as.numeric(x[8])+1
          x[10]<-as.numeric(x[10])+1
          return(x)
        })))

        RI[,1]<-paste0(RI[,3],"_","IR","_",RI[,4],"_",RI[,5],"_",RI[,6],"_",RI[,7],"_",RI[,8],"_",RI[,9],"_",RI[,10],"_",RI[,11])
        write.table(RI[,c(1,16)],file = paste0(path_rmats,"/PSI_JCEC.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      if ("A3SS.MATS.JCEC.txt" %in% file_list){
        A3SS<-read.delim(paste0(path_rmats,"/A3SS.MATS.JCEC.txt"),header = T)
        A3SS<-A3SS[,c(1:11,17,18,13,14,21)]
        A3SS<-as.data.frame(t(apply(A3SS, 1, function(x){
          x[6]<-as.numeric(x[6])+1
          x[8]<-as.numeric(x[8])+1
          x[10]<-as.numeric(x[10])+1
          return(x)
        })))

        A3SS[,1]<-paste0(A3SS[,3],"_","A3SS","_",A3SS[,4],"_",A3SS[,5],"_",A3SS[,6],"_",A3SS[,7],"_",A3SS[,8],"_",A3SS[,9],"_",A3SS[,10],"_",A3SS[,11])
        write.table(A3SS[,c(1,16)],file = paste0(path_rmats,"/PSI_JCEC.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      if ("A5SS.MATS.JCEC.txt" %in% file_list){
        A5SS<-read.delim(paste0(path_rmats,"/A5SS.MATS.JCEC.txt"),header = T)
        A5SS<-A5SS[,c(1:11,17,18,13,14,21)]
        A5SS<-as.data.frame(t(apply(A5SS, 1, function(x){
          x[6]<-as.numeric(x[6])+1
          x[8]<-as.numeric(x[8])+1
          x[10]<-as.numeric(x[10])+1
          return(x)
        })))

        A5SS[,1]<-paste0(A5SS[,3],"_","A5SS","_",A5SS[,4],"_",A5SS[,5],"_",A5SS[,6],"_",A5SS[,7],"_",A5SS[,8],"_",A5SS[,9],"_",A5SS[,10],"_",A5SS[,11])
        write.table(A5SS[,c(1,16)],file = paste0(path_rmats,"/PSI_JCEC.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      if ("MXE.MATS.JCEC.txt" %in% file_list){
        MXE<-read.delim(paste0(path_rmats,"/MXE.MATS.JCEC.txt"),header = T)
        MXE<-MXE[,c(1:13,19:20,15:16,23)]
        MXE<-as.data.frame(t(apply(MXE, 1, function(x){
          x[6]<-as.numeric(x[6])+1
          x[8]<-as.numeric(x[8])+1
          x[10]<-as.numeric(x[10])+1
          x[12]<-as.numeric(x[12])+1
          return(x)
        })))

        MXE[,1]<-paste0(MXE[,3],"_","MEX","_",MXE[,4],"_",MXE[,5],"_",MXE[,10],"_",MXE[,11],"_",MXE[,6],"_",MXE[,7],"_",MXE[,8],"_",MXE[,9],"_",MXE[,12],"_",MXE[,13])
        write.table(MXE[,c(1,18)],file = paste0(path_rmats,"/PSI_JCEC.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)

      }
      data<-data.table::fread(paste0(path_rmats,"/PSI_JCEC.txt"),sep = ",",header = F)
      data[,1]<-apply(data,1,function(x){
        return(gsub(" ","",x[1]))
      })
      data.table::fwrite(data,file = paste0(path_rmats,"/PSI_JCEC.txt"),
                         row.names = F,col.names = F,
                         sep = ",",quote = F)
      data<-data.table::fread(paste0(path_rmats,"/PSI_JCEC.txt"),sep = ",",header = F)
      data[is.na(data)]<-0
      return(data)
    }
    if (subtype == "JC"){

    }
  }


  ######JUM#########
  if (software == "JUM"){
    path_jum<-file_path
    file_list<-list.files(path_jum)
    if (subtype == "simplified"){
      if (file.exists(paste0(path_jum,"/PSI_simplified.txt"))) file.remove(paste0(path_jum,"/PSI_simplified.txt"))
      if ("AS_differential_JUM_output_cassette_exon_events_pvalue_0.05_final_simplified.txt" %in% file_list){
        SE<-read.delim(paste0(path_jum,"/AS_differential_JUM_output_cassette_exon_events_pvalue_0.05_final_simplified.txt"),header = T)
        SE<-SE[complete.cases(SE),]
        SE[,3]<-apply(SE,1,function(x){
          return(strsplit(x[2],"_")[[1]][1])
        })
        SE[,3]<-paste0("chr",SE[,3])
        SE[,4]<-apply(SE,1,function(x){
          return(strsplit(x[2],"_")[[1]][2])
        })
        SE[,5]<-apply(SE,1,function(x){
          return(strsplit(x[2],"_")[[1]][3])
        })
        SE[,6]<-apply(SE,1,function(x){
          return(strsplit(x[2],"_")[[1]][4])
        })
        SE[,6]<-as.numeric(SE[,6])+1
        SE[,7]<-apply(SE,1,function(x){
          return(strsplit(x[2],"_")[[1]][5])
        })
        SE[,8]<-apply(SE,1,function(x){
          return(strsplit(x[2],"_")[[1]][6])
        })
        SE[,8]<-as.numeric(SE[,8])+2
        SE[,2]<-paste0(SE[,1],"_","ES","_",SE[,3],"_",SE[,4],"_",SE[,6],"_",SE[,7],"_","x","_",SE[,5],"_",SE[,8],"_","x")
        write.table(SE[,2],file = paste0(path_jum,"/PSI_simplified.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      if ("AS_differential_JUM_output_intron_retention_pvalue_0.05_final_simplified.txt" %in% file_list){
        RI<-read.delim(paste0(path_jum,"/AS_differential_JUM_output_intron_retention_pvalue_0.05_final_simplified.txt"),header = T)
        RI<-RI[complete.cases(RI),]
        RI[,3]<-apply(RI,1,function(x){
          return(strsplit(x[2],"_")[[1]][1])
        })
        RI[,3]<-paste0("chr",RI[,3])
        RI[,4]<-apply(RI,1,function(x){
          return(strsplit(x[2],"_")[[1]][2])
        })
        RI[,5]<-apply(RI,1,function(x){
          return(strsplit(x[2],"_")[[1]][3])
        })
        RI[,6]<-apply(RI,1,function(x){
          return(strsplit(x[2],"_")[[1]][4])
        })
        RI[,6]<-as.numeric(RI[,6])+2

        RI[,2]<-paste0(RI[,1],"_","IR","_",RI[,3],"_",RI[,4],"_","x","_","x","_","x","_",RI[,5],"_",RI[,6],"_","x")
        write.table(RI[,2],file = paste0(path_jum,"/PSI_simplified.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      if ("AS_differential_JUM_output_A3SS_events_pvalue_0.05_final_simplified.txt" %in% file_list){
        A3SS<-read.delim(paste0(path_jum,"/AS_differential_JUM_output_A3SS_events_pvalue_0.05_final_simplified.txt"),header = T)
        A3SS<-A3SS[complete.cases(A3SS),]
        A3SS[,3]<-apply(A3SS,1,function(x){
          return(strsplit(x[2],"_")[[1]][1])
        })
        A3SS[,3]<-paste0("chr",A3SS[,3])
        A3SS[,4]<-apply(A3SS,1,function(x){
          return(strsplit(x[2],"_")[[1]][2])
        })
        A3SS[,5]<-apply(A3SS,1,function(x){
          return(strsplit(x[2],"_")[[1]][3])
        })
        A3SS[,6]<-apply(A3SS,1,function(x){
          return(strsplit(x[2],"_")[[1]][4])
        })
        A3SS[,6]<-ifelse(A3SS[,4] == "+",as.numeric(A3SS[,6])+2,A3SS[,6])
        A3SS[,7]<-apply(A3SS,1,function(x){
          return(strsplit(x[2],"_")[[1]][5])
        })
        A3SS[,7]<-as.numeric(A3SS[,7])+2
        A3SS[,2]<-ifelse(A3SS[,4] == "+",paste0(A3SS[,1],"_","A3SS","_",A3SS[,3],"_",A3SS[,4],"_",A3SS[,6],"_","x","_",A3SS[,7],"_","x","_","x","_",A3SS[,5]),
                         paste0(A3SS[,1],"_","A3SS","_",A3SS[,3],"_",A3SS[,4],"_","x","_",A3SS[,6],"_","x","_",A3SS[,5],"_",A3SS[,7],"_","x"))

        write.table(A3SS[,2],file = paste0(path_jum,"/PSI_simplified.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      if ("AS_differential_JUM_output_A5SS_events_pvalue_0.05_final_simplified.txt" %in% file_list){
        A5SS<-read.delim(paste0(path_jum,"/AS_differential_JUM_output_A5SS_events_pvalue_0.05_final_simplified.txt"),header = T)
        A5SS<-A5SS[complete.cases(A5SS),]
        A5SS[,3]<-apply(A5SS,1,function(x){
          return(strsplit(x[2],"_")[[1]][1])
        })
        A5SS[,3]<-paste0("chr",A5SS[,3])
        A5SS[,4]<-apply(A5SS,1,function(x){
          return(strsplit(x[2],"_")[[1]][2])
        })
        A5SS[,5]<-apply(A5SS,1,function(x){
          return(strsplit(x[2],"_")[[1]][3])
        })
        A5SS[,5]<-A5SS[,5]
        A5SS[,6]<-apply(A5SS,1,function(x){
          return(strsplit(x[2],"_")[[1]][4])
        })
        A5SS[,6]<-ifelse(A5SS[,4] == "-",as.numeric(A5SS[,6])+2,A5SS[,6])
        A5SS[,7]<-apply(A5SS,1,function(x){
          return(strsplit(x[2],"_")[[1]][5])
        })
        A5SS[,7]<-ifelse(A5SS[,4] == "-",as.numeric(A5SS[,7])+2,A5SS[,7])
        A5SS[,2]<-ifelse(A5SS[,4] == "+",paste0(A5SS[,1],"_","A5SS","_",A5SS[,3],"_",A5SS[,4],"_","x","_",A5SS[,6],"_","x","_",A5SS[,5],"_",A5SS[,7],"_","x"),
                         paste0(A5SS[,1],"_","A5SS","_",A5SS[,3],"_",A5SS[,4],"_",A5SS[,6],"_","x","_",A5SS[,7],"_","x","_","x","_",A5SS[,5]))
        write.table(A5SS[,2],file = paste0(path_jum,"/PSI_simplified.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      if ("AS_differential_JUM_output_MXE_events_pvalue_0.05_final_simplified.txt" %in% file_list){
        MEX<-read.delim(paste0(path_jum,"/AS_differential_JUM_output_MXE_events_pvalue_0.05_final_simplified.txt"),header = T)
        MEX<-MEX[complete.cases(MEX),]
        MEX[,3]<-apply(MEX,1,function(x){
          return(strsplit(x[2],"_")[[1]][1])
        })
        MEX[,3]<-paste0("chr",MEX[,3])
        MEX[,4]<-apply(MEX,1,function(x){
          return(strsplit(x[2],"_")[[1]][2])
        })
        MEX[,5]<-apply(MEX,1,function(x){
          return(strsplit(x[2],"_")[[1]][3])
        })
        MEX[,6]<-apply(MEX,1,function(x){
          return(strsplit(x[2],"_")[[1]][4])
        })
        MEX[,7]<-apply(MEX,1,function(x){
          return(strsplit(x[6],"-")[[1]][2])
        })
        MEX[,6]<-apply(MEX,1,function(x){
          return(strsplit(x[6],"-")[[1]][1])
        })
        MEX[,6]<-as.numeric(MEX[,6])+1
        MEX[,7]<-as.numeric(MEX[,7])
        MEX[,8]<-apply(MEX,1,function(x){
          return(strsplit(x[2],"_")[[1]][5])
        })
        MEX[,9]<-apply(MEX,1,function(x){
          return(strsplit(x[8],"-")[[1]][2])
        })
        MEX[,8]<-apply(MEX,1,function(x){
          return(strsplit(x[8],"-")[[1]][1])
        })
        MEX[,8]<-as.numeric(MEX[,8])+1
        MEX[,9]<-as.numeric(MEX[,9])
        MEX[,10]<-apply(MEX,1,function(x){
          return(strsplit(x[2],"_")[[1]][6])
        })
        MEX[,10]<-as.numeric(MEX[,10])+2
        MEX[,2]<-paste0(MEX[,1],"_","MEX","_",MEX[,3],"_",MEX[,4],"_","x","_",MEX[,5],"_",MEX[,6],"_",MEX[,7],"_",MEX[,8],"_",MEX[,9],"_",MEX[,10],"_","x")
        write.table(MEX[,2],file = paste0(path_jum,"/PSI_simplified.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      data<-read.delim(paste0(path_jum,"/PSI_simplified.txt"),sep = ",",header = F)
      data[,1]<-gsub(" ","",data[,1])
      write.table(data,file = paste0(path_jum,"/PSI_simplified.txt"),
                  row.names = F,col.names = F,
                  sep = ",",quote = F)
      data<-fread(paste0(path_jum,"/PSI_simplified.txt"),sep = ",",header = F)
      data[is.na(data)]<-0
      return(data)
    }
    if (subtype == "detailed"){
      if (file.exists(paste0(path_jum,"/PSI_detailed.txt"))) file.remove(paste0(path_jum,"/PSI_detailed.txt"))
      if ("AS_differential_JUM_output_cassette_exon_events_pvalue_0.05_final_detailed.txt" %in% file_list){
        SE_pre<-read.delim(paste0(path_jum,"/AS_differential_JUM_output_cassette_exon_events_pvalue_0.05_final_detailed.txt"),header = T)
        SE_pre<-SE_pre[which(SE_pre$sub_junction_strand %in% c("+","-")),]
        SE_pre[,1]<-paste0(SE_pre[,1],"_",SE_pre[,2])
        SE<-matrix(NA,nrow = length(unique(SE_pre[,1])),ncol = ncol(SE_pre))
        SE[,1]<-unique(SE_pre[,1])
        SE[,2]<-apply(SE,1,function(x){
          return(strsplit(x[1],"_")[[1]][1])
        })
        SE[,3]<-apply(SE,1,function(x){
          return(strsplit(x[1],"_")[[1]][2])
        })
        SE[,3]<-paste0("chr",SE[,3])
        SE[,4]<-apply(SE,1,function(x){
          return(strsplit(x[1],"_")[[1]][3])
        })
        SE[,5]<-apply(SE,1,function(x){
          return(strsplit(x[1],"_")[[1]][4])
        })
        SE[,6]<-apply(SE,1,function(x){
          return(strsplit(x[1],"_")[[1]][5])
        })
        SE[,6]<-as.numeric(SE[,6])+1
        SE[,7]<-apply(SE,1,function(x){
          return(strsplit(x[1],"_")[[1]][6])
        })
        SE[,8]<-apply(SE,1,function(x){
          return(strsplit(x[1],"_")[[1]][7])
        })
        SE[,8]<-as.numeric(SE[,8])+2
        SE[,9]<-paste0(SE[,2],"_","ES","_",SE[,3],"_",SE[,4],"_",SE[,6],"_",SE[,7],"_","x","_",SE[,5],"_",SE[,8],"_","x")

        sample_num<-(ncol(SE_pre)-17)/2
        for (id in 1:nrow(SE)) {
          event<-SE[id,1]
          event_mat<-SE_pre[which(SE_pre[,1] == event),1:(16+sample_num)]
          SE_type<-as.data.frame(table(event_mat$sub_junction_size))
          inclusion<-which(event_mat$sub_junction_size %in% SE_type$Var1[which(SE_type$Freq == 1)])
          skipping<-which(event_mat$sub_junction_size %in% SE_type$Var1[which(SE_type$Freq == 2)])
          event_mat_deal<-as.data.frame(event_mat[,17:ncol(event_mat)])
          event_mat_deal<-apply(event_mat_deal, 2, as.numeric)
          psi<-colSums(event_mat_deal[inclusion,])/colSums(event_mat_deal)
          SE[id,10:(9+sample_num)]<-round(psi,2)
        }
        write.table(SE[,9:(9+sample_num)],file = paste0(path_jum,"/PSI_detailed.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      if ("AS_differential_JUM_output_intron_retention_pvalue_0.05_final_detailed.txt" %in% file_list){
        RI_pre<-read.delim(paste0(path_jum,"/AS_differential_JUM_output_intron_retention_pvalue_0.05_final_detailed.txt"),header = T)
        RI_pre[,1]<-paste0(RI_pre[,1],"_",RI_pre[,2])
        RI_pre<-RI_pre[which(RI_pre$sub_junction_strand %in% c("+","-")),]
        RI<-matrix(NA,nrow = length(unique(RI_pre[,1])),ncol = ncol(RI_pre))
        RI[,1]<-unique(RI_pre[,1])
        RI[,2]<-apply(RI,1,function(x){
          return(strsplit(x[1],"_")[[1]][1])
        })
        RI[,3]<-apply(RI,1,function(x){
          return(strsplit(x[1],"_")[[1]][2])
        })
        RI[,3]<-paste0("chr",RI[,3])
        RI[,4]<-apply(RI,1,function(x){
          return(strsplit(x[1],"_")[[1]][3])
        })
        RI[,5]<-apply(RI,1,function(x){
          return(strsplit(x[1],"_")[[1]][4])
        })
        RI[,6]<-apply(RI,1,function(x){
          return(strsplit(x[1],"_")[[1]][5])
        })
        RI[,6]<-as.numeric(RI[,6])+2

        RI[,7]<-paste0(RI[,2],"_","IR","_",RI[,3],"_",RI[,4],"_","x","_","x","_","x","_",RI[,5],"_",RI[,6],"_","x")
        sample_num<-(ncol(RI_pre)-17)/2
        for (id in 1:nrow(RI)) {
          event<-RI[id,1]
          event_mat<-RI_pre[which(RI_pre[,1] == event),1:(16+sample_num)]
          RI_type<-as.data.frame(table(event_mat$sub_junction_size))
          total<-sort(as.numeric(event_mat$sub_junction_size))[(0.5*nrow(event_mat)+1):nrow(event_mat)]
          inclusion<-which(event_mat$sub_junction_size %in% total)
          skipping<-which(!(event_mat$sub_junction_size %in% total))
          event_mat_deal<-as.data.frame(event_mat[,17:ncol(event_mat)])
          event_mat_deal<-apply(event_mat_deal, 2, as.numeric)
          psi<-colSums(event_mat_deal[inclusion,])/colSums(event_mat_deal)
          RI[id,8:(7+sample_num)]<-round(psi,2)
        }
        write.table(RI[,7:(7+sample_num)],file = paste0(path_jum,"/PSI_detailed.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      if ("AS_differential_JUM_output_A3SS_events_pvalue_0.05_final_detailed.txt" %in% file_list){
        A3SS_pre<-read.delim(paste0(path_jum,"/AS_differential_JUM_output_A3SS_events_pvalue_0.05_final_detailed.txt"),header = T)
        id_deal<-unique(A3SS_pre$AS_event_ID[which(A3SS_pre$sub_junction_ID == "J003")])
        A3SS_pre_deal<-matrix(NA,1,1)
        for (ii in unique(A3SS_pre$AS_event_ID)) {
          if (ii %in% id_deal){
            mat<-as.data.frame((A3SS_pre[which(A3SS_pre$AS_event_ID == ii),]))
            id_all<-strsplit(ii,"_")[[1]]
            id_loc<-ifelse(id_all[2] == "+",id_all[3],id_all[length(id_all)])
            id_change<-setdiff(id_all[3:length(id_all)],id_loc)
            com_mat<-t(as.data.frame(combn(1:length(unique(mat$sub_junction_ID)), 2)))
            deal_mat<-matrix(NA,nrow = nrow(com_mat)*2,ncol = ncol(A3SS_pre))
            for (dd in 1:nrow(com_mat)) {
              id1<-as.numeric(com_mat[dd,1])
              id2<-as.numeric(com_mat[dd,2])
              deal_mat[dd*2-1,]<-unlist(mat[which(mat$sub_junction_ID == mat$sub_junction_ID[id1]),])
              deal_mat[dd*2,]<-unlist(mat[which(mat$sub_junction_ID == mat$sub_junction_ID[id2]),])
              deal_mat[c(dd*2-1,dd*2),2]<-ifelse(id_all[2] == "+",
                                                 tj(c(id_all[1],id_all[2],id_loc,id_change[id1],id_change[id2]),type = "_"),
                                                 tj(c(id_all[1],id_all[2],id_change[id1],id_change[id2],id_loc),type = "_"))
            }
            }else{
            deal_mat<-as.matrix(A3SS_pre[which(A3SS_pre$AS_event_ID == ii),])
          }
          if (is.na(A3SS_pre_deal[1,1])){
            A3SS_pre_deal<-deal_mat
          }else{
            A3SS_pre_deal<-rbind(A3SS_pre_deal,deal_mat)
          }
        }
        # cat(dim(A3SS_pre),"\t")
        # cat(dim(A3SS_pre_deal),"\t")
        # cat(dim(A3SS_pre_deal)[1]-dim(A3SS_pre)[1],"\t")
        A3SS_pre_deal<-as.data.frame(A3SS_pre_deal)
        A3SS_pre_deal[,1]<-paste0(A3SS_pre_deal[,1],"_",A3SS_pre_deal[,2])
        A3SS_pre_deal<-A3SS_pre_deal[which(A3SS_pre_deal$sub_junction_strand %in% c("+","-")),]
        A3SS<-matrix(NA,nrow = length(unique(A3SS_pre_deal[,1])),ncol = ncol(A3SS_pre_deal))
        A3SS[,1]<-unique(A3SS_pre_deal[,1])
        A3SS[,2]<-apply(A3SS,1,function(x){
          return(strsplit(x[1],"_")[[1]][1])
        })
        A3SS[,3]<-apply(A3SS,1,function(x){
          return(strsplit(x[1],"_")[[1]][2])
        })
        A3SS[,3]<-paste0("chr",A3SS[,3])
        A3SS[,4]<-apply(A3SS,1,function(x){
          return(strsplit(x[1],"_")[[1]][3])
        })
        A3SS[,5]<-apply(A3SS,1,function(x){
          return(strsplit(x[1],"_")[[1]][4])
        })
        A3SS[,6]<-apply(A3SS,1,function(x){
          return(strsplit(x[1],"_")[[1]][5])
        })
        A3SS[,6]<-ifelse(A3SS[,4] == "+",as.numeric(A3SS[,6])+2,A3SS[,6])
        A3SS[,7]<-apply(A3SS,1,function(x){
          return(strsplit(x[1],"_")[[1]][6])
        })
        A3SS[,7]<-as.numeric(A3SS[,7])+2
        A3SS[,8]<-ifelse(A3SS[,4] == "+",paste0(A3SS[,2],"_","A3SS","_",A3SS[,3],"_",A3SS[,4],"_",A3SS[,6],"_","x","_",A3SS[,7],"_","x","_","x","_",A3SS[,5]),
               paste0(A3SS[,2],"_","A3SS","_",A3SS[,3],"_",A3SS[,4],"_","x","_",A3SS[,6],"_","x","_",A3SS[,5],"_",A3SS[,7],"_","x"))
        sample_num<-(ncol(A3SS_pre_deal)-17)/2
        for (id in 1:nrow(A3SS)) {
          event<-A3SS[id,1]
          event_mat<-A3SS_pre_deal[which(A3SS_pre_deal[,1] == event),1:(16+sample_num)]
          inclusion<-which(event_mat$sub_junction_size %in% min(event_mat$sub_junction_size))
          skipping<-which(event_mat$sub_junction_size %in% max(event_mat$sub_junction_size))
          event_mat_deal<-as.data.frame(event_mat[,17:ncol(event_mat)])
          event_mat_deal<-apply(event_mat_deal, 2, as.numeric)
          psi<-event_mat_deal[inclusion,]/colSums(event_mat_deal)
          A3SS[id,9:(8+sample_num)]<-round(psi,2)
        }
        write.table(A3SS[,8:(8+sample_num)],file = paste0(path_jum,"/PSI_detailed.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      if ("AS_differential_JUM_output_A5SS_events_pvalue_0.05_final_detailed.txt" %in% file_list){
        A5SS_pre<-read.delim(paste0(path_jum,"/AS_differential_JUM_output_A5SS_events_pvalue_0.05_final_detailed.txt"),header = T)
        id_deal<-unique(A5SS_pre$AS_event_ID[which(A5SS_pre$sub_junction_ID == "J003")])
        A5SS_pre_deal<-matrix(NA,1,1)
        for (ii in unique(A5SS_pre$AS_event_ID)) {
          if (ii %in% id_deal){
            mat<-as.data.frame((A5SS_pre[which(A5SS_pre$AS_event_ID == ii),]))
            id_all<-strsplit(ii,"_")[[1]]
            id_loc<-ifelse(id_all[2] == "-",id_all[3],id_all[length(id_all)])
            id_change<-setdiff(id_all[3:length(id_all)],id_loc)
            com_mat<-t(as.data.frame(combn(1:length(unique(mat$sub_junction_ID)), 2)))
            deal_mat<-matrix(NA,nrow = nrow(com_mat)*2,ncol = ncol(A5SS_pre))
            for (dd in 1:nrow(com_mat)) {
              id1<-as.numeric(com_mat[dd,1])
              id2<-as.numeric(com_mat[dd,2])
              deal_mat[dd*2-1,]<-unlist(mat[which(mat$sub_junction_ID == mat$sub_junction_ID[id1]),])
              deal_mat[dd*2,]<-unlist(mat[which(mat$sub_junction_ID == mat$sub_junction_ID[id2]),])
              deal_mat[c(dd*2-1,dd*2),2]<-ifelse(id_all[2] == "+",
                                                 tj(c(id_all[1],id_all[2],id_loc,id_change[id1],id_change[id2]),type = "_"),
                                                 tj(c(id_all[1],id_all[2],id_change[id1],id_change[id2],id_loc),type = "_"))
            }
          }else{
            deal_mat<-as.matrix(A5SS_pre[which(A5SS_pre$AS_event_ID == ii),])
          }
          if (is.na(A5SS_pre_deal[1,1])){
            A5SS_pre_deal<-deal_mat
          }else{
            A5SS_pre_deal<-rbind(A5SS_pre_deal,deal_mat)
          }
        }
        # cat(dim(A5SS_pre),"\t")
        # cat(dim(A5SS_pre_deal),"\t")
        # cat(dim(A5SS_pre_deal)[1]-dim(A5SS_pre)[1],"\t")
        A5SS_pre_deal<-as.data.frame(A5SS_pre_deal)
        A5SS_pre_deal[,1]<-paste0(A5SS_pre_deal[,1],"_",A5SS_pre_deal[,2])
        A5SS_pre_deal<-A5SS_pre_deal[which(A5SS_pre_deal$sub_junction_strand %in% c("+","-")),]
        A5SS<-matrix(NA,nrow = length(unique(A5SS_pre_deal[,1])),ncol = ncol(A5SS_pre_deal))
        A5SS[,1]<-unique(A5SS_pre_deal[,1])
        A5SS[,2]<-apply(A5SS,1,function(x){
          return(strsplit(x[1],"_")[[1]][1])
        })
        A5SS[,3]<-apply(A5SS,1,function(x){
          return(strsplit(x[1],"_")[[1]][2])
        })
        A5SS[,3]<-paste0("chr",A5SS[,3])
        A5SS[,4]<-apply(A5SS,1,function(x){
          return(strsplit(x[1],"_")[[1]][3])
        })
        A5SS[,5]<-apply(A5SS,1,function(x){
          return(strsplit(x[1],"_")[[1]][4])
        })
        A5SS[,5]<-as.numeric(A5SS[,5])
        A5SS[,6]<-apply(A5SS,1,function(x){
          return(strsplit(x[1],"_")[[1]][5])
        })
        A5SS[,6]<-ifelse(A5SS[,4] == "-",as.numeric(A5SS[,6])+2,A5SS[,6])
        A5SS[,7]<-apply(A5SS,1,function(x){
          return(strsplit(x[1],"_")[[1]][6])
        })
        A5SS[,7]<-ifelse(A5SS[,4] == "-",as.numeric(A5SS[,7])+2,A5SS[,7])
        A5SS[,8]<-ifelse(A5SS[,4] == "+",paste0(A5SS[,2],"_","A5SS","_",A5SS[,3],"_",A5SS[,4],"_","x","_",A5SS[,6],"_","x","_",A5SS[,5],"_",A5SS[,7],"_","x"),
                                         paste0(A5SS[,2],"_","A5SS","_",A5SS[,3],"_",A5SS[,4],"_",A5SS[,6],"_","x","_",A5SS[,7],"_","x","_","x","_",A5SS[,5]))
        sample_num<-(ncol(A5SS_pre_deal)-17)/2
        for (id in 1:nrow(A5SS)) {
          event<-A5SS[id,1]
          event_mat<-A5SS_pre_deal[which(A5SS_pre_deal[,1] == event),1:(16+sample_num)]
          inclusion<-which(event_mat$sub_junction_size %in% min(event_mat$sub_junction_size))
          skipping<-which(event_mat$sub_junction_size %in% max(event_mat$sub_junction_size))
          event_mat_deal<-as.data.frame(event_mat[,17:ncol(event_mat)])
          event_mat_deal<-apply(event_mat_deal, 2, as.numeric)
          psi<-event_mat_deal[inclusion,]/colSums(event_mat_deal)
          A5SS[id,9:(8+sample_num)]<-round(psi,2)
        }
        write.table(A5SS[,8:(8+sample_num)],file = paste0(path_jum,"/PSI_detailed.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      if ("AS_differential_JUM_output_MXE_events_pvalue_0.05_final_detailed.txt" %in% file_list){
        MEX_pre<-read.delim(paste0(path_jum,"/AS_differential_JUM_output_MXE_events_pvalue_0.05_final_detailed.txt"),header = T)
        MEX_pre<-MEX_pre[which(MEX_pre$sub_junction_strand %in% c("+","-")),]
        MEX_pre[,1]<-paste0(MEX_pre[,1],"_",MEX_pre[,2])
        MEX<-matrix(NA,nrow = length(unique(MEX_pre[,1])),ncol = ncol(MEX_pre))
        MEX[,1]<-unique(MEX_pre[,1])
        MEX[,2]<-apply(MEX,1,function(x){
          return(strsplit(x[1],"_")[[1]][1])
        })
        MEX[,3]<-apply(MEX,1,function(x){
          return(strsplit(x[1],"_")[[1]][2])
        })
        MEX[,3]<-paste0("chr",MEX[,3])
        MEX[,4]<-apply(MEX,1,function(x){
          return(strsplit(x[1],"_")[[1]][3])
        })
        MEX[,5]<-apply(MEX,1,function(x){
          return(strsplit(x[1],"_")[[1]][4])
        })
        MEX[,6]<-apply(MEX,1,function(x){
          return(strsplit(x[1],"_")[[1]][5])
        })
        MEX[,7]<-apply(MEX,1,function(x){
          return(strsplit(x[6],"-")[[1]][2])
        })
        MEX[,6]<-apply(MEX,1,function(x){
          return(strsplit(x[6],"-")[[1]][1])
        })
        MEX[,6]<-as.numeric(MEX[,6])+1
        MEX[,7]<-as.numeric(MEX[,7])
        MEX[,8]<-apply(MEX,1,function(x){
          return(strsplit(x[1],"_")[[1]][6])
        })
        MEX[,9]<-apply(MEX,1,function(x){
          return(strsplit(x[8],"-")[[1]][2])
        })
        MEX[,8]<-apply(MEX,1,function(x){
          return(strsplit(x[8],"-")[[1]][1])
        })
        MEX[,8]<-as.numeric(MEX[,8])+1
        MEX[,9]<-as.numeric(MEX[,9])
        MEX[,10]<-apply(MEX,1,function(x){
          return(strsplit(x[1],"_")[[1]][7])
        })
        MEX[,10]<-as.numeric(MEX[,10])+2
        MEX[,11]<-paste0(MEX[,2],"_","MEX","_",MEX[,3],"_",MEX[,4],"_","x","_",MEX[,5],"_",MEX[,6],"_",MEX[,7],"_",MEX[,8],"_",MEX[,9],"_",MEX[,10],"_","x")
        sample_num<-(ncol(MEX_pre)-17)/2
        for (id in 1:nrow(MEX)) {
          event<-MEX[id,1]
          event_mat<-MEX_pre[which(MEX_pre[,1] == event),1:(16+sample_num)]
          id_pre<-strsplit(tj(strsplit(event,"_")[[1]][c(4,5,7)],"-"),"-")[[1]]
          id_pre2<-which((event_mat$sub_junction_start_coor %in% id_pre)&(event_mat$sub_junction_end_coor %in% id_pre))
          inclusion<-which(event_mat$sub_junction_ID == event_mat$sub_junction_ID[id_pre2])
          event_mat_deal<-as.data.frame(event_mat[,17:ncol(event_mat)])
          event_mat_deal<-apply(event_mat_deal, 2, as.numeric)
          psi<-colSums(event_mat_deal[inclusion,])/colSums(event_mat_deal)
          MEX[id,12:(11+sample_num)]<-round(psi,2)
        }

        write.table(MEX[,11:(11+sample_num)],file = paste0(path_jum,"/PSI_detailed.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }
      data<-read.delim(paste0(path_jum,"/PSI_detailed.txt"),sep = ",",header = F)
      data[,1]<-gsub(" ","",data[,1])
      write.table(data,file = paste0(path_jum,"/PSI_detailed.txt"),
                  row.names = F,col.names = F,
                  sep = ",",quote = F)
      data<-fread(paste0(path_jum,"/PSI_detailed.txt"),sep = ",",header = F)
      data[is.na(data)]<-0
      return(data)
    }
  }


  ######SUPPA#########
  if (software == "SUPPA"){
    path_suppa<-file_path
    file_list<-list.files(path_suppa)
    data<-read.delim(paste0(path_suppa,list.files(path_suppa)[grep("*.psi",list.files(path_suppa))]))
    name<-as.data.frame(rownames(data))
    name[,2]<-apply(name, 1, function(x){
      return(strsplit(x[1],";")[[1]][1])
    })
    name[,3]<-apply(name, 1, function(x){
      a<-strsplit(x[1],";")[[1]][2]
      return(strsplit(a,":")[[1]][1])
    })
    name[,4]<-apply(name, 1, function(x){
      return(strsplit(x[1],":")[[1]][2])
    })
    name[,4]<-paste0("chr",name[,4])
    name[,5]<-name[,1]
    name[,6]<-apply(name, 1, function(x){
      return(rev(strsplit(x[1],":")[[1]])[1])
    })
    SE_id<-which(name[,3] == "SE")
    name[SE_id,5]<-apply(name[SE_id,],1,function(x){
      e1s2<-strsplit(x[1],":")[[1]][3]
      e2s3<-strsplit(x[1],":")[[1]][4]
      a1<-strsplit(e1s2,"-")[[1]][2]
      a2<-strsplit(e2s3,"-")[[1]][1]
      a3<-"x"
      a4<-strsplit(e1s2,"-")[[1]][1]
      a5<-strsplit(e2s3,"-")[[1]][2]
      a6<-"x"
      return(paste0(a1,"_",a2,"_",a3,"_",a4,"_",a5,"_",a6))
    })
    RI_id<-which(name[,3] == "RI")
    name[RI_id,5]<-apply(name[RI_id,],1,function(x){
      e1s2<-strsplit(x[1],":")[[1]][4]
      a1<-strsplit(x[1],":")[[1]][3]
      a2<-strsplit(x[1],":")[[1]][5]
      a3<-strsplit(x[1],":")[[1]][3]
      a4<-strsplit(e1s2,"-")[[1]][1]
      a5<-strsplit(e1s2,"-")[[1]][2]
      a6<-strsplit(x[1],":")[[1]][5]
      return(paste0(a1,"_",a2,"_",a3,"_",a4,"_",a5,"_",a6))
    })
    name[RI_id,3]<-"IR"
    A3_id<-which(name[,3] == "A3")
    name[A3_id,5]<-apply(name[A3_id,],1,function(x){
      # browser()
      s2<-strsplit(x[1],":")[[1]][3]
      s3<-strsplit(x[1],":")[[1]][4]
      a1<-strsplit(s2,"-")[[1]][1]
      a2<-strsplit(s2,"-")[[1]][2]
      a3<-strsplit(s3,"-")[[1]][1]
      a4<-strsplit(s3,"-")[[1]][2]
      if (x[6] == "+"){
        return(paste0(a2,"_","x","_",a4,"_","x","_","x","_",a1))
      }else{
        return(paste0("x","_",a1,"_","x","_",a3,"_",a2,"_","x"))
      }

    })
    name[A3_id,3]<-"A3SS"
    A5_id<-which(name[,3] == "A5")
    name[A5_id,5]<-apply(name[A5_id,],1,function(x){
      # browser()
      s2<-strsplit(x[1],":")[[1]][3]
      s3<-strsplit(x[1],":")[[1]][4]
      a1<-strsplit(s2,"-")[[1]][1]
      a2<-strsplit(s2,"-")[[1]][2]
      a3<-strsplit(s3,"-")[[1]][1]
      a4<-strsplit(s3,"-")[[1]][2]
      if (x[6] == "+"){
        return(paste0("x","_",a1,"_","x","_",a3,"_",a2,"_","x"))
      }else{
        return(paste0(a2,"_","x","_",a4,"_","x","_","x","_",a1))
      }

    })
    name[A5_id,3]<-"A5SS"
    MX_id<-which(name[,3] == "MX")
    name[MX_id,5]<-apply(name[MX_id,],1,function(x){
      e1s2<-strsplit(x[1],":")[[1]][3]
      e2s4<-strsplit(x[1],":")[[1]][4]
      e1s3<-strsplit(x[1],":")[[1]][5]
      e3s4<-strsplit(x[1],":")[[1]][6]
      a1<-"x"
      a2<-strsplit(e1s2,"-")[[1]][1]
      a3<-strsplit(e1s2,"-")[[1]][2]
      a4<-strsplit(e2s4,"-")[[1]][1]
      a5<-strsplit(e1s3,"-")[[1]][2]
      a6<-strsplit(e3s4,"-")[[1]][1]
      a7<-strsplit(e3s4,"-")[[1]][2]
      a8<-"x"
      return(paste0(a1,"_",a2,"_",a3,"_",a4,"_",a5,"_",a6,"_",a7,"_",a8))
    })
    name[MX_id,3]<-"MEX"
    name[,7]<-paste0(name[,2],"_",name[,3],"_",name[,4],"_",name[,6],"_",name[,5])
    rownames(data)<-name[,7]
    return(data)
  }
}

#' id_normalization
#'
#' @param symbol Gene symbol.
#' @param type Splicing event type
#' @param chr Chromosome number.
#' @param chain Positive and negative chains.
#' @param loc1,loc2,loc3,loc4,loc5,loc6,loc7,loc8 Splicing event coordinates. The unknown position is "x". The user needs to provide key locations or all locations in standard format.
#'
#' @return Standard splicing event name.
#' @export
#'
#'
id_normalization<-function(symbol,type="x",chr="x",chain="x",loc1="x",loc2="x",loc3="x",loc4="x",loc5="x",loc6="x",loc7="x",loc8="x"){
  variables <- list(type, chr, chain, loc1, loc2, loc3, loc4, loc5, loc6, loc7, loc8)
  for (var in which(variables != "x")) {
    if (length(symbol) != length(variables[[var]])){
      stop("You should input the same length with symbol!")
    }
  }
  id_dat<-data.frame(symbol = symbol,
                     type = type,
                     chr = chr,
                     chain = chain,
                     loc1 = loc1,
                     loc2 = loc2,
                     loc3 = loc3,
                     loc4 = loc4,
                     loc5 = loc5,
                     loc6 = loc6,
                     loc7 = loc7,
                     loc8 = loc8)
  id_nor<-unlist(apply(id_dat,1,function(x){
    if (x[2]=="MEX") {
      return(tj(x[1:12],type = "_"))
    }else{
      return(tj(x[1:10],type = "_"))
    }
  }))
  return(id_nor)
}

#' id_change
#'
#' @param id1 Source name needs to be converted.
#' @param id2 Target name needs to be converted to.
#' @param err_len Mismatch length.
#'
#' @return Match the name in id1 to id2.
#' @export
#'
#'
id_change<-function(id1,id2,err_len=2){
  del<-function(xx){
    if (!("x" %in% xx)){
      xx<-sort(unique(as.numeric(xx)))
      return(xx[-which(xx == max(xx) | xx == min(xx))])
    }else{
      xx<-sort(xx[which(xx != "x")])
      return(xx)
    }
  }

  id1_mat<-do.call(rbind,strsplit(id1,"_"))
  id1_mat<-cbind(id1,id1_mat)
  id1_mat<-as.data.frame(id1_mat)
  if (!('MEX' %in% id1_mat[,3])) {
    id1_mat$loc7<-NA
    id1_mat$loc8<-NA
    if (!('ES' %in% id1_mat[,3])){
      id2_mat$loc6<-NA
      if (unique(id1_mat[,3]) == "IR"){
        id2_mat$loc5<-NA
      }
    }
  }
  colnames(id1_mat)<-c("id1","symbol","type","chr","chain","loc1","loc2",
                       "loc3","loc4","loc5","loc6","loc7","loc8")
  id1_mat<-id1_mat[,1:13]
  id1_mat$id_new<-unlist(apply(id1_mat,1,function(x){
    if (x[5]=="MEX"){
      a<-tryCatch(tj(del(x[6:13]),"_"),warning = function(w){
        message(paste0(x[1],"_",w))
        return(NA)
      })

    }else{
      if (x[5]=="IR"){
        a<-tryCatch(tj(del(x[6:11]),"_"),warning = function(w){
          # browser()
          message(paste0(x[1],"_",w))
          return(NA)
        })
      }else{
        a<-tryCatch(tj(del(x[6:11]),"_"),warning = function(w){
          message(paste0(x[1],"_",w))
          return(NA)
        })
      }
    }
    return(a)
  }))
  id1_mat$n1<-unlist(apply(id1_mat,1,function(x){
    return(strsplit(x[14],"_")[[1]][1])
  }))
  id1_mat$n2<-unlist(apply(id1_mat,1,function(x){
    return(strsplit(x[14],"_")[[1]][2])
  }))
  id1_mat$n3<-unlist(apply(id1_mat,1,function(x){
    return(strsplit(x[14],"_")[[1]][3])
  }))
  id1_mat$n4<-unlist(apply(id1_mat,1,function(x){
    return(strsplit(x[14],"_")[[1]][4])
  }))
  id1_mat$n5<-unlist(apply(id1_mat,1,function(x){
    return(strsplit(x[14],"_")[[1]][5])
  }))
  id1_mat$n6<-unlist(apply(id1_mat,1,function(x){
    return(strsplit(x[14],"_")[[1]][6])
  }))
  id2_mat<-do.call(rbind,strsplit(id2,"_"))
  id2_mat<-cbind(id2,id2_mat)
  id2_mat<-as.data.frame(id2_mat)
  if (!('MEX' %in% id2_mat[,3])) {
    id2_mat$loc7<-NA
    id2_mat$loc8<-NA
    if (!('ES' %in% id2_mat[,3])){
      id2_mat$loc6<-NA
      if (unique(id2_mat[,3]) == "IR"){
        id2_mat$loc5<-NA
      }
    }
  }
  colnames(id2_mat)<-c("id2","symbol","type","chr","chain","loc1","loc2",
                       "loc3","loc4","loc5","loc6","loc7","loc8")
  id2_mat$id_new<-unlist(apply(id2_mat,1,function(x){
    # browser()
    if (x[3]=="MEX"){
      a<-tryCatch(tj(del(x[6:13]),"_"),warning = function(w){
        message(paste0(x[1],"_",w))
        return(NA)
      })

    }else{
      if (x[3]=="IR"){
        tmp<-x[6:11][!is.na(x[6:11])]
        a<-tryCatch(tj(del(tmp),"_"),warning = function(w){
          # browser()
          message(paste0(x[1],"_",w))
          return(NA)
        })
      }else{
        a<-tryCatch(tj(del(x[6:11]),"_"),warning = function(w){
          message(paste0(x[1],"_",w))
          return(NA)
        })
      }
    }
    return(a)
  }))
  id2_mat$n1<-unlist(apply(id2_mat,1,function(x){
    return(strsplit(x[14],"_")[[1]][1])
  }))
  id2_mat$n2<-unlist(apply(id2_mat,1,function(x){
    return(strsplit(x[14],"_")[[1]][2])
  }))
  id2_mat$n3<-unlist(apply(id2_mat,1,function(x){
    return(strsplit(x[14],"_")[[1]][3])
  }))
  id2_mat$n4<-unlist(apply(id2_mat,1,function(x){
    return(strsplit(x[14],"_")[[1]][4])
  }))
  id2_mat$n5<-unlist(apply(id2_mat,1,function(x){
    return(strsplit(x[14],"_")[[1]][5])
  }))
  id2_mat$n6<-unlist(apply(id2_mat,1,function(x){
    return(strsplit(x[14],"_")[[1]][6])
  }))

  # id1_mat<-id1_mat[which(id1_mat$symbol %in% id2_mat$symbol),]
  id_mat<-matrix(NA,1,1)
  event_type<-unique(id1_mat$type)
  if ("ES" %in% event_type){
    id1_mat_SE<-id1_mat[which(id1_mat$type=="ES"),]
    id2_mat_SE<-id2_mat[which(id2_mat$type=="ES"),]
    for (chr in unique(id1_mat_SE$chr)){
      id1_mat_SE_part<-id1_mat_SE[which(id1_mat_SE$chr == chr),c("id1","n1","n2","n3","n4")]
      id2_mat_SE_part<-id2_mat_SE[which(id2_mat_SE$chr == chr),c("id2","n1","n2","n3","n4")]
      id1_mat_SE_part$trans<-id1_mat_SE_part$id1
      id1_mat_SE_part$label<-"No"
      if (nrow(id1_mat_SE_part)>0 && (nrow(id2_mat_SE_part)>0)){
        id1_n1<-t(sapply(as.numeric(id1_mat_SE_part[,2]),rep,nrow(id2_mat_SE_part)))
        id1_n2<-t(sapply(as.numeric(id1_mat_SE_part[,3]),rep,nrow(id2_mat_SE_part)))
        id1_n3<-t(sapply(as.numeric(id1_mat_SE_part[,4]),rep,nrow(id2_mat_SE_part)))
        id1_n4<-t(sapply(as.numeric(id1_mat_SE_part[,5]),rep,nrow(id2_mat_SE_part)))
        id2_n1<-sapply(as.numeric(id2_mat_SE_part[,2]),rep,nrow(id1_mat_SE_part))
        id2_n2<-sapply(as.numeric(id2_mat_SE_part[,3]),rep,nrow(id1_mat_SE_part))
        id2_n3<-sapply(as.numeric(id2_mat_SE_part[,4]),rep,nrow(id1_mat_SE_part))
        id2_n4<-sapply(as.numeric(id2_mat_SE_part[,5]),rep,nrow(id1_mat_SE_part))
        n1<-abs(id1_n1-id2_n1)
        n2<-abs(id1_n2-id2_n2)
        n3<-abs(id1_n3-id2_n3)
        n4<-abs(id1_n4-id2_n4)
        n_all<-cbind(cbind(n1,n2),cbind(n3,n4))
        n_all_loc<-unlist(apply(n_all, 1, function(x){
          n11<-intersect(which(x[1:(length(x)/4)]<=err_len),which(x[(length(x)/4+1):(length(x)/2)]<=err_len))
          n12<-intersect(which(x[(length(x)/2+1):(length(x)*3/4)]<=err_len),which(x[(length(x)*3/4+1):(length(x))]<=err_len))
          n13<-intersect(n11,n12)
          if (length(n13)==0) {
            return(NA)
          }else{
            return(n13[1])
          }}))
        id1_mat_SE_part$trans[which(!is.na(n_all_loc))]<-id2_mat_SE_part$id2[n_all_loc[which(!is.na(n_all_loc))]]
        id1_mat_SE_part$label[which(!is.na(n_all_loc))]<-"Yes"
      }
      if (is.na(id_mat[1,1])) {
        id_mat<-id1_mat_SE_part[,c("id1","trans","label")]
      }else{
        use<-id1_mat_SE_part[,c("id1","trans","label")]
        id_mat<-rbind(id_mat,use)
      }
      # cat(chr,"\n")
    }
  }
  if ("IR" %in% event_type){
    id1_mat_IR<-id1_mat[which(id1_mat$type=="IR"),]
    id2_mat_IR<-id2_mat[which(id2_mat$type=="IR"),]
    for (chr in unique(id1_mat_IR$chr)){
      id1_mat_IR_part<-id1_mat_IR[which(id1_mat_IR$chr == chr),c("id1","n1","n2")]
      id2_mat_IR_part<-id2_mat_IR[which(id2_mat_IR$chr == chr),c("id2","n1","n2")]
      id1_mat_IR_part$trans<-id1_mat_IR_part$id1
      id1_mat_IR_part$label<-"No"
      if (nrow(id1_mat_IR_part)>0 && (nrow(id2_mat_IR_part)>0)){
        id1_n1<-t(sapply(as.numeric(id1_mat_IR_part[,2]),rep,nrow(id2_mat_IR_part)))
        id1_n2<-t(sapply(as.numeric(id1_mat_IR_part[,3]),rep,nrow(id2_mat_IR_part)))
        id2_n1<-sapply(as.numeric(id2_mat_IR_part[,2]),rep,nrow(id1_mat_IR_part))
        id2_n2<-sapply(as.numeric(id2_mat_IR_part[,3]),rep,nrow(id1_mat_IR_part))
        n1<-abs(id1_n1-id2_n1)
        n2<-abs(id1_n2-id2_n2)
        n_all<-cbind(n1,n2)
        n_all_loc<-unlist(apply(n_all, 1, function(x){
          n11<-intersect(which(x[1:(length(x)/2)]<=err_len),which(x[(length(x)/2+1):(length(x))]<=err_len))
          if (length(n11)==0) {
            return(NA)
          }else{
            return(n11[1])
          }}))
        id1_mat_IR_part$trans[which(!is.na(n_all_loc))]<-id2_mat_IR_part$id2[n_all_loc[which(!is.na(n_all_loc))]]
        id1_mat_IR_part$label[which(!is.na(n_all_loc))]<-"Yes"
      }
      if (is.na(id_mat[1,1])) {
        id_mat<-id1_mat_IR_part[,c("id1","trans","label")]
      }else{
        use<-id1_mat_IR_part[,c("id1","trans","label")]
        id_mat<-rbind(id_mat,use)
      }
      # cat(chr,"\n")
    }
  }
  if ("A3SS" %in% event_type){
    id1_mat_A3SS<-id1_mat[which(id1_mat$type=="A3SS"),]
    id2_mat_A3SS<-id2_mat[which(id2_mat$type=="A3SS"),]
    for (chr in unique(id1_mat_A3SS$chr)){
      id1_mat_A3SS_part<-id1_mat_A3SS[which(id1_mat_A3SS$chr == chr),c("id1","n1","n2","n3")]
      id2_mat_A3SS_part<-id2_mat_A3SS[which(id2_mat_A3SS$chr == chr),c("id2","n1","n2","n3")]
      id1_mat_A3SS_part$trans<-id1_mat_A3SS_part$id1
      id1_mat_A3SS_part$label<-"No"
      if (nrow(id1_mat_A3SS_part)>0 && (nrow(id2_mat_A3SS_part)>0)){
        id1_n1<-t(sapply(as.numeric(id1_mat_A3SS_part[,2]),rep,nrow(id2_mat_A3SS_part)))
        id1_n2<-t(sapply(as.numeric(id1_mat_A3SS_part[,3]),rep,nrow(id2_mat_A3SS_part)))
        id1_n3<-t(sapply(as.numeric(id1_mat_A3SS_part[,4]),rep,nrow(id2_mat_A3SS_part)))
        id2_n1<-sapply(as.numeric(id2_mat_A3SS_part[,2]),rep,nrow(id1_mat_A3SS_part))
        id2_n2<-sapply(as.numeric(id2_mat_A3SS_part[,3]),rep,nrow(id1_mat_A3SS_part))
        id2_n3<-sapply(as.numeric(id2_mat_A3SS_part[,4]),rep,nrow(id1_mat_A3SS_part))
        n1<-abs(id1_n1-id2_n1)
        n2<-abs(id1_n2-id2_n2)
        n3<-abs(id1_n3-id2_n3)
        n_all<-cbind(cbind(n1,n2),n3)
        n_all_loc<-unlist(apply(n_all, 1, function(x){
          n11<-intersect(which(x[1:(length(x)/3)]<=err_len),which(x[(length(x)/3+1):(length(x)*2/3)]<=err_len))
          n12<-intersect(n11,which(x[(length(x)*2/3+1):(length(x))]<=err_len))
          if (length(n12)==0) {
            return(NA)
          }else{
            return(n12[1])
          }}))
        id1_mat_A3SS_part$trans[which(!is.na(n_all_loc))]<-id2_mat_A3SS_part$id2[n_all_loc[which(!is.na(n_all_loc))]]
        id1_mat_A3SS_part$label[which(!is.na(n_all_loc))]<-"Yes"
      }
      if (is.na(id_mat[1,1])) {
        id_mat<-id1_mat_A3SS_part[,c("id1","trans","label")]
      }else{
        use<-id1_mat_A3SS_part[,c("id1","trans","label")]
        id_mat<-rbind(id_mat,use)
      }
      # cat(chr,"\n")
    }
  }
  if ("A5SS" %in% event_type){
    id1_mat_A5SS<-id1_mat[which(id1_mat$type=="A5SS"),]
    id2_mat_A5SS<-id2_mat[which(id2_mat$type=="A5SS"),]
    for (chr in unique(id1_mat_A5SS$chr)){
      id1_mat_A5SS_part<-id1_mat_A5SS[which(id1_mat_A5SS$chr == chr),c("id1","n1","n2","n3")]
      id2_mat_A5SS_part<-id2_mat_A5SS[which(id2_mat_A5SS$chr == chr),c("id2","n1","n2","n3")]
      id1_mat_A5SS_part$trans<-id1_mat_A5SS_part$id1
      id1_mat_A5SS_part$label<-"No"
      if (nrow(id1_mat_A5SS_part)>0 && (nrow(id2_mat_A5SS_part)>0)){
        id1_n1<-t(sapply(as.numeric(id1_mat_A5SS_part[,2]),rep,nrow(id2_mat_A5SS_part)))
        id1_n2<-t(sapply(as.numeric(id1_mat_A5SS_part[,3]),rep,nrow(id2_mat_A5SS_part)))
        id1_n3<-t(sapply(as.numeric(id1_mat_A5SS_part[,4]),rep,nrow(id2_mat_A5SS_part)))
        id2_n1<-sapply(as.numeric(id2_mat_A5SS_part[,2]),rep,nrow(id1_mat_A5SS_part))
        id2_n2<-sapply(as.numeric(id2_mat_A5SS_part[,3]),rep,nrow(id1_mat_A5SS_part))
        id2_n3<-sapply(as.numeric(id2_mat_A5SS_part[,4]),rep,nrow(id1_mat_A5SS_part))
        n1<-abs(id1_n1-id2_n1)
        n2<-abs(id1_n2-id2_n2)
        n3<-abs(id1_n3-id2_n3)
        n_all<-cbind(cbind(n1,n2),n3)
        n_all_loc<-unlist(apply(n_all, 1, function(x){
          n11<-intersect(which(x[1:(length(x)/3)]<=err_len),which(x[(length(x)/3+1):(length(x)*2/3)]<=err_len))
          n12<-intersect(n11,which(x[(length(x)*2/3+1):(length(x))]<=err_len))
          if (length(n12)==0) {
            return(NA)
          }else{
            return(n12[1])
          }}))
        id1_mat_A5SS_part$trans[which(!is.na(n_all_loc))]<-id2_mat_A5SS_part$id2[n_all_loc[which(!is.na(n_all_loc))]]
        id1_mat_A5SS_part$label[which(!is.na(n_all_loc))]<-"Yes"
      }
      if (is.na(id_mat[1,1])) {
        id_mat<-id1_mat_A5SS_part[,c("id1","trans","label")]
      }else{
        use<-id1_mat_A5SS_part[,c("id1","trans","label")]
        id_mat<-rbind(id_mat,use)
      }
      # cat(chr,"\n")
    }
  }
  if ("MEX" %in% event_type){
    id1_mat_MEX<-id1_mat[which(id1_mat$type=="MEX"),]
    id2_mat_MEX<-id2_mat[which(id2_mat$type=="MEX"),]
    for (chr in unique(id1_mat_MEX$chr)){
      id1_mat_MEX_part<-id1_mat_MEX[which(id1_mat_MEX$chr == chr),c("id1","n1","n2","n3","n4","n5","n6")]
      id2_mat_MEX_part<-id2_mat_MEX[which(id2_mat_MEX$chr == chr),c("id2","n1","n2","n3","n4","n5","n6")]
      id1_mat_MEX_part$trans<-id1_mat_MEX_part$id1
      id1_mat_MEX_part$label<-"No"
      if (nrow(id1_mat_MEX_part)>0 && (nrow(id2_mat_MEX_part)>0)){
        id1_n1<-t(sapply(as.numeric(id1_mat_MEX_part[,2]),rep,nrow(id2_mat_MEX_part)))
        id1_n2<-t(sapply(as.numeric(id1_mat_MEX_part[,3]),rep,nrow(id2_mat_MEX_part)))
        id1_n3<-t(sapply(as.numeric(id1_mat_MEX_part[,4]),rep,nrow(id2_mat_MEX_part)))
        id1_n4<-t(sapply(as.numeric(id1_mat_MEX_part[,5]),rep,nrow(id2_mat_MEX_part)))
        id1_n5<-t(sapply(as.numeric(id1_mat_MEX_part[,6]),rep,nrow(id2_mat_MEX_part)))
        id1_n6<-t(sapply(as.numeric(id1_mat_MEX_part[,7]),rep,nrow(id2_mat_MEX_part)))
        id2_n1<-sapply(as.numeric(id2_mat_MEX_part[,2]),rep,nrow(id1_mat_MEX_part))
        id2_n2<-sapply(as.numeric(id2_mat_MEX_part[,3]),rep,nrow(id1_mat_MEX_part))
        id2_n3<-sapply(as.numeric(id2_mat_MEX_part[,4]),rep,nrow(id1_mat_MEX_part))
        id2_n4<-sapply(as.numeric(id2_mat_MEX_part[,5]),rep,nrow(id1_mat_MEX_part))
        id2_n5<-t(sapply(as.numeric(id2_mat_MEX_part[,6]),rep,nrow(id1_mat_MEX_part)))
        id2_n6<-t(sapply(as.numeric(id2_mat_MEX_part[,7]),rep,nrow(id1_mat_MEX_part)))
        n1<-abs(id1_n1-id2_n1)
        n2<-abs(id1_n2-id2_n2)
        n3<-abs(id1_n3-id2_n3)
        n4<-abs(id1_n4-id2_n4)
        n5<-abs(id1_n5-id2_n5)
        n6<-abs(id1_n6-id2_n6)
        n_all<-cbind(cbind(cbind(n1,n2),cbind(n3,n4)),cbind(n5,n6))
        n_all_loc<-unlist(apply(n_all, 1, function(x){
          n11<-intersect(which(x[1:(length(x)/6)]<=err_len),which(x[(length(x)/6+1):(length(x)/3)]<=err_len))
          n12<-intersect(which(x[(length(x)/3+1):(length(x)/3)]<=err_len),which(x[(length(x)/2+1):(length(x)*2/3)]<=err_len))
          n13<-intersect(which(x[(length(x)*2/3+1):(length(x)*5/6)]<=err_len),which(x[(length(x)*5/6+1):(length(x))]<=err_len))
          n14<-intersect(intersect(n11,n12),n13)
          if (length(n14)==0) {
            return(NA)
          }else{
            return(n14[1])
          }}))
        id1_mat_MEX_part$trans[which(!is.na(n_all_loc))]<-id2_mat_MEX_part$id2[n_all_loc[which(!is.na(n_all_loc))]]
        id1_mat_MEX_part$label[which(!is.na(n_all_loc))]<-"Yes"
      }
      if (is.na(id_mat[1,1])) {
        id_mat<-id1_mat_MEX_part[,c("id1","trans","label")]
      }else{
        use<-id1_mat_MEX_part[,c("id1","trans","label")]
        id_mat<-rbind(id_mat,use)
      }
      # cat(chr,"\n")
    }
  }
  rownames(id_mat)<-id_mat[,1]
  id_mat<-id_mat[id1_mat[,1],]
  countOccurrences <- function(x) {
    tab <- table(x)
    occurrences <- tab[x]
    return(occurrences)
  }
  # 使用ave()函数将计算结果应用于第四列
  id_mat$count <- ave(id_mat[, 2], id_mat[, 2], FUN = countOccurrences)

  return(id_mat)
}

