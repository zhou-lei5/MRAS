#' data_integration_preprocess_rmats
#'
#' @param path_rmats path
#' @param type type
#' @param design design
#' @param sample_num number
#' @param psi_cutoff psi_cutoff
#' @param reads_cutoff reads_cutoff
#' @param reads_num reads_num
#' @param colname_ann colname_ann
#' @param gsub_pattern_up gsub_pattern_up
#' @param gsub_pattern_down gsub_pattern_down
#'
#' @return result
#' @export
#'
data_integration_preprocess_rmats<-function(path_rmats,type="default",design=NULL,
                                            sample_num,psi_cutoff,
                                            reads_cutoff = NULL,reads_num = 10,colname_ann = FALSE,
                                            gsub_pattern_up = NULL, gsub_pattern_down = NULL){
  data<-data_integration_rmats(path_rmats = path_rmats,type = type,design = design)
  psi_pair<-data_preprocess_rmats(data = data,sample_num = sample_num,psi_cutoff = psi_cutoff,
                                  reads_cutoff = reads_cutoff,reads_num = reads_num,colname_ann = colname_ann,path_rmats = path_rmats,
                                  gsub_pattern_up = gsub_pattern_up,gsub_pattern_down = gsub_pattern_down)
  return(psi_pair)
}


#' data_integration_rmats
#'
#' @param path_rmats path_rmats
#' @param type type
#' @param design design
#'
#' @return result
#' @export
#'

data_integration_rmats<-function(path_rmats,type="default",design=NULL){
  if (is.null(design)){
    type<-match.arg(type,choices = c("default","all"))
    if (type == "default"){
      SE<-read.delim(paste0(path_rmats,"/SE.MATS.JCEC.txt"),header = T)
      RI<-read.delim(paste0(path_rmats,"/RI.MATS.JCEC.txt"),header = T)
      A3SS<-read.delim(paste0(path_rmats,"/A3SS.MATS.JCEC.txt"),header = T)
      A5SS<-read.delim(paste0(path_rmats,"/A5SS.MATS.JCEC.txt"),header = T)

      SE<-SE[,c(1:11,17,18,13,14,21)]
      RI<-RI[,c(1:11,17,18,13,14,21)]
      A3SS<-A3SS[,c(1:11,17,18,13,14,21)]
      A5SS<-A5SS[,c(1:11,17,18,13,14,21)]

      SE<-as.data.frame(t(apply(SE, 1, function(x){
        x[6]<-as.numeric(x[6])+1
        x[8]<-as.numeric(x[8])+1
        x[10]<-as.numeric(x[10])+1
        return(x)
      })))
      RI<-as.data.frame(t(apply(RI, 1, function(x){
        x[6]<-as.numeric(x[6])+1
        x[8]<-as.numeric(x[8])+1
        x[10]<-as.numeric(x[10])+1
        return(x)
      })))
      A3SS<-as.data.frame(t(apply(A3SS, 1, function(x){
        x[6]<-as.numeric(x[6])+1
        x[8]<-as.numeric(x[8])+1
        x[10]<-as.numeric(x[10])+1
        return(x)
      })))
      A5SS<-as.data.frame(t(apply(A5SS, 1, function(x){
        x[6]<-as.numeric(x[6])+1
        x[8]<-as.numeric(x[8])+1
        x[10]<-as.numeric(x[10])+1
        return(x)
      })))

      SE[,1]<-paste0(SE[,3],"_","ES","_",SE[,4],"_",SE[,5],"_",SE[,6],"_",SE[,7],"_",SE[,8],"_",SE[,9],"_",SE[,10],"_",SE[,11])
      RI[,1]<-paste0(RI[,3],"_","IR","_",RI[,4],"_",RI[,5],"_",RI[,6],"_",RI[,7],"_",RI[,8],"_",RI[,9],"_",RI[,10],"_",RI[,11])
      A3SS[,1]<-paste0(A3SS[,3],"_","A3SS","_",A3SS[,4],"_",A3SS[,5],"_",A3SS[,6],"_",A3SS[,7],"_",A3SS[,8],"_",A3SS[,9],"_",A3SS[,10],"_",A3SS[,11])
      A5SS[,1]<-paste0(A5SS[,3],"_","A5SS","_",A5SS[,4],"_",A5SS[,5],"_",A5SS[,6],"_",A5SS[,7],"_",A5SS[,8],"_",A5SS[,9],"_",A5SS[,10],"_",A5SS[,11])

      write.table(SE[,c(1,12:16)],file = paste0(path_rmats,"/PSI.txt"),row.names = F,col.names = F,sep = ",",quote = F)
      write.table(RI[,c(1,12:16)],file = paste0(path_rmats,"/PSI.txt"),row.names = F,col.names = F,sep = ",",quote = F,append = T)
      write.table(A3SS[,c(1,12:16)],file = paste0(path_rmats,"/PSI.txt"),row.names = F,col.names = F,sep = ",",quote = F,append = T)
      write.table(A5SS[,c(1,12:16)],file = paste0(path_rmats,"/PSI.txt"),row.names = F,col.names = F,sep = ",",quote = F,append = T)
       }else{
        SE<-read.delim(paste0(path_rmats,"/SE.MATS.JCEC.txt"),header = T)
        RI<-read.delim(paste0(path_rmats,"/RI.MATS.JCEC.txt"),header = T)
        A3SS<-read.delim(paste0(path_rmats,"/A3SS.MATS.JCEC.txt"),header = T)
        A5SS<-read.delim(paste0(path_rmats,"/A5SS.MATS.JCEC.txt"),header = T)
        MXE<-read.delim(paste0(path_rmats,"/MXE.MATS.JCEC.txt"),header = T)

        SE<-SE[,c(1:11,17,18,13,14,21)]
        RI<-RI[,c(1:11,17,18,13,14,21)]
        A3SS<-A3SS[,c(1:11,17,18,13,14,21)]
        A5SS<-A5SS[,c(1:11,17,18,13,14,21)]
        MXE<-MXE[,c(1:13,19:20,15:16,23)]

        SE<-as.data.frame(t(apply(SE, 1, function(x){
          x[6]<-as.numeric(x[6])+1
          x[8]<-as.numeric(x[8])+1
          x[10]<-as.numeric(x[10])+1
          return(x)
        })))
        RI<-as.data.frame(t(apply(RI, 1, function(x){
          x[6]<-as.numeric(x[6])+1
          x[8]<-as.numeric(x[8])+1
          x[10]<-as.numeric(x[10])+1
          return(x)
        })))
        A3SS<-as.data.frame(t(apply(A3SS, 1, function(x){
          x[6]<-as.numeric(x[6])+1
          x[8]<-as.numeric(x[8])+1
          x[10]<-as.numeric(x[10])+1
          return(x)
        })))
        A5SS<-as.data.frame(t(apply(A5SS, 1, function(x){
          x[6]<-as.numeric(x[6])+1
          x[8]<-as.numeric(x[8])+1
          x[10]<-as.numeric(x[10])+1
          return(x)
        })))
        MXE<-as.data.frame(t(apply(MXE, 1, function(x){
            x[6]<-as.numeric(x[6])+1
            x[8]<-as.numeric(x[8])+1
            x[10]<-as.numeric(x[10])+1
            x[12]<-as.numeric(x[12])+1
          return(x)
        })))

        SE[,1]<-paste0(SE[,3],"_","ES","_",SE[,4],"_",SE[,5],"_",SE[,6],"_",SE[,7],"_",SE[,8],"_",SE[,9],"_",SE[,10],"_",SE[,11])
        RI[,1]<-paste0(RI[,3],"_","IR","_",RI[,4],"_",RI[,5],"_",RI[,6],"_",RI[,7],"_",RI[,8],"_",RI[,9],"_",RI[,10],"_",RI[,11])
        A3SS[,1]<-paste0(A3SS[,3],"_","A3SS","_",A3SS[,4],"_",A3SS[,5],"_",A3SS[,6],"_",A3SS[,7],"_",A3SS[,8],"_",A3SS[,9],"_",A3SS[,10],"_",A3SS[,11])
        A5SS[,1]<-paste0(A5SS[,3],"_","A5SS","_",A5SS[,4],"_",A5SS[,5],"_",A5SS[,6],"_",A5SS[,7],"_",A5SS[,8],"_",A5SS[,9],"_",A5SS[,10],"_",A5SS[,11])
        MXE[,1]<-paste0(MXE[,3],"_","MEX","_",MXE[,4],"_",MXE[,5],"_",MXE[,10],"_",MXE[,11],"_",MXE[,6],"_",MXE[,7],"_",MXE[,8],"_",MXE[,9],"_",MXE[,12],"_",MXE[,13])

        write.table(SE[,c(1,12:16)],file = paste0(path_rmats,"/PSI.txt"),row.names = F,col.names = F,sep = ",",quote = F)
        write.table(RI[,c(1,12:16)],file = paste0(path_rmats,"/PSI.txt"),row.names = F,col.names = F,sep = ",",quote = F,append = T)
        write.table(A3SS[,c(1,12:16)],file = paste0(path_rmats,"/PSI.txt"),row.names = F,col.names = F,sep = ",",quote = F,append = T)
        write.table(A5SS[,c(1,12:16)],file = paste0(path_rmats,"/PSI.txt"),row.names = F,col.names = F,sep = ",",quote = F,append = T)
        write.table(MXE[,c(1,14:18)],file = paste0(path_rmats,"/PSI.txt"),row.names = F,col.names = F,sep = ",",quote = F,append = T)
      }
  }else{
    if ("SE" %in% design){
      SE<-read.delim(paste0(path_rmats,"/SE.MATS.JCEC.txt"),header = T)
      SE<-SE[,c(1:11,17,18,13,14,21)]
      SE<-as.data.frame(t(apply(SE, 1, function(x){
        x[6]<-as.numeric(x[6])+1
        x[8]<-as.numeric(x[8])+1
        x[10]<-as.numeric(x[10])+1
        return(x)
      })))

      SE[,1]<-paste0(SE[,3],"_","ES","_",SE[,4],"_",SE[,5],"_",SE[,6],"_",SE[,7],"_",SE[,8],"_",SE[,9],"_",SE[,10],"_",SE[,11])
      write.table(SE[,c(1,12:16)],file = paste0(path_rmats,"/PSI.txt"),
                  row.names = F,col.names = F,
                  sep = ",",quote = F,append = T)
     }else{
      if ("RI" %in% design){
        RI<-read.delim(paste0(path_rmats,"/RI.MATS.JCEC.txt"),header = T)
        RI<-RI[,c(1:11,17,18,13,14,21)]
        RI<-as.data.frame(t(apply(RI, 1, function(x){
          x[6]<-as.numeric(x[6])+1
          x[8]<-as.numeric(x[8])+1
          x[10]<-as.numeric(x[10])+1
          return(x)
        })))

        RI[,1]<-paste0(RI[,3],"_","IR","_",RI[,4],"_",RI[,5],"_",RI[,6],"_",RI[,7],"_",RI[,8],"_",RI[,9],"_",RI[,10],"_",RI[,11])
        write.table(RI[,c(1,12:16)],file = paste0(path_rmats,"/PSI.txt"),
                    row.names = F,col.names = F,
                    sep = ",",quote = F,append = T)
      }else{
        if ("A3SS" %in% design){
          A3SS<-read.delim(paste0(path_rmats,"/A3SS.MATS.JCEC.txt"),header = T)
          A3SS<-A3SS[,c(1:11,17,18,13,14,21)]
          A3SS<-as.data.frame(t(apply(A3SS, 1, function(x){
            x[6]<-as.numeric(x[6])+1
            x[8]<-as.numeric(x[8])+1
            x[10]<-as.numeric(x[10])+1
            return(x)
          })))

          A3SS[,1]<-paste0(A3SS[,3],"_","A3SS","_",A3SS[,4],"_",A3SS[,5],"_",A3SS[,6],"_",A3SS[,7],"_",A3SS[,8],"_",A3SS[,9],"_",A3SS[,10],"_",A3SS[,11])
          write.table(A3SS[,c(1,12:16)],file = paste0(path_rmats,"/PSI.txt"),
                      row.names = F,col.names = F,
                      sep = ",",quote = F,append = T)
        }else{
          if ("A5SS" %in% design){
            A5SS<-read.delim(paste0(path_rmats,"/A5SS.MATS.JCEC.txt"),header = T)
            A5SS<-A5SS[,c(1:11,17,18,13,14,21)]
            A5SS<-as.data.frame(t(apply(A5SS, 1, function(x){
              x[6]<-as.numeric(x[6])+1
              x[8]<-as.numeric(x[8])+1
              x[10]<-as.numeric(x[10])+1
              return(x)
            })))

            A5SS[,1]<-paste0(A5SS[,3],"_","A5SS","_",A5SS[,4],"_",A5SS[,5],"_",A5SS[,6],"_",A5SS[,7],"_",A5SS[,8],"_",A5SS[,9],"_",A5SS[,10],"_",A5SS[,11])
            write.table(A5SS[,c(1,12:16)],file = paste0(path_rmats,"/PSI.txt"),
                        row.names = F,col.names = F,
                        sep = ",",quote = F,append = T)
          }else{
            if ("MXE" %in% design){
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
            write.table(MXE[,c(1,14:18)],file = paste0(path_rmats,"/PSI.txt"),
                        row.names = F,col.names = F,
                        sep = ",",quote = F,append = T)

            }
          }
        }
      }
    }

  }
  data<-fread(paste0(path_rmats,"/PSI.txt"),sep = ",",header = F)
  data[is.na(data)]<-0
  return(data)
}



#' data_preprocess_rmats
#'
#' @param data data
#' @param sample_num sample_num
#' @param psi_cutoff psi_cutoff
#' @param reads_cutoff reads_cutoff
#' @param reads_num reads_num
#' @param colname_ann colname_ann
#' @param path_rmats path_rmats
#' @param gsub_pattern_up gsub_pattern_up
#' @param gsub_pattern_down gsub_pattern_down
#'
#' @return result
#' @export
#'

data_preprocess_rmats<-function(data,sample_num,psi_cutoff,
                                reads_cutoff = NULL,reads_num = 10,colname_ann = FALSE,path_rmats = NULL,
                                gsub_pattern_up = NULL, gsub_pattern_down = NULL){

  pc1<-data[,4:(sample_num+3)]
  pc2<-data[,(sample_num+4):(2*sample_num+3)]
  pc3<-data[,(2*sample_num+4):(3*sample_num+3)]

  f3<-filter_rmats_psi(pc3,psi_cutoff)
  if (is.null(reads_cutoff) | is.null(reads_num)){
    psi_pair<-data[f3,c(1,(2*sample_num+4):(3*sample_num+3))]
    psi_pair<-as.data.frame(psi_pair)
  }else{
    f1<-filter_ramts_reads(pc1,reads_cutoff,reads_num)
    f2<-filter_ramts_reads(pc2,reads_cutoff,reads_num)
    psi_pair<-data[f1 & f2 & f3,c(1,(2*sample_num+4):(3*sample_num+3))]
    psi_pair<-as.data.frame(psi_pair)
  }
  if (colname_ann){
    colname<-read.table(paste0(path_rmats,"/bam_path.txt"))
    colname<-t(apply(colname, 1, function(x){
      if (!(is.null(gsub_pattern_up))) {
        a<-gsub(gsub_pattern_up,"",x[1])
        if (!(is.null(gsub_pattern_down))) {
          b<-gsub(gsub_pattern_down,"",a)
        }
      }else{
        b<-gsub(gsub_pattern_down,"",x[1])
      }
      return(b)
    }))
    colnames(psi_pair)<-colname
    }
  write.table(psi_pair,file = paste0(path_rmats,"/PSI_final.txt"),
                row.names = T,col.names = T,
                quote = F,sep = "\t")
  return(psi_pair)
}


filter_rmats_psi<-function(rmats_matrix,n){
  results<-matrix(TRUE,nrow = nrow(rmats_matrix),ncol = 1)
  results<-apply(rmats_matrix,1,function(x){
    length0<-length(which(as.numeric(x)==0))
    length1<-length(which(as.numeric(x)==1))
    result<-ifelse((length0+length1)>n,FALSE,TRUE)
    return(result)
  })
  return(results)
}
filter_ramts_reads<-function(reads_matrix,n,reads_num=10){
  results<-matrix(TRUE,nrow = nrow(reads_matrix),ncol = 1)
  results<-apply(reads_matrix,1,function(x){
    length10<-length(which(as.numeric(x)<reads_num))
    result<-ifelse((length10)>n,FALSE,TRUE)
    return(result)
  })
  return(results)
}
