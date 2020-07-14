###source function for output preperations

library(stringr)
library(dplyr)
library(RColorBrewer)
library(pals)
library(tidyr)
library(qpcR)


##prepares output for publication.
output_prep<- function(filepath, output_name){
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9, 17, 21, 22)] 
  
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  b<-c()
  for (state in levels(a$Protein.State)){##
    temp1<-a[which(a$Protein.State ==state ),] ##creates temporary df, temp1, with Protein states going through all unique protein states
    for (time in unique(temp1$Deut.Time)){
      temp2<-temp1[which(temp1$Deut.Time ==time),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(n_tmp[1:2],paste(n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , paste(n_tmp[9],"_",nb,sep=""), paste(n_tmp[10],"_",nb,sep=""),
               paste(n_tmp[11],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                                     'Charge')), df_List)}
      b=rbind(b, bp)}} ## b has all information bound together again to have all information df
  
  ###order columns
  #### names(b)[grep("Exp.Cent", colnames(b))] gives position of all variable with pattern and return position
  ord1<-c('Protein.State', 'Deut.Time', 'Start', 'End', 'Sequence', 'Search.RT', 'Charge',
          names(b)[grep("Experiment", colnames(b))],
          names(b)[grep("Exp.Cent", colnames(b))], names(b)[grep("X..Deut", colnames(b))], names(b)[grep("Deut.._", colnames(b))] )
  b<-b[,ord1]
  ###calculate means +sd and write to out dataframe. Later assign names, choose only important
  out<-data.frame(b[,1:8], rowMeans(b[,grep("Exp.Cent", colnames(b))]), apply(b[,grep("Exp.Cent", colnames(b))],1,sd), 
                  rowMeans(b[,grep("X..Deut", colnames(b))]),apply(b[,grep("X..Deut", colnames(b))],1,sd))
  
  ###prepare Full deuteration data.frame to be bound with out data.frame. 
  FD2<-data.frame(FD[,c(1,2,4:6,8,7,3,9)], rep(0,times=dim(FD)[1]), FD[,10], rep(0,times=dim(FD)[1]))
  undeut2<-data.frame(undeut[,c(1,2,4:6,8,7,3,9)], rep(0,times=dim(undeut)[1]),  rep(0,times=dim(undeut)[1]), rep(0,times=dim(undeut)[1]))
  
  nm_final<-c(names(out)[1:5], "Retention_Time_[min]","Charge","Experiment",
              "Mean.Peptide_Mass_[Da]", "st.dev_Peptide_Mass",
              "Mean.Deut.Uptake_[Da]", "st.dev_Deut_Uptake")
  
  colnames(FD2)<-nm_final
  colnames(undeut2)<-nm_final
  colnames(out)<-nm_final
  
  
  out<-rbind(out, FD2)
  out<-rbind(out, undeut2)
  out<-out[,c(1:6,9,11,12)]
  out<-(arrange(out, Start, End, Protein.State))
  ##find a way to order correctly assign negative time to zero and
  ##very long time for FD and order based of these columns
  ##if bp does not work use out
  ord=out$Deut.Time
  ord=type.convert(ord, as.is = TRUE)
  ord[ord=="0s"]<- c("-5.0")
  ord[ord== "FD"]<- c("10000000000000.0s")
  ord=as.numeric(str_sub(ord, end=-2))
  
  bp<-data.frame(out, ord)
  bp<-arrange(bp, Start, End, Protein.State, ord)
  bp<-bp[,-dim(bp)[2]]
  ###write output
  write.csv(bp, output_name)
  return()}

####################
##prepares output for timepoints analysis. It gives raw uptake data in columns. 
output_tp<- function(filepath){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (time in unique(a$Deut.Time)){
    temp1<-a[which(a$Deut.Time ==time),]
    st_l<-c()
    nbs=0
    for (state in unique(temp1$Protein.State)){##
      temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      st<-gsub(c(' '),'',state)
      df_nm_st<-paste(st, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , 
               paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Deut.Time, Start, End, Charge)
  ord=b$Deut.Time
  ord=as.numeric(str_sub(ord, end=-2))
  
  bp1<-data.frame(b, ord)
  bp1<-arrange(bp1, ord, Start, End)
  b<-bp1[,-dim(bp1)[2]]
  
  uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(b[,1:6], b[,grep("Deut.._", colnames(b))])
  return(uptake.tp)}

##########
##takes raw data from HDXexaminer and organize data into columns with % deuteration
output_tp_proc<- function(filepath){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (time in unique(a$Deut.Time)){
    temp1<-a[which(a$Deut.Time ==time),]
    st_l<-c()
    nbs=0
    for (state in unique(temp1$Protein.State)){##
      temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      st<-gsub(c(' '),'',state)
      df_nm_st<-paste(st, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , 
               paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
    
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Deut.Time, Start, End, Charge)
  ord=b$Deut.Time
  ord=as.numeric(str_sub(ord, end=-2))
  
  bp1<-data.frame(b, ord)
  bp1<-arrange(bp1, ord, Start, End)
  b<-bp1[,-dim(bp1)[2]]
  uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(b[,1:6], b[,grep("Deut.._", colnames(b))])
  return(procent.tp)}


############
####

initial_analysis_timepoint<-function(df, nb_timepoint) {
  ### calculate standard deviation of sample 1 for all data points
  nb_sets=(dim(df)[2]-6)/nb_timepoint
  sd1<-c(); for ( j in 1:nb_timepoint) {
    for (i in 1:dim(df)[1]) {sd1<-c(sd1, sd(df[i,(6+(j-1)*nb_sets+1):(6+(j-1)*nb_sets+3)]))} 
  }
  #print(sd1)
  sd2<-data.frame(matrix(sd1, ncol =nb_timepoint , byrow = FALSE))
  names(sd2)<-c()
  return(sd2)
}
####


###calculates standard devation for timepoint
sd_timepoint<-function(df, replicates=3) {
  ### calculate standard deviation of sample 1 for all data points
  nb_sets=(dim(df)[2]-6)/replicates
  nm_root<-colnames(df[,c(6+((1:nb_sets)-1)*replicates+1)])
  nm_root<-str_sub(nm_root, end = -3)
  sd_nm<-paste("sd_", nm_root, sep="")
  
  #print(nb_sets)
  sd1<-c(); for ( j in 1:nb_sets) {
    for (i in 1:dim(df)[1]) {sd1<-c(sd1, sd(df[i,(6+(j-1)*replicates+1):(6+(j-1)*replicates+replicates)]))} 
  }
  #print(sd1)
  sd2<-data.frame(matrix(sd1, ncol=nb_sets , byrow = FALSE))
  colnames(sd2)<-sd_nm
  sd2<-data.frame(df[,1:6], sd2)
  return(sd2)
} 
###average 
ave_timepoint<-function(df, replicates=3) {
  ### calculate standard deviation of sample 1 for all data points
  nb_sets=(dim(df)[2]-6)/replicates
  nm_root<-colnames(df[,c(6+((1:nb_sets)-1)*replicates+1)])
  nm_root<-str_sub(nm_root, end = -3)
  ave_nm<-paste("av_", nm_root, sep="")
  
  #print(nb_sets)
  ave1<-c(); for ( j in 1:nb_sets) {
    for (i in 1:dim(df)[1]) {ave1<-c(ave1, rowMeans(df[i,(6+(j-1)*replicates+1):(6+(j-1)*replicates+replicates)]))} 
  }
  #print(sd1)
  ave2<-data.frame(matrix(ave1, ncol=nb_sets , byrow = FALSE))
  colnames(ave2)<-ave_nm
  ave2<-data.frame(df[,1:6], ave2)
  return(ave2)
} 

pv_timepoint<-function(df, replicates=3) {
  ### calculate standard deviation of sample 1 for all data points
  nb_sets=(dim(df)[2]-6)/replicates
  nm_root<-colnames(df[,c(6+((1:nb_sets)-1)*replicates+1)])
  nm_root<-str_sub(nm_root, end = -3)
  pv_nm<-paste("pv_", nm_root, sep="")
  
  
  pv1<-c(); for ( j in 1:nb_sets) {
    for (i in 1:dim(df)[1]) {x1<-df[i,7:(7+replicates-1)]; x2<-df[i,(6+(j-1)*replicates+1):(6+(j-1)*replicates+replicates)]
    tt<-c(); tt<-t.test(x1, x2)
    pv1<-c(pv1, tt$p.val)} 
  }
  #print(sd1)
  pv2<-data.frame(matrix(pv1, ncol=nb_sets , byrow = FALSE))
  colnames(pv2)<-pv_nm
  pv2<-data.frame(df[,1:6], pv2)
  return(pv2)
}

dif_ave<-function(df){
  da1<-data.frame(df[,1:6],df[,7:dim(df)[2]]-df[,7])
  return(da1)
}


verbose_timepoint_output<-function(filepath, output_name){
  df<-output_tp(filepath)
  pv1<-pv_timepoint(df)
  s1<-sd_timepoint(df)
  av1<-ave_timepoint(df)
  df_List<-list(av1, s1, pv1) 
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  ord=as.numeric(str_sub(bp$Deut.Time, end=-2))
  bp<-data.frame(bp, ord)
  bp<-arrange(bp, ord, Start, End, Charge)
  bp<-bp[,-dim(bp)[2]]
  write.csv(bp, output_name)
  return(bp)
}

verbose_timepoint_output_proc<-function(filepath, output_name){
  df<-output_tp_proc(filepath)
  pv1<-pv_timepoint(df)
  s1<-sd_timepoint(df)
  av1<-ave_timepoint(df)
  df_List<-list(av1, s1, pv1) 
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  ord=as.numeric(str_sub(bp$Deut.Time, end=-2))
  bp<-data.frame(bp, ord)
  bp<-arrange(bp, ord, Start, End, Charge)
  bp<-bp[,-dim(bp)[2]]
  write.csv(bp, output_name)
  return(bp)
}

######

output_tcourse<- function(filepath){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  
  
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(a$Protein.State)){
    temp1<-a[which(a$Protein.State ==state),]
    st_l<-c()
    nbs=0
    for (time in unique(temp1$Deut.Time)){##
      temp2<-temp1[which(temp1$Deut.Time ==time ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      df_nm_st<-paste(time, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(n_tmp[1],paste("t",time,n_tmp[2],"_",nb,sep=""),paste("t",time,n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , 
               paste("t",time, "_",n_tmp[9],"_",nb,sep=""), paste("t",time, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
    
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Protein.State,Start, End,  Charge)
  uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(b[,1:6], b[,grep("Deut.._", colnames(b))])
  return(uptake.tp)}


output_tcourse_proc<- function(filepath, output_file){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  b<-c()
  
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(a$Protein.State)){
    temp1<-a[which(a$Protein.State ==state),]
    st_l<-c()
    nbs=0
    for (time in unique(temp1$Deut.Time)){##
      temp2<-temp1[which(temp1$Deut.Time ==time ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      df_nm_st<-paste(time, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(n_tmp[1],paste("t",time,n_tmp[2],"_",nb,sep=""),paste("t",time,n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , 
               paste("t",time, "_",n_tmp[9],"_",nb,sep=""), paste("t",time, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
    
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Protein.State, Start, End,  Charge)
  uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(b[,1:6], b[,grep("Deut.._", colnames(b))])
  return(procent.tp)}

output_tcourse_proc_csv<-function(filepath, output){
  write.csv(output_tcourse_proc(filepath), output)}

output_tcourse_csv<-function(filepath, output){
  write.csv(output_tcourse(filepath), output)}

output_tcourse_proc_csv<-function(filepath, output){
  write.csv(output_tcourse_proc(filepath), output)}

output_tp_csv<-function(filepath, output){
  write.csv(output_tp(filepath), output)}

output_tp_proc_csv<-function(filepath, output){
  write.csv(output_tp_proc(filepath), output)}


output_tp_states<- function(filepath, states){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (time in unique(a$Deut.Time)){
    temp1<-a[which(a$Deut.Time ==time),]
    st_l<-c()
    nbs=0
    for (state in states){##
      temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      st<-gsub(c(' '),'',state)
      df_nm_st<-paste(st, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , 
               paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Deut.Time, Start, End, Charge)
  ord=b$Deut.Time
  ord=as.numeric(str_sub(ord, end=-2))
  
  bp1<-data.frame(b, ord)
  bp1<-arrange(bp1, ord, Start, End)
  b<-bp1[,-dim(bp1)[2]]
  
  uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(b[,1:6], b[,grep("Deut.._", colnames(b))])
  return(uptake.tp)}

verbose_timepoint_output_subset<-function(filepath, states, output_name){
  df<-output_tp_states(filepath, states)
  pv1<-pv_timepoint(df)
  s1<-sd_timepoint(df)
  av1<-ave_timepoint(df)
  df_List<-list(av1, s1, pv1) 
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  ord=as.numeric(str_sub(bp$Deut.Time, end=-2))
  bp<-data.frame(bp, ord)
  bp<-arrange(bp, ord, Start, End, Charge)
  bp<-bp[,-dim(bp)[2]]
  write.csv(bp, output_name)
  return(bp)
}

timepoint_output_subset<-function(filepath, states){
  df<-output_tp_states(filepath, states)
  pv1<-pv_timepoint(df)
  s1<-sd_timepoint(df)
  av1<-ave_timepoint(df)
  df_List<-list(av1, s1, pv1) 
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  ord=as.numeric(str_sub(bp$Deut.Time, end=-2))
  bp<-data.frame(bp, ord)
  bp<-arrange(bp, ord, Start, End, Charge)
  bp<-bp[,-dim(bp)[2]]
  return(bp)
}




output_tp_proc_states<- function(filepath, states){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (time in unique(a$Deut.Time)){
    temp1<-a[which(a$Deut.Time ==time),]
    st_l<-c()
    nbs=0
    for (state in states){##
      temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      st<-gsub(c(' '),'',state)
      df_nm_st<-paste(st, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , 
               paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Deut.Time, Start, End, Charge)
  ord=b$Deut.Time
  ord=as.numeric(str_sub(ord, end=-2))
  
  bp1<-data.frame(b, ord)
  bp1<-arrange(bp1, ord, Start, End)
  b<-bp1[,-dim(bp1)[2]]
  
  uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(b[,1:6], b[,grep("Deut.._", colnames(b))])
  return(procent.tp)}

verbose_timepoint_output_proc_subset<-function(filepath, states, output_name){
  df<-output_tp_proc_states(filepath, states)
  pv1<-pv_timepoint(df)
  s1<-sd_timepoint(df)
  av1<-ave_timepoint(df)
  df_List<-list(av1, s1, pv1) 
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  ord=as.numeric(str_sub(bp$Deut.Time, end=-2))
  bp<-data.frame(bp, ord)
  bp<-arrange(bp, ord, Start, End, Charge)
  bp<-bp[,-dim(bp)[2]]
  write.csv(bp, output_name)
  return(bp)
}

timepoint_output_proc_subset<-function(filepath, states){
  df<-output_tp_proc_states(filepath, states)
  pv1<-pv_timepoint(df)
  s1<-sd_timepoint(df)
  av1<-ave_timepoint(df)
  df_List<-list(av1, s1, pv1) 
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  ord=as.numeric(str_sub(bp$Deut.Time, end=-2))
  bp<-data.frame(bp, ord)
  bp<-arrange(bp, ord, Start, End, Charge)
  bp<-bp[,-dim(bp)[2]]
  return(bp)
}



output_tp_noSeq<- function(filepath){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (time in unique(a$Deut.Time)){
    temp1<-a[which(a$Deut.Time ==time),]
    st_l<-c()
    nbs=0
    for (state in unique(temp1$Protein.State)){##
      temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      st<-gsub(c(' '),'',state)
      df_nm_st<-paste(st, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:5] ,
               paste(st,n_tmp[6],"_",nb,sep=""),n_tmp[7:8] ,
               paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Deut.Time, Start, End, Charge)
  ord=b$Deut.Time
  ord=as.numeric(str_sub(ord, end=-2))
  
  bp1<-data.frame(b, ord)
  bp1<-arrange(bp1, ord, Start, End)
  b<-bp1[,-dim(bp1)[2]]
  
  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )
  
  
  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  return(uptake.tp)}

output_tp_noSeq_proc<- function(filepath){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (time in unique(a$Deut.Time)){
    temp1<-a[which(a$Deut.Time ==time),]
    st_l<-c()
    nbs=0
    for (state in unique(temp1$Protein.State)){##
      temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      st<-gsub(c(' '),'',state)
      df_nm_st<-paste(st, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:5] ,
               paste(st,n_tmp[6],"_",nb,sep=""),n_tmp[7:8] ,
               paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Deut.Time, Start, End, Charge)
  ord=b$Deut.Time
  ord=as.numeric(str_sub(ord, end=-2))
  
  bp1<-data.frame(b, ord)
  bp1<-arrange(bp1, ord, Start, End)
  b<-bp1[,-dim(bp1)[2]]
  
  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )
  
  
  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  return(procent.tp)}


output_tp_noSeq_states<- function(filepath, states=unique(a$Protein.State)){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (time in unique(a$Deut.Time)){
    temp1<-a[which(a$Deut.Time ==time),]
    st_l<-c()
    nbs=0
    for (state in states){##
      temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      st<-gsub(c(' '),'',state)
      df_nm_st<-paste(st, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:5] ,
               paste(st,n_tmp[6],"_",nb,sep=""),n_tmp[7:8] ,
               paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Deut.Time, Start, End, Charge)
  ord=b$Deut.Time
  ord=as.numeric(str_sub(ord, end=-2))
  
  bp1<-data.frame(b, ord)
  bp1<-arrange(bp1, ord, Start, End)
  b<-bp1[,-dim(bp1)[2]]
  
  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )
  
  
  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  return(uptake.tp)}

output_tp_noSeq_states_proc<- function(filepath, states=unique(a$Protein.State)){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (time in unique(a$Deut.Time)){
    temp1<-a[which(a$Deut.Time ==time),]
    st_l<-c()
    nbs=0
    for (state in states){##
      temp2<-temp1[which(temp1$Protein.State ==state ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      st<-gsub(c(' '),'',state)
      df_nm_st<-paste(st, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(paste(st,n_tmp[1],"_",nb,sep=""), n_tmp[2],paste(st,n_tmp[3],"_",nb,sep=""),n_tmp[4:5] ,
               paste(st,n_tmp[6],"_",nb,sep=""),n_tmp[7:8] ,
               paste(st, "_",n_tmp[9],"_",nb,sep=""), paste(st, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Deut.Time', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Deut.Time, Start, End, Charge)
  ord=b$Deut.Time
  ord=as.numeric(str_sub(ord, end=-2))
  
  bp1<-data.frame(b, ord)
  bp1<-arrange(bp1, ord, Start, End)
  b<-bp1[,-dim(bp1)[2]]
  
  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )
  
  
  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  return(procent.tp)}

#####input timecourses
######
####



output_FD<- function(filepath, state=unique(a$Protein.State)){
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  #### Full deuteration to analysis. 
  
  fd<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(FD$Protein.State)){
    temp1<-FD[which(FD$Protein.State ==state),]
    st_l<-c()
    nb=0
    fs<-c()
    for (exp in unique(temp1$Experiment)){
      nb=nb+1
      df_nm<-paste("FD_",nb,sep="")
      temp2<-temp1[which(temp1$Experiment == exp),]
      n_tmp<-names(temp2)
      nms<-c(n_tmp[1],paste("t_FD",n_tmp[2],"_",nb,sep=""),paste("t_FD",n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , 
             paste("t_FD", "_",n_tmp[9],"_",nb,sep=""), paste("t_FD", "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
      colnames(temp2)<-nms
      assign(df_nm, temp2) 
      fs<-c(fs, df_nm) ##creates number of data.frames that equals to number of replicates
      df_List<-mget(fs) 
      ##will merge all the replicates dataframes to bp dataframe
      fp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                   'Charge')), df_List)}
    assign(df_nm_st, fp) }
  
  nb_X..deut=grep("X..Deut", colnames(fp))
  uptake.fp<-data.frame(fp[,c(1, 4:8, nb_X..deut)])
  nb_deut=grep("Deut.._", colnames(fp))
  procent.fp<-data.frame(fp[,c(1, 4:8, nb_deut)])
  
  return(uptake.fp)}

output_UD<- function(filepath, state=unique(a$Protein.State)){
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9, 21, 22)] ##choose only useful columns
  rownames(a)<-1:dim(a)[1] ##name rows
  a<-na.omit(a) ##remove missing values, removes non-deuterated state
  zero<-a[which(a$Deut.Time == '0.00s'),] ##this line can be used to take advantage of control sample
  
  ud<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(zero$Protein.State)){
    temp1<-zero[which(zero$Protein.State ==state),]
    st_l<-c()
    nb=0
    us<-c()
    for (exp in unique(temp1$Experiment)){
      nb=nb+1
      df_nm<-paste("zero_",nb,sep="")
      temp2<-temp1[which(temp1$Experiment == exp),]
      n_tmp<-names(temp2)
      nms<-c(n_tmp[1],paste("t_zero",n_tmp[2],"_",nb,sep=""),paste("t_zero",n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , 
             paste("t_zero", "_",n_tmp[9],"_",nb,sep=""), paste("t_zero", "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
      colnames(temp2)<-nms
      assign(df_nm, temp2) 
      us<-c(us, df_nm) ##creates number of data.frames that equals to number of replicates
      df_List<-mget(us) 
      ##will merge all the replicates dataframes to bp dataframe
      ud<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                   'Charge')), df_List)}
    assign(df_nm_st, ud) }
  nb_X..deut=grep("X..Deut", colnames(ud))
  uptake.ud<-data.frame(ud[,c(1, 4:8, nb_X..deut)])
  nb_deut=grep("Deut.._", colnames(ud))
  procent.ud<-data.frame(ud[,c(1, 4:8, nb_deut)])
  
  
  return(uptake.ud)} 


output_FD_proc<- function(filepath, state=unique(a$Protein.State)){
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  #### Full deuteration to analysis. 
  
  fd<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(FD$Protein.State)){
    temp1<-FD[which(FD$Protein.State ==state),]
    st_l<-c()
    nb=0
    fs<-c()
    for (exp in unique(temp1$Experiment)){
      nb=nb+1
      df_nm<-paste("FD_",nb,sep="")
      temp2<-temp1[which(temp1$Experiment == exp),]
      n_tmp<-names(temp2)
      nms<-c(n_tmp[1],paste("t_FD",n_tmp[2],"_",nb,sep=""),paste("t_FD",n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , 
             paste("t_FD", "_",n_tmp[9],"_",nb,sep=""), paste("t_FD", "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
      colnames(temp2)<-nms
      assign(df_nm, temp2) 
      fs<-c(fs, df_nm) ##creates number of data.frames that equals to number of replicates
      df_List<-mget(fs) 
      ##will merge all the replicates dataframes to bp dataframe
      fp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                   'Charge')), df_List)}
    assign(df_nm_st, fp) }
  
  nb_X..deut=grep("X..Deut", colnames(fp))
  uptake.fp<-data.frame(fp[,c(1, 4:8, nb_X..deut)])
  nb_deut=grep("Deut.._", colnames(fp))
  procent.fp<-data.frame(fp[,c(1, 4:8, nb_deut)])
  
  return(procent.fp)}

output_UD_proc<- function(filepath, state=unique(a$Protein.State)){
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  ud<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(zero$Protein.State)){
    temp1<-zero[which(zero$Protein.State ==state),]
    st_l<-c()
    nb=0
    us<-c()
    for (exp in unique(temp1$Experiment)){
      nb=nb+1
      df_nm<-paste("zero_",nb,sep="")
      temp2<-temp1[which(temp1$Experiment == exp),]
      n_tmp<-names(temp2)
      nms<-c(n_tmp[1],paste("t_zero",n_tmp[2],"_",nb,sep=""),paste("t_zero",n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , 
             paste("t_zero", "_",n_tmp[9],"_",nb,sep=""), paste("t_zero", "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
      colnames(temp2)<-nms
      assign(df_nm, temp2) 
      us<-c(us, df_nm) ##creates number of data.frames that equals to number of replicates
      df_List<-mget(us) 
      ##will merge all the replicates dataframes to bp dataframe
      ud<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                   'Charge')), df_List)}
    assign(df_nm_st, ud) }
  nb_X..deut=grep("X..Deut", colnames(ud))
  uptake.ud<-data.frame(ud[,c(1, 4:8, nb_X..deut)])
  nb_deut=grep("Deut.._", colnames(ud))
  procent.ud<-data.frame(ud[,c(1, 4:8, nb_deut)])
  
  return(procent.ud)}

output_tcourse<- function(filepath, state=unique(a$Protein.State)){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  
  
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(a$Protein.State)){
    temp1<-a[which(a$Protein.State ==state),]
    st_l<-c()
    nbs=0
    for (time in unique(temp1$Deut.Time)){##
      temp2<-temp1[which(temp1$Deut.Time ==time ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      df_nm_st<-paste(time, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(n_tmp[1],paste("t",time,n_tmp[2],"_",nb,sep=""),paste("t",time,n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , 
               paste("t",time, "_",n_tmp[9],"_",nb,sep=""), paste("t",time, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
    
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Start, End,  Charge)
  
  uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(b[,1:6], b[,grep("Deut.._", colnames(b))])
  return(uptake.tp)}




output_tcourse_proc<- function(filepath, output_file, state=unique(a$Protein.State)){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  b<-c()
  
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(a$Protein.State)){
    temp1<-a[which(a$Protein.State ==state),]
    st_l<-c()
    nbs=0
    for (time in unique(temp1$Deut.Time)){##
      temp2<-temp1[which(temp1$Deut.Time ==time ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      df_nm_st<-paste(time, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(n_tmp[1],paste("t",time,n_tmp[2],"_",nb,sep=""),paste("t",time,n_tmp[3],"_",nb,sep=""),n_tmp[4:8] , 
               paste("t",time, "_",n_tmp[9],"_",nb,sep=""), paste("t",time, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
    
  } ## b has all information bound together again to have all information df
  b<-arrange(b,  Start, End,  Charge)
  uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(b[,1:6], b[,grep("Deut.._", colnames(b))])
  return(procent.tp)}


#####

verbose_timecourse_output_proc<-function(filepath, output_name, state){
  df<-output_tcourse_proc(filepath, state)
  fp<-output_FD_proc(filepath, state)
  ud<-output_UD_proc(filepath, state)
  
  ##average for timepoints
  s1<-sd_timepoint(df)
  av1<-ave_timepoint(df)
  
  
  ### get averages or non-averages values for FD and undeuterated. 
  ###if function to check for number of replicates. if replicates =1 or smaller, 
  if ((dim(fp)[2]-6) >1){
    av_fp<-ave_timepoint(fp, replicates=(dim(fp)[2]-6))
    sd_fp<-sd_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) ==1) {
    av_fp=fp
    sd_fp<-sd_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) == 0){
    print("Full deuteration sample not provided")
  }
  
  if ((dim(ud)[2]-6) >1){
    av_ud<-ave_timepoint(ud, replicates=(dim(ud)[2]-6))
    sd_ud<-sd_timepoint(ud, replicates=(dim(ud)[2]-6))
    
  } else if ((dim(ud)[2]-6) ==1) {
    av_ud=ud
    sd_ud<-sd_timepoint(ud, replicates=(dim(ud)[2]-6))
  }  else if ((dim(ud)[2]-6) == 0){
    print("Non-deuterated sample not provided")
  }
  
  if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list(av_ud,av1, av_fp, sd_ud, s1, sd_fp) 
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list(av1, av_fp, s1, sd_fp)
  } else if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list(av_ud,av1, sd_ud, s1)
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list(av1, s1)
  }
  
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  bp<-arrange(bp, Start, End, Charge)
  write.csv(bp, output_name)
  return(bp)
}



verbose_timecourse_output<-function(filepath, output_name, state){
  df<-output_tcourse(filepath, state)
  fp<-output_FD(filepath, state)
  ud<-output_UD(filepath, state)
  
  ##average for timepoints
  s1<-sd_timepoint(df)
  av1<-ave_timepoint(df)
  
  
  ### get averages or non-averages values for FD and undeuterated. 
  ###if function to check for number of replicates. if replicates =1 or smaller, 
  if ((dim(fp)[2]-6) >1){
    av_fp<-ave_timepoint(fp, replicates=(dim(fp)[2]-6))
    sd_fp<-sd_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) ==1) {
    av_fp=fp
    sd_fp<-sd_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) == 0){
    print("Full deuteration sample not provided")
  }
  
  if ((dim(ud)[2]-6) >1){
    av_ud<-ave_timepoint(ud, replicates=(dim(ud)[2]-6))
    sd_ud<-sd_timepoint(ud, replicates=(dim(ud)[2]-6))
    
  } else if ((dim(ud)[2]-6) ==1) {
    av_ud=ud
    sd_ud<-sd_timepoint(ud, replicates=(dim(ud)[2]-6))
  }  else if ((dim(ud)[2]-6) == 0){
    print("Non-deuterated sample not provided")
  }
  
  if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list(av_ud,av1, av_fp, sd_ud, s1, sd_fp) 
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list(av1, av_fp, s1, sd_fp)
  } else if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list(av_ud,av1, sd_ud, s1)
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list(av1, s1)
  }
  
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  bp<-arrange(bp, Start, End, Charge)
  write.csv(bp, output_name)
  return(bp)
}

average_timecourse<-function(filepath, state){
  df<-output_tcourse(filepath, state)
  fp<-output_FD(filepath, state)
  ud<-output_UD(filepath, state)
  
  ##average for timepoints
  av1<-ave_timepoint(df)
  
  ### get averages or non-averages values for FD and undeuterated. 
  ###if function to check for number of replicates. if replicates =1 or smaller, 
  if ((dim(fp)[2]-6) >1){
    av_fp<-ave_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) ==1) {
    av_fp=fp
  } else if ((dim(fp)[2]-6) == 0){
    print("Full deuteration sample not provided")
  }
  
  if ((dim(ud)[2]-6) >1){
    av_ud<-ave_timepoint(ud, replicates=(dim(ud)[2]-6))
  } else if ((dim(ud)[2]-6) ==1) {
    av_ud=ud
  }  else if ((dim(ud)[2]-6) == 0){
    print("Non-deuterated sample not provided")
  }
  
  if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list(av_ud,av1, av_fp) 
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list(av1, av_fp)
  } else if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list(av_ud,av1)
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list(av1)
  }
  
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  bp<-arrange(bp, Start, End, Charge)
  return(bp)
}

sd_timecourse<-function(filepath, output_name, state){
  df<-output_tcourse(filepath, state)
  fp<-output_FD(filepath, state)
  ud<-output_UD(filepath, state)
  
  ##average for timepoints
  s1<-sd_timepoint(df)
  
  
  ### get averages or non-averages values for FD and undeuterated. 
  ###if function to check for number of replicates. if replicates =1 or smaller, 
  if ((dim(fp)[2]-6) >1){
    sd_fp<-sd_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) ==1) {
    sd_fp<-sd_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) == 0){
    print("Full deuteration sample not provided")
  }
  
  if ((dim(ud)[2]-6) >1){
    sd_ud<-sd_timepoint(ud, replicates=(dim(ud)[2]-6))
    
  } else if ((dim(ud)[2]-6) ==1) {
    sd_ud<-sd_timepoint(ud, replicates=(dim(ud)[2]-6))
  }  else if ((dim(ud)[2]-6) == 0){
    print("Non-deuterated sample not provided")
  }
  
  if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list( sd_ud, s1, sd_fp) 
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list( s1, sd_fp)
  } else if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list( sd_ud, s1)
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list( s1)
  }
  
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  bp<-arrange(bp, Start, End, Charge)
  write.csv(bp, output_name)
  return(bp)
}


average_timecourse_proc<-function(filepath, state){
  df<-output_tcourse_proc(filepath, state)
  fp<-output_FD_proc(filepath, state)
  ud<-output_UD_proc(filepath, state)
  
  ##average for timepoints
  av1<-ave_timepoint(df)
  
  ### get averages or non-averages values for FD and undeuterated. 
  ###if function to check for number of replicates. if replicates =1 or smaller, 
  if ((dim(fp)[2]-6) >1){
    av_fp<-ave_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) ==1) {
    av_fp=fp
  } else if ((dim(fp)[2]-6) == 0){
    print("Full deuteration sample not provided")
  }
  
  if ((dim(ud)[2]-6) >1){
    av_ud<-ave_timepoint(ud, replicates=(dim(ud)[2]-6))
  } else if ((dim(ud)[2]-6) ==1) {
    av_ud=ud
  }  else if ((dim(ud)[2]-6) == 0){
    print("Non-deuterated sample not provided")
  }
  
  if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list(av_ud,av1, av_fp) 
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list(av1, av_fp)
  } else if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list(av_ud,av1)
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list(av1)
  }
  
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  bp<-arrange(bp, Start, End, Charge)
  return(bp)
}

sd_timecourse_proc<-function(filepath, output_name, state){
  df<-output_tcourse_proc(filepath, state)
  fp<-output_FD_proc(filepath, state)
  ud<-output_UD_proc(filepath, state)
  
  ##average for timepoints
  s1<-sd_timepoint(df)
  
  
  ### get averages or non-averages values for FD and undeuterated. 
  ###if function to check for number of replicates. if replicates =1 or smaller, 
  if ((dim(fp)[2]-6) >1){
    sd_fp<-sd_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) ==1) {
    sd_fp<-sd_timepoint(fp, replicates=(dim(fp)[2]-6))
  } else if ((dim(fp)[2]-6) == 0){
    print("Full deuteration sample not provided")
  }
  
  if ((dim(ud)[2]-6) >1){
    sd_ud<-sd_timepoint(ud, replicates=(dim(ud)[2]-6))
    
  } else if ((dim(ud)[2]-6) ==1) {
    sd_ud<-sd_timepoint(ud, replicates=(dim(ud)[2]-6))
  }  else if ((dim(ud)[2]-6) == 0){
    print("Non-deuterated sample not provided")
  }
  
  if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list( sd_ud, s1, sd_fp) 
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) != 0  ){
    df_List<-list( s1, sd_fp)
  } else if ((dim(ud)[2]-6) != 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list( sd_ud, s1)
  } else if ((dim(ud)[2]-6) == 0 & (dim(fp)[2]-6) == 0  ){
    df_List<-list( s1)
  }
  
  ##will merge all the replicates dataframes to bp dataframe
  bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Sequence', 'Search.RT',
                                               'Charge')), df_List)
  bp<-arrange(bp, Start, End, Charge)
  write.csv(bp, output_name)
  return(bp)
}


output_tcourse_noSeq<- function(filepath, state=unique(a$Protein.State)){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  
  
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(a$Protein.State)){
    temp1<-a[which(a$Protein.State ==state),]
    st_l<-c()
    nbs=0
    for (time in unique(temp1$Deut.Time)){##
      temp2<-temp1[which(temp1$Deut.Time ==time ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      df_nm_st<-paste(time, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(n_tmp[1],paste("t",time,n_tmp[2],"_",nb,sep=""),paste("t",time,n_tmp[3],"_",nb,sep=""),
               n_tmp[4:5] , paste("t",time,n_tmp[6],"_",nb,sep=""),n_tmp[7:8], 
               paste("t",time, "_",n_tmp[9],"_",nb,sep=""), paste("t",time, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
    
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Start, End,  Charge)
  
  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )
  
  
  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  
  return(uptake.tp)}




output_tcourse_noSeq_proc<- function(filepath, state=unique(a$Protein.State)){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  
  
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(a$Protein.State)){
    temp1<-a[which(a$Protein.State ==state),]
    st_l<-c()
    nbs=0
    for (time in unique(temp1$Deut.Time)){##
      temp2<-temp1[which(temp1$Deut.Time ==time ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      df_nm_st<-paste(time, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(n_tmp[1],paste("t",time,n_tmp[2],"_",nb,sep=""),paste("t",time,n_tmp[3],"_",nb,sep=""),
               n_tmp[4:5] , paste("t",time,n_tmp[6],"_",nb,sep=""),n_tmp[7:8], 
               paste("t",time, "_",n_tmp[9],"_",nb,sep=""), paste("t",time, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
    
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Start, End,  Charge)
  
  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )
  
  
  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  
  return(procent.tp)}



output_tcourse_noSeq_states<- function(filepath, states=unique(a$Protein.State)){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  
  
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in states){
    temp1<-a[which(a$Protein.State ==state),]
    st_l<-c()
    nbs=0
    for (time in unique(temp1$Deut.Time)){##
      temp2<-temp1[which(temp1$Deut.Time ==time ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      df_nm_st<-paste(time, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(n_tmp[1],paste("t",time,n_tmp[2],"_",nb,sep=""),paste("t",time,n_tmp[3],"_",nb,sep=""),
               n_tmp[4:5] , paste("t",time,n_tmp[6],"_",nb,sep=""),n_tmp[7:8], 
               paste("t",time, "_",n_tmp[9],"_",nb,sep=""), paste("t",time, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
    
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Start, End,  Charge)
  
  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )
  
  
  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  
  return(uptake.tp)}




output_tcourse_noSeq_states_proc<- function(filepath, states=unique(a$Protein.State)){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  
  
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in states){
    temp1<-a[which(a$Protein.State ==state),]
    st_l<-c()
    nbs=0
    for (time in unique(temp1$Deut.Time)){##
      temp2<-temp1[which(temp1$Deut.Time ==time ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      df_nm_st<-paste(time, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(n_tmp[1],paste("t",time,n_tmp[2],"_",nb,sep=""),paste("t",time,n_tmp[3],"_",nb,sep=""),
               n_tmp[4:5] , paste("t",time,n_tmp[6],"_",nb,sep=""),n_tmp[7:8], 
               paste("t",time, "_",n_tmp[9],"_",nb,sep=""), paste("t",time, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
    
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Start, End,  Charge)
  
  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )
  
  
  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  
  return(procent.tp)}
output_tcourse_noSeq_proc<- function(filepath, state=unique(a$Protein.State)){
  
  a<-read.csv(file=filepath,  header = F, skip = 1)### load the file witout headers
  nm<-read.csv(file=filepath, header = T, row.names = NULL, nrows = 1)###load the header names
  nm1<-colnames(nm)
  colnames(a)<-c(nm1) ##assign names to the columns
  a<-a[,c(1:6,8,9,  21, 22)] 
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    undeut<-a[which(a$Deut.Time == '0s'),]}
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    FD<-a[which(a$Deut.Time == 'FD'),] 
    a<-a[-which(a$Deut.Time == c('FD')),]}
  
  a<-na.omit(a)
  rownames(a)<-1:dim(a)[1] ##name rows
  
  ##loop below will go through Protein states, timepoints and Experiments to get replicates 
  ##it will save a dataframe in wide format instead of long format, result of this loop is dataframe named "b"
  
  
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  
  b<-c()
  ##creates temporary df, temp1, with Protein states going through all unique protein states
  for (state in unique(a$Protein.State)){
    temp1<-a[which(a$Protein.State ==state),]
    st_l<-c()
    nbs=0
    for (time in unique(temp1$Deut.Time)){##
      temp2<-temp1[which(temp1$Deut.Time ==time ),]##creates temporary df, temp2 from one state of protein with the same timepoints
      nb=0
      nbs=nbs+1
      #print(c(time, state))
      df_nm_st<-paste(time, "_", nbs,sep="")
      st_l<-c(st_l, df_nm_st)
      bs<-c()
      for (exp in unique(temp2$Experiment)){
        nb=nb+1
        df_nm<-paste("b",nb,sep="")
        temp3<-temp2[which(temp2$Experiment == exp),]
        n_tmp<-names(temp3)
        nms<-c(n_tmp[1],paste("t",time,n_tmp[2],"_",nb,sep=""),paste("t",time,n_tmp[3],"_",nb,sep=""),
               n_tmp[4:5] , paste("t",time,n_tmp[6],"_",nb,sep=""),n_tmp[7:8], 
               paste("t",time, "_",n_tmp[9],"_",nb,sep=""), paste("t",time, "_", n_tmp[10],"_",nb,sep="")) ## creates names for the dataframe
        
        colnames(temp3)<-nms
        assign(df_nm, temp3) 
        bs<-c(bs, df_nm) ##creates number of data.frames that equals to number of replicates
        df_List<-mget(bs) 
        ##will merge all the replicates dataframes to bp dataframe
        bp<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Search.RT',
                                                     'Charge')), df_List)}
      assign(df_nm_st, bp) 
    }
    df_List2<-mget(st_l)
    bp2<-Reduce(function(x, y) merge(x, y, by = c('Protein.State', 'Start','End', 'Search.RT',
                                                  'Charge')), df_List2)
    b=rbind(b, bp2)
    
  } ## b has all information bound together again to have all information df
  b<-arrange(b, Start, End,  Charge)
  
  df_description<-data.frame(b[,1:3],  b[,grep("Sequence", colnames(b))][1], b[,4:5])
  colnames(df_description)<-c(names(b[,1:3]),"Sequence",names(b[,4:5]) )
  
  
  #uptake.tp<-data.frame(b[,1:6], b[,grep("X..Deut", colnames(b))])
  uptake.tp<-data.frame(df_description, b[,grep("X..Deut", colnames(b))])
  procent.tp<-data.frame(df_description, b[,grep("Deut.._", colnames(b))])
  
  return(procent.tp)}

#####
#####
#####




ppar<-function(mfcol2){
  par(mfrow=mfcol2,mar = c(1, 1, 1, 1), oma=c(3,2.4,1,1), cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  return()
}
ppar_bottom_legend<-function(mfcol2){
  par(mfrow=mfcol2,mar = c(1, 1, 1, 1), oma=c(4, 2.4,1.5,1), cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  return()}  ### function I use to make fast changes in numer of plots in plot area



ppar_wider<-function(mfcol2){
  par(mfrow=mfcol2,mar = c(1.5, 1.5, 1.5, 1.5), oma=c(2,2,2,2), cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  return()}  ### function I use to make fast changes in numer of plots in plot area

##function by Farid Cheraghi, 
###https://stackoverflow.com/questions/9292563/reset-the-graphical-parameters-back-to-default-values-without-use-of-dev-off
reset_par <- function(){
  op <- structure(list(xlog = FALSE, ylog = FALSE, adj = 0.5, ann = TRUE,
                       ask = FALSE, bg = "transparent", bty = "o", cex = 1, 
                       cex.axis = 1, cex.lab = 1, cex.main = 1.2, cex.sub = 1, 
                       col = "black", col.axis = "black", col.lab = "black", 
                       col.main = "black", col.sub = "black", crt = 0, err = 0L, 
                       family = "", fg = "black", fig = c(0, 1, 0, 1), 
                       fin = c(6.99999895833333, 6.99999895833333), font = 1L, 
                       font.axis = 1L, font.lab = 1L, font.main = 2L, 
                       font.sub = 1L, lab = c(5L, 5L, 7L), las = 0L, 
                       lend = "round", lheight = 1, ljoin = "round", lmitre = 10, 
                       lty = "solid", lwd = 1, mai = c(1.02, 0.82, 0.82, 0.42), 
                       mar = c(5.1, 4.1, 4.1, 2.1), mex = 1, mfcol = c(1L, 1L), 
                       mfg = c(1L, 1L, 1L,1L), mfrow = c(1L, 1L), 
                       mgp = c(3, 1, 0), mkh = 0.001, new = FALSE, 
                       oma = c(0, 0, 0, 0), omd = c(0, 1, 0, 1), 
                       omi = c(0, 0, 0,0), pch = 1L, 
                       pin = c(5.75999895833333, 5.15999895833333),
                       plt = c(0.117142874574832, 0.939999991071427, 
                               0.145714307397962, 0.882857125425167), 
                       ps = 12L, pty = "m", smo = 1, srt = 0, tck = NA_real_, 
                       tcl = -0.5, usr = c(0.568, 1.432, 0.568, 1.432), 
                       xaxp = c(0.6, 1.4, 4), xaxs = "r", xaxt = "s", 
                       xpd = FALSE, yaxp = c(0.6, 1.4, 4), yaxs = "r", 
                       yaxt = "s", ylbias = 0.2), 
                  .Names = c("xlog", "ylog", "adj", "ann", "ask", "bg", 
                             "bty", "cex", "cex.axis", "cex.lab", "cex.main", "cex.sub", 
                             "col", "col.axis", "col.lab", "col.main", "col.sub", "crt", 
                             "err", "family", "fg", "fig", "fin", "font", "font.axis", 
                             "font.lab", "font.main", "font.sub", "lab", "las", "lend", 
                             "lheight", "ljoin", "lmitre", "lty", "lwd", "mai", "mar", 
                             "mex", "mfcol", "mfg", "mfrow", "mgp", "mkh", "new", "oma",
                             "omd", "omi", "pch", "pin", "plt", "ps", "pty", "smo", 
                             "srt", "tck", "tcl", "usr", "xaxp", "xaxs", "xaxt", "xpd", 
                             "yaxp", "yaxs", "yaxt", "ylbias"))
  par(op)
}




CI_2pts<-function(s1, s2, replicates=3){ 
  tvalue=abs(qt(0.01/2, replicates*2-2))
  sp1<-sqrt(sum(s1^2*(replicates-1))/((replicates-1)*length(s1)))
  sp2<-sqrt(sum(s2^2*(replicates-1))/((replicates-1)*length(s2)[1]))
  spa<-sqrt((sp1^2+sp2^2)/replicates)
  CI<-spa*tvalue ###if number of replicates changes this number has to change!!! Loook in the table!!!
  return(CI)
}

CI_tp<-function(df, replicates=3, pv_cutoff=0.01 ){
  CI_all<-c()
  tvalue=abs(qt(pv_cutoff/2, replicates*2-2))
  for ( i in 8:dim(df)[2]){
    sp1<-sqrt(sum(df[,7]^2*(replicates-1))/((replicates-1)*length(df[,7])))
    sp2<-sqrt(sum(df[,i]^2*(replicates-1))/((replicates-1)*length(df[,i])[1]))
    spa<-sqrt((sp1^2+sp2^2)/replicates)
    CI<-spa*tvalue 
    CI_all<-c(CI_all, CI)}
  return(CI_all)
}




color_ranges_Blue_Red_heat_map<-function(ranges, colors_initial) {
  
  blue1<-brewer.pal(n = 9, name = "Blues")
  red1<-brewer.pal(n = 9, name = "Reds")
  cbr1<-c(colors_initial, rev(colorRampPalette(red1)(floor(length(ranges)/2))),
          colorRampPalette(blue1)(floor(length(ranges)/2) ))
  return(cbr1)
}



coverage_residue<-function(df1, start_col, end_col){
  
  ##prep of coverage
  vc<-c()
  vc1<-c()
  for ( i in 1:dim(df1)[1]){
    vc<-rep(0, length=max(df1[,end_col]))
    vc[df1[i,start_col]:df1[i,end_col]]=1##make multiple vectors which have 1 at position which peptide covers
    vc1<-c(vc1, vc)}
  cov.pos=data.frame(matrix(vc1, ncol =max(df1[,end_col]) , byrow = TRUE))  ##make coverage dataframe
  coverage<-colSums(cov.pos) ##sumarized occurances of peptides
  return(coverage)
}


ranges_function<-function(df_ave, values_df){
  print(paste("values range from ", round(range(values_df)[1],2), "to",
              round(range(values_df)[2],2), "%" ))
  
  lbs<-str_sub(colnames(df_ave[7:dim(df_ave)[2]]), start=4, end=-9)
  for ( i in 1:dim(values_df)[2]){
    lbs1=paste(lbs[1], lbs[i+1])
    print(paste(lbs1,"range", round(range(values_df[,i])[1],2), 
                round(range(values_df[,i])[2],2)))
  }
}





###################

av_tp<-function(df, cola=c(1, brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"))) {
  n1=max((dim(df)[2]-7), 3)
  #cola<-c(1, brewer.pal(n = n1, name = "Paired"))
  plot(df[,7], type="n", xlab="", ylab="", lwd=2, col=cola[1], ylim=c(0, max(df[,7:dim(df)[2]])))
  
  for ( i in 7:dim(df)[2]){
    points(df[,i], type="l", xlab="", ylab="", col=cola[i-6])}
  axis(1, at=seq(0, 1000, by=10), cex.axis=1, labels=F,tcl=-0.2)
  axis(2, at=seq(0, 1000, by=5), cex.axis=1, labels=F,tcl=-0.2)}




legend_raw_ave<-function(df, cola=c(1, brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"))){
  n1=max((dim(df)[2]-7), 3)
  #cola<-c(1, brewer.pal(n = n1, name = "Paired"))
  ##draw boxplots ave and sd1
  exp_du<-expression('D'[2]*'O uptake')
  mtext(c("Index"),  c(SOUTH<-1),line=0.7, outer=TRUE)
  mtext(exp_du,  c(WEST<-2),line=0.7, outer=TRUE)
  nm1<-str_sub(colnames(df[7:dim(df)[2]]), start = 4, end = -9)
  legend(c("right"), nm1,
         fill=cola,  bty="n", cex=0.6, inset=c(-0.2,0), xpd = TRUE )
}
plots_av_tp<-function(df,replicates=3, cola=c(1, brewer.pal(n = max((dim(ave_timepoint(df))[2]-7), 3), name = "Paired"))){
  par(mar = c(1,1,1,6), mfrow=c(length(unique(df$Deut.Time)), 1))
  av1<-ave_timepoint(df, replicates)
  for ( i in(unique(df$Deut.Time))){
    av_tp(av1[av1$Deut.Time==i,], cola)
    legend_raw_ave(av1, cola)
    mtext(i,  c(NORTH<-3),line=-1, outer=FALSE, cex=0.5)}
}


dif_tp<-function(df, cola=c(brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"))) {
  n1=max((dim(df)[2]-7), 3)
  #cola<-c(brewer.pal(n = n1, name = "Paired"))
  plot(df[,7], type="n", xlab="", ylab="", lwd=2, col=cola[1], 
       ylim=c(min(df[,7:dim(df)[2]])-0.5, max(df[,7:dim(df)[2]])))
  abline(h=0)
  for ( i in 8:dim(df)[2]){
    points(df[,i], type="l", xlab="", ylab="", col=cola[i-7])}
  axis(1, at=seq(0, 1000, by=10), cex.axis=1, labels=F,tcl=-0.2)
  axis(2, at=seq(0, 1000, by=2.5), cex.axis=1, labels=F,tcl=-0.2)
  
}


lab_dif<-function(df, cola=c(brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"))){
  n1=max((dim(df)[2]-7), 3)
  #cola<-c(brewer.pal(n = n1, name = "Paired"))
  exp_ddu<-expression(Delta*' D'[2]*'O uptake [Da]')
  mtext(c("Index"),  c(SOUTH<-1),line=0.7, outer=TRUE)
  mtext(exp_ddu,  c(WEST<-2),line=0.7, outer=TRUE)
  nm1<-str_sub(colnames(df[8:dim(df)[2]]), start = 4, end=-9)
  legend(c("right"), nm1,
         fill=cola,  bty="n", cex=0.6, inset=c(-0.2,0), xpd = TRUE )
  
}

lab_vol<-function(df, cola=c(brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"))){
  n1=max((dim(df)[2]-7), 3)
  #cola<-c(brewer.pal(n = n1, name = "Paired"))
  nm1<-str_sub(colnames(df[8:dim(df)[2]]), start = 4, end=-9)
  legend(c("right"), nm1,
         fill=cola,  bty="n", cex=0.6, inset=c(-0.2,0), xpd = TRUE )
  
}


plots_diff_tp<-function(df, replicates=3, cola=c(brewer.pal(n = max((dim(ave_timepoint(df))[2]-7), 3), name = "Paired"))){
  reset_par()
  ppar(c(1,1))
  par(mar = c(1,1,1,6), mfrow=c(length(unique(df$Deut.Time)), 1))
  av1<-ave_timepoint(df, replicates)
  da1<-dif_ave(av1)
  for ( i in(unique(df$Deut.Time))){
    dif_tp(da1[da1$Deut.Time==i,], cola)
    lab_dif(da1, cola)
    mtext(i,  c(NORTH<-3),line=-1, outer=FALSE, cex=0.5)}
}


vol_tp<-function(df1, pv, CI, pv_cutoff=0.01, cola=c(brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"))){ 
  abs.a<-abs(df1[,8:dim(df1)[2]])>CI
  cl1<-pv[,8:dim(df1)[2]]< pv_cutoff 
  col_res<-data.frame(abs.a*cl1)
  n1=max((dim(df1)[2]-7), 3)###either 3 or length of the variants compared to chosen state
  #cola<-c(brewer.pal(n = n1, name = "Paired"))
  plot(df1[,8], pv[,8], log="y", ylim=c(1, 10^-7), xlab="", ylab="", col=1,
       pch=20, yaxt="n", type="n", xlim=c(min(df1[,8:dim(df1)[2]]),max(df1[,8:dim(df1)[2]]) ))
  rect(xleft=c(min(df1[,8:dim(df1)[2]])-5),xright = -CI, ybottom = pv_cutoff, ytop=10^-50, col="grey95", lty=3)
  rect(xleft=CI,xright = max(df1[,8:dim(df1)[2]] )+5, ybottom = pv_cutoff, ytop=10^-50, col="grey95", lty=3)
  for ( i in 8:dim(df1)[2]){
    cola1<-c("black", cola[i-7] )
    c1<-cola1[(col_res[,i-7])+1]
    points(df1[,i], pv[,i],col=c1, pch=20, type="p") }
  box(lwd=2)
  exp_ddu<-expression(Delta*' D'[2]*'O uptake [Da]')
  mtext(c(exp_ddu, "p-value"),  c(SOUTH<-1, WEST<-2),line=0.7, outer=TRUE)
  axis(2, at=c(10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10), cex.axis=1, labels=T)
  at.y <- outer(1:9, 10^(2:(-7)))
  axis(2, at=at.y, labels=F,tcl=-0.2)}

plots_vol_tp<-function(df,replicates=3, pv_cutoff=0.01, 
                       cola=c(brewer.pal(n = max((dim(ave_timepoint(df))[2]-7), 3), name = "Paired"))){
  reset_par()
  par(mar = c(1,1,1,6), mfrow=c(length(unique(df$Deut.Time)), 1))
  av1<-ave_timepoint(df, replicates)
  da1<-dif_ave(av1)
  pv1<-pv_timepoint(df, replicates)
  s1<-sd_timepoint(df,  replicates)
  
  for ( i in(unique(df$Deut.Time))){
    CI<-max(CI_tp(s1[s1$Deut.Time==i,], replicates, pv_cutoff))
    vol_tp(da1[da1$Deut.Time==i,], pv1[pv1$Deut.Time==i,], CI, pv_cutoff, cola)
    lab_vol(da1, cola)
    mtext(i,  c(NORTH<-3),line=-1, outer=FALSE, cex=0.5)
  }}


significant_peptide_uptake<-function(df_av, pv, sd, pv_cutoff=0.01, replicates=3){
  CI=CI_tp(sd, replicates,pv_cutoff)
  abs.a<-c()
  for ( i in 8:dim(df_av)[2]){
    abs1<-abs(df_av[,i])>CI[i-7]
    abs.a<-c(abs.a, abs1)
  }
  abs.a=data.frame(matrix(abs.a, ncol =length(CI) , byrow = FALSE))
  
  cl1a<-pv[,8:dim(df_av)[2]]<pv_cutoff
  cl1<-abs.a*cl1a
  return(cl1)
}

heat_map_tp<-function(df, pv, sd, ranges=c(-Inf, seq(-30, 30, by=10), Inf), pv_cutoff=0.01, replicates=3){
  #####
  #preparation significance per residue & coverage
  
  cl1<-significant_peptide_uptake(df, pv, sd, pv_cutoff, replicates)
  
  start_col<-which(colnames(df)=='Start')
  end_col<-which(colnames(df)=='End')
  fc.d<-data.frame((df[,7]-df[,8:dim(df)[2]])/(df[,7])*10*cl1)###vector which has significant average
  fc.d[fc.d=="-Inf"]<-0
  
  ac<-data.frame(matrix(ncol=dim(fc.d)[2] , nrow =max(df[,end_col]), rep(0, length.out=max(df[,end_col])*dim(fc.d)[2])))##make mock matrix
  #ac<-data.frame(matrix(ncol=dim(fc.d)[2] , nrow =max(df[,end_col]), rep(0, length.out=dim(max(df[,end_col]))*dim(fc.d)[2])))##make mock matrix
  for ( j in 1:dim(ac)[2]){
    for ( i in 1:dim(df)[1]){
      
      sum.nm<-(ac[df[i,start_col]:df[i,end_col],j]+fc.d[i,j])
      ###create a vector of which is a sum of significance at position
      ac[df[i,start_col]:df[i,end_col],j]=sum.nm[1] ## assign new value to position
    }}
  
  ##prep of coverage
  coverage<-coverage_residue(df, start_col, end_col)
  
  ave.p.cov<-ac/coverage ## sums of the significant avererages divided by coverage. 
  ave.p.cov[ave.p.cov=="NaN"]<-0 ##remove NAN divisions introduced by division by no coverage
  
  ranges_function(df,ave.p.cov ) ###print ranges of data
  
  #####preparation of average per residue data.frame, which will have 
  
  
  xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
  
  cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c("white", "grey45"))
  
  
  for ( i in 1:(length(xli)-1)){
    ave.p.cov[xli[i]< ave.p.cov & ave.p.cov < xli[i+1]] <- num_ass[i]
  }
  ave.p.cov[ave.p.cov==0]<- (-10000)
  
  si_apc<-abs(ave.p.cov)-9999
  cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)###define lack of coverage
  si_apc<-si_apc*cv_mis+1
  
  ###define missing coverage
  
  
  plot(c(1,1), type="n", xlab="", ylab="", lwd=2, col=4, xlim=c(min(df[,start_col]), max(df[,end_col])-5), 
       ylim=c(0, (dim(df)[2]-7)), yaxt="n") ## mock plot, just to have it drawn correct limits set up
  xl <- 1; yb <- (0); xr <- max(df[,end_col]); yt <- (1)
  for ( i in 0:(dim(df)[2]-8)){
    yb=i; yt=i+1 ##loop to have initial values for y postions in loop to use multiple postion
    rect(head(seq(xl,xr,(xr-xl)/xr),-1),yb,
         tail(seq(xl,xr,(xr-xl)/xr),-1), yt,col=cbr1[si_apc[,i+1]], border = NA)
    
    ###coverage
    # rect(head(seq(xl,xr,(xr-xl)/a[dim(a)[1],4]),-1)*cc,yb,
    #      tail(seq(xl,xr,(xr-xl)/a[dim(a)[1],4]),-1)*cc, yt,col=cc[col_cv_mis+1], border = NA)
  }
  axis(1, at=seq(0, 700, by=25),  tcl=-0.2, labels = F)
  abline(h=0:7, lwd=0.5, col="grey30")
  box(lwd=2)
  return()
}



legend_heat_map_tp<-function(df){
  nm1<-str_sub(colnames(df[8:dim(df)[2]]), start=4, end=-9)
  axis(2, at=0.5:(dim(df)[2]-7.5), labels=nm1, las=2, line = , cex.axis=0.7)
}

pallette_legend<-function(col_pallette){
  xl <- 1; yb <- 1; xr <- 1.5; yt <- 2 ## seq values
  plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")# mock plot
  rect(xl,head(seq(yb,yt,(yt-yb)/length(col_pallette)),-1),xr,tail(seq(yb,yt,(yt-yb)/length(col_pallette)),-1),col=col_pallette)# make rect
  return()}

pallette_ll<-function(pallette, lab){
  ppar(c(1,1))
  pallette_legend(pallette)
  yb=1; yt=2
  ypos=seq(yb,yt,(yt-yb)/length(pallette))
  ypos.a<-c()
  
  for ( i in 1:(length(ypos)-1)){
    ypos.a<-c(ypos.a, (ypos[i]+ypos[i+1])/2)}
  axis(4, at=ypos.a,labels =lab,las=2, cex.axis=0.6, tcl=0, pos=1.5)
}

plot_heat_map_tp<-function(df, lim1=3.5, ranges=c(-Inf, -3,-2,-1, 0,1, 2,3, Inf),
                           pv_cutoff=0.01, replicates=3){
  pv1<-pv_timepoint(df, replicates)
  s1<-sd_timepoint(df, replicates)
  av1<-ave_timepoint(df, replicates)
  par(mfrow=c(length(unique(av1$Deut.Time)),1),
      mar = c(1.5, lim1, 1, 1.1), oma=c(3,2.4,1,1), 
      cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  for ( i in(unique(av1$Deut.Time))){
    print(paste("For time point", i))
    a1=av1[av1$Deut.Time==i,]
    p1=pv1[pv1$Deut.Time==i,]
    sd1=s1[s1$Deut.Time==i,]
    colmp<-heat_map_tp(a1, p1, sd1, ranges, pv_cutoff, replicates=3)
    legend_heat_map_tp(av1)
    mtext(i, side=3, outer=FALSE, line=0, cex=0.65)}
  mtext(c("Residues"),  c(NORTH<-1),line=0.7, outer=TRUE, cex=0.8)
  return()}


peptide_pv_tp<-function(df, pv, sd, nb_row, ranges=c(-Inf, seq(-30, 30, by=10), Inf), 
                        pv_cutoff=0.01, replicates=3){
  #####
  #preparation significance per residue & coverage
  cl1<-significant_peptide_uptake(df, pv, sd, pv_cutoff, replicates)
  
  start_col<-which(colnames(df)=='Start')
  end_col<-which(colnames(df)=='End')
  
  fc.d<-data.frame((df[,7]-df[,8:dim(df)[2]])/df[,7]*10*cl1)
  
  # fc.d<-data.frame((df[,7]-df[,8:dim(df)[2]])*10*cl1)###vector which has significant average
  # 
  # for ( i in 1:dim(fc.d)[2]){
  #   fc.d[,i]<-(fc.d[,i]/((df[,7]+df[,i+7])*0.5)) }
  
  
  
  si.fv<-(fc.d)
  
  ranges_function(df,si.fv ) ###print ranges of data
  
  xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
  for ( i in 1:(length(xli)-1)){
    si.fv[xli[i]< si.fv & si.fv < xli[i+1]] <- num_ass[i]}
  si.fv[si.fv==0]<- (-10000)
  si.fv<-abs(si.fv)-9999
  ######
  ##color set up
  
  cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c("black"))
  
  
  for ( j in 1:dim(si.fv)[2]){
    y1=c(rep(1:nb_row, times=floor(dim(df)[1]/nb_row)), 1:(dim(df)[1]%%nb_row)) ##y values on the plot corresponding to peptides index
    plot(c(1,1), type="n", xlab="", ylab="", lwd=2, col=4, xlim=c(min(df[,start_col]), max(df[,end_col])), 
         ylim=c(nb_row, 0), yaxt="n", xaxt="n") ## mock plot, just to have it drawn correct limits set up
    
    for ( i in 1:dim(df)[1]){
      #print(c(df[i,start_col], df[i,end_col],y1[i], cbr1[si.fv[i,j]] ))
      points(c(df[i,start_col], df[i,end_col]), c(y1[i], y1[i]), type="l",
             col=cbr1[si.fv[i,j]])}
    
    main_nm=str_sub(colnames(df)[j+7], end = -9, start=4)
    mtext(main_nm, side=1, outer=FALSE, line=0, cex=0.6)
    axis(3, at=seq(0, max(df[, end_col])+25, by=25),  tcl=-0.2, labels = F)
    axis(3, at=seq(0, max(df[, end_col])+25, by=50), labels = T, cex.axis=0.65)
    box(lwd=2)
  }}

plot_peptide_sig_tp<-function(df1, replicates=3, nb_pep_row=100, ranges=c(-Inf, seq(-30, 30, by=10), Inf), pv_cutoff=0.01){
  pv1<-pv_timepoint(df1,  replicates)
  s1<-sd_timepoint(df1, replicates)
  av1<-ave_timepoint(df1, replicates)
  for ( k in(unique(av1$Deut.Time))){
    print(paste("For timepoint", k))
    a1=av1[av1$Deut.Time==k,]
    p1=pv1[pv1$Deut.Time==k,]
    sd1=s1[s1$Deut.Time==k,]
    peptide_pv_tp(a1, p1, sd1, nb_pep_row,ranges, pv_cutoff, replicates)
    #legend_heat_map_tp(av1)
    mtext(k, side=4, outer=FALSE, line=0.2, cex=0.7)}
  mtext(c("Residues"),  c(NORTH<-3),line=0, outer=TRUE, cex=0.7)}

###########


boxplot_tp<-function(df, replicates=3){
  ppar(c(length(unique(df$Deut.Time)),1))
  av1<-ave_timepoint(df, replicates)
  s1<-sd_timepoint(df, replicates)
  nm1<-str_sub(colnames(av1[c(6+(1:(dim(s1)[2]-6)))]), start = 4, end=-9)
  for ( i in(unique(df$Deut.Time))){
    boxplot(av1[av1$Deut.Time==i,7:dim(av1)[2]], names = nm1, col=1)
    mtext(i, side=3, outer=FALSE, line=-0.9, cex=0.5 )}
  mtext("Average", side=2, outer=TRUE, line=0.8, cex=1 )
}



heat_map_tp_maxuptake<-function(df, pv, sd, ranges=c(-Inf, seq(-30, 30, by=10), Inf),
                                pv_cutoff=0.01, replicates=3){
  #####
  #preparation significance per residue & coverage
  cl1<-significant_peptide_uptake(df, pv, sd, pv_cutoff, replicates)
  
  start_col<-which(colnames(df)=='Start')
  end_col<-which(colnames(df)=='End')
  fc.d<-data.frame((df[,7]-df[,8:dim(df)[2]])/(df[,7])*10*cl1)###vector which has significant average
  fc.d[fc.d=="-Inf"]<-0
  max.ac1<-c() 
  for ( j in 1:dim(fc.d)[2]){
    ac<-c()
    ac1<-c()
    ac2<-c()
    for ( i in 1:dim(df)[1]){
      ac<-rep(0, length=max(df[,end_col]))
      ac[df[i,start_col]:df[i,end_col]]=fc.d[i,j]
      ##make multiple vectors which have 1 at position which peptide covers
      ac1<-c(ac1, ac)}
    ac2=data.frame(matrix(ac1, nrow =dim(df)[1], byrow=T))
    max.a<-c()
    for ( k in 1:dim(ac2)[2]){ 
      ind1<-which.max(abs(ac2[,k]))
      nb1<-(ac2[ind1,k])
      max.a<-c(max.a, nb1)}
    max.ac1<-c(max.ac1, max.a)}
  max.ac2=data.frame(matrix(max.ac1, ncol = dim(fc.d)[2]))
  
  
  
  ##prep of coverage
  coverage<-coverage_residue(df, start_col, end_col) 
  
  ranges_function(df,max.ac2 ) ###print ranges of data
  
  #####preparation of average per residue data.frame, which will have 
  xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
  for ( i in 1:(length(xli)-1)){
    max.ac2[xli[i]< max.ac2 & max.ac2 < xli[i+1]] <- num_ass[i]
  }
  max.ac2[max.ac2==0]<- (-10000)
  cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)
  si_apc<-abs(max.ac2)-9999
  ###define lack of coverage
  si_apc<-si_apc*cv_mis+1
  
  ###define missing coverage
  ##pallette definition
  ##color set up
  
  cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c("white", "grey45"))
  
  
  plot(c(1,1), type="n", xlab="", ylab="", lwd=2, col=4, xlim=c(min(df[,start_col]), max(df[,end_col])-5), 
       ylim=c(0, (dim(df)[2]-7)), yaxt="n") ## mock plot, just to have it drawn correct limits set up
  xl <- 1; yb <- (0); xr <- max(df[,end_col]); yt <- (1)
  for ( i in 0:(dim(df)[2]-8)){
    yb=i; yt=i+1 ##loop to have initial values for y postions in loop to use multiple postion
    rect(head(seq(xl,xr,(xr-xl)/xr),-1),yb,
         tail(seq(xl,xr,(xr-xl)/xr),-1), yt,col=cbr1[si_apc[,i+1]], border = NA)
    
    ###coverage
    # rect(head(seq(xl,xr,(xr-xl)/a[dim(a)[1],4]),-1)*cc,yb,
    #      tail(seq(xl,xr,(xr-xl)/a[dim(a)[1],4]),-1)*cc, yt,col=cc[col_cv_mis+1], border = NA)
  }
  axis(1, at=seq(0, 700, by=25),  tcl=-0.2, labels = F)
  abline(h=0:7, lwd=0.5, col="grey30")
  box(lwd=2)
  return()
}

plot_heat_map_max_uptake_tp<-function(df, replicates=3, lim1=3.5, ranges=c(-Inf, seq(-30, 30, by=10), Inf), pv_cutoff=0.01){
  pv1<-pv_timepoint(df, replicates)
  s1<-sd_timepoint(df, replicates)
  av1<-ave_timepoint(df, replicates)
  par(mfrow=c(length(unique(av1$Deut.Time)),1),
      mar = c(1.5, lim1, 1, 1.1), oma=c(3,2.4,1,1), 
      cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  for ( i in(unique(av1$Deut.Time))){
    a1=av1[av1$Deut.Time==i,]
    p1=pv1[pv1$Deut.Time==i,]
    sd1=s1[s1$Deut.Time==i,]
    print(paste("For timepoint", i))
    colmp<-heat_map_tp_maxuptake(a1, p1, sd1, ranges, pv_cutoff, replicates)
    legend_heat_map_tp(av1)
    mtext(i, side=3, outer=FALSE, line=0, cex=0.65)}
  mtext(c("Residues"),  c(NORTH<-1),line=0.7, outer=TRUE, cex=0.8)
  return()}

#heat_map_tp_maxuptake(av1, pv1, sd,ranges=c(-Inf, -3,-2,-1, 0,1, 2,3, Inf) )
#plot_heat_map_max_uptake_tp(h)

pymol_str<-function(ind1) {ind2<-paste(as.character(ind1), sep="' '", collapse="+")
return(ind2)}



pymol_script_significant_residue<-function(df, ranges=c(-Inf, seq(-30, 30, by=10), Inf), 
                                           pv_cutoff=0.01,  replicates=3){
  #####from HDX get data and 
  
  for ( deut_time in(unique(df$Deut.Time))){
    pv<-pv_timepoint(df[df$Deut.Time==deut_time,])
    sd<-sd_timepoint(df[df$Deut.Time==deut_time,])
    df1<-ave_timepoint(df[df$Deut.Time==deut_time,])
    
    
    #preparation significance per residue & coverage
    cl1<-significant_peptide_uptake(df1, pv, sd, pv_cutoff, replicates)
    
    
    start_col<-which(colnames(df1)=='Start')
    end_col<-which(colnames(df1)=='End')
    
    ###preparation of the coloring ranges. Difference of average. 
    fc.d<-data.frame((df1[,7]-df1[,8:dim(df1)[2]])/(df1[,7])*10*cl1)###vector which has significant average
    ###per residue maximum uptake value is chosen. -> 
    max.ac1<-c() 
    for ( j in 1:dim(fc.d)[2]){
      ac<-c()
      ac1<-c()
      ac2<-c()
      for ( i in 1:dim(df1)[1]){
        ac<-rep(0, length=max(df1[,end_col]))
        ac[df1[i,start_col]:df1[i,end_col]]=fc.d[i,j]
        ##make multiple vectors which have 1 at position which peptide covers
        ac1<-c(ac1, ac)}
      ac2=data.frame(matrix(ac1, nrow =dim(df1)[1], byrow=T))
      max.a<-c()
      for ( k in 1:dim(ac2)[2]){ 
        ind1<-which.max(abs(ac2[,k]))
        nb1<-(ac2[ind1,k])
        max.a<-c(max.a, nb1)}
      max.ac1<-c(max.ac1, max.a)}
    max.ac2=data.frame(matrix(max.ac1, ncol = dim(fc.d)[2]))
    
    
    ##prep of coverage
    coverage<-coverage_residue(df1, start_col, end_col)
    
    #####preparation of average per residue data.frame, which will have 
    xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
    for ( i in 1:(length(xli)-1)){
      max.ac2[xli[i]< max.ac2 & max.ac2 < xli[i+1]] <- num_ass[i]}
    max.ac2[max.ac2==0]<- (-10000)
    
    si_apc<-abs(max.ac2)-9999
    cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)###define lack of coverage
    si_apc<-si_apc*cv_mis+1
    
    
    cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c("white", "grey45"))
    
    
    ##assign name to colors in pallettes
    col_nm<-c("no_cov", "NSig")
    for ( i in 1:(length(ranges)-1)){
      col_nm<-c(col_nm, paste("col_", ranges[i],"_", ranges[i+1], sep=""))}
    
    rgb_col<-col2rgb(cbr1, alpha = FALSE) ## function to return rgb values for colors in pallette
    
    ##set_color 0%, [0 , 0 , 120], make command for pymol, to set colors
    set_colors<-c()
    for ( i in 1:length(col_nm)){
      set_colors<-c(set_colors, paste("set_color ", col_nm[i],", [",  rgb_col[1,i], 
                                      ",", rgb_col[2,i], ",", rgb_col[3,i], "]", sep=""))}
    ###write outputs per each state in the 
    nm1<-str_sub(colnames(df1[8:dim(df1)[2]]), start=4, end=-9)
    for (j in 1:dim(si_apc)[2]){
      output_name<-paste("pymol_maxuptake_", nm1[j],"_", deut_time, ".txt", sep="")
      res.txt<-c()
      for ( i in 1:length(col_nm)){
        if (length(which(si_apc[,j]==i)) !=0){
          res.txt<-c(res.txt, (paste(c("color ", col_nm[i],", resi ",
                                       pymol_str(which(si_apc[,j] ==i))), sep="", collapse="")))}}
      fileConn<-file(output_name)
      writeLines(c("hide","show cartoon","color grey", "bg white",set_colors,res.txt ), fileConn)
      close(fileConn)}}
  
  leg_nm<-c("No coverage", "Not Sig")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],":", ranges[i+1], "%", sep=""))}
  
  pallette_ll(cbr1, leg_nm)
  return()}


pymol_script_significant_peptide<-function(df, ranges=c(-Inf, seq(-30, 30, by=10), Inf), 
                                           pv_cutoff=0.01, replicates=3){
  #####from HDX get data and 
  for ( deut_time in(unique(df$Deut.Time))){
    
    pv<-pv_timepoint(df[df$Deut.Time==deut_time,], replicates)
    sd<-sd_timepoint(df[df$Deut.Time==deut_time,], replicates)
    df1<-ave_timepoint(df[df$Deut.Time==deut_time,], replicates)
    
    
    #preparation significance per residue & coverage
    cl1<-significant_peptide_uptake(df1, pv, sd, pv_cutoff, replicates)
    
    
    start_col<-which(colnames(df1)=='Start')
    end_col<-which(colnames(df1)=='End')
    
    ###preparation of the coloring ranges. Difference of average. 
    fc.d<-data.frame((df1[,7]-df1[,8:dim(df1)[2]])/(df1[,7])*10*cl1)###vector which has significant average
    ##sumarized occurances of peptides
    si.f=fc.d
    #####preparation of average per residue data.frame, which will have 
    xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
    for ( i in 1:(length(xli)-1)){
      si.f[xli[i]<  si.f &  si.f < xli[i+1]] <- num_ass[i]}
    si.f[ si.f==0]<- (-10000)
    si_apc<-abs(si.f)-9999
    
    
    ##pallette definition
    ##color set up
    
    cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c("grey45"))
    
    
    ##assign name to colors in pallettes
    col_nm<-c("NSig")
    for ( i in 1:(length(ranges)-1)){
      col_nm<-c(col_nm, paste("col_", ranges[i],"_", ranges[i+1], sep=""))}
    
    
    rgb_col<-col2rgb(cbr1, alpha = FALSE) ## function to return rgb values for colors in pallette
    
    ##set_color 0%, [0 , 0 , 120], make command for pymol, to set colors
    set_colors<-c()
    for ( i in 1:length(col_nm)){
      set_colors<-c(set_colors, paste("set_color ", col_nm[i],", [",  rgb_col[1,i], 
                                      ",", rgb_col[2,i], ",", rgb_col[3,i], "]", sep=""))}
    
    
    ###write outputs per each state in the 
    nm1<-str_sub(colnames(df1[8:dim(df1)[2]]), start=4, end=-9)
    for (j in 1:dim(si_apc)[2]){
      output_name<-paste("pymol_all_peptides_", nm1[j],"_", deut_time, ".txt", sep="")
      
      res.txt<-c()
      for ( i in 1:length(col_nm)){
        if (length(which(si_apc[,j]==i)) !=0 ) {
          pep_nb<-which(si_apc[,j] == i)
          for ( k in pep_nb){line<-c()
          line<- paste(c("color ", col_nm[i],", resi ", df1[k,start_col], "-", df1[k,end_col] , sep=""))
          res.txt<-c(res.txt, 
                     paste(line, sep="' '", collapse=""))}
        }}
      fileConn<-file(output_name)
      writeLines(c("hide","show cartoon","color black", "bg white",
                   set_colors, res.txt ), fileConn)
      close(fileConn)}}
  leg_nm<-c("Not Sig")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],":", ranges[i+1], "%", sep=""))}
  pallette_ll(cbr1, leg_nm)
  return()}

legend_heat_map<-function(ranges=c(-Inf,seq(-30, 30, by=10), Inf )){
  ##color set up
  cbr1<-color_ranges_Blue_Red_heat_map(ranges, c("white", "grey45"))
  
  leg_nm<-c("no coverage", "not-significant")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],"%:", ranges[i+1], "%", sep=""))}
  
  pallette_ll(cbr1, leg_nm)}


pymol_script_average_residue<-function(df, ranges=c(-Inf, seq(-30, 30, by=10), Inf), 
                                       pv_cutoff=0.01, replicates=3){
  #####from HDX get data and 
  
  for ( deut_time in(unique(df$Deut.Time))){
    pv<-pv_timepoint(df[df$Deut.Time==deut_time,],  replicates)
    sd<-sd_timepoint(df[df$Deut.Time==deut_time,], replicates)
    df1<-ave_timepoint(df[df$Deut.Time==deut_time,], replicates)
    
    
    #preparation significance per residue & coverage
    cl1<-significant_peptide_uptake(df1, pv, sd, pv_cutoff, replicates)
    
    start_col<-which(colnames(df1)=='Start')
    end_col<-which(colnames(df1)=='End')
    
    fc.d<-data.frame((df1[,7]-df1[,8:dim(df1)[2]])/(df1[,7])*10*cl1)###vector which has significant average
    
    ac<-data.frame(matrix(ncol=dim(fc.d)[2] , nrow =max(df1[,end_col]), rep(0, length.out=max(df1[,end_col])*dim(fc.d)[2])))##make mock matrix
    #ac<-data.frame(matrix(ncol=dim(fc.d)[2] , nrow =max(df1[,end_col]), rep(0, length.out=dim(max(df1[,end_col]))*dim(fc.d)[2])))##make mock matrix
    for ( j in 1:dim(ac)[2]){
      for ( i in 1:dim(df1)[1]){
        
        sum.nm<-(ac[df1[i,start_col]:df1[i,end_col],j]+fc.d[i,j])
        ###create a vector of which is a sum of significance at position
        ac[df1[i,start_col]:df1[i,end_col],j]=sum.nm[1] ## assign new value to position
      }}
    
    
    ##prep of coverage
    coverage<-coverage_residue(df1, start_col, end_col)
    
    ave.p.cov<-ac/coverage ## sums of the significant avererages divided by coverage. 
    ave.p.cov[ave.p.cov=="NaN"]<-0 ##remove NAN divisions introduced by division by no coverage
    
    
    #####preparation of average per residue data.frame, which will have 
    
    xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
    cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c("white", "grey45"))
    
    
    for ( i in 1:(length(xli)-1)){
      ave.p.cov[xli[i]< ave.p.cov & ave.p.cov < xli[i+1]] <- num_ass[i]}
    ave.p.cov[ave.p.cov==0]<- (-10000)
    
    si_apc<-abs(ave.p.cov)-9999
    cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)###define lack of coverage
    si_apc<-si_apc*cv_mis+1
    
    ##assign name to colors in pallettes
    col_nm<-c("no_cov", "NSig")
    for ( i in 1:(length(ranges)-1)){
      col_nm<-c(col_nm, paste("col_", ranges[i],"_", ranges[i+1], sep=""))}
    
    rgb_col<-col2rgb(cbr1, alpha = FALSE) ## function to return rgb values for colors in pallette
    
    ##set_color 0%, [0 , 0 , 120], make command for pymol, to set colors
    set_colors<-c()
    for ( i in 1:length(col_nm)){
      set_colors<-c(set_colors, paste("set_color ", col_nm[i],", [",  rgb_col[1,i], 
                                      ",", rgb_col[2,i], ",", rgb_col[3,i], "]", sep=""))}
    ###write outputs per each state in the 
    nm1<-str_sub(colnames(df1[8:dim(df1)[2]]), start=4, end=-9)
    for (j in 1:dim(si_apc)[2]){
      output_name<-paste("pymol_average_", nm1[j],"_", deut_time, ".txt", sep="")
      res.txt<-c()
      for ( i in 1:length(col_nm)){
        if (length(which(si_apc[,j]==i)) !=0){
          res.txt<-c(res.txt, (paste(c("color ", col_nm[i],", resi ",
                                       pymol_str(which(si_apc[,j] ==i))), sep="", collapse="")))}}
      fileConn<-file(output_name)
      writeLines(c("hide","show cartoon","color grey", "bg white",set_colors,res.txt ), fileConn)
      close(fileConn)}}
  
  leg_nm<-c("No coverage", "Not Sig")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],":", ranges[i+1], "%", sep=""))}
  
  pallette_ll(cbr1, leg_nm)
  return()}
#####

ranges_function_tc<-function(df_ave, values_df){
  print(paste("values range from ", round(range(values_df)[1],2), "to",
              round(range(values_df)[2],2), "%" ))
  
  lbs<-str_sub(colnames(df_ave[7:dim(df_ave)[2]]), start=4, end=-9)
  for ( i in 1:dim(values_df)[2]){
    
    print(paste(lbs[i],"range", round(range(values_df[,i])[1],2), 
                round(range(values_df[,i])[2],2)))
  }
}

color_ranges_Spectral<-function(ranges, colors_initial) {
  if (length(ranges)-1 <11){
    n1=max(length(ranges), 3)
    
    cbr1<-c(colors_initial,rev(brewer.pal(n = n1, name = "Spectral")))
  } else {
    Spectral<-brewer.pal(n = 11, name = "Spectral")
    cbr1<-c(colors_initial, rev(colorRampPalette(Spectral)(floor(length(ranges))-1)))}
  return(cbr1)
}

color_ranges_Spectral<-function(ranges, colors_initial) {
  Spectral<-brewer.pal(n = 11, name = "Spectral")
  cbr1<-c(colors_initial, rev(colorRampPalette(Spectral)(floor(length(ranges))-1)))
  return(cbr1)
}






####
heat_map_tc<-function(df, ranges=c(seq(0, 100, by=10), Inf)){
  #####
  #preparation significance per residue & coverage
  
  start_col<-which(colnames(df)=='Start')
  end_col<-which(colnames(df)=='End')
  fc.d<-df[,7:dim(df)[2]]###vector which has significant average
  ac<-data.frame(matrix(ncol= dim(df)[1], nrow =max(df[,end_col]), rep(0, length.out=max(df[,end_col])*dim(df)[1])))##make mock matrix
  
  sum.nm=0
  ac1<-c()
  #ac<-data.frame(matrix(ncol=dim(fc.d)[2] , nrow =max(df[,end_col]), rep(0, length.out=dim(max(df[,end_col]))*dim(fc.d)[2])))##make mock matrix
  for ( j in 1:dim(fc.d)[2]){
    for ( i in 1:dim(df)[1]){
      
      sum.nm<-fc.d[i,j]
      ac[df[i,start_col]:df[i,end_col],i]<-sum.nm[1]
      ###create a vector of which is a sum of significance at position #assign new value to position
    } 
    ac1<-c(ac1, rowSums(ac))}
  ac2=data.frame(matrix(ac1, ncol =dim(fc.d)[2] , byrow = FALSE))
  ##prep of coverage
  coverage<-coverage_residue(df, start_col, end_col)
  ave.p.cov<-ac2/coverage ## sums of the significant avererages divided by coverage. 
  ave.p.cov[ave.p.cov=="NaN"]<-0 ##remove NAN divisions introduced by division by no coverage
  ranges_function_tc(df,ave.p.cov)
  #####preparation of average per residue data.frame, which will have 
  xli=ranges; num_ass<-c(-10000001:(-10000001-(length(xli)-2)))
  for ( i in 1:(length(xli)-1)){
    ave.p.cov[xli[i]< ave.p.cov & ave.p.cov < xli[i+1]] <- num_ass[i]
  }
  
  ave.p.cov[ave.p.cov==0]<- (-10000000)
  
  si_apc<-abs(ave.p.cov)-9999999
  cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)###define lack of coverage
  si_apc<-si_apc*cv_mis+1
  
  ###define missing coverage
  cbr1<-color_ranges_Spectral(ranges=xli, c("white", 1))
  
  plot(c(1,1), type="n", xlab="", ylab="", lwd=2, col=4, xlim=c(min(df[,start_col]), max(df[,end_col])-5), 
       ylim=c(0, (dim(si_apc)[2])), yaxt="n") ## mock plot, just to have it drawn correct limits set up
  xl <- 1; yb <- (0); xr <- max(df[,end_col]); yt <- (1)
  for ( i in 0:(dim(si_apc)[2]-1)){
    yb=i; yt=i+1 ##loop to have initial values for y postions in loop to use multiple postion
    rect(head(seq(xl,xr,(xr-xl)/xr),-1),yb,
         tail(seq(xl,xr,(xr-xl)/xr),-1), yt,col=cbr1[si_apc[,i+1]], border = NA)
    
    ###coverage
    # rect(head(seq(xl,xr,(xr-xl)/a[dim(a)[1],4]),-1)*cc,yb,
    #      tail(seq(xl,xr,(xr-xl)/a[dim(a)[1],4]),-1)*cc, yt,col=cc[col_cv_mis+1], border = NA)
  }
  axis(1, at=seq(0, 700, by=25),  tcl=-0.2, labels = F)
  abline(h=0:7, lwd=0.5, col="grey30")
  box(lwd=2)
  return()
}

legend_heat_map_tc<-function(df){
  nm1<-str_sub(colnames(df[7:dim(df)[2]]), start=4, end=-9)
  axis(2, at=0.5:(dim(df)[2]-6.5), labels=nm1, las=2, line = , cex.axis=0.7)
}



legend_heat_map_timecourse<-function(ranges=c(-Inf,seq(-30, 30, by=10), Inf )){
  cbr1<-color_ranges_Spectral(ranges, c("white"))
  
  leg_nm<-c("no coverage")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],"%:", ranges[i+1], "%", sep=""))}
  
  pallette_ll(cbr1, leg_nm)}


plot_heat_map_tc<-function(df, replicates=3, lim1=3.5, ranges=c(-Inf, seq(0, 100, by=10), Inf)){
  av1<-ave_timepoint(df, replicates)
  par(mfrow=c(length(unique(av1$Protein.State)),1),
      mar = c(1.5, lim1, 1, 1.1), oma=c(3,2.4,1,1), 
      cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  for ( i in(unique(av1$Protein.State))){
    a1=av1[av1$Protein.State==i,]
    colmp<-heat_map_tc(a1, ranges)
    legend_heat_map_tc(a1)
    mtext(i, side=3, outer=FALSE, line=0, cex=0.65)}
  mtext(c("Residues"),  c(NORTH<-1),line=0.7, outer=TRUE, cex=0.8)
  return()}




pv_timecourse<-function(df_c,df_v, replicates=3) {
  ### calculate standard deviation of sample 1 for all data points
  nb_sets=(dim(df_c)[2]-6)/replicates
  nm_root<-colnames(df_c[,c(6+((1:nb_sets)-1)*replicates+1)])
  nm_root<-str_sub(nm_root, end = -3)
  pv_nm<-paste("pv_", nm_root, sep="")
  
  pv1<-c(); for ( j in 1:nb_sets) {
    combined_df<-merge(df_c, df_v, by = c('Start','End', 'Sequence', 'Search.RT',
                                          'Charge'))
    for (i in 1:dim(combined_df)[1]) {x1<-combined_df[i,(6+(j-1)*replicates+1):(6+(j-1)*replicates+replicates)]; 
    x2<-combined_df[i,(19+(j-1)*replicates+1):(19+(j-1)*replicates+replicates)]
    tt<-c(); tt<-t.test(x1, x2)
    pv1<-c(pv1, tt$p.val)} 
  }
  #print(sd1)
  pv2<-data.frame(matrix(pv1, ncol=nb_sets , byrow = FALSE))
  colnames(pv2)<-pv_nm
  pv2<-data.frame(combined_df[,1:6], pv2)
  return(pv2)
}


CI_tc<-function(sd_c, sd_v, replicates=3, pv_cutoff=0.01 ){
  CI_all<-c()
  tvalue=abs(qt(pv_cutoff/2, replicates*2-2))
  
  for ( i in 7:dim(sd_c)[2]){
    sp1<-sqrt(sum(sd_c[,i]^2*(replicates-1))/((replicates-1)*length(sd_c[,i])))
    sp2<-sqrt(sum(sd_v[,i]^2*(replicates-1))/((replicates-1)*length(sd_v[,i])[1]))
    spa<-sqrt((sp1^2+sp2^2)/replicates)
    CI<-spa*tvalue 
    CI_all<-c(CI_all, CI)}
  return(CI_all)
}

prep_timecourse_plot_ave<-function(control_df, variant_df, replicates=3){
  av_c<-  ave_timepoint(control_df, replicates)
  av_v<-  ave_timepoint(variant_df, replicates)
  comb_av<-merge(av_c, av_v, by = c('Start','End', 'Sequence', 'Search.RT','Charge'))
  comb_av<-arrange(comb_av, c(Start))
  sh_avc<-comb_av[,c(1:dim(av_c)[2])]
  sh_avv<-comb_av[,c(1:5, (dim(av_c)[2]+1):(dim(av_c)[2]*2-5))]
  return(list(sh_avc, sh_avv))
}

prep_timecourse_plot_sd<-function(control_df_up, variant_df_up, replicates=3, pv_cutoff=0.01){
  sd_c<-  sd_timepoint(control_df_up, replicates)
  sd_v<-  sd_timepoint(variant_df_up, replicates)
  CI_all<-CI_tc(sd_c, sd_v, replicates, pv_cutoff)
  return(CI_all)}

legend_tc_bottom<-function(df, cols=cola){
  nm1<-str_sub(colnames(df[7:dim(df)[2]]), start = 5, end = -10)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", nm1, xpd = TRUE, horiz = TRUE, inset = c(0,   0), bty = "n",  fill = cols, cex = 0.8)
}




robot_plot_All<-function(thP, th, replicates=3, 
                         pvalue=0.01, states=unique(thP$Protein.State), CI_factor=1){
  for ( state in states[2:length(states)]) {
    
    control_df<- thP[thP$Protein.State==states[1],]
    variant_df<- thP[thP$Protein.State==state,]
    
    control_df_up<- th[th$Protein.State==states[1],]
    variant_df_up<- th[th$Protein.State==state,]
    
    pv1<-pv_timecourse(df_c = control_df_up, df_v=variant_df_up, replicates=3)
    lav.proc<-prep_timecourse_plot_ave(control_df, variant_df, replicates=3)
    lav.proc_up<-prep_timecourse_plot_ave(control_df, variant_df, replicates=3)
    
    sh_avc<-lav.proc[[1]]
    sh_avv<-lav.proc[[2]]
    sh_avc_up<-lav.proc_up[[1]]
    sh_avv_up<-lav.proc_up[[2]]
    
    CI_all<-prep_timecourse_plot_sd(control_df = control_df_up, variant_df = variant_df_up, replicates=3, pv_cutoff=pvalue)
    CI_all<-CI_all*CI_factor
    
    cola<-(brewer.pal(n = length(7:dim(sh_avc)[2])+1, name = "Oranges"))
    
    plot(x=1, type = "n", ylim=c(-115, 115), xlim=c(min(thP$Start), max(thP$End)), ylab="", 
         xlab="", yaxt="n")
    axis(1, at=seq(0, 1000, by=10), cex.axis=1, labels=F,tcl=-0.2)
    axis(2, at=seq(-1000, 1000, by=50), cex.axis=1, labels=c(rev(seq(50,1000, by=50)), seq(0,1000, by=50)))
    axis(2, at=seq(-1000, 1000, by=10), cex.axis=1, labels=F,tcl=-0.2)
    exp_ddu<-expression('% Deuteration')
    mtext(c("Residue"),  c(SOUTH<-1),line=0.3, outer=TRUE, cex=0.8)
    mtext(exp_ddu,  c(WEST<-2),line=0.7, outer=TRUE, cex=0.85)
    
    nb1=1
    peptide_all<-c()
    for (i in 7:dim(sh_avc)[2]) {
      peptide_all<-c(peptide_all, which(pv1[, i]<pvalue & 
                                          abs(sh_avc_up[, i]-sh_avv_up[, i]) > CI_all[i-6]))
    }
    peptide_all<-sort(unique(peptide_all))
    
    colg<-(brewer.pal(n = length(7:dim(sh_avc)[2])+2, name = "Blues"))
    
    for ( i in dim(sh_avc)[2]:7){
      xpoly<-c((sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, 
               rev((sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2))
      ypoly<-c(sh_avc[peptide_all,i], rev(sh_avv[peptide_all,i]*(-1)))
      polygon(x =xpoly,                           # X-Coordinates of polygon
              y = ypoly,                             # Y-Coordinates of polygon
              col = colg[i-6])}
    abline(h=0)
    
    for ( j in 7:dim(sh_avc)[2]){
      peptide_index<-which(pv1[,j]<pvalue & abs(sh_avc_up[,j]-sh_avv_up[,j]) > CI_all[j-6])
      nb1=nb1+1
      for ( i in peptide_all){
        points(c(sh_avc$Start[i], sh_avc$End[i]), c(sh_avc[i,j],sh_avc[i,j] ), type="l", col="grey45")
        points(c(sh_avv$Start[i], sh_avv$End[i]), c(sh_avv[i,j],sh_avv[i,j])*(-1), type="l", col="grey45")
      }
      for ( pinx in peptide_index){
        points(c(sh_avc$Start[pinx], sh_avc$End[pinx]), c(sh_avc[pinx,j],sh_avc[pinx,j] ), type="l", col=cola[nb1],  lwd=2)
        points(c(sh_avv$Start[pinx], sh_avv$End[pinx]), c(sh_avv[pinx,j],sh_avv[pinx,j])*(-1), type="l", col=cola[nb1], lwd=2)
      }
      
      points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avc[peptide_all,j] ), type="p", col="grey45", pch=20, lwd=2)
      points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avv[peptide_all,j])*(-1), type="p", col="grey45",pch=20, lwd=2)
      points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avc[peptide_all,j] ), type="l", col=cola[nb1], pch=20)
      points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avv[peptide_all,j])*(-1), type="l", col=cola[nb1],pch=20)
      
      points(c(sh_avc$Start[peptide_index]+sh_avc$End[peptide_index])/2, c(sh_avc[peptide_index,j] ), type="p", col=cola[nb1], pch=20, lwd=2)
      points(c(sh_avc$Start[peptide_index]+sh_avc$End[peptide_index])/2, c(sh_avv[peptide_index,j])*(-1), type="p", col=cola[nb1],pch=20, lwd=2)
    }  
    text(x=(min(thP$Start)+max(thP$End))/2,y=106, states[1], cex=0.7)
    text(x=(min(thP$Start)+max(thP$End))/2,y=-106, state, cex=0.7)
  }
  legend_tc_bottom(sh_avc, cola[2:length(cola)])
  reset_par()
  ppar(c(1,1))}

robot_indexes_df<-function(thP, th, replicates=3, 
                           pvalue=0.01, states=unique(thP$Protein.State), CI_factor=1){
  
  for ( state in states[2:length(states)]) {
    
    control_df<- thP[thP$Protein.State==states[1],]
    variant_df<- thP[thP$Protein.State==state,]
    
    control_df_up<- th[th$Protein.State==states[1],]
    variant_df_up<- th[th$Protein.State==state,]
    
    pv1<-pv_timecourse(df_c = control_df_up, df_v=variant_df_up, replicates=3)
    lav.proc<-prep_timecourse_plot_ave(control_df, variant_df, replicates=3)
    lav.proc_up<-prep_timecourse_plot_ave(control_df, variant_df, replicates=3)
    
    sh_avc<-lav.proc[[1]]
    sh_avv<-lav.proc[[2]]
    sh_avc_up<-lav.proc_up[[1]]
    sh_avv_up<-lav.proc_up[[2]]
    
    CI_all<-prep_timecourse_plot_sd(control_df = control_df_up, variant_df = variant_df_up, replicates=3, pv_cutoff = pvalue)
    CI_all<-CI_all*CI_factor
    peptide_all<-c()
    for (i in 7:dim(sh_avc)[2]) {
      peptide_all<-c(peptide_all, which(pv1[, i]<pvalue & 
                                          abs(sh_avc_up[, i]-sh_avv_up[, i]) > CI_all[i-6]))
    }
    peptide_all<-sort(unique(peptide_all))
  }
  
  sig_pep<-data.frame(1:length(peptide_all),sh_avc[peptide_all,c(1:3,5)])
  colnames(sig_pep)[1]<-"nb"
  return(sig_pep)}



robot_indexes<-function(thP, th, replicates=3, 
                        pvalue=0.01, states=unique(thP$Protein.State), CI_factor=1){
  
  for ( state in states[2:length(states)]) {
    print(paste(state, states[1]))
    control_df<- thP[thP$Protein.State==states[1],]
    variant_df<- thP[thP$Protein.State==state,]
    
    control_df_up<- th[th$Protein.State==states[1],]
    variant_df_up<- th[th$Protein.State==state,]
    
    pv1<-pv_timecourse(df_c = control_df_up, df_v=variant_df_up, replicates=3)
    lav.proc<-prep_timecourse_plot_ave(control_df, variant_df, replicates=3)
    lav.proc_up<-prep_timecourse_plot_ave(control_df, variant_df, replicates=3)
    
    sh_avc<-lav.proc[[1]]
    sh_avv<-lav.proc[[2]]
    sh_avc_up<-lav.proc_up[[1]]
    sh_avv_up<-lav.proc_up[[2]]
    
    CI_all<-prep_timecourse_plot_sd(control_df = control_df_up, variant_df = variant_df_up, replicates=3, pv_cutoff=pvalue)
    CI_all<-CI_all*CI_factor
    peptide_all<-c()
    for (i in 7:dim(sh_avc)[2]) {
      peptide_all<-c(peptide_all, which(pv1[, i]<pvalue & 
                                          abs(sh_avc_up[, i]-sh_avv_up[, i]) > CI_all[i-6]))
    }
    peptide_all<-sort(unique(peptide_all))
    
  }
  
  return(peptide_all)}


robot_2states_indexes<-function(thP, th,indexes, states, replicates=3, 
                                pvalue=0.01, 
                                ylim=c(-110, 120), xlim=c(min(thP$Start), max(thP$End)),
                                CI_factor=1){
  
  control_df<- thP[thP$Protein.State==states[1],]
  variant_df<- thP[thP$Protein.State==states[2],]
  
  control_df_up<- th[th$Protein.State==states[1],]
  variant_df_up<- th[th$Protein.State==states[2],]
  
  pv1<-pv_timecourse(df_c = control_df_up, df_v=variant_df_up, replicates=3)
  lav.proc<-prep_timecourse_plot_ave(control_df, variant_df, replicates=3)
  lav.proc_up<-prep_timecourse_plot_ave(control_df, variant_df, replicates=3)
  
  sh_avc<-lav.proc[[1]]
  sh_avv<-lav.proc[[2]]
  sh_avc_up<-lav.proc_up[[1]]
  sh_avv_up<-lav.proc_up[[2]]
  
  CI_all<-prep_timecourse_plot_sd(control_df = control_df_up, variant_df = variant_df_up, replicates=3, pv_cutoff=pvalue)
  CI_all=CI_all*CI_factor
  
  cola<-(brewer.pal(n = length(7:dim(sh_avc)[2])+1, name = "Oranges"))
  
  plot(x=1, type = "n", ylim=ylim, xlim=xlim, ylab="", 
       xlab="", yaxt="n")
  axis(1, at=seq(0, 1000, by=10), cex.axis=1, labels=F,tcl=-0.2)
  axis(2, at=seq(-1000, 1000, by=50), cex.axis=1, labels=c(rev(seq(50,1000, by=50)), seq(0,1000, by=50)))
  axis(2, at=seq(-1000, 1000, by=10), cex.axis=1, labels=F,tcl=-0.2)
  exp_ddu<-expression('% Deuteration')
  mtext(c("Residue"),  c(SOUTH<-1),line=0.7, outer=TRUE, cex=0.8)
  mtext(exp_ddu,  c(WEST<-2),line=0.7, outer=TRUE, cex=0.8)
  cov_nb<-c()
  
  
  nb1=1
  peptide_all<-indexes
  
  
  colg<-(brewer.pal(n = length(7:dim(sh_avc)[2])+2, name = "Blues"))
  
  for ( i in dim(sh_avc)[2]:7){
    xpoly<-c((sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, 
             rev((sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2))
    ypoly<-c(sh_avc[peptide_all,i], rev(sh_avv[peptide_all,i]*(-1)))
    polygon(x =xpoly,                           # X-Coordinates of polygon
            y = ypoly,                             # Y-Coordinates of polygon
            col = colg[i-6])}
  abline(h=0)
  
  for ( j in 7:dim(sh_avc)[2]){
    peptide_index<-which(pv1[peptide_all,j]<pvalue & abs(sh_avc_up[peptide_all,j]-sh_avv_up[peptide_all,j]) > CI_all[j-6])
    peptide_index<-peptide_all[peptide_index]
    nb1=nb1+1
    for ( i in peptide_all){
      points(c(sh_avc$Start[i], sh_avc$End[i]), c(sh_avc[i,j],sh_avc[i,j] ), type="l", col="grey45")
      points(c(sh_avv$Start[i], sh_avv$End[i]), c(sh_avv[i,j],sh_avv[i,j])*(-1), type="l", col="grey45")
    }
    for ( pinx in peptide_index){
      points(c(sh_avc$Start[pinx], sh_avc$End[pinx]), c(sh_avc[pinx,j],sh_avc[pinx,j] ), type="l", col=cola[nb1],  lwd=2)
      points(c(sh_avv$Start[pinx], sh_avv$End[pinx]), c(sh_avv[pinx,j],sh_avv[pinx,j])*(-1), type="l", col=cola[nb1], lwd=2)
    }
    
    points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avc[peptide_all,j] ), type="p", col="grey45", pch=20, lwd=2)
    points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avv[peptide_all,j])*(-1), type="p", col="grey45",pch=20, lwd=2)
    points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avc[peptide_all,j] ), type="l", col=cola[nb1], pch=20)
    points(c(sh_avc$Start[peptide_all]+sh_avc$End[peptide_all])/2, c(sh_avv[peptide_all,j])*(-1), type="l", col=cola[nb1],pch=20)
    
    points(c(sh_avc$Start[peptide_index]+sh_avc$End[peptide_index])/2, c(sh_avc[peptide_index,j] ), type="p", col=cola[nb1], pch=20, lwd=2)
    points(c(sh_avc$Start[peptide_index]+sh_avc$End[peptide_index])/2, c(sh_avv[peptide_index,j])*(-1), type="p", col=cola[nb1],pch=20, lwd=2)
    text(x=(min(thP$Start)+max(thP$End))/2,y=ylim[2]-10, states[1], cex=0.7)
    text(x=(min(thP$Start)+max(thP$End))/2,y=ylim[1], states[2], cex=0.7)
  }
  
  for ( indx in indexes){
    cov_nb<-c(cov_nb, sh_avc$Start[indx]:sh_avc$End[indx])}
  cov_nb<-unique(cov_nb)
  col_index<- as.numeric(min(sh_avc$Start):max(sh_avc$End) %in% cov_nb)
  xl <- min(sh_avc$Start) ; yb <- (-1); xr <- max(sh_avc$End); yt <- (1)
  ##loop to have initial values for y postions in loop to use multiple postion
  rect(head(seq(xl-0.5,xr+0.5,1),-1),ylim[2],
       tail(seq(xl-0.5,xr+0.5,1),-1), ylim[2]+5,col=c("grey35", cola[4])[col_index+1], border = NA)
  
  legend_tc_bottom(sh_avc, cola[2:length(cola)])
  
  
  reset_par()
  ppar(c(1,1))
}



# 
# plots_av_tcourse<-function(df, 
#                            cola=c(1, brewer.pal(n = max((dim(ave_timepoint(df))[2]-7), 3), name = "Paired"))){
#   par(mar = c(1,1,1,6), mfrow=c(length(unique(df$Protein.State)), 1))
#   av1<-ave_timepoint(df)
#   for ( i in(unique(df$Protein.State))){
#     av_tc(av1[av1$Protein.State==i,], cola)
#     mtext(i,  c(North<-3), line=0, outer=FALSE, cex=0.4)}
#   legend_raw_ave_tc(av1, cola)
# }

plots_av_tcourse<-function(df, replicates=3,
                           cola=c(1, brewer.pal(n = max((dim(ave_timepoint(df))[2]-7), 3), name = "Paired"))){
  par(mar = c(1,1,1,1), mfrow=c(length(unique(df$Protein.State)), 1), oma = c(3.5, 2.75, 1, 1))
  av1<-ave_timepoint(df, replicates)
  for ( i in(unique(df$Protein.State))){
    av_tc(av1[av1$Protein.State==i,], cola)
    mtext(i,  c(North<-3), line=0, outer=FALSE, cex=0.45)}
  legend_raw_ave_tc(av1, cola)
}



av_tc<-function(df, cola=c(1, brewer.pal(n = max((dim(ave_timepoint(df))[2]-7), 3), name = "Paired"))) {
  plot(df[,7], type="n", xlab="", ylab="", lwd=2, cex.axis=0.8,
       col=cola[1], ylim=c(-5, max(df[,7:dim(df)[2]])+10))
  for ( i in 7:dim(df)[2]){
    points(df[,i], type="l", xlab="", ylab="", col=cola[i-6])}
  axis(1, at=seq(0, 1000, by=10), cex.axis=0.75, labels=F,tcl=-0.2)
  axis(2, at=seq(0, 1000, by=5), cex.axis=0.75, labels=F,tcl=-0.2)}


legend_raw_ave_tc<-function(df,cola=c(1, brewer.pal(n = max((dim(ave_timepoint(df))[2]-7), 3), name = "Paired"))){
  ##draw boxplots ave and sd1
  mtext(c("Index"),  c(SOUTH<-1),line=0.7, outer=TRUE, cex=0.7)
  mtext("% Deuteration",  c(WEST<-2),line=0.7, outer=TRUE, cex=0.7)
  nm1<-str_sub(colnames(df[7:dim(df)[2]]), start = 4, end = -9)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", nm1, xpd = TRUE, horiz = TRUE, inset = c(0,   0), bty = "n",  fill = cola, cex = 0.8)
}


#### EXTREme


extreme_input<-function(hm_dir, replicates,timepoints ){
  setwd(hm_dir)
  dirs<-list.dirs(full.names = TRUE, recursive = TRUE)
  
  
  for (directories in dirs){
    print(directories)
    nd<-list.files(directories, pattern="Non-D")
    fd<-list.files(directories, pattern="Full-D")
    tp<-c(list.files(directories,pattern="s-"), list.files(directories,pattern="m-"), list.files(directories, pattern="h-"))
    if (all(is.na(c(nd, fd, tp)))==TRUE) {
      next         ## << Only line that differs
    } else {
      
      tp1<-as.data.frame(tp) %>% separate(tp, into = paste("V", 1:3, sep =""),  sep="([\\-])", extra="merge")
      
      
      hms=str_sub(tp1[,1], start=-1)
      nbs.t=as.numeric(str_sub(tp1[,1], end=-2))
      
      for ( i in 1:length(hms)){
        if (hms[i]=="h"){
          nbs.t[i]<-nbs.t[i]*3600
        } else if (hms[i]=="m"){
          nbs.t[i]<-nbs.t[i]*60
        } else (nbs.t[i]=nbs.t[i])
      }
      
      tp1[,1]<-nbs.t
      
      for (charge in unique(tp1[,3])){
        ind_one_charge=which(tp1[,3]==charge)
        one_rep<-tp1[which(tp1[,3]==charge),]
        one_nbs.t<-nbs.t[ind_one_charge]
        
        nd1<-as.data.frame(nd) %>% separate(nd, into = paste("V", 1:4, sep = c("-")), sep="([\\-])", extra="merge")
        fd1<-as.data.frame(fd) %>% separate(fd, into = paste("V", 1:4, sep = c("-")), sep="([\\-])", extra="merge")
        
        
        replicates_each<-as.data.frame(table(one_nbs.t))
        #sort_tp<-sort(tp1[,1])
        
        
        mock.df<-data.frame(c(NA, NA), c(NA, NA))
        colnames(mock.df)<-c("V1", "V2")
        ##preparation of files to read
        
        ###
        
        nm_nd<-c(paste("t-", str_sub(c(nd[which(nd1[,4]==charge)[1]]), end=-5), sep=""))
        nm_nd<-gsub("-", "_", nm_nd)
        if(nm_nd=="t_NA"){
          nm_nd<-"mock.df"
        } else {my.files <- c(nd[which(nd1[,4]==charge)[1]])
        for (k in 1:(length(my.files))){
          full_path=paste(directories,"/", my.files[k], sep="")
          # import the file
          cur.file <- read.csv(file = full_path, header=F)
          my.name <- nm_nd[k]
          # assign the name to the object
          assign(paste(my.name), cur.file)}
        } 
        
        nm_fd<-c(paste("t-", str_sub(c(fd[which(fd1[,4]==charge)[1]]), end=-5), sep=""))
        nm_fd<-gsub("-", "_", nm_fd)
        if(nm_fd=="t_NA"){
          nm_fd<-"mock.df"
        } else {my.files <- c(fd[which(fd1[,4]==charge)[1]])
        
        for (k in 1:(length(my.files))){
          full_path=paste(directories,"/", my.files[k], sep="")
          # import the file
          cur.file <- read.csv(file = full_path, header=F)
          my.name <- nm_fd[k]
          # assign the name to the object
          assign(paste(my.name), cur.file)}
        } 
        
        
        
        
        nm1<-c(paste("t-", str_sub(c(tp), end=-5), sep="")[ind_one_charge] )
        nm1<-gsub("-", "_", nm1)
        my.files <- c( tp[ind_one_charge])
        for (k in 1:(length(my.files))){
          full_path=paste(directories,"/", my.files[k], sep="")
          # import the file
          cur.file <- read.csv(file = full_path, header=F)
          my.name <- nm1[k]
          # assign the name to the object
          assign(paste(my.name), cur.file)}
        #print(my.files)
        #my.files<-my.files[!my.files %in% NA]
        
        
        
        
        #### order the files correctly + add mock dif to lists of files that need to bound together. 
        
        order_files<-order(one_rep)
        nm1_rb<-c()
        for ( ti in 1:length(timepoints)){
          nb_rep=length(which(one_rep[,1]==timepoints[ti]))
          nm_cb=which(one_rep[,1]==timepoints[ti])
          if (nb_rep == replicates ){
            nm1_rb<-c(nm1_rb, nm1[nm_cb])
          } else if (nb_rep<replicates & nb_rep >0){
            dif_mock=replicates-nb_rep
            nm1_rb<-c(nm1_rb, nm1[nm_cb])
            nm1_rb<-c(nm1_rb, rep("mock.df", dif_mock))
          } else if (nb_rep==0){
            nm1_rb<-c(nm1_rb, rep("mock.df", replicates))
          }
        }
        nm1_rb<-c(nm_nd,nm1_rb, nm_fd)
        
        dfB<-c()
        dfB<-mget(nm1_rb[1])[[1]]
        for ( rbi in 2:length(nm1_rb)){
          dfB<-qpcR:::cbind.na(dfB, mget(nm1_rb[rbi])[[1]])}
        
        
        labs<-c(paste("undeut", sep=""))
        
        for ( ti in 1:length(timepoints)){
          for (tj in seq(replicates)){
            labs<-c(labs,paste(timepoints[ti],".0", tj, " sec", sep=""))}}
        labs<-c(labs,paste("TD",  sep=""))
        
        labsy<-rep(" ", length.out=length(labs))
        
        nms_df<-c()
        for ( ti in 1: length(labsy)){
          nms_df<-c(nms_df, labs[ti], labsy[ti])}
        
        colnames(dfB)<-nms_df
        
        dir_nm<-gsub("/", "", directories)
        
        output_name<-paste(hm_dir,"/", dir_nm,"_", charge,  ".csv", sep="")
        
        write.table(dfB, output_name,na = "",row.names = FALSE, sep = ",")}}
    
  }}




#################
####

legend_raw_ave_proc<-function(df, cola=c(1, brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"))) {
  #cola<-c(1, brewer.pal(n = n1, name = "Paired")
  ##draw boxplots ave and sd1
  exp_du<-expression('% D'[2]*'O ')
  mtext(c("Index"),  c(SOUTH<-1),line=0.7, outer=TRUE)
  mtext(exp_du,  c(WEST<-2),line=0.7, outer=TRUE)
  nm1<-str_sub(colnames(df[7:dim(df)[2]]), start = 4, end = -9)
  legend(c("right"), nm1,
         fill=cola,  bty="n", cex=0.6, inset=c(-0.2,0), xpd = TRUE )
}
plots_av_tp_proc<-function(df,replicates=3, cola=c(1, brewer.pal(n = max((dim(ave_timepoint(df))[2]-7), 3), name = "Paired"))) {
  #cola<-c(1, brewer.pal(n = n1, name = "Paired")){
  par(mar = c(1,1,1,6), mfrow=c(length(unique(df$Deut.Time)), 1))
  av1<-ave_timepoint(df, replicates)
  for ( i in(unique(df$Deut.Time))){
    av_tp(av1[av1$Deut.Time==i,], cola)
    legend_raw_ave_proc(av1, cola)
    mtext(i,  c(NORTH<-3),line=-1, outer=FALSE, cex=0.5)}
}




lab_dif_proc<-function(df, cola=c(brewer.pal(n = max((dim(df)[2]-7), 3), name = "Paired"))){
  
  exp_ddu<-expression(Delta*' %D'[2]*'O')
  mtext(c("Index"),  c(SOUTH<-1),line=0.7, outer=TRUE)
  mtext(exp_ddu,  c(WEST<-2),line=0.7, outer=TRUE)
  nm1<-str_sub(colnames(df[8:dim(df)[2]]), start = 4, end=-9)
  legend(c("right"), nm1,
         fill=cola,  bty="n", cex=0.6, inset=c(-0.2,0), xpd = TRUE )
  
}



plots_diff_tp_proc<-function(df,replicates=3, cola=c(brewer.pal(n = max((dim(ave_timepoint(df))[2]-7), 3), name = "Paired"))){
  reset_par()
  ppar(c(1,1))
  par(mar = c(1,1,1,6), mfrow=c(length(unique(df$Deut.Time)), 1))
  av1<-ave_timepoint(df, replicates)
  da1<-dif_ave(av1)
  for ( i in(unique(df$Deut.Time))){
    dif_tp(da1[da1$Deut.Time==i,], cola)
    lab_dif_proc(da1, cola)
    mtext(i,  c(NORTH<-3),line=-1, outer=FALSE, cex=0.5)}
}






heat_map_tp_proc<-function(df,dfup, pv, sd, ranges=c(-Inf, seq(-30, 30, by=10), Inf), 
                           pv_cutoff=0.01, replicates=3){
  #####
  #preparation significance per residue & coverage
  
  cl1<-significant_peptide_uptake(dfup, pv, sd, pv_cutoff, replicates)
  
  start_col<-which(colnames(df)=='Start')
  end_col<-which(colnames(df)=='End')
  fc.d<-data.frame((df[,7]-df[,8:dim(df)[2]])*cl1)###vector which has significant average
  
  ac<-data.frame(matrix(ncol=dim(fc.d)[2] , nrow =max(df[,end_col]), rep(0, length.out=max(df[,end_col])*dim(fc.d)[2])))##make mock matrix
  #ac<-data.frame(matrix(ncol=dim(fc.d)[2] , nrow =max(df[,end_col]), rep(0, length.out=dim(max(df[,end_col]))*dim(fc.d)[2])))##make mock matrix
  for ( j in 1:dim(ac)[2]){
    for ( i in 1:dim(df)[1]){
      
      sum.nm<-(ac[df[i,start_col]:df[i,end_col],j]+fc.d[i,j])
      ###create a vector of which is a sum of significance at position
      ac[df[i,start_col]:df[i,end_col],j]=sum.nm[1] ## assign new value to position
    }}
  ##prep of coverage
  coverage<-coverage_residue(df, start_col, end_col)
  
  ave.p.cov<-ac/coverage ## sums of the significant avererages divided by coverage. 
  ave.p.cov[ave.p.cov=="NaN"]<-0 ##remove NAN divisions introduced by division by no coverage
  
  ranges_function(df,ave.p.cov )
  
  #####preparation of average per residue data.frame, which will have 
  xli=ranges; num_ass<-c(-10001:(-10001-(length(xli)-2)))
  for ( i in 1:(length(xli)-1)){
    ave.p.cov[xli[i]< ave.p.cov & ave.p.cov < xli[i+1]] <- num_ass[i]
  }
  ave.p.cov[ave.p.cov==0]<- (-10000)
  
  si_apc<-abs(ave.p.cov)-9999
  cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)###define lack of coverage
  si_apc<-si_apc*cv_mis+1
  
  ###define missing coverage
  cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c("white", "grey45"))
  
  
  
  plot(c(1,1), type="n", xlab="", ylab="", lwd=2, col=4, xlim=c(min(df[,start_col]), max(df[,end_col])-5), 
       ylim=c(0, (dim(df)[2]-7)), yaxt="n") ## mock plot, just to have it drawn correct limits set up
  xl <- 1; yb <- (0); xr <- max(df[,end_col]); yt <- (1)
  for ( i in 0:(dim(df)[2]-8)){
    yb=i; yt=i+1 ##loop to have initial values for y postions in loop to use multiple postion
    rect(head(seq(xl,xr,(xr-xl)/xr),-1),yb,
         tail(seq(xl,xr,(xr-xl)/xr),-1), yt,col=cbr1[si_apc[,i+1]], border = NA)
    
    ###coverage
    # rect(head(seq(xl,xr,(xr-xl)/a[dim(a)[1],4]),-1)*cc,yb,
    #      tail(seq(xl,xr,(xr-xl)/a[dim(a)[1],4]),-1)*cc, yt,col=cc[col_cv_mis+1], border = NA)
  }
  axis(1, at=seq(0, 700, by=25),  tcl=-0.2, labels = F)
  abline(h=0:7, lwd=0.5, col="grey30")
  box(lwd=2)
  return()
}


plot_heat_map_tp_proc<-function(input_proc, input_up, lim1=3.5, ranges=c(-Inf, -3,-2,-1, 0,1, 2,3, Inf), 
                                pv_cutoff=0.01, replicates=3){
  pv1<-pv_timepoint(input_up, replicates)
  s1<-sd_timepoint(input_up, replicates)
  av1<-ave_timepoint(input_proc, replicates)
  avu<-ave_timepoint(input_up, replicates)
  par(mfrow=c(length(unique(av1$Deut.Time)),1),
      mar = c(1.5, lim1, 1, 1.1), oma=c(3,2.4,1,1), 
      cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  for ( i in(unique(av1$Deut.Time))){
    print(paste("For time point", i))
    a1=av1[av1$Deut.Time==i,]
    au=avu[avu$Deut.Time==i,]
    p1=pv1[pv1$Deut.Time==i,]
    sd1=s1[s1$Deut.Time==i,]
    colmp<-heat_map_tp_proc(a1,au, p1, sd1, ranges, pv_cutoff, replicates)
    legend_heat_map_tp(av1)
    mtext(i, side=3, outer=FALSE, line=0, cex=0.65)}
  mtext(c("Residues"),  c(NORTH<-1),line=0.7, outer=TRUE, cex=0.8)
  return()}

heat_map_tp_maxuptake_proc<-function(df, dfup, pv, sd, ranges=c(-Inf, seq(-30, 30, by=10), Inf), 
                                     pv_cutoff=0.01, replicates=3){
  #####
  #preparation significance per residue & coverage
  cl1<-significant_peptide_uptake(dfup, pv, sd, pv_cutoff, replicates)
  
  start_col<-which(colnames(df)=='Start')
  end_col<-which(colnames(df)=='End')
  fc.d<-data.frame((df[,7]-df[,8:dim(df)[2]])*cl1)###vector which has significant average
  
  max.ac1<-c() 
  for ( j in 1:dim(fc.d)[2]){
    ac<-c()
    ac1<-c()
    ac2<-c()
    for ( i in 1:dim(df)[1]){
      ac<-rep(0, length=max(df[,end_col]))
      ac[df[i,start_col]:df[i,end_col]]=fc.d[i,j]
      ##make multiple vectors which have 1 at position which peptide covers
      ac1<-c(ac1, ac)}
    ac2=data.frame(matrix(ac1, nrow =dim(df)[1], byrow=T))
    max.a<-c()
    for ( k in 1:dim(ac2)[2]){ 
      ind1<-which.max(abs(ac2[,k]))
      nb1<-(ac2[ind1,k])
      max.a<-c(max.a, nb1)}
    max.ac1<-c(max.ac1, max.a)}
  max.ac2=data.frame(matrix(max.ac1, ncol = dim(fc.d)[2]))
  
  
  ##prep of coverage
  coverage<-coverage_residue(df, start_col, end_col)
  
  ranges_function(df,max.ac2 )
  
  #ave.p.cov<-ac/coverage ## sums of the significant avererages divided by coverage. 
  #ave.p.cov[ave.p.cov=="NaN"]<-0 ##remove NAN divisions introduced by division by no coverage
  
  
  #####preparation of average per residue data.frame, which will have 
  xli=ranges; num_ass<-c(-10001:(-10001-(length(xli)-2)))
  for ( i in 1:(length(xli)-1)){
    max.ac2[xli[i]< max.ac2 & max.ac2 < xli[i+1]] <- num_ass[i]
  }
  max.ac2[max.ac2==0]<- (-10000)
  
  si_apc<-abs(max.ac2)-9999
  cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)###define lack of coverage
  si_apc<-si_apc*cv_mis+1
  
  ###define missing coverage
  ##pallette definition
  cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c("white", "grey45"))
  
  
  
  
  plot(c(1,1), type="n", xlab="", ylab="", lwd=2, col=4, xlim=c(min(df[,start_col]), max(df[,end_col])-5), 
       ylim=c(0, (dim(df)[2]-7)), yaxt="n") ## mock plot, just to have it drawn correct limits set up
  xl <- 1; yb <- (0); xr <- max(df[,end_col]); yt <- (1)
  for ( i in 0:(dim(df)[2]-8)){
    yb=i; yt=i+1 ##loop to have initial values for y postions in loop to use multiple postion
    rect(head(seq(xl,xr,(xr-xl)/xr),-1),yb,
         tail(seq(xl,xr,(xr-xl)/xr),-1), yt,col=cbr1[si_apc[,i+1]], border = NA)
    
    ###coverage
    # rect(head(seq(xl,xr,(xr-xl)/a[dim(a)[1],4]),-1)*cc,yb,
    #      tail(seq(xl,xr,(xr-xl)/a[dim(a)[1],4]),-1)*cc, yt,col=cc[col_cv_mis+1], border = NA)
  }
  axis(1, at=seq(0, 700, by=25),  tcl=-0.2, labels = F)
  abline(h=0:7, lwd=0.5, col="grey30")
  box(lwd=2)
  return()
}

plot_heat_map_max_uptake_tp_proc<-function(input_proc, input_up, lim1=3.5, ranges=c(-Inf, -3,-2,-1, 0,1, 2,3, Inf), 
                                           pv_cutoff=0.01, replicates=3){
  pv1<-pv_timepoint(input_up, replicates)
  s1<-sd_timepoint(input_up, replicates)
  av1<-ave_timepoint(input_proc, replicates)
  avu<-ave_timepoint(input_up, replicates)
  par(mfrow=c(length(unique(av1$Deut.Time)),1),
      mar = c(1.5, lim1, 1, 1.1), oma=c(3,2.4,1,1), 
      cex.axis=1, cex.main=1, cex.lab=1.1, mgp=c(0.1, 0.4, 0), ps=14, font=2, bg="white", font.lab=2, font.axis=2)
  for ( i in(unique(av1$Deut.Time))){
    print(paste("For time point", i))
    a1=av1[av1$Deut.Time==i,]
    au=avu[avu$Deut.Time==i,]
    p1=pv1[pv1$Deut.Time==i,]
    sd=s1[s1$Deut.Time==i,]
    colmp<-heat_map_tp_maxuptake_proc(a1, au, p1, sd, ranges, pv_cutoff, replicates)
    legend_heat_map_tp(av1)
    mtext(i, side=3, outer=FALSE, line=0, cex=0.65)}
  mtext(c("Residues"),  c(NORTH<-1),line=0.7, outer=TRUE, cex=0.8)
  return()}



peptide_pv_tp_proc<-function(df, dfup, pv, sd,nb_row=100, ranges=c(-Inf, seq(-30, 30, by=10), Inf), 
                             pv_cutoff=0.01, replicates=3){
  #preparation significance per residue & coverage
  cl1<-significant_peptide_uptake(dfup, pv, sd, pv_cutoff, replicates)
  
  start_col<-which(colnames(df)=='Start')
  end_col<-which(colnames(df)=='End')
  fc.d<-data.frame((df[,7]-df[,8:dim(df)[2]])/10*cl1)###vector which has significant average
  
  si.fv<-(fc.d)
  ranges_function(df,si.fv )
  
  
  xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
  for ( i in 1:(length(xli)-1)){
    si.fv[xli[i]< si.fv & si.fv < xli[i+1]] <- num_ass[i] }
  si.fv[si.fv==0]<- (-10000)
  
  si.fv<-abs(si.fv)-9999
  
  cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c("black"))
  
  
  
  for ( j in 1:dim(si.fv)[2]){
    y1=c(rep(1:nb_row, times=floor(dim(df)[1]/nb_row)), 1:(dim(df)[1]%%nb_row)) ##y values on the plot corresponding to peptides index
    plot(c(1,1), type="n", xlab="", ylab="", lwd=2, col=4, xlim=c(min(df[,start_col]), max(df[,end_col])), 
         ylim=c(nb_row, 0), yaxt="n", xaxt="n") ## mock plot, just to have it drawn correct limits set up
    
    for ( i in 1:dim(df)[1]){
      #print(c(df[i,start_col], df[i,end_col],y1[i], cbr1[si.fv[i,j]] ))
      points(c(df[i,start_col], df[i,end_col]), c(y1[i], y1[i]), type="l",
             col=cbr1[si.fv[i,j]])
      
    }
    
    main_nm=str_sub(colnames(df)[j+7], end = -9, start=4)
    mtext(main_nm, side=1, outer=FALSE, line=0, cex=0.6)
    axis(3, at=seq(0, max(df[, end_col])+25, by=25),  tcl=-0.2, labels = F)
    axis(3, at=seq(0, max(df[, end_col])+25, by=50), labels = T, cex.axis=0.65)
    box(lwd=2)
  }}

plot_peptide_sig_tp_proc<-function(input_proc, input_up, nb_pep_row=100, lim1=2, ranges=c(-Inf, seq(-30, 30, by=10), Inf), 
                                   pv_cutoff=0.01, replicates=3){
  pv1<-pv_timepoint(input_up, replicates)
  s1<-sd_timepoint(input_up, replicates)
  av1<-ave_timepoint(input_proc, replicates)
  avu<-ave_timepoint(input_up, replicates)
  for ( i in(unique(av1$Deut.Time))){
    print(paste("For time point", i))
    a1=av1[av1$Deut.Time==i,]
    au=avu[avu$Deut.Time==i,]
    p1=pv1[pv1$Deut.Time==i,]
    sd=s1[s1$Deut.Time==i,]
    peptide_pv_tp_proc(a1, au, p1, sd, nb_pep_row,ranges, pv_cutoff, replicates)
    #legend_heat_map_tp(av1)
    mtext(i, side=4, outer=FALSE, line=0.2, cex=0.7)}
  mtext(c("Residues"),  c(NORTH<-3),line=0, outer=TRUE, cex=0.7)}

###########





pymol_script_significant_residue_proc<-function(df,dfup, ranges=c(-Inf, seq(-30, 30, by=10), Inf), 
                                                pv_cutoff=0.01, replicates=3){
  #####from HDX get data and 
  
  for ( deut.time in(unique(df$Deut.Time))){ 
    pv<-pv_timepoint(dfup[dfup$Deut.Time==deut.time,], replicates = 3)
    sd<-sd_timepoint(dfup[dfup$Deut.Time==deut.time,], replicates = 3)
    dfu<-ave_timepoint(dfup[dfup$Deut.Time==deut.time,], replicates = 3)
    df1<-ave_timepoint(df[df$Deut.Time==deut.time,], replicates = 3)
    
    #preparation significance per residue & coverage
    cl1<-significant_peptide_uptake(dfu, pv, sd, pv_cutoff,replicates)
    
    start_col<-which(colnames(df1)=='Start')
    end_col<-which(colnames(df1)=='End')
    
    ###preparation of the coloring ranges. Difference of average. 
    fc.d<-data.frame((df1[,7]-df1[,8:dim(df1)[2]])/10*cl1)###vector which has significant average
    ###per residue maximum uptake value is chosen. -> 
    max.ac1<-c() 
    for ( j in 1:dim(fc.d)[2]){
      ac<-c()
      ac1<-c()
      ac2<-c()
      for ( i in 1:dim(df1)[1]){
        ac<-rep(0, length=max(df1[,end_col]))
        ac[df1[i,start_col]:df1[i,end_col]]=fc.d[i,j]
        ##make multiple vectors which have 1 at position which peptide covers
        ac1<-c(ac1, ac)}
      ac2=data.frame(matrix(ac1, nrow =dim(df1)[1], byrow=T))
      max.a<-c()
      for ( k in 1:dim(ac2)[2]){ 
        ind1<-which.max(abs(ac2[,k]))
        nb1<-(ac2[ind1,k])
        max.a<-c(max.a, nb1)}
      max.ac1<-c(max.ac1, max.a)}
    max.ac2=data.frame(matrix(max.ac1, ncol = dim(fc.d)[2]))
    
    
    ##prep of coverage
    coverage<-coverage_residue(df1, start_col, end_col)
    
    #####preparation of average per residue data.frame, which will have 
    xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
    for ( i in 1:(length(xli)-1)){
      max.ac2[xli[i]< max.ac2 & max.ac2 < xli[i+1]] <- num_ass[i]}
    max.ac2[max.ac2==0]<- (-10000)
    
    si_apc<-abs(max.ac2)-9999
    cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)###define lack of coverage
    si_apc<-si_apc*cv_mis+1
    
    cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c("white", "grey45"))
    
    
    ##assign name to colors in pallettes
    col_nm<-c("no_cov", "NSig")
    for ( i in 1:(length(ranges)-1)){
      col_nm<-c(col_nm, paste("col_", ranges[i],"_", ranges[i+1], sep=""))}
    
    rgb_col<-col2rgb(cbr1, alpha = FALSE) ## function to return rgb values for colors in pallette
    
    ##set_color 0%, [0 , 0 , 120], make command for pymol, to set colors
    set_colors<-c()
    for ( i in 1:length(col_nm)){
      set_colors<-c(set_colors, paste("set_color ", col_nm[i],", [",  rgb_col[1,i], 
                                      ",", rgb_col[2,i], ",", rgb_col[3,i], "]", sep=""))}
    ###write outputs per each state in the 
    nm1<-str_sub(colnames(df1[8:dim(df1)[2]]), start=4, end=-9)
    for (j in 1:dim(si_apc)[2]){
      output_name<-paste("pymol_maxuptake_proc_", nm1[j],"_", deut.time, ".txt", sep="")
      res.txt<-c()
      for ( i in 1:length(col_nm)){
        if (length(which(si_apc[,j]==i)) !=0){
          res.txt<-c(res.txt, (paste(c("color ", col_nm[i],", resi ",
                                       pymol_str(which(si_apc[,j] ==i))), sep="", collapse="")))}}
      fileConn<-file(output_name)
      writeLines(c("hide","show cartoon","color grey", "bg white",set_colors,res.txt ), fileConn)
      close(fileConn)}}
  
  leg_nm<-c("No coverage", "Not Sig")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],":", ranges[i+1], "%", sep=""))}
  
  pallette_ll(cbr1, leg_nm)
  return()}


pymol_script_significant_peptide_proc<-function(df,dfup, ranges=c(-Inf, seq(-30, 30, by=10), Inf), 
                                                pv_cutoff=0.01, replicates=3){
  #####from HDX get data and 
  for ( deut.time in(unique(df$Deut.Time))){
    pv<-pv_timepoint(dfup[dfup$Deut.Time==deut.time,], replicates)
    sd<-sd_timepoint(dfup[dfup$Deut.Time==deut.time,], replicates)
    dfu<-ave_timepoint(dfup[dfup$Deut.Time==deut.time,], replicates)
    df1<-ave_timepoint(df[df$Deut.Time==deut.time,], replicates)

    #preparation significance per residue & coverage
    cl1<-significant_peptide_uptake(dfu, pv, sd, pv_cutoff, replicates)
    
    start_col<-which(colnames(df1)=='Start')
    end_col<-which(colnames(df1)=='End')

    ###preparation of the coloring ranges. Difference of average. 
    fc.d<-data.frame((df1[,7]-df1[,8:dim(df1)[2]])/10*cl1)###vector which has significant average
    ##sumarized occurances of peptides
    si.f=fc.d
    #####preparation of average per residue data.frame, which will have 
    xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
    for ( i in 1:(length(xli)-1)){
      si.f[xli[i]<  si.f &  si.f < xli[i+1]] <- num_ass[i]}
    si.f[ si.f==0]<- (-10000)
    si_apc<-abs(si.f)-9999
    
    
    cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c( "grey45"))
    
    ##assign name to colors in pallettes
    col_nm<-c("NSig")
    for ( i in 1:(length(ranges)-1)){
      col_nm<-c(col_nm, paste("col_", ranges[i],"_", ranges[i+1], sep=""))}
    
    
    rgb_col<-col2rgb(cbr1, alpha = FALSE) ## function to return rgb values for colors in pallette
    
    ##set_color 0%, [0 , 0 , 120], make command for pymol, to set colors
    set_colors<-c()
    for ( i in 1:length(col_nm)){
      set_colors<-c(set_colors, paste("set_color ", col_nm[i],", [",  rgb_col[1,i], 
                                      ",", rgb_col[2,i], ",", rgb_col[3,i], "]", sep=""))}
    
    
    ###write outputs per each state in the 
    nm1<-str_sub(colnames(df1[8:dim(df1)[2]]), start=4, end=-9)
    for (j in 1:dim(si_apc)[2]){
      output_name<-paste("pymol_all_peptides_proc_", nm1[j],"_", deut.time, ".txt", sep="")
      
      res.txt<-c()
      for ( i in 1:length(col_nm)){
        if (length(which(si_apc[,j]==i)) !=0 ) {
          pep_nb<-which(si_apc[,j] == i)
          for ( k in pep_nb){line<-c()
          line<- paste(c("color ", col_nm[i],", resi ", df1[k,start_col], "-", df1[k,end_col] , sep=""))
          res.txt<-c(res.txt, 
                     paste(line, sep="' '", collapse=""))}
        }}
      fileConn<-file(output_name)
      writeLines(c("hide","show cartoon","color black", "bg white",
                   set_colors, res.txt ), fileConn)
      close(fileConn)}}
  leg_nm<-c("Not Sig")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],":", ranges[i+1], "%", sep=""))}
  
  pallette_ll(cbr1, leg_nm)
  return()}


####
