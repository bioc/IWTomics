IWTomicsTest <- function(regionsFeatures,
                         id_region1=idRegions(regionsFeatures)[1],id_region2=NULL, 
                         id_features_subset=idFeatures(regionsFeatures),
                         mu=0,statistics="mean",probs=0.5,max_scale=NULL,paired=FALSE,B=1000){
  if(class(regionsFeatures)!='IWTomicsData')
    stop('invalid regionsFeatures. IWTomicsData object expected.')
  if(!validObject(regionsFeatures))
    stop('invalid regionsFeatures.')
  if((length(id_region2)==1)&&(id_region2==''))
    id_region2=NULL
  regionsFeatures_tot=initialize(regionsFeatures,
                                 test=list(input=list(id_region1=id_region1,id_region2=id_region2,
                                                      id_features_subset=id_features_subset,mu=mu,
                                                      statistics=statistics,probs=probs,max_scale=max_scale,
                                                      paired=paired,B=B),
                                           result=list()))
  regionsFeatures=regionsFeatures[setdiff(union(id_region1,id_region2),''),id_features_subset]
  regionsFeatures=initialize(regionsFeatures,
                             test=list(input=list(id_region1=id_region1,id_region2=id_region2,
                                                  id_features_subset=id_features_subset,mu=mu,
                                                  statistics=statistics,probs=probs,max_scale=max_scale,
                                                  paired=paired,B=B),
                                       result=list()))
  if(!is.null(id_region2)){
    if(paired){
      if(sum(lengthRegions(regionsFeatures)[id_region1]!=lengthRegions(regionsFeatures)[id_region2],na.rm=TRUE))
        stop('cannot perform a paired test with different sample size in the two region datasets. ')
    }
  }
  
  regionsFeatures@test$result=vector("list",length(id_region1))
  for(i in seq_along(id_region1)){
    id_region1_i=id_region1[i]
    id_region2_i=id_region2[i]
    if((!is.null(id_region2_i))&&(id_region2_i==''))
      id_region2_i=NULL
    if(is.null(id_region2_i)){
      write2='...'
    }else{
      write2=paste0(' vs. \'',nameRegions(regionsFeatures[id_region2_i,]),'\'...')
    }
    message('Performing IWT for \'',nameRegions(regionsFeatures)[id_region1_i],'\'',write2)
    features_test_i=lapply(features(regionsFeatures),function(feature) feature[c(id_region1_i,id_region2_i)])
    features_test_i_not_NA=lapply(features_test_i,function(feature_test) rowSums(!is.na(Reduce(cbind,feature_test)))!=0)
    features_test_i=mapply(function(feature_test,features_test_i_not_NA) lapply(feature_test,function(feature) as.matrix(feature[features_test_i_not_NA,])),
                           features_test_i,features_test_i_not_NA,SIMPLIFY=FALSE)
    if(alignment(regionsFeatures)=='scale'){
      if(sum(unlist(lapply(features_test_i,function(feature_test) lapply(feature_test,function(feature) sum(is.na(feature[nrow(feature),])))))))
        stop('IWTomicsTest is incompatible with \'scale\' alignment and regions of different length. Smooth data first.')
    }
    if(sum((unlist(lapply(mu,length))!=1)&(unlist(lapply(mu,length))!=unlist(lapply(features_test_i,function(feature) nrow(feature[[1]]))))))
      stop('invalid mu. Grids different from feature grids.')
    mu_i=as.list(rep(mu,length.out=length(id_features_subset)))
    names(mu_i)=id_features_subset
    if(is.null(max_scale)){
      max_scale_i=lapply(features_test_i,function(feature_test) nrow(feature_test[[id_region1_i]]))
    }else{
      if(is.list(max_scale)){
        max_scale_i=as.list(rep(max_scale[[i]],length.out=length(id_features_subset)))
      }else{
        max_scale_i=as.list(rep(max_scale,length.out=length(id_features_subset)))
      }
      names(max_scale_i)=id_features_subset
    }
    regionsFeatures@test$result[[i]]=list()
    for(id_feature in id_features_subset){
      message('  Performing IWT for feature \'',nameFeatures(regionsFeatures)[id_feature],'\'...')
      p=nrow(features_test_i[[id_feature]][[id_region1_i]])
      if((max_scale_i[[id_feature]]<1)||(max_scale_i[[id_feature]]>p))
        warning('invalid max_scale. Setting it to the default value.',call.=FALSE,immediate.=TRUE)
      max_scale_i[[id_feature]]=ifelse(max_scale_i[[id_feature]]<1,p,min(max_scale_i[[id_feature]],p))
      
      if(is.null(id_region2_i)){
        regionsFeatures@test$result[[i]][[id_feature]]=.IWTomics.1pop(features_test_i[[id_feature]][[id_region1_i]],
                                                                      mu=mu_i[[id_feature]],statistics=statistics,probs=probs,
                                                                      max_scale=max_scale_i[[id_feature]],B=B)
        if(regionsFeatures@test$result[[i]][[id_feature]]$exact)
          warning('number of iteration B greater than number of possible permutations. Exact p-values computed.',call.=FALSE,immediate.=TRUE)
      }else{
        tmp=.IWTomics.2pop(features_test_i[[id_feature]][[id_region1_i]],features_test_i[[id_feature]][[id_region2_i]],
                           mu=mu_i[[id_feature]],statistics=statistics,probs=probs,max_scale=max_scale_i[[id_feature]],paired=paired,B=B)
        if(TRUE %in% tmp$no.pval)
          warning('p-value not fully computable in some points, because of too many NAs present.',call.=FALSE,immediate.=TRUE)
        regionsFeatures@test$result[[i]][[id_feature]]=tmp$result
        if(regionsFeatures@test$result[[i]][[id_feature]]$exact)
          warning('number of iteration B greater than number of possible permutations. Exact p-values computed.',call.=FALSE,immediate.=TRUE)
      }
      regionsFeatures@test$result[[i]][[id_feature]]$notNA=features_test_i_not_NA[[id_feature]]
    }
  }
  
  return(regionsFeatures)
}

.pval.correct <- function(pval.matrix,maxrow){
  p=ncol(pval.matrix)
  pval.matrix.2x=cbind(pval.matrix[,p:1],pval.matrix[,p:1])
  corrected.pval.matrix=matrix(nrow=p,ncol=p)
  corrected.pval.matrix[p,]=rev(pval.matrix[p,])
  for(var in 1:p){
    pval_var=pval.matrix.2x[p,var]
    start=var
    end=var
    for(row in (p-1):maxrow){
      end=end+1
      pval_cone=pval.matrix.2x[row,start:end]
      pval_var=suppressWarnings(max(pval_var,pval_cone,na.rm=TRUE))
      corrected.pval.matrix[row,var]=pval_var
    }
  }
  corrected.pval.matrix=corrected.pval.matrix[,p:1]
  return(corrected.pval.matrix)
}

.IWTomics.1pop <- function(data,mu=0,statistics='mean',probs=0.5,max_scale=nrow(data),B=1000,recycle=FALSE){
  # recycle     if TRUE, edges are recycled.
  # max_scale   max interval length for the adjustment.
  
  if(statistics=='median'){
    statistics='quantile'
    probs=0.5
  }
  max_scale=min(max_scale,nrow(data))
  
  result=list(test='1pop',mu=mu,max_scale=max_scale)
  n=ncol(data)
  p=nrow(data)
  data=data-mu
  exact=(B>=(2^n))
  
  # Univariate permutations
  message('    Point-wise tests...')
  if(statistics=='mean'){
    T0_plot=rowMeans(data,na.rm=TRUE)
    T0=(T0_plot)^2
  }
  if(statistics=='quantile'){
    T0_plot=apply(data,1,quantile,probs=probs,na.rm=TRUE)
    T0=(T0_plot)^2
    if(is.matrix(T0_plot)){
      T0_plot=colSums(T0_plot)
      T0=colSums(T0)
    }
  }
  if(statistics=='variance'){
    T0_plot=apply(data,1,var,na.rm=TRUE)
    T0=T0_plot
  }
  if(exact){
    T_perm=do.call(cbind,lapply(seq_len(n+1)-1,
                                function(m){
                                  signs_change=combn(n,m)
                                  T_perm=apply(signs_change,2,
                                               function(change){
                                                 data_perm=data
                                                 data_perm[,change]=-data_perm[,change]
                                                 if(statistics=='mean')
                                                   return((rowMeans(data_perm,na.rm=TRUE))^2)
                                                 if(statistics=='quantile'){
                                                   T_perm=apply(data_perm,1,quantile,probs=probs,na.rm=TRUE)^2
                                                   if(is.matrix(T_perm))
                                                     return(colSums(T_perm))
                                                   return(T_perm)
                                                   }
                                                if(statistics=='variance')
                                                  return(apply(data_perm,1,var,na.rm=TRUE))
                                               })
                                  return(T_perm)
                                }))
  }else{
    T_perm=do.call(cbind,lapply(seq.int(B-1),
                                function(perm){
                                  signs=rbinom(n,1,0.5)*2-1
                                  data_perm=data*rep(signs,each=nrow(data))
                                  if(statistics=='mean')
                                    return((rowMeans(data_perm,na.rm=TRUE))^2)
                                  if(statistics=='quantile'){
                                    T_perm=apply(data_perm,1,quantile,probs=probs,na.rm=TRUE)^2
                                    if(is.matrix(T_perm))
                                      return(colSums(T_perm))
                                    return(T_perm)
                                  }
                                  if(statistics=='variance')
                                    return(apply(data_perm,1,var,na.rm=TRUE))
                                }))
    T_perm=cbind(T_perm,T0)
  }
  if(statistics!='variance'){
    pval=rowSums(T_perm>=T0)/B
  }else{
    pval=pmin(2*rowSums(T_perm>=T0)/B,2*rowSums(T_perm<=T0)/B)
  }
  
  # Combination
  message('    Interval-wise tests...')
  # Asymmetric combination matrix:
  matrix_pval_asymm=matrix(nrow=p,ncol=p)
  matrix_pval_asymm[p,]=pval
  T0_2x=c(T0,T0)
  T_perm_2x=rbind(T_perm,T_perm)
  maxrow=p-max_scale+1
  #message('      Creating the p-value matrix:')
  if(recycle==TRUE){
    for(i in (p-1):maxrow){ # rows
      for(j in seq.int(p)){ # columns
        inf=j
        sup=(p-i)+j
        T0_temp=sum(T0_2x[inf:sup])
        T_temp=colSums(T_perm_2x[inf:sup,])
        matrix_pval_asymm[i,j]=ifelse(statistics!='variance',sum(T_temp>=T0_temp)/B,min(2*sum(T_temp>=T0_temp)/B,2*sum(T_temp<=T0_temp)/B))
      }
      #message('               end of row ',p-i+1,' out of ',p,'...')
    }
  }else{
    for(i in (p-1):maxrow){ # rows
      for(j in seq.int(i)){ # columns
        inf=j
        sup=(p-i)+j
        T0_temp=sum(T0_2x[inf:sup])
        T_temp=colSums(T_perm_2x[inf:sup,])
        matrix_pval_asymm[i,j]=ifelse(statistics!='variance',sum(T_temp>=T0_temp)/B,min(2*sum(T_temp>=T0_temp)/B,2*sum(T_temp<=T0_temp)/B))
      }
      #message('               end of row ',p-i+1,' out of ',p,'...')
    }
  }
  corrected.pval.matrix=.pval.correct(matrix_pval_asymm,maxrow)
  corrected.pval=corrected.pval.matrix[maxrow,]
  
  result$T0_plot=T0_plot
  result$adjusted_pval=corrected.pval
  result$adjusted_pval_matrix=corrected.pval.matrix
  result$unadjusted_pval=pval
  result$pval_matrix=matrix_pval_asymm
  result$exact=exact
  class(result)='ITWomics.1pop'
  return(result)
}

.IWTomics.2pop <- function(data1,data2,mu=0,statistics='mean',probs=0.5,max_scale=nrow(data1),paired=FALSE,B=1000,recycle=FALSE){
  # recycle     if TRUE, edges are recycled.
  # max_scale   max interval length for the adjustment.
  
  if(statistics=='median'){
    statistics='quantile'
    probs=0.5
  }
  n1=ncol(data1)
  n2=ncol(data2)
  n=n1+n2
  data1=data1-mu
  if(paired){
    exact=(B>=(2^n1))
  }else{
    exact=(B>=choose(n,n1))
  }
  # Non computable p-values (all NAs in one group)
  allNA=(rowSums(!is.na(data1))==0)|(rowSums(!is.na(data2))==0)
  data1=data1[!allNA,]
  data2=data2[!allNA,]
  max_scale=min(max_scale,nrow(data1))
  
  p=nrow(data1)
  result=list(test='2pop',mu=mu,max_scale=max_scale)
  data=cbind(data1,data2)
  
  # Univariate permutations
  message('    Point-wise tests...')
  if(statistics=='mean'){
    T0_plot=rowMeans(data1,na.rm=TRUE)-rowMeans(data2,na.rm=TRUE)
    T0=(T0_plot)^2
  }
  if(statistics=='quantile'){
    T0_plot=apply(data1,1,quantile,probs=probs,na.rm=TRUE)-apply(data2,1,quantile,probs=probs,na.rm=TRUE)
    T0=(T0_plot)^2
    if(is.matrix(T0_plot)){
      T0_plot=colSums(T0_plot)
      T0=colSums(T0)
    }
  }
  if(statistics=='variance'){
    T0_plot=apply(data1,1,var,na.rm=TRUE)/apply(data2,1,var,na.rm=TRUE)
    T0=T0_plot
  }
  if(exact){
    if(paired){
      T_perm=do.call(cbind,lapply(seq_len(n1+1)-1,
                                  function(m){
                                    group_change=combn(n1,m)
                                    T_perm=apply(group_change,2,
                                                 function(change){
                                                   data_perm=data
                                                   data_perm[,c(change,n1+change)]=data_perm[,c(n1+change,change)]
                                                   if(statistics=='mean')
                                                     return((rowMeans(as.matrix(data_perm[,seq.int(n1)]),na.rm=TRUE)-rowMeans(as.matrix(data_perm[,n1+seq.int(n2)]),na.rm=TRUE))^2)
                                                   if(statistics=='quantile'){
                                                     T_perm=(apply(as.matrix(data_perm[,seq.int(n1)]),1,quantile,probs=probs,na.rm=TRUE)-apply(as.matrix(data_perm[,n1+seq.int(n2)]),1,quantile,probs=probs,na.rm=TRUE))^2
                                                     if(is.matrix(T_perm))
                                                       return(colSums(T_perm))
                                                     return(T_perm)
                                                   }
                                                   if(statistics=='variance')
                                                     return((apply(as.matrix(data_perm[,seq.int(n1)]),1,var,na.rm=TRUE)/apply(as.matrix(data_perm[,n1+seq.int(n2)]),1,var,na.rm=TRUE)))
                                                 })
                                    return(T_perm)
                                  }))
    }else{
      first_group=combn(n,n1)
      T_perm=apply(first_group,2,
                   function(group){
                     data_perm=data[,c(group,setdiff(seq_len(n1+n2),group))]
                     if(statistics=='mean')
                       return((rowMeans(as.matrix(data_perm[,seq.int(n1)]),na.rm=TRUE)-rowMeans(as.matrix(data_perm[,n1+seq.int(n2)]),na.rm=TRUE))^2)
                     if(statistics=='quantile'){
                       T_perm=(apply(as.matrix(data_perm[,seq.int(n1)]),1,quantile,probs=probs,na.rm=TRUE)-apply(as.matrix(data_perm[,n1+seq.int(n2)]),1,quantile,probs=probs,na.rm=TRUE))^2
                       if(is.matrix(T_perm))
                         return(colSums(T_perm))
                       return(T_perm)
                     }
                     if(statistics=='variance')
                       return((apply(as.matrix(data_perm[,seq.int(n1)]),1,var,na.rm=TRUE)/apply(as.matrix(data_perm[,n1+seq.int(n2)]),1,var,na.rm=TRUE)))
                   })
    }
  }else{
    T_perm=do.call(cbind,lapply(seq.int(B-1),
                                function(perm){
                                  if(paired){
                                    couple.perm=rbinom(n1,1,0.5)
                                    data_perm=data[,c(n1*couple.perm,-n1*couple.perm)+seq.int(2*n1)]
                                  }else{
                                    permutation=sample(n,n1)
                                    data_perm=data[,c(permutation,setdiff(seq_len(n),permutation))]
                                  }
                                  if(statistics=='mean')
                                    return((rowMeans(as.matrix(data_perm[,seq.int(n1)]),na.rm=TRUE)-rowMeans(as.matrix(data_perm[,n1+seq.int(n2)]),na.rm=TRUE))^2)
                                  if(statistics=='quantile'){
                                    T_perm=(apply(as.matrix(data_perm[,seq.int(n1)]),1,quantile,probs=probs,na.rm=TRUE)-apply(as.matrix(data_perm[,n1+seq.int(n2)]),1,quantile,probs=probs,na.rm=TRUE))^2
                                    if(is.matrix(T_perm))
                                      return(colSums(T_perm))
                                    return(T_perm)
                                  }
                                  if(statistics=='variance')
                                    return((apply(as.matrix(data_perm[,seq.int(n1)]),1,var,na.rm=TRUE)/apply(as.matrix(data_perm[,n1+seq.int(n2)]),1,var,na.rm=TRUE)))
                                }))
    T_perm=cbind(T_perm,T0)
  }
  
  # Not fully computable p-values (some permutations do not produce any test statistics because of the NAs)
  if(statistics=='variance'){
    no.pval=rowSums(is.na(data))>=(min(n1,n2)-1)
  }else{
    no.pval=rowSums(is.na(data))>=min(n1,n2)
  }
  #T_perm[no.pval,]=NaN # do not compute any p-value when it is not fully computable
  # do not compute any p-value if NaN is in T_perm
  #if(statistics!='variance'){
  #  pval=rowSums(T_perm>=T0)/B
  #}else{
  #  pval=pmin(2*rowSums(T_perm>=T0)/B,2*rowSums(T_perm<=T0)/B)
  #}
  # compute p-value omitting NaN
  if(statistics!='variance'){
    pval=rowSums(T_perm>=T0,na.rm=TRUE)/rowSums(!is.nan(T_perm))
  }else{
    pval=pmin(2*rowSums(T_perm>=T0,na.rm=TRUE)/rowSums(!is.nan(T_perm)),2*rowSums(T_perm<=T0,na.rm=TRUE)/rowSums(!is.nan(T_perm)))
  }
  
  # Combination
  message('    Interval-wise tests...')
  # Asymmetric combination matrix:
  matrix_pval_asymm=matrix(nrow=p,ncol=p)
  matrix_pval_asymm[p,]=pval
  T0_2x=c(T0,T0)
  T_perm_2x=rbind(T_perm,T_perm)
  maxrow=p-max_scale+1
  #message('      Creating the p-value matrix:')
  if(recycle==TRUE){
    for(i in (p-1):maxrow){ # rows
      for(j in seq.int(p)){ # columns
        inf=j
        sup=(p-i)+j
        T0_temp=sum(T0_2x[inf:sup])
        T_temp=colSums(T_perm_2x[inf:sup,])
        # do not compute any p-value if NaN is in T_temp
        #matrix_pval_asymm[i,j]=ifelse(statistics!='variance',sum(T_temp>=T0_temp)/B,min(2*sum(T_temp>=T0_temp)/B,2*sum(T_temp<=T0_temp)/B))
        # compute p-value omitting NaN
        matrix_pval_asymm[i,j]=ifelse(statistics!='variance',sum(T_temp>=T0_temp,na.rm=TRUE)/sum(!is.nan(T_temp)),min(2*sum(T_temp>=T0_temp,na.rm=TRUE)/sum(!is.nan(T_temp)),2*sum(T_temp<=T0_temp,na.rm=TRUE)/sum(!is.nan(T_temp))))
      }
      #message('               end of row ',p-i+1,' out of ',p,'...')
    }
  }else{
    for(i in (p-1):maxrow){ # rows
      for(j in seq.int(i)){ # columns
        inf=j
        sup=(p-i)+j
        T0_temp=sum(T0_2x[inf:sup])
        T_temp=colSums(T_perm_2x[inf:sup,])
        # do not compute any p-value if NaN is in T_temp
        #matrix_pval_asymm[i,j]=ifelse(statistics!='variance',sum(T_temp>=T0_temp)/B,min(2*sum(T_temp>=T0_temp)/B,2*sum(T_temp<=T0_temp)/B))
        # compute p-value omitting NaN
        matrix_pval_asymm[i,j]=ifelse(statistics!='variance',sum(T_temp>=T0_temp,na.rm=TRUE)/sum(!is.nan(T_temp)),min(2*sum(T_temp>=T0_temp,na.rm=TRUE)/sum(!is.nan(T_temp)),2*sum(T_temp<=T0_temp,na.rm=TRUE)/sum(!is.nan(T_temp))))
      }
      #message('               end of row ',p-i+1,' out of ',p,'...')
    }
  }
  corrected.pval.matrix=.pval.correct(matrix_pval_asymm,maxrow)
  corrected.pval=corrected.pval.matrix[maxrow,]
  
  result$T0_plot=rep(NA,length(allNA))
  result$T0_plot[!allNA]=T0_plot
  result$adjusted_pval=rep(NA,length(allNA))
  result$adjusted_pval[!allNA]=corrected.pval
  result$adjusted_pval_matrix=matrix(NA,nrow=length(allNA),ncol=length(allNA))
  result$adjusted_pval_matrix[(sum(allNA)+1):length(allNA),!allNA]=corrected.pval.matrix
  result$unadjusted_pval=rep(NA,length(allNA))
  result$unadjusted_pval[!allNA]=pval
  result$pval_matrix=matrix(NA,nrow=length(allNA),ncol=length(allNA))
  result$pval_matrix[(sum(allNA)+1):length(allNA),!allNA]=matrix_pval_asymm
  result$exact=exact
  class(result)='ITWomics.2pop'
  return(list(result=result,no.pval=no.pval))
}
