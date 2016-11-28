setGeneric("smooth",function(x,...) standardGeneric("smooth"))
setMethod("smooth","IWTomicsData",
          function(x,type='locpoly',
                   id_regions_subset=idRegions(x),id_features_subset=idFeatures(x),
                   resolution_new=resolution(x)[id_features_subset],
                   scale_grid=unlist(lapply(lengthFeatures(x)[id_features_subset],function(length) max(unlist(length[id_regions_subset])))),
                   fill_gaps=TRUE,bandwidth=5,degree=3,dist_knots=10,parallel=FALSE){
  if(!(type %in% c('locpoly','kernel','splines')))
    stop('invalid smoothing type \'',type,'\'. Available types are \'locpoly\', \'kernel\' and \'splines\'.')
  if(sum(!(id_regions_subset %in% idRegions(x))))
    stop('invalid id_regions_subset. The region datasets provided are not listed in x.')
  if(sum(!(id_features_subset %in% idFeatures(x))))
    stop('invalid id_features_subset. The features provided are not listed in x.')
  
  if(parallel){
    core.number <- detectCores()
    worker.number <- core.number-1
  }
  
  if(alignment(x)=='scale'){
    scale_grid=rep(scale_grid,length.out=length(id_features_subset))
    names(scale_grid)=id_features_subset
    if(length(setdiff(idRegions(x),id_regions_subset))>0){
      length=lapply(lengthFeatures(x)[id_features_subset],function(length) unique(unlist(length)))
      if(TRUE %in% (unlist(lapply(length,length)>1))){
        stop('with alignment=\'scale\' all the regions should be scaled to the same length. Run the function over all 
             region datasets to scale regions.')
      }else{
        if(TRUE %in% (scale_grid!=unlist(length)))
          stop('each feature should be measured on the same grid in all the regions. Run the function over all 
               region datasets to change feature grid.')
      }
    }
    length.old=lapply(lengthFeatures(x)[id_features_subset],function(length) unique(unlist(length)))
    if(parallel){
      cl <- makeCluster(worker.number)
      for(id_feature in id_features_subset){
        for(id_region in id_regions_subset){
          feature_NA=is.na(features(x)[[id_feature]][[id_region]])
          if(sum(feature_NA)==0){
            abscissa=seq_len(nrow(feature_NA))
            abscissa_new=seq(0.5,nrow(feature_NA)+0.5,length.out=scale_grid[id_feature]+1)
            abscissa_new=(abscissa_new[-length(abscissa_new)]+abscissa_new[-1])/2
            length_smooth=rep(length(abscissa_new),ncol(feature_NA))
            if(type=='locpoly')
              feature_smooth=parApply(cl,features(x)[[id_feature]][[id_region]],2,.locpoly_single,
                                      abscissa=abscissa,abscissa_new=abscissa_new,degree=degree,bandwidth=bandwidth)
            if(type=='kernel')
              feature_smooth=parApply(cl,features(x)[[id_feature]][[id_region]],2,.ksmooth_single,
                                      abscissa=abscissa,abscissa_new=abscissa_new,bandwidth=bandwidth)
            if(type=='splines')
              feature_smooth=.splines_single(features(x)[[id_feature]][[id_region]],
                                             abscissa=abscissa,abscissa_new=abscissa_new,degree=degree,dist_knots=dist_knots)
          } else {
            abscissa=lapply(x@length_features[[id_feature]][[id_region]],function(length) seq_len(length))
            data=mapply(function(feature,abscissa,scale_grid){
                          y=feature[abscissa]
                          abscissa_new=seq(0.5,length(abscissa)+0.5,length.out=scale_grid+1)
                          abscissa_new=(abscissa_new[-length(abscissa_new)]+abscissa_new[-1])/2
                          abscissa=abscissa[!is.na(y)]
                          y=feature[abscissa]
                          return(list(abscissa=abscissa,y=y,abscissa_new=abscissa_new))
                        },split(features(x)[[id_feature]][[id_region]],col(features(x)[[id_feature]][[id_region]])),abscissa,
                        MoreArgs=list(scale_grid=scale_grid[id_feature]),SIMPLIFY=FALSE)
            abscissa=lapply(data,function(data) data$abscissa)
            abscissa_new=lapply(data,function(data) data$abscissa_new)
            y=lapply(data,function(data) data$y)
            length_smooth=unlist(lapply(abscissa_new,length))
            if(type=='locpoly')
              feature_smooth=Reduce(cbind,clusterMap(cl,.locpoly_single,abscissa,y,abscissa_new,degree=degree,bandwidth=bandwidth,SIMPLIFY=FALSE))
            if(type=='kernel')
              feature_smooth=Reduce(cbind,clusterMap(cl,.ksmooth_single,abscissa,y,abscissa_new,bandwidth=bandwidth,SIMPLIFY=FALSE))
            if(type=='splines')
              feature_smooth=Reduce(cbind,clusterMap(cl,.splines_single,abscissa,y,abscissa_new,degree=degree,dist_knots=dist_knots,SIMPLIFY=FALSE))
          }
          colnames(feature_smooth)=NULL
          x@features[[id_feature]][[id_region]]=feature_smooth
          names(length_smooth)=NULL
          x@length_features[[id_feature]][[id_region]]=length_smooth
        }
      }
      stopCluster(cl)
    } else {
      for(id_feature in id_features_subset){
        for(id_region in id_regions_subset){
          feature_NA=is.na(features(x)[[id_feature]][[id_region]])
          if(sum(feature_NA)==0){
            abscissa=seq_len(nrow(feature_NA))
            abscissa_new=seq(0.5,nrow(feature_NA)+0.5,length.out=scale_grid[id_feature]+1)
            abscissa_new=(abscissa_new[-length(abscissa_new)]+abscissa_new[-1])/2
            length_smooth=rep(length(abscissa_new),ncol(feature_NA))
            if(type=='locpoly')
              feature_smooth=apply(features(x)[[id_feature]][[id_region]],2,.locpoly_single,
                                   abscissa=abscissa,abscissa_new=abscissa_new,degree=degree,bandwidth=bandwidth)
            if(type=='kernel')
              feature_smooth=apply(features(x)[[id_feature]][[id_region]],2,.ksmooth_single,
                                   abscissa=abscissa,abscissa_new=abscissa_new,bandwidth=bandwidth)
            if(type=='splines')
              feature_smooth=.splines_single(features(x)[[id_feature]][[id_region]],
                                             abscissa=abscissa,abscissa_new=abscissa_new,degree=degree,dist_knots=dist_knots)
          } else {
            abscissa=lapply(x@length_features[[id_feature]][[id_region]],function(length) seq_len(length))
            data=mapply(function(feature,abscissa,scale_grid){
                          y=feature[abscissa]
                          abscissa_new=seq(0.5,length(abscissa)+0.5,length.out=scale_grid+1)
                          abscissa_new=(abscissa_new[-length(abscissa_new)]+abscissa_new[-1])/2
                          abscissa=abscissa[!is.na(y)]
                          y=feature[abscissa]
                          return(list(abscissa=abscissa,y=y,abscissa_new=abscissa_new))
                        },split(features(x)[[id_feature]][[id_region]],col(features(x)[[id_feature]][[id_region]])),abscissa,
                        MoreArgs=list(scale_grid=scale_grid[id_feature]),SIMPLIFY=FALSE)
            abscissa=lapply(data,function(data) data$abscissa)
            abscissa_new=lapply(data,function(data) data$abscissa_new)
            y=lapply(data,function(data) data$y)
            length_smooth=unlist(lapply(abscissa_new,length))
            if(type=='locpoly')
              feature_smooth=Reduce(cbind,mapply(.locpoly_single,abscissa,y,abscissa_new,degree=degree,bandwidth=bandwidth,SIMPLIFY=FALSE))
            if(type=='kernel')
              feature_smooth=Reduce(cbind,mapply(.ksmooth_single,abscissa,y,abscissa_new,bandwidth=bandwidth,SIMPLIFY=FALSE))
            if(type=='splines')
              feature_smooth=Reduce(cbind,mapply(.splines_single,abscissa,y,abscissa_new,degree=degree,dist_knots=dist_knots,SIMPLIFY=FALSE))
          }
          colnames(feature_smooth)=NULL
          x@features[[id_feature]][[id_region]]=feature_smooth
          names(length_smooth)=NULL
          x@length_features[[id_feature]][[id_region]]=length_smooth
        }
      }
    }
    index=unlist(lapply(length.old,length))==1
    x@metadata$feature_datasets[id_features_subset[index],'resolution']=resolution(x)[id_features_subset[index]]*
      unlist(length.old[index])/scale_grid[index]
    x@metadata$feature_datasets[id_features_subset[!index],'resolution']=NA
  } else {
    resolution_new=rep(resolution_new,length.out=length(id_features_subset))
    names(resolution_new)=id_features_subset
    if(length(setdiff(idRegions(x),id_regions_subset))>0){
      if(TRUE %in% (resolution_new!=resolution(x)[id_features_subset]))
        stop('each feature should be measured at the same resolution in all the regions. Run the function over all 
             region datasets to change feature resolutions.')
    }
    if(parallel){
      cl <- makeCluster(worker.number)
      for(id_feature in id_features_subset){
        by=resolution_new[id_feature]/(resolution(x)[id_feature])
        fill_gaps_feature=fill_gaps|(by!=1)
        for(id_region in id_regions_subset){
          feature_NA=is.na(features(x)[[id_feature]][[id_region]])
          if(sum(feature_NA)==0){
            abscissa=seq_len(nrow(feature_NA))
            abscissa_new=seq(0.5,nrow(feature_NA)+0.5,by)
            abscissa_new=(abscissa_new[-length(abscissa_new)]+abscissa_new[-1])/2
            length_smooth=rep(length(abscissa_new),ncol(feature_NA))
            if(type=='locpoly')
              feature_smooth=parApply(cl,features(x)[[id_feature]][[id_region]],2,.locpoly_single,
                                      abscissa=abscissa,abscissa_new=abscissa_new,degree=degree,bandwidth=bandwidth)
            if(type=='kernel')
              feature_smooth=parApply(cl,features(x)[[id_feature]][[id_region]],2,.ksmooth_single,
                                      abscissa=abscissa,abscissa_new=abscissa_new,bandwidth=bandwidth)
            if(type=='splines')
              feature_smooth=.splines_single(features(x)[[id_feature]][[id_region]],
                                             abscissa=abscissa,abscissa_new=abscissa_new,degree=degree,dist_knots=dist_knots)
          } else {
            if(alignment(x)=='left')
              abscissa=lapply(x@length_features[[id_feature]][[id_region]],function(length) seq_len(length))
            if(alignment(x)=='right')
              abscissa=lapply(x@length_features[[id_feature]][[id_region]],function(length,max.length) seq(max.length-length+1,max.length),nrow(features(x)[[id_feature]][[id_region]]))
            if(alignment(x)=='center')
              abscissa=lapply(x@length_features[[id_feature]][[id_region]],function(length,max.length) max.length%/%2-length%/%2+seq_len(length),nrow(features(x)[[id_feature]][[id_region]]))
            data=mapply(function(feature,abscissa,by){
                          y=feature[abscissa]
                          abscissa_new=seq(abscissa[1]-0.5,abscissa[length(abscissa)]+0.5,by)
                          abscissa_new=(abscissa_new[-length(abscissa_new)]+abscissa_new[-1])/2
                          abscissa=abscissa[!is.na(y)]
                          y=feature[abscissa]
                          return(list(abscissa=abscissa,y=y,abscissa_new=abscissa_new))
                        },split(features(x)[[id_feature]][[id_region]],col(features(x)[[id_feature]][[id_region]])),abscissa,
                        MoreArgs=list(by),SIMPLIFY=FALSE)
            abscissa=lapply(data,function(data) data$abscissa)
            abscissa_new=lapply(data,function(data) data$abscissa_new)
            y=lapply(data,function(data) data$y)
            length_smooth=unlist(lapply(abscissa_new,length))
            if(type=='locpoly')
              feature_smooth=clusterMap(cl,.locpoly_single,abscissa,y,abscissa_new,degree=degree,bandwidth=bandwidth,SIMPLIFY=FALSE)
            if(type=='kernel')
              feature_smooth=clusterMap(cl,.ksmooth_single,abscissa,y,abscissa_new,bandwidth=bandwidth,SIMPLIFY=FALSE)
            if(type=='splines'){
              feature_smooth=clusterMap(cl,.splines_single,abscissa,y,abscissa_new,degree=degree,dist_knots=dist_knots,SIMPLIFY=FALSE)
              feature_smooth=lapply(feature_smooth,as.vector)
            }
            if(!fill_gaps_feature)
              feature_smooth=mapply(function(feature_smooth,abscissa_new,abscissa){feature_smooth[abscissa_new %in% setdiff(abscissa_new,abscissa)]=NA;return(feature_smooth)},
                                    feature_smooth,abscissa_new,abscissa,SIMPLIFY=FALSE)
          }
          names(feature_smooth)=NULL
          x@features[[id_feature]][[id_region]]=feature_smooth
          names(length_smooth)=NULL
          x@length_features[[id_feature]][[id_region]]=length_smooth
        }
        locations_list_id=id_regions_subset[!unlist(lapply(features(x)[[id_feature]][id_regions_subset],is.matrix))]
        if(length(locations_list_id)>0){
          length=x@length_features[[id_feature]][locations_list_id]
          if(length(unique(unlist(length)))==1){
            x@features[[id_feature]][locations_list_id]=lapply(x@features[[id_feature]][locations_list_id],function(feature) do.call(cbind,feature))
          }else{
            length.max=max(unlist(length))
            if(alignment(x) %in% c('left','scale')){
              length.NA.right=lapply(length,function(length) length.max-length)
              for(id_region in locations_list_id)
                x@features[[id_feature]][[id_region]]=mapply(function(feature,length.NA.right) c(feature,rep.int(NA,length.NA.right)),
                                                             features(x)[[id_feature]][[id_region]],length.NA.right[[id_region]])
            }
            if(alignment(x)=='right'){
              length.NA.left=lapply(length,function(length) length.max-length)
              for(id_region in locations_list_id)
                x@features[[id_feature]][[id_region]]=mapply(function(feature,length.NA.left) c(rep.int(NA,length.NA.left),feature),
                                                             features(x)[[id_feature]][[id_region]],length.NA.left[[id_region]])
            }
            if(alignment(x)=='center'){
              center=lapply(length,'%/%',2) 
              # the alignment is approximate if there are regions with an odd number of windows and regions with an even number of regions
              length.NA.right=lapply(mapply('-',length,center,SIMPLIFY=FALSE),function(lenght_center) (length.max-length.max%/%2)-lenght_center)
              length.NA.left=lapply(center,function(center) length.max%/%2-center)
              for(id_region in locations_list_id)
                x@features[[id_feature]][[id_region]]=mapply(function(feature,length.NA.left,length.NA.right) as.matrix(c(rep.int(NA,length.NA.left),feature,rep.int(NA,length.NA.right))),
                                                             features(x)[[id_feature]][[id_region]],length.NA.left[[id_region]],length.NA.right[[id_region]])
            }
          }
        }
      }
      stopCluster(cl)
    } else {
      for(id_feature in id_features_subset){
        by=resolution_new[id_feature]/(resolution(x)[id_feature])
        fill_gaps_feature=fill_gaps|(by!=1)
        for(id_region in id_regions_subset){
          feature_NA=is.na(features(x)[[id_feature]][[id_region]])
          if(sum(feature_NA)==0){
            abscissa=seq_len(nrow(feature_NA))
            abscissa_new=seq(0.5,nrow(feature_NA)+0.5,by)
            abscissa_new=(abscissa_new[-length(abscissa_new)]+abscissa_new[-1])/2
            length_smooth=rep(length(abscissa_new),ncol(feature_NA))
            if(type=='locpoly')
              feature_smooth=apply(features(x)[[id_feature]][[id_region]],2,.locpoly_single,
                                   abscissa=abscissa,abscissa_new=abscissa_new,degree=degree,bandwidth=bandwidth)
            if(type=='kernel')
              feature_smooth=apply(features(x)[[id_feature]][[id_region]],2,.ksmooth_single,
                                   abscissa=abscissa,abscissa_new=abscissa_new,bandwidth=bandwidth)
            if(type=='splines')
              feature_smooth=.splines_single(features(x)[[id_feature]][[id_region]],
                                             abscissa=abscissa,abscissa_new=abscissa_new,degree=degree,dist_knots=dist_knots)
          } else {
            if(alignment(x)=='left')
              abscissa=lapply(x@length_features[[id_feature]][[id_region]],function(length) seq_len(length))
            if(alignment(x)=='right')
              abscissa=lapply(x@length_features[[id_feature]][[id_region]],function(length,max.length) seq(max.length-length+1,max.length),nrow(features(x)[[id_feature]][[id_region]]))
            if(alignment(x)=='center')
              abscissa=lapply(x@length_features[[id_feature]][[id_region]],function(length,max.length) max.length%/%2-length%/%2+seq_len(length),nrow(features(x)[[id_feature]][[id_region]]))
            data=mapply(function(feature,abscissa,by){
                          y=feature[abscissa]
                          abscissa_new=seq(abscissa[1]-0.5,abscissa[length(abscissa)]+0.5,by)
                          abscissa_new=(abscissa_new[-length(abscissa_new)]+abscissa_new[-1])/2
                          abscissa=abscissa[!is.na(y)]
                          y=feature[abscissa]
                          return(list(abscissa=abscissa,y=y,abscissa_new=abscissa_new))
                        },split(features(x)[[id_feature]][[id_region]],col(features(x)[[id_feature]][[id_region]])),abscissa,
                        MoreArgs=list(by),SIMPLIFY=FALSE)
            abscissa=lapply(data,function(data) data$abscissa)
            abscissa_new=lapply(data,function(data) data$abscissa_new)
            y=lapply(data,function(data) data$y)
            length_smooth=unlist(lapply(abscissa_new,length))
            if(type=='locpoly')
              feature_smooth=mapply(.locpoly_single,abscissa,y,abscissa_new,degree=degree,bandwidth=bandwidth,SIMPLIFY=FALSE)
            if(type=='kernel')
              feature_smooth=mapply(.ksmooth_single,abscissa,y,abscissa_new,bandwidth=bandwidth,SIMPLIFY=FALSE)
            if(type=='splines'){
              feature_smooth=mapply(.splines_single,abscissa,y,abscissa_new,degree=degree,dist_knots=dist_knots,SIMPLIFY=FALSE)
              feature_smooth=lapply(feature_smooth,as.vector)
            }
            if(!fill_gaps_feature)
              feature_smooth=mapply(function(feature_smooth,abscissa_new,abscissa){feature_smooth[abscissa_new %in% setdiff(abscissa_new,abscissa)]=NA;return(feature_smooth)},
                                    feature_smooth,abscissa_new,abscissa,SIMPLIFY=FALSE)
          }
          names(feature_smooth)=NULL
          x@features[[id_feature]][[id_region]]=feature_smooth
          names(length_smooth)=NULL
          x@length_features[[id_feature]][[id_region]]=length_smooth
        }
        locations_list_id=id_regions_subset[!unlist(lapply(features(x)[[id_feature]][id_regions_subset],is.matrix))]
        if(length(locations_list_id)>0){
          length=x@length_features[[id_feature]][locations_list_id]
          if(length(unique(unlist(length)))==1){
            x@features[[id_feature]][locations_list_id]=lapply(x@features[[id_feature]][locations_list_id],function(feature) do.call(cbind,feature))
          }else{
            length.max=max(unlist(length))
            if(alignment(x) %in% c('left','scale')){
              length.NA.right=lapply(length,function(length) length.max-length)
              for(id_region in locations_list_id)
                x@features[[id_feature]][[id_region]]=mapply(function(feature,length.NA.right) c(feature,rep.int(NA,length.NA.right)),
                                                             features(x)[[id_feature]][[id_region]],length.NA.right[[id_region]])
            }
            if(alignment(x)=='right'){
              length.NA.left=lapply(length,function(length) length.max-length)
              for(id_region in locations_list_id)
                x@features[[id_feature]][[id_region]]=mapply(function(feature,length.NA.left) c(rep.int(NA,length.NA.left),feature),
                                                             features(x)[[id_feature]][[id_region]],length.NA.left[[id_region]])
            }
            if(alignment(x)=='center'){
              center=lapply(length,'%/%',2) 
              # the alignment is approximate if there are regions with an odd number of windows and regions with an even number of regions
              length.NA.right=lapply(mapply('-',length,center,SIMPLIFY=FALSE),function(lenght_center) (length.max-length.max%/%2)-lenght_center)
              length.NA.left=lapply(center,function(center) length.max%/%2-center)
              for(id_region in locations_list_id)
                x@features[[id_feature]][[id_region]]=mapply(function(feature,length.NA.left,length.NA.right) as.matrix(c(rep.int(NA,length.NA.left),feature,rep.int(NA,length.NA.right))),
                                                             features(x)[[id_feature]][[id_region]],length.NA.left[[id_region]],length.NA.right[[id_region]])
            }
          }
        }
      }
    }
    x@metadata$feature_datasets[id_features_subset,'resolution']=resolution_new
  }
  x@test=NULL
  return(x)
})

.ksmooth_single <- function(abscissa,y,abscissa_new,bandwidth=5){
  y.new=ksmooth(abscissa,y,kernel='normal',bandwidth=bandwidth,x.points=abscissa_new,range.x=range(abscissa_new))$y
  return(y.new)
}

.locpoly_single <- function(abscissa,y,abscissa_new,degree=3,bandwidth=5){
  y.new=locpoly(abscissa,y,drv=0,degree=degree,bandwidth=bandwidth,gridsize=length(abscissa_new),range.x=range(abscissa_new))$y
  return(y.new)
}

.splines_single <- function(abscissa,y,abscissa_new,degree=3,dist_knots=10){
  breaks=seq(min(c(abscissa,abscissa_new)),max(c(abscissa,abscissa_new)),length.out=floor(max(c(abscissa,abscissa_new))/dist_knots)+1)
  basis=create.bspline.basis(range(c(abscissa,abscissa_new)),norder=degree+1,breaks=breaks)
  basis_fitted=smooth.basis(argvals=abscissa,y=y,basis)
  y.new=eval.fd(abscissa_new,basis_fitted$fd)
  return(y.new)
}