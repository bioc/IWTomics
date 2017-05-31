plotTest <- function(regionsFeatures,alpha=0.05,scale_threshold=NULL,nlevel=100,type='boxplot',
                     N_regions=pmin(lengthRegions(regionsFeatures),10),
                     probs=c(0.25,0.5,0.75),average=TRUE,size=TRUE,
                     id_features_subset=idFeatures(regionsFeatures),
                     col=1+seq_len(nRegions(regionsFeatures)),
                     ask=TRUE,xlab='Windows',ylim=NULL,...){
  if(class(regionsFeatures)!='IWTomicsData')
    stop('invalid regionsFeatures. IWTomicsData object expected.')
  if(!validObject(regionsFeatures))
    stop('invalid regionsFeatures.')
  if(is.null(regionsFeatures@test))
    stop('No test results present in regionFeatures.')
  if(sum(!(id_features_subset %in% idFeatures(regionsFeatures))))
    stop('invalid id_features_subset. The features provided are not listed in regionFeatures.')
  if(!(type %in% c('curves','boxplot')))
    stop('invalid plot type \'',type,'\'. Available types are \'curves\' and \'boxplot\'.')
  ylim_ext=ylim
  col=rep(col,length.out=nRegions(regionsFeatures))
  names(col)=idRegions(regionsFeatures)
  N_regions=rep(N_regions,length.out=nRegions(regionsFeatures))
  names(N_regions)=idRegions(regionsFeatures)
  
  devAskNewPage(ask)
  for(i in seq_len(nTests(regionsFeatures))){
    if(is.null(scale_threshold)){
      scale_threshold_i=lapply(.testResults(regionsFeatures)[[i]][id_features_subset],function(feature) feature$max_scale)
    }else{
      if(is.list(scale_threshold)){
        scale_threshold_i=as.list(rep(scale_threshold[[i]],length.out=length(id_features_subset)))
      }else{
        scale_threshold_i=as.list(rep(scale_threshold,length.out=length(id_features_subset)))
      }
      names(scale_threshold_i)=id_features_subset
    }
    id_region1_i=idRegionsTest(regionsFeatures,i)[[1]][1]
    id_region2_i=idRegionsTest(regionsFeatures,i)[[1]][2]
    if((!is.null(id_region2_i))&&(id_region2_i==''))
      id_region2_i=NULL
    for(id_feature in id_features_subset){
      if(.testResults(regionsFeatures)[[i]][[id_feature]]$test=='1pop'){
        plot_data=plot(regionsFeatures,type=type,
                       N_regions=N_regions[id_region1_i],probs=probs,average=average,size=size,
                       id_regions_subset=id_region1_i,id_features_subset=id_feature,
                       col=col[id_region1_i],plot=FALSE)
      }else{
        plot_data=plot(regionsFeatures,type=type,
                       N_regions=N_regions[c(id_region1_i,id_region2_i)],probs=probs,average=average,size=size,
                       id_regions_subset=c(id_region1_i,id_region2_i),id_features_subset=id_feature,
                       col=col[c(id_region1_i,id_region2_i)],plot=FALSE)
      }
      if(size){
        layout(rbind(1:2,cbind(3:6,0)),widths=c(8,2),heights=c(4,2,2,1,1))
        mar.left=7
      }else{
        layout(rbind(1:2,cbind(3:4,0)),widths=c(8,2),heights=c(4,2,2))
        mar.left=4
      }
      p=ncol(.testResults(regionsFeatures)[[i]][[id_feature]]$adjusted_pval_matrix)
      max_scale=.testResults(regionsFeatures)[[i]][[id_feature]]$max_scale
      if((scale_threshold_i[[id_feature]]<1)||(scale_threshold_i[[id_feature]]>p)){
        warning('invalid scale_threshold. Setting it to the default value.',call.=FALSE,immediate.=TRUE)
        scale_threshold_i[[id_feature]]=max_scale
      }
      
      # Adjusted p-value heatmap
      par(mar=c(4,mar.left,4,1.5),mgp=c(2.5,1,0),xpd=FALSE,las=1)
      col_heatmap=rainbow(nlevel,start=0.15,end=0.67)[nlevel:1]
      x_plot=plot_data$x_plot[[id_feature]]
      image(x_plot,1:max_scale-0.5,t(.testResults(regionsFeatures)[[i]][[id_feature]]$adjusted_pval_matrix[p:(p-max_scale+1),]), 
            col=col_heatmap,ylab="Maximum interval length",main=paste(nameFeatures(regionsFeatures)[id_feature],"\nAdjusted p-value heatmap",sep='\n'), 
            xlab=xlab,zlim=c(0,1),...)
      abline(h=scale_threshold_i[[id_feature]])
      abline(a=-2*max_scale/(max(x_plot)-min(x_plot)+diff(x_plot[1:2]))*(x_plot[1]-diff(x_plot[1:2])/2),
             b=2*max_scale/(max(x_plot)-min(x_plot)+diff(x_plot[1:2])))
      abline(a=2*max_scale/(max(x_plot)-min(x_plot)+diff(x_plot[1:2]))*(rev(x_plot)[1]+diff(x_plot[1:2])/2),
             b=-2*max_scale/(max(x_plot)-min(x_plot)+diff(x_plot[1:2])))
      box()
      par(mar=c(4,0.5,4,3))
      image(1,seq.int(nlevel)-(nlevel+1)/2,t(as.matrix(seq.int(nlevel))),axes=FALSE,xlab='',ylab='',col=col_heatmap)
      axis(side=4,at=seq(0,nlevel,length.out=6)-nlevel/2,labels=format(seq(0,1,length.out=6),2),padj=0.4,...)
      box()
      
      # Adjusted p-value plot at scale_threshold
      par(mar=c(4,mar.left,1.5,1.5))
      pval_scale_threshold=adjusted_pval(regionsFeatures,i,id_feature,scale_threshold_i[[id_feature]])[[1]][[1]]
      plot(1,type="n",xlim=c(x_plot[1]-diff(x_plot[1:2])/2,rev(x_plot)[1]+diff(x_plot[1:2])/2),ylim=c(0,1),ylab="p-value",xlab=xlab,xaxs="i",
           main=paste0("Adjusted p-values - Threshold ",scale_threshold_i[[id_feature]]),...)
      low.p.value=which((pval_scale_threshold<=alpha)&(!is.infinite(pval_scale_threshold)))
      if(length(low.p.value)){
        for(j in seq_along(low.p.value))
          rect(x_plot[low.p.value[j]]-diff(x_plot[1:2])/2,par("usr")[3],
               x_plot[low.p.value[j]]+diff(x_plot[1:2])/2,par("usr")[4],col="gray90",density=-2,border=NA)
      }
      box()
      for (j in 0:10)
        abline(h=j/10,col="lightgray",lty="dotted")
      points(x_plot,pval_scale_threshold,pch=16)
      
      # Data plot
      if(is.null(ylim_ext)){
        if(average){
          ylim=pmin(range(c(unlist(plot_data$features_plot[[1]]),unlist(plot_data$features_average[[1]])),na.rm=TRUE),c(0,+Inf))
        }else{
          ylim=pmin(range(unlist(plot_data$features_plot[[1]]),na.rm=TRUE),c(0,+Inf))
        }
      }else{
        ylim=ylim_ext
      }
      par(mar=c(4,mar.left,1.5,1.5),las=0)
      if(type=='curves'){
        matplot(x_plot,plot_data$features_plot[[id_feature]],type='l',col=plot_data$col_plot,xlim=par('usr')[1:2],ylim=ylim,xaxs="i",
                main=nameFeatures(regionsFeatures)[id_feature],xlab=xlab,
                ylab=nameFeatures(regionsFeatures)[id_feature],...)
      }
      if(type=='boxplot'){
        plot(1,type="n",xlim=par('usr')[1:2],ylim=ylim,xaxs="i",
             main=nameFeatures(regionsFeatures)[id_feature],xlab=xlab,
             ylab=nameFeatures(regionsFeatures)[id_feature],...)
        for(id_region in c(id_region1_i,id_region2_i)){
          polygon(c(x_plot,rev(x_plot)),c(plot_data$features_plot[[id_feature]][[id_region]][,1],rev(plot_data$features_plot[[id_feature]][[id_region]][,length(probs)])),
                  col=plot_data$col_plot[[id_region]][1],border=FALSE)
          for(j in seq_along(probs))
            lines(x_plot,plot_data$features_plot[[id_feature]][[id_region]][,j],col=plot_data$col_plot[[id_region]][2],lty=2,lwd=1.5)
        }
      }
      if(average)
        matplot(x_plot,plot_data$features_average[[id_feature]],type='l',col=col[c(id_region1_i,id_region2_i)],lty=1,lwd=2,add=TRUE)
      if(.testResults(regionsFeatures)[[i]][[id_feature]]$test=='1pop')
        lines(x_plot,rep(.testResults(regionsFeatures)[[i]][[id_feature]]$mu,p),col='black',type='l',lwd=2)
      args=as.list(match.call())
      if(is.null(args$cex)){
        cex=ifelse(is.null(args$cex.lab),1,args$cex.lab)
      }else{
        cex=args$cex
      }
      legend(par('usr')[2],mean(par('usr')[3:4]),legend=nameRegions(regionsFeatures)[c(id_region1_i,id_region2_i)],xpd=NA,bty='n',lty=1,lwd=2,col=col[c(id_region1_i,id_region2_i)],yjust=0.5,cex=cex)
      if(size){
        par(mar=c(1,mar.left,2,1.5))
        image(x_plot,seq_along(c(id_region1_i,id_region2_i)),plot_data$features_position_size[[id_feature]],
              col=cm.colors(101),xlim=par('usr')[1:2],ylim=range(seq_along(c(id_region1_i,id_region2_i)))+c(-0.5,0.5),xaxs="i",axes=FALSE,xlab='',ylab='',...)
        axis(side=3,...)
        axis(side=2,at=seq_along(c(id_region1_i,id_region2_i)),labels=rev(nameRegions(regionsFeatures)[c(id_region1_i,id_region2_i)]),tick=FALSE,las=1,line=-0.5,...)
        x.rect=range(x_plot)+c(-1,1)*diff(x_plot)[1]/2
        rect(x.rect[1],seq_along(c(id_region1_i,id_region2_i))-0.5,x.rect[2],seq_along(c(id_region1_i,id_region2_i))+0.5,border='black')
        par(mar=c(4,mar.left,0,1.5))
        image(-50:50,1,as.matrix(seq.int(101)),xlim=c(-100,100),axes=FALSE,xlab='Sample size',ylab='',col=cm.colors(101),...)
        labels=seq(min(plot_data$features_position_size[[id_feature]]),max(plot_data$features_position_size[[id_feature]]),length.out=5)
        if(length(unique(labels))==1)
          labels=-2:2+labels
        axis(side=1,at=c(-50,-25,0,25,50),labels=labels,...)
        rect(-50.5,par('usr')[3],50.5,par('usr')[4],border='black')
      }
    }
  }
}
