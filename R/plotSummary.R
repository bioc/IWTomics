plotSummary <- function(regionsFeatures,alpha=0.05,only_significant=FALSE,scale_threshold=NULL,nlevel=100,
                        groupby='test',test=1:nTests(regionsFeatures),gaps_tests=NULL,
                        id_features_subset=idFeaturesTest(regionsFeatures),gaps_features=NULL,
                        ask=TRUE,filenames=NA,align_lab=NA,cellwidth=NA,cellheight=NA,
                        xlab='Windows',ylab=ifelse(groupby=='test','Features','Tests'),...){

  if(class(regionsFeatures)!='IWTomicsData')
    stop('invalid regionsFeatures. IWTomicsData object expected.')
  if(!validObject(regionsFeatures))
    stop('invalid regionsFeatures.')
  if(is.null(regionsFeatures@test))
    stop('No test results present in regionFeatures.')
  if(!(groupby %in% c('test','feature')))
    stop('invalid groupby \'',groupby,'\'. Available groupby are \'test\' and \'feature\'.')
  if(sum(!(test %in% 1:nTests(regionsFeatures))))
    stop('invalid test number.')
  if(sum(!(id_features_subset %in% idFeaturesTest(regionsFeatures))))
    stop('invalid id_features_subset. The features provided are not listed in regionsFeatures.')
  if(is.null(scale_threshold)){
    scale_threshold=lapply(.testResults(regionsFeatures)[test],
                           function(test){
                             unlist(lapply(test[id_features_subset],function(feature) feature$max_scale))
                           })
  }else{
    if(is.list(scale_threshold)){
      scale_threshold=as.list(rep(scale_threshold,length.out=length(test)))
      scale_threshold=lapply(scale_threshold,function(scale_threshold) rep(scale_threshold,length.out=length(id_features_subset)))
    }else{
      scale_threshold=lapply(seq_along(test),function(test) return(rep(scale_threshold,length.out=length(id_features_subset))))
    }
    scale_threshold=lapply(scale_threshold,
                           function(scale_threshold){
                             names(scale_threshold)=id_features_subset
                             return(scale_threshold)})
  }
  if(alignment(regionsFeatures)=='scale')
    align_lab=NA
  
  devAskNewPage(ask)
  if(groupby=='test'){
    if(length(unique(resolution(regionsFeatures)[id_features_subset]))!=1)
      stop('groupby \'test\' but selected features with different resolution. Smooth data first to have the same resolution.')
    if(!is.na(filenames)[1])
      if(length(filenames)<length(test))
        stop('filenames should contain one file path for each test as defined by groupby.')
    if(!is.null(gaps_features))
       if(FALSE %in% (gaps_features %in% seq_along(id_features_subset)))
         stop('invalid gap_features.')
    for(i in seq_along(test)){
      is.first=ifelse(i==1,TRUE,FALSE)
      filename_i=filenames[i]
      scale_threshold_i=scale_threshold[[i]]
      id_region1_i=idRegionsTest(regionsFeatures,test[i])[[1]][1]
      id_region2_i=idRegionsTest(regionsFeatures,test[i])[[1]][2]
      if((!is.null(id_region2_i))&&(id_region2_i==''))
        id_region2_i=NULL
      id_features_subset_i=id_features_subset
      gaps_features_i=gaps_features
      result_i=.testResults(regionsFeatures,test[i],id_features_subset_i)[[1]]
      pval_scale_threshold=Reduce(cbind,adjusted_pval(regionsFeatures,test[i],id_features_subset_i,scale_threshold_i)[[1]])
      # Significant p-value
      significant=pval_scale_threshold
      significant[significant<1e-5]=1e-5
      significant[significant>alpha]=1
      if(only_significant){
        if(!is.null(gaps_features_i)){
          with_gaps=rep(NA,length(id_features_subset_i)+length(gaps_features_i))
          with_gaps[-(gaps_features_i+seq_along(gaps_features_i))]=id_features_subset_i
        }
        id_features_subset_i=id_features_subset_i[colSums(significant<=alpha,na.rm=TRUE)>0]
        significant=significant[,id_features_subset_i]
        result_i=result_i[id_features_subset_i]
        scale_threshold_i=scale_threshold_i[id_features_subset_i]
        if(!is.null(gaps_features_i)){
          with_gaps=with_gaps[with_gaps %in% c(id_features_subset_i,NA)]
          gaps_features_i=which(is.na(with_gaps))-seq_len(sum(is.na(with_gaps)))
          gaps_features_i=gaps_features_i[gaps_features_i>0]
          if(length(gaps_features_i)==0)
            gaps_features_i=NULL
        }
      }
      log_significant=as.matrix(-log10(significant))
      colnames(log_significant)=nameFeatures(regionsFeatures)[id_features_subset_i]
      T0_plot=as.matrix(Reduce(cbind,lapply(result_i,function(feature) feature$T0_plot)))
      colnames(T0_plot)=id_features_subset_i
      # Recode the T0_plot matrix in 1 and -1
      if(testInput(regionsFeatures)$statistics=='variance'){
        T0_plot[T0_plot<1]=-1
        T0_plot[T0_plot>1]=1
      }else{
        T0_plot[T0_plot<0]=-1
        T0_plot[T0_plot>=0]=1
      }
      # Put a sign to pvalues based on T0_plot
      log_significant=log_significant*T0_plot
      
      if(alignment(regionsFeatures)=='left')
        row.names(log_significant)=seq_len(nrow(log_significant))
      if(alignment(regionsFeatures)=='right')
        row.names(log_significant)=-(nrow(log_significant):1)
      if(alignment(regionsFeatures)=='center')
        if(nrow(log_significant)%%2!=0){
          row.names(log_significant)=seq_len(nrow(log_significant))-(nrow(log_significant)+1)/2
        }else{
          row.names(log_significant)=setdiff(seq(-nrow(log_significant)/2,nrow(log_significant)/2),0)
        }
      if(alignment(regionsFeatures)=='scale')
        row.names(log_significant)=round(seq(0,1,length.out=nrow(log_significant)),2)
      
      if(result_i[[1]]$test=="1pop"){
        main=paste0(nameRegions(regionsFeatures)[id_region1_i])
      }else{
        main=paste0(nameRegions(regionsFeatures)[id_region1_i],' vs ',nameRegions(regionsFeatures)[id_region2_i])
      }
      .pheatmap(t(log_significant),gaps_row=gaps_features_i,border_color="grey60",filename=filename_i,is.first=is.first,
                breaks=seq(-5,5,length.out=nlevel),color=colorRampPalette(c("navy","white","red"))(n=nlevel-1),
                legend.main='-log10(p-value)',main=main,cex.main=2,xlab=xlab,ylab=ylab,thresholds=scale_threshold_i,
                zero_lab=align_lab,zero_lab_pos=alignment(regionsFeatures),cellwidth=cellwidth,cellheight=cellheight,...)
    }
  } else {
    if(!is.na(filenames[1]))
      if(length(filenames)<length(id_features_subset))
        stop('filenames should contain one file path for each test as defined by groupby.')
    if(!is.na(filenames)[1])
      if(length(filenames)<length(id_features_subset))
        stop('filenames should contain one file path for each test as defined by groupby.')
    if(!is.null(gaps_tests))
      if(FALSE %in% (gaps_tests %in% 1:nTests(regionsFeatures)))
        stop('invalid gap_tests.')
    for(i in seq_along(id_features_subset)){
      is.first=ifelse(i==1,TRUE,FALSE)
      filename_i=filenames[i]
      scale_threshold_i=lapply(scale_threshold,function(scale_threshold) scale_threshold[id_features_subset[i]])
      id_region1_i=testInput(regionsFeatures)$id_region1[test]
      id_region2_i=testInput(regionsFeatures)$id_region2[test]
      if(is.null(id_region2_i))
        id_region2_i=rep('',length(id_region1_i))
      id_features_subset_i=id_features_subset[i]
      gaps_tests_i=gaps_tests
      gaps_tests_i=match(gaps_tests_i,test)
      gaps_tests_i=gaps_tests_i[!is.na(gaps_tests_i)]
      result_i=.testResults(regionsFeatures,test,id_features_subset_i)
      pval_scale_threshold=mapply(function(test,pval){
                             pval_scale_threshold=rep(NA,length(test[[1]]$notNA))
                             pval_scale_threshold[test[[1]]$notNA]=pval[[1]]
                             return(pval_scale_threshold)
                           },result_i,adjusted_pval(regionsFeatures,test,id_features_subset_i,scale_threshold_i))
      
      # Significant p-value
      significant=pval_scale_threshold
      significant[significant<1e-5]=1e-5
      significant[significant>alpha]=1
      if(only_significant){
        if(!is.null(gaps_tests_i)){
          with_gaps=rep(NA,length(id_region1_i)+length(gaps_tests_i))
          with_gaps[-(gaps_tests_i+seq_along(gaps_tests_i))]=seq_along(id_region1_i)
        }
        index=which(colSums(significant<=alpha,na.rm=TRUE)>0)
        id_region1_i=id_region1_i[index]
        id_region2_i=id_region2_i[index]
        significant=significant[,index]
        result_i=result_i[index]
        scale_threshold_i=scale_threshold_i[index]
        if(!is.null(gaps_tests_i)){
          with_gaps=with_gaps[with_gaps %in% c(index,NA)]
          gaps_tests_i=which(is.na(with_gaps))-seq_len(sum(is.na(with_gaps)))
          gaps_tests_i=gaps_tests_i[gaps_tests_i>0]
          if(length(gaps_tests_i)==0)
            gaps_tests_i=NULL
        }
      }
      log_significant=as.matrix(-log10(significant))
      colnames(log_significant)[id_region2_i!='']=paste0(nameRegions(regionsFeatures)[id_region1_i[id_region2_i!='']],' vs ',nameRegions(regionsFeatures)[id_region2_i[id_region2_i!='']])
      colnames(log_significant)[id_region2_i=='']=nameRegions(regionsFeatures)[id_region1_i[id_region2_i=='']]
      T0_plot=as.matrix(Reduce(cbind,lapply(result_i,
                                            function(test){
                                              T0_plot=rep(NA,length(test[[1]]$notNA))
                                              T0_plot[test[[1]]$notNA]=test[[1]]$T0_plot
                                              return(T0_plot)
                                            })))
      colnames(T0_plot)[id_region2_i!='']=paste0(nameRegions(regionsFeatures)[id_region1_i[id_region2_i!='']],' vs ',nameRegions(regionsFeatures)[id_region2_i[id_region2_i!='']])
      colnames(T0_plot)[id_region2_i=='']=nameRegions(regionsFeatures)[id_region1_i[id_region2_i=='']]
      # Recode the T0_plot matrix in 1 and -1
      if(testInput(regionsFeatures)$statistics=='variance'){
        T0_plot[T0_plot<1]=-1
        T0_plot[T0_plot>1]=1
      }else{
        T0_plot[T0_plot<0]=-1
        T0_plot[T0_plot>=0]=1
      }
      # Put a sign to pvalues based on T0_plot
      log_significant=log_significant*T0_plot
      
      if(alignment(regionsFeatures)=='left')
        row.names(log_significant)=seq_len(nrow(log_significant))
      if(alignment(regionsFeatures)=='right')
        row.names(log_significant)=-(nrow(log_significant):1)
      if(alignment(regionsFeatures)=='center')
        if(nrow(log_significant)%%2!=0){
          row.names(log_significant)=seq_len(nrow(log_significant))-(nrow(log_significant)+1)/2
        }else{
          row.names(log_significant)=setdiff(seq(-nrow(log_significant)/2,nrow(log_significant)/2),0)
        }
      if(alignment(regionsFeatures)=='scale')
        row.names(log_significant)=round(seq(0,1,length.out=nrow(log_significant)),2)
      
      main=paste0(nameFeatures(regionsFeatures)[id_features_subset_i])
      .pheatmap(t(log_significant),gaps_row=gaps_tests_i,border_color="grey60",filename=filename_i,is.first=is.first,
                breaks=seq(-5,5,length.out=nlevel),color=colorRampPalette(c("navy","white","red"))(n=nlevel-1),
                legend.main='-log10(p-value)',main=main,cex.main=2,xlab=xlab,ylab=ylab,thresholds=unlist(scale_threshold_i),
                zero_lab=align_lab,zero_lab_pos=alignment(regionsFeatures),cellwidth=cellwidth,cellheight=cellheight,...)
    }
  }
}

.pheatmap <- function (mat, color = colorRampPalette(c("navy","white","red"))(100), 
                       breaks = NA, border_color = "grey60", 
                       cellwidth = NA, cellheight = NA, scale = "none", 
                       legend = TRUE, legend.main = '',legend_breaks = NA, legend_labels = NA, 
                       show_rownames = TRUE, show_colnames = TRUE, main = NA, fontsize = 10, 
                       fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = FALSE, 
                       number_format = "%.2f", number_color = "grey30", fontsize_number = 0.8 * fontsize, 
                       gaps_row = NULL, gaps_col = NULL, labels_row = NULL, 
                       labels_col = NULL, filename = NA, is.first=TRUE, width = NA, height = NA, 
                       silent = FALSE, xlab = NA, ylab = NA, cex.main = 1.3, thresholds = NA, zero_lab = NA, zero_lab_pos = "center", ...) {
  # modified from pheatmap package
  
  # Set labels
  if (is.null(labels_row)) {
    labels_row = rownames(mat)
  }
  if (is.null(labels_col)) {
    labels_col = colnames(mat)
  }
  
  # Preprocess matrix
  mat = as.matrix(mat)
  if (scale != "none") {
    mat = .scale_mat(mat, scale)
    if (.is.na2(breaks)) {
      breaks = .generate_breaks(mat, length(color), center = TRUE)
    }
  }
  
  # Format numbers to be displayed in cells
  if (is.matrix(display_numbers) | is.data.frame(display_numbers)) {
    if (nrow(display_numbers) != nrow(mat) | ncol(display_numbers) != ncol(mat)) {
      stop("If display_numbers provided as matrix, its dimensions have to match with mat")
    }
    display_numbers = as.matrix(display_numbers)
    fmat = matrix(as.character(display_numbers), nrow = nrow(display_numbers), ncol = ncol(display_numbers))
    fmat_draw = TRUE
  } else {
    if (display_numbers) {
      fmat = matrix(sprintf(number_format, mat), nrow = nrow(mat), ncol = ncol(mat))
      fmat_draw = TRUE
    } else {
      fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
      fmat_draw = FALSE
    }
  }
  
  attr(fmat, "draw") = fmat_draw
  
  # Colors and scales
  if (!.is.na2(legend_breaks) & !.is.na2(legend_labels)) {
    if (length(legend_breaks) != length(legend_labels)) {
      stop("Lengths of legend_breaks and legend_labels must be the same")
    }
  }
  if (.is.na2(breaks)) {
    breaks = .generate_breaks(as.vector(mat), length(color))
  }
  if (legend & .is.na2(legend_breaks)) {
    legend = grid.pretty(range(as.vector(breaks)))
    names(legend) = legend
  } else if (legend & !.is.na2(legend_breaks)) {
    legend = legend_breaks[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
    if (!.is.na2(legend_labels)) {
      legend_labels = legend_labels[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
      names(legend) = legend_labels
    } else {
      names(legend) = legend
    }
  } else {
    legend = NA
  }
  mat = .scale_colours(mat, col = color, breaks = breaks)
  
  if (!show_rownames) {
    labels_row = NULL
  }
  
  if (!show_colnames) {
    labels_col = NULL
  }
  
  # Draw heatmap
  gt = .heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, 
                      cellheight = cellheight, filename = filename, width = width, 
                      height = height, breaks = breaks, color = color, legend = legend, legend.main = legend.main,
                      main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                      fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, 
                      number_color = number_color, gaps_row = gaps_row, gaps_col = gaps_col, thresholds = thresholds, zero_lab = zero_lab, zero_lab_pos = zero_lab_pos, 
                      labels_row = labels_row, labels_col = labels_col, xlab = xlab, ylab = ylab, cex.main = cex.main,...)
  
  if (is.na(filename) & !silent) {
    if(!is.first)
      grid.newpage()
    grid.draw(gt)
  }
  
  invisible(list(gtable = gt))
}

.lo <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, 
                 legend, legend.main, main, fontsize, fontsize_row, thresholds, zero_lab, 
                 fontsize_col, gaps_row, gaps_col, xlab, ylab, cex.main = cex.main, ...) {
  # modified from pheatmap package
  
  # Get height of colnames and length of rownames
  if (!is.null(coln[1])) {
    longest_coln = which.max(strwidth(coln, units = "in"))
    gp = list(fontsize = fontsize_col, ...)
    coln_height = unit(1, "grobheight", textGrob(coln[longest_coln], rot = 90, gp = do.call(gpar, gp))) + unit(10, "bigpts")
  } else {
    coln_height = unit(5, "bigpts")
  }
  
  if (!is.null(rown[1])) {
    longest_rown = which.max(strwidth(rown, units = "in"))
    gp = list(fontsize = fontsize_row, ...)
    rown_width = unit(1, "grobwidth", textGrob(rown[longest_rown], gp = do.call(gpar, gp))) + unit(5, "bigpts")
  } else {
    rown_width = unit(5, "bigpts")
  }
  
  gp = list(fontsize = fontsize, ...)
  # Axis label positions
  if (!.is.na2(xlab)) {
    xlab_height = unit(1, "grobheight", textGrob(xlab, gp = do.call(gpar, gp))) + unit(5, "bigpts")
  } else {
    xlab_height = unit(1, "bigpts")
  }
  if (!.is.na2(ylab)) {
    ylab_width = unit(1, "grobheight", textGrob(ylab, rot=270, gp = do.call(gpar, gp))) + unit(1, "bigpts")
  } else {
    ylab_width = unit(5, "bigpts")
  }
  
  # Threshold position
  if (!.is.na2(thresholds)){
    thresholds_width = unit(1.1, "grobwidth", textGrob("Threshold", gp = gpar(...)))+unit(2, "bigpts")
    thresholds_lab_height = unit(1, "grobheight", textGrob("Threshold", just = "bottom", gp = do.call(gpar, gp))) + unit(5, "bigpts")
  } else {
    thresholds_width = unit(0, "bigpts")
    thresholds_lab_height = unit(0, "bigpts")
  }
  
  # Zero label position
  if(!.is.na2(zero_lab)){
    thresholds_lab_height = unit(1, "grobheight", textGrob("Threshold", just = "bottom", gp = do.call(gpar, gp))) + unit(5, "bigpts")
  } 
  
  # Legend position
  if (!.is.na2(legend)) {
    longest_break = which.max(nchar(names(legend)))
    longest_break = unit(1.1, "grobwidth", textGrob(as.character(names(legend))[longest_break], gp = do.call(gpar, gp)))
    legend_width = unit(25, "bigpts") + longest_break * 1.2
  } else {
    legend_width = unit(0, "bigpts")
  }
  if (!.is.na2(legend.main)) {
    legend_main_width = unit(1.1, "grobwidth", textGrob(legend.main, rot=90, gp = gpar(...)))+unit(6, "bigpts")
  } else {
    legend_main_width = unit(0, "bigpts")
  }
  
  # Set main title height
  if (is.na(main)) {
    main_height = unit(0, "npc")
  } else {
    main_height = unit(2.2, "grobheight", textGrob(main, gp = gpar(fontsize = cex.main * fontsize, ...)))
  }
  textheight = unit(fontsize, "bigpts")
  
  # Set cell sizes
  if (is.na(cellwidth)) {
    mat_width = unit(1, "npc") - thresholds_width - rown_width - ylab_width - legend_width - legend_main_width
  } else {
    mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * unit(6, "bigpts")
  }
  
  if (is.na(cellheight)) {
    mat_height = unit(1, "npc") - main_height - thresholds_lab_height - coln_height - xlab_height 
  } else {
    mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * unit(6, "bigpts")
  }
  
  # Produce gtable
  gt = gtable(widths = unit.c(unit(0, "bigpts"), unit(0, "bigpts"), thresholds_width, mat_width, rown_width, ylab_width, legend_width, legend_main_width, unit(0, "bigpts")), 
              heights = unit.c(main_height, unit(0, "bigpts"), unit(0, "bigpts"), thresholds_lab_height, mat_height, coln_height, xlab_height), vp = viewport(gp = do.call(gpar,  gp)))
  cw = convertWidth(mat_width - (length(gaps_col) * unit(6, "bigpts")), "bigpts", valueOnly = TRUE)/ncol
  ch = convertHeight(mat_height - (length(gaps_row) * unit(6, "bigpts")), "bigpts", valueOnly = TRUE)/nrow
  
  # Return minimal cell dimension in bigpts to decide if borders are drawn
  mindim = min(cw, ch)
  
  res = list(gt = gt, mindim = mindim)
  
  return(res)
}

.find_coordinates <- function (n, gaps, m = 1:n) {
  if (length(gaps) == 0) {
    return(list(coord = unit(m/n, "npc"), size = unit(1/n, "npc")))
  }
  
  if (max(gaps) > n) {
    stop("Gaps do not match with matrix size")
  }
  
  size = (1/n) * (unit(1, "npc") - length(gaps) * unit("6", "bigpts"))
  gaps2 = apply(sapply(gaps, function(gap, x) { x > gap }, m), 1, sum)
  
  coord = m * size + (gaps2 * unit("6", "bigpts"))
  return(list(coord = coord, size = size))
}

.draw_matrix=function (matrix, border_color, gaps_rows, gaps_cols, fmat, fontsize_number, number_color, zero_lab, zero_lab_pos) {
  # modified from pheatmap package
  
  n = nrow(matrix)
  m = ncol(matrix)
  
  coord_x = .find_coordinates(m, gaps_cols)
  coord_y = .find_coordinates(n, gaps_rows)
  
  x = coord_x$coord - 0.5 * coord_x$size
  y = unit(1, "npc") - (coord_y$coord - 0.5 * coord_y$size)
  
  res = gList()
  for(i in seq_len(nrow(matrix))){
    index_i=which(!is.na(matrix[i,]))
    x_i = (coord_x$coord - 0.5 * coord_x$size)[index_i]
    y_i = y[i]
    coord = expand.grid(y = y_i, x = x_i)
    res[[paste0("rect",i)]] = rectGrob(x = coord$x, y = coord$y, width = coord_x$size, 
                                       height = coord_y$size, gp = gpar(fill = matrix[i,][index_i], col = border_color))
  }
  
  if (attr(fmat, "draw")) {
    res[["text"]] = textGrob(x = coord$x, y = coord$y, label = fmat, gp = gpar(col = number_color, fontsize = fontsize_number))
  }
  
  if(!.is.na2(zero_lab)) {
    if(zero_lab_pos=='center')
      x = rep(unit(0.5, "npc"),2)
    if(zero_lab_pos=='right')
      x = rep(unit(1, "npc"),2)
    if(zero_lab_pos=='left')
      x = rep(unit(0, "npc"),2)
    
    #res[["zero_line"]] = linesGrob(x = x, y = unit(c(0,1), "npc")+unit(c(0,3), "bigpts"),gp=gpar(lwd = 2))
    res[["zero_line"]] = linesGrob(x = x, y = y[c(length(y),1)]+c(-0.5,0.5)*coord_y$size+unit(c(0,3), "bigpts"),gp=gpar(lwd = 2))
  }
  
  res = gTree(children = res)
  
  return(res)
}

.draw_colnames <- function (coln, gaps, ...) {
  # modified from pheatmap package
  
  coord = .find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), just='right', rot = 90, gp = gpar(...))
  
  return(res)
}

.draw_rownames <- function (rown, gaps, ...) {
  # modified from pheatmap package
  
  coord = .find_coordinates(length(rown), gaps)
  y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
  
  res = textGrob(rown, x = unit(3, "bigpts"), y = y, just='left', gp = gpar(...))
  
  return(res)
}

.draw_thresholds <- function (thresholds, gaps, ...) {
  
  coord = .find_coordinates(length(thresholds), gaps)
  y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
  
  res = textGrob(thresholds, x =unit(0.9, "npc") ,y = y, just='right', gp = gpar(...))
  
  return(res)
}

.draw_thresholds_lab <- function(thresholds_lab, ...) {
  
  res = textGrob(thresholds_lab, x = unit(0.9, "npc"), y = unit(15, "bigpts"), just=c('right','top'), gp = gpar(...))
  
  return(res)
}

.draw_zero_lab <- function(zero_lab, zero_lab_pos, ...) {
  
  if(zero_lab_pos=='center')
    x = unit(0.5, "npc")
  if(zero_lab_pos=='right')
    x = unit(1, "npc")
  if(zero_lab_pos=='left')
    x = unit(0, "npc")
  
  res = textGrob(zero_lab, x = x, y = unit(15, "bigpts"), just=c(zero_lab_pos,'top'), gp = gpar(...))
  
  return(res)
}

.draw_legend <- function (color, breaks, legend,...) {
  # modified from pheatmap package
  
  height = min(unit(1, "npc"), unit(150, "bigpts"))
  
  legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
  legend_pos = height * legend_pos + (unit(1, "npc") - height)
  
  breaks = (breaks - min(breaks))/(max(breaks) - min(breaks))
  breaks = height * breaks + (unit(1, "npc") - height)
  
  h = breaks[-1] - breaks[-length(breaks)]
  
  rect = rectGrob(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, 
                  hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
  text = textGrob(names(legend), x = unit(14, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))
  
  res = grobTree(rect, text)
  return(res)
}

.draw_legend_main <- function (breaks, legend, legend.main,...) {
  height = min(unit(1, "npc"), unit(150, "bigpts"))
  
  legend_main_pos=(sum(range(legend))/2-min(breaks))/(max(breaks) - min(breaks))
  legend_main_pos=height * legend_main_pos + (unit(1, "npc") - height)
  
  label = textGrob(legend.main, x = 0, y = legend_main_pos, just='centre', rot=270, gp = gpar(...))
  
  res = grobTree(label)
  return(res)
}

.draw_main = function(text, ...){
  # modified from pheatmap package
  
  res = textGrob(text, gp = gpar(fontface = "bold", ...))
  
  return(res)
}

.draw_xlab = function(text, ...){
  
  res = textGrob(text, just = "bottom", gp = gpar(...))
  
  return(res)
}

.draw_ylab = function(text, ...){
  
  res = textGrob(text, just = "center", rot = 270, gp = gpar(...))
  
  return(res)
}

.heatmap_motor <- function (matrix, border_color, cellwidth, cellheight, filename, width, height, 
                            breaks, color, legend, legend.main, main, fontsize, fontsize_row, 
                            fontsize_col, fmat, fontsize_number, number_color, gaps_col, thresholds, zero_lab, zero_lab_pos,
                            gaps_row, labels_row, labels_col, xlab, ylab, cex.main,...) {
  # modified from pheatmap package
  
  # Set layout
  lo = .lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), 
           ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, 
           legend = legend, legend.main=legend.main, 
           main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
           thresholds = thresholds, zero_lab = zero_lab,
           fontsize_col = fontsize_col, gaps_row = gaps_row, gaps_col = gaps_col, 
           xlab = xlab, ylab = ylab, cex.main = cex.main, ...)
  
  res = lo$gt
  mindim = lo$mindim
  
  if (!is.na(filename)) {
    if (is.na(height)) {
      height = convertHeight(gtable_height(res), "inches", valueOnly = TRUE)
    }
    if (is.na(width)) {
      width = convertWidth(gtable_width(res), "inches", valueOnly = TRUE)
    }
    # Get file type
    r = regexpr("\\.[a-zA-Z]*$", filename)
    if (r == -1) 
      stop("Improper filename")
    ending = substr(filename, r + 1, r + attr(r, "match.length"))
    
    f = switch(ending, pdf = function(x, ...) pdf(x, ...), 
               png = function(x, ...) png(x, units = "in", res = 300, ...), 
               jpeg = function(x, ...) jpeg(x, units = "in", res = 300, ...), 
               jpg = function(x, ...) jpeg(x, units = "in", res = 300, ...), 
               tiff = function(x, ...) tiff(x, units = "in", res = 300, compression = "lzw", ...), 
               bmp = function(x, ...) bmp(x, units = "in", res = 300, ...), 
               stop("File type should be: pdf, png, bmp, jpg, tiff"))
    
    f(filename, height = height, width = width)
    gt = .heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, 
                        border_color = border_color, breaks = breaks, 
                        color = color, legend = legend, legend.main=legend.main, filename = NA, 
                        main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                        fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, thresholds = thresholds, zero_lab = zero_lab, zero_lab_pos = zero_lab_pos,
                        number_color = number_color, labels_row = labels_row, xlab = xlab, ylab = ylab, cex.main = cex.main, 
                        labels_col = labels_col, gaps_col = gaps_col, gaps_row = gaps_row, ...)
    grid.draw(gt)
    dev.off()
    
    return(gt)
  }
  
  # Omit border color if cell size is too small 
  if (mindim < 3) 
    border_color = NA
  
  # Draw title
  if (!is.na(main)) {
    elem = .draw_main(main, fontsize = cex.main * fontsize, ...)
    res = gtable_add_grob(res, elem, t = 1, l = 4, name = "main")
  }
  
  # Draw thresholds
  if (!.is.na2(thresholds)) {
    elem = .draw_thresholds(thresholds, gaps_row, fontsize = fontsize)
    res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", name = "thresholds")
    elem = .draw_thresholds_lab("Threshold", fontsize = fontsize)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", name = "thresholds_lab")
  }
  
  # Draw matrix
  elem = .draw_matrix(matrix, border_color, gaps_row, gaps_col, fmat, fontsize_number, number_color, zero_lab, zero_lab_pos)
  res = gtable_add_grob(res, elem, t = 5, l = 4, clip = "off", name = "matrix")
  
  # Draw zero label
  if (!.is.na2(zero_lab)) {
    elem = .draw_zero_lab(zero_lab, zero_lab_pos, fontsize = fontsize)
    res = gtable_add_grob(res, elem, t = 4, l = 4, clip = "off", name = "zero_lab")
  }
  
  # Draw colnames
  if (length(labels_col) != 0) {
    pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, ...)
    elem = do.call(.draw_colnames, pars)
    res = gtable_add_grob(res, elem, t = 6, l = 4, clip = "off", name = "col_names")
  }
  
  # Draw rownames
  if (length(labels_row) != 0) {
    pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row,  ...)
    elem = do.call(.draw_rownames, pars)
    res = gtable_add_grob(res, elem, t = 5, l = 5, clip = "off", name = "row_names")
  }
  
  # Draw xlab
  if (!.is.na2(xlab)){
    elem = .draw_xlab(xlab,fontsize = fontsize)
    res = gtable_add_grob(res, elem, t = 7, l = 4, clip = "off", name = "xlab")
  }
  
  # Draw ylab
  if (!.is.na2(ylab)){
    elem = .draw_ylab(ylab,fontsize = fontsize)
    res = gtable_add_grob(res, elem, t = 5, l = 6, clip = "off", name = "ylab")
  }
  
  # Draw legend
  if (!.is.na2(legend)) {
    elem = .draw_legend(color, breaks, legend, fontsize = fontsize, ...)
    #t = ifelse(is.null(labels_row), 5, 4)
    res = gtable_add_grob(res, elem, t = 5, l = 7, b = 5, clip = "off", name = "legend")
  }
  if (!.is.na2(legend.main)) {
    elem = .draw_legend_main(breaks, legend, legend.main, fontsize = fontsize, ...)
    #t = ifelse(is.null(labels_row), 4, 3)
    res = gtable_add_grob(res, elem, t = 5, l = 8, b = 5, clip = "off", name = "legend_main")
  }
  
  return(res)
}

.generate_breaks = function(x, n, center = FALSE){
  # modified from pheatmap package
  
  if(center){
    m = max(abs(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))))
    res = seq(-m, m, length.out = n + 1)
  }
  else{
    res = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1)
  }
  
  return(res)
}

.scale_vec_colours <- function (x, col, breaks = NA){
  # modified from pheatmap package
  
  return(col[as.numeric(cut(x, breaks = breaks, include.lowest = TRUE))])
}

.scale_colours <- function (mat, col, breaks = NA) {
  # modified from pheatmap package
  
  mat = as.matrix(mat)
  return(matrix(.scale_vec_colours(as.vector(mat), col = col, breaks = breaks), nrow(mat), 
                ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

.scale_rows = function(x){
  # modified from pheatmap package
  
  m = apply(x, 1, mean, na.rm = TRUE)
  s = apply(x, 1, sd, na.rm = TRUE)
  return((x - m) / s)
}

.scale_mat = function(mat, scale){
  # modified from pheatmap package
  
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = .scale_rows(mat), column = t(.scale_rows(t(mat))))
  return(mat)
}

.is.na2 = function(x){
  # modified from pheatmap package
  
  if(is.list(x) | length(x) > 1){
    return(FALSE)
  }
  if(length(x) == 0){
    return(TRUE)
  }
  
  return(is.na(x))
}

.identity2 = function(x, ...){
  # modified from pheatmap package
  
  return(x)
}

