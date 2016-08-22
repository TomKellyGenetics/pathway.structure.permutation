##' @name Plots for Tests of Pathway Structure in iGraph
##' @rdname pathway.structure.plot
##'
##' @title Extensions to iGraph for Plotting Pathway Structure
##'
##' @description Functions to plot outcomes for testing pathway structure with permutation analysis.
##'
##'
##' @param table_obj List. Output of \code{\link[permutation.structure.permutation]{pathway.structure}}
##' @param main,sub,xlab,ylab,cex.lab,col,cex,pch Arguments to pass on to \code{\link[graphics]{plot}} or \code{\link[stats]{plot.density}}
##' @param guide,legend,annotate Logical. Whether guide lines, legend, or p-value annotation is included on the density plot. Default to TRUE.
##' @param cex.legend,cex.annotation numeric. Relative expansion of legend and annotation test in density plot.
##' @import igraph
NULL

##' @rdname pathway.structure.plot
##' @examples
##'
##' #generate example graphs
##' library("igraph")
##' graph <- make_graph(unlist(lapply(letters, function(x) rep(x, 2)))[2:51], directed = TRUE)
##' plot(graph)
##' 
##' #Example states
##' test.structure(graph, "a", "c")
##' source <- sample(letters, 10)
##' target <- sample(letters, 10)
##' V(graph)$color <- ifelse(names(V(graph)) %in% source, "lightblue", "grey75")  
##' V(graph)$color <- ifelse(names(V(graph)) %in% target, "palevioletred", V(graph)$color)  
##' V(graph)$color <- ifelse(names(V(graph)) %in% intersect(source, target), "mediumpurple2", V(graph)$color)  
##' plot(graph, layout = layout.fruchterman.reingold, vertex.color= V(graph)$color, vertex.label.family = "mono", vertex.size = 10, vertex.label.color = "black", vertex.frame.color= "grey50", main = g$name)
##' matrix.structure(graph, source, target)
##' table(matrix.structure(graph, source, target))
##' perm_table <- permutation.structure(graph, source, target, letters)
##' 
##' #Plot Raw Output of Permutation Test
##' pathway_perm_plot_raw(perm_table)
##' @export
pathway_perm_plot_raw <- function(table_obj, main = NULL, sub = NULL, xlab="up events", ylab="down events", cex.lab = 1, col = "black", cex = 1, pch = 1){
  plot(table_obj$exp$up, table_obj$exp$down, main = main, sub = sub, xlab = xlab, ylab = ylab, cex.lab = cex.lab, col = col, cex = cex, pch = pch)
  points(table_obj$obs$up, table_obj$obs$down, col="red")
  legend("topright", fill=c("red", "black"), legend=c("observed", "expected"))
}

##' @rdname pathway.structure.plot
##' @examples
##'
##' #Plot Density Output of Permutation Test
##' pathway_perm_plot_density(perm_table)
##' @export
pathway_perm_plot_density <- function(table_obj, main = NULL, xlab = "up - down events", ylab = "density", cex.lab = 1, guide = TRUE, legend = TRUE, annotate = TRUE, cex.legend = 1, cex.annotation = 1, col="black"){
  if(is.null(main)) main <- ""
  plot(density(table_obj$exp$up - table_obj$exp$down), main = main, xlab = xlab, ylab = ylab, cex.lab = cex.lab, col = col)
  if(guide){
    abline(v=table_obj$obs$up - table_obj$obs$down, col="red")
    abline(v=quantile(table_obj$exp$up - table_obj$exp$down, 0.025), col="grey50")
    abline(v=quantile(table_obj$exp$up - table_obj$exp$down, 0.05), col="grey75")
    abline(v=quantile(table_obj$exp$up - table_obj$exp$down, 0.95), col="grey75")
    abline(v=quantile(table_obj$exp$up - table_obj$exp$down, 0.975), col="grey50")
    if(legend){
      legend("right", fill=c("red", "grey75", "grey50"), legend=c("observed", "90% interval", "95% interval"), xjust = 1, cex = cex.legend)
    }
  }
  if(annotate){
    text(min(table_obj$exp$up - table_obj$exp$down), max(density(table_obj$exp$up - table_obj$exp$down)$y), labels = paste("\n\nemp p-val\ndownstream\n", table_obj$down), cex = cex.annotation)
    text(max(table_obj$exp$up - table_obj$exp$down), max(density(table_obj$exp$up - table_obj$exp$down)$y), labels = paste("\n\nemp p-val\nupstream\n", table_obj$up), cex = cex.annotation)
  }
}

##' @rdname pathway.structure.plot
##' @examples
##'
##' #Plot Density Output of Permutation Test
##' perm_table <- permutation.structure(graph, source, target, letters)
##' pathway_perm_plot_density(perm_table)
##' perm_table2 <- permutation.structure(graph, source, target, letters, fixed_interval = T)
##' pathway_perm_lines_density(perm_table2, col="lightblue")
##' @export
pathway_perm_lines_density <- function(table_obj, col="black"){
  lines(density(table_obj$exp$up - table_obj$exp$down), col = col)
}

