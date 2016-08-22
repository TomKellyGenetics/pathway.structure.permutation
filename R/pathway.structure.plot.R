##' @name Plots for Tests of Pathway Structure in iGraph
##' @rdname pathway.structure.plot
##'
##' @title Extensions to iGraph for Plotting Pathway Structure
##'
##' @description Functions to plot outcomes for testing pathway structure with permutation analysis.
##'
##'
##' @param table_obj List. Output of \code{\link[permutation.structure]{pathway.structure.permuation}}
##' @param main,sub,xlab,ylab,col,cex,pch Arguments to passs on to \code{\link[plot]{graphics}} or \code{\link[plot.density]{stats}}
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
pathway_perm_plot_raw <- function(table_obj, main = NULL, sub = NULL, xlab="up events", ylab="down events", col = "black", cex = 1, pch = 19){
  plot(test_perm$exp$up, test_perm$exp$down, main = main, sub = sub, xlab = xlab, ylab = ylab, col = col, cex = cex, pch = pch)
  points(test_perm$obs$up, test_perm$obs$down, col="red")
  legend("topright", fill=c("red", "black"), legend=c("observed", "expected"))
}

##' @rdname pathway.structure.plot
##' @examples
##'
##' #Plot Density Output of Permutation Test
##' pathway_perm_plot_density(perm_table)
##' @export
pathway_perm_plot_density <- function(table_obj, main = NULL, xlab = "up - down events", ylab = "density"){
  plot(density(test_perm$exp$up - test_perm$exp$down), main = main, xlab = xlab, ylab = ylab)
  abline(v=test_perm$obs$up - test_perm$obs$down, col="red")
  abline(v=quantile(test_perm$exp$up - test_perm$exp$down, 0.025), col="grey50")
  abline(v=quantile(test_perm$exp$up - test_perm$exp$down, 0.05), col="grey75")
  abline(v=quantile(test_perm$exp$up - test_perm$exp$down, 0.95), col="grey75")
  abline(v=quantile(test_perm$exp$up - test_perm$exp$down, 0.975), col="grey50")
  text(min(test_perm$exp$up - test_perm$exp$down), max(density(test_perm$exp$up - test_perm$exp$down)$y), labels = paste("emp p-val\ndownstream\n", test_perm$down))
  text(max(test_perm$exp$up - test_perm$exp$down), max(density(test_perm$exp$up - test_perm$exp$down)$y), labels = paste("emp p-val\nupstream\n", test_perm$up))
  legend("right", fill=c("red", "grey75", "grey50"), legend=c("observed", "90% interval", "95% interval"))
}

