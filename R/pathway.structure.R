##' @name Tests of Pathway Structure in iGraph
##' @rdname pathway_test
##'
##' @title Extensions to iGraph for Pathway Structure
##'
##' @description Functions to compute the pathway structure of different states within a network with testing pathway structure with permutation analysis. Uses shortest paths to compute the structure between nodes of differing states. Shortest paths are computed in advance where possible to reduce computational redundancy.
##'
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted as long as a shortest path can be computed.
##' @param target string: vector of target states for testing pathway stucture. Must be a single string for target_node. These are cross referenced against V(graph)$name.
##' @param source string: vector of source states for testing pathway stucture. Must be a single string for source_node. These are cross referenced against V(graph)$name.
##' @param universe string vector of potential nodes to be assigned target and\/or source states. This may be V(graph)$name, and subset thereof, or a larger pool of nodes to assign states for permutation analysis.
##' @param shortest.paths.in Defaults to NULL leading to computing the shortest paths from the input graph where necessary, these may be given as computed in advance (or passed from higher functions) to reduce computational redundancy.
##' @param reps scalar numeric. Number of permutations to statistically test the structure of the network.
##' @import igraph
NULL

##' @rdname pathway_test
##' @examples
##'
##' #generate example graphs
##' library("igraph")
##' g1 <- make_ring(10)
##' V(g1)$name <- letters[1:10]
##' g2 <- make_star(10)
##' V(g2)$name <- letters[1:10]
##'
##' #test pathway structure between two points
##' pathway_test(g1, "a", "c")
##' pathway_test(g2, "a", "c")
##' @export
pathway_test <- function(graph, target_node, source_node, shortest.paths.in = NULL){
  if(is.null(shortest.paths.in)) shortest.paths.in <- shortest.paths(graph, mode="in")
  target_to_source_vec <- shortest.paths.in[match(source_node, rownames(shortest.paths.in)), match(target_node, colnames(shortest.paths.in))]
  source_to_target_vec <- shortest.paths.in[match(target_node, rownames(shortest.paths.in)), match(source_node, colnames(shortest.paths.in))]
  if(target_to_source_vec==0) target_to_source_vec <- Inf
  if(source_to_target_vec==0) source_to_target_vec <- Inf
  source_status <- ifelse(target_to_source_vec > source_to_target_vec, "down", ifelse(target_to_source_vec < source_to_target_vec, "up", "loop"))
  return(source_status)
}

##' @rdname pathway_test_matrix
##' @examples
##'
##' #test pathway structure between two vectors of points
##' pathway_test_matrix(g1, letters[1:5], letters[5:7])
##' pathway_test_matrix(g2, letters[1:5], letters[5:7])
##' @export
pathway_test_matrix <- function(graph, target_vec, source_vec, pathway_nodes=NULL, shortest.paths.in = NULL){
  if(is.null(shortest.paths.in)) shortest.paths.in <- shortest.paths(graph, mode="in")
  if(is.null(pathway_nodes)) pathway_nodes <- names(V(graph))
  overlap <- intersect(target_vec, source_vec)
  target_vec <- setdiff(target_vec, overlap)
  source_vec <- setdiff(source_vec, overlap)
  test_matrix <- matrix(NA, nrow=length(target_vec[target_vec %in% pathway_nodes]),
                        ncol=length(source_vec[source_vec %in% pathway_nodes]))
  dim(test_matrix)
  for(ii in 1:length(target_vec[target_vec %in% pathway_nodes])){
    for(jj in 1:length(source_vec[source_vec %in% pathway_nodes])){
      test_matrix[ii,jj] <- pathway_test(graph, target_vec[target_vec %in% pathway_nodes][ii],
                                          source_vec[source_vec %in% pathway_nodes][jj],
                                         shortest.paths.in = shortest.paths.in)
    }
    print(ii)
  }
  rownames(test_matrix) <- target_vec[target_vec %in% pathway_nodes]
  colnames(test_matrix) <- source_vec[source_vec %in% pathway_nodes]
  print(table(test_matrix))
  return(test_matrix)
}

##' @rdname test_permutation
##' @examples
##'
##' #test pathway structure between two vectors of points with permutations
##' pathway_test_permutation(g1, letters[1:5], letters[5:7], letters)
##' pathway_test_permutation(g2, letters[1:5], letters[5:7], letters)
##' @export
pathway_test_permutation <- function(graph, target, source, universe, shortest.paths.in = NULL, reps=1000){
  if(is.null(shortest.paths.in)) shortest.paths.in <- shortest.paths(graph, mode="in")
  hits_up_minus_down <- rep(NA, reps)
  for(ii in 1:reps){
    target_sim <- sample(universe, length(target))
    source_sim <- sample(universe, length(source))
    test_matrix <- pathway_test_matrix(graph, target_sim, source_sim, shortest.paths.in = shortest.paths.in)
    hits <- as.list(table(test_matrix))
    hits_up_minus_down[ii] <- ifelse(is.null(hits$up), 0, hits$up) - ifelse(is.null(hits$down), 0, hits$down)
    print(ii)
  }
  test_matrix <- pathway_test_matrix(graph, target, source, shortest.paths.in = shortest.paths.in)
  hits <- as.list(table(test_matrix))
  hits$up - hits$down
  abline(v=hits$up - hits$down)
  print(paste(sum(hits$up - hits$down < hits_up_minus_down) / length(hits_up_minus_down), "target downstream"))
  print(paste(sum(hits_up_minus_down < hits$up - hits$down) / length(hits_up_minus_down), "target upstream"))
  table_obj <- list()
  table_obj$down <- sum(hits$up - hits$down < hits_up_minus_down) / length(hits_up_minus_down)
  table_obj$up <- sum(hits_up_minus_down < hits$up - hits$down) / length(hits_up_minus_down)
  return(table_obj)
}
