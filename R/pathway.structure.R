##' @name Tests of Pathway Structure in iGraph
##' @rdname pathway.structure
##'
##' @title Extensions to iGraph for Pathway Structure
##'
##' @description Functions to compute the pathway structure of different states within a network with testing pathway structure with permutation analysis. Uses shortest paths to compute the structure between nodes of differing states. Shortest paths are computed in advance where possible to reduce computational redundancy.
##'
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted as long as a shortest path can be computed.
##' @param target string: vector of target states for testing pathway structure. Must be a single string for target_node. These are cross referenced against V(graph)$name.
##' @param source string: vector of source states for testing pathway structure. Must be a single string for source_node. These are cross referenced against V(graph)$name.
##' @param universe string vector of potential nodes to be assigned target and\/or source states. This may be V(graph)$name, and subset thereof, or a larger pool of nodes to assign states for permutation analysis.
##' @param shortest.paths.in Defaults to NULL leading to computing the shortest paths from the input graph where necessary, these may be given as computed in advance (or passed from higher functions) to reduce computational redundancy.
##' @param reps scalar numeric. Number of permutations to statistically test the structure of the network.
##' @param fixed_intersect logical. Defaults to FALSE. Whether number of intersecting states is fixed to the same in permutations as the input states. 
##' @import igraph
NULL

##' @rdname pathway.structure
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
##' test.structure(g1, "a", "c")
##' test.structure(g2, "a", "c")
##' @export
test.structure <- function(graph, target_node, source_node, shortest.paths.in = NULL){
  if(is.null(shortest.paths.in)) shortest.paths.in <- shortest.paths(graph, mode="in")
  target_to_source_vec <- shortest.paths.in[match(source_node, rownames(shortest.paths.in)), match(target_node, colnames(shortest.paths.in))]
  source_to_target_vec <- shortest.paths.in[match(target_node, rownames(shortest.paths.in)), match(source_node, colnames(shortest.paths.in))]
  if(is.na(target_to_source_vec)) target_to_source_vec <- 0
  if(is.na(source_to_target_vec)) source_to_target_vec <- 0
  if(target_to_source_vec==0) target_to_source_vec <- Inf
  if(source_to_target_vec==0) source_to_target_vec <- Inf
  source_status <- ifelse(target_to_source_vec > source_to_target_vec, "down", ifelse(target_to_source_vec < source_to_target_vec, "up", "loop"))
  return(source_status)
}

##' @rdname pathway.structure
##' @examples
##'
##' #test pathway structure between two vectors of points
##' matrix.structure(g1, letters[1:5], letters[5:7])
##' matrix.structure(g2, letters[1:5], letters[5:7])
##' @export
matrix.structure <- function(graph, target_vec, source_vec, pathway_structure_nodes=NULL, shortest.paths.in = NULL){
  if(is.null(shortest.paths.in)) shortest.paths.in <- shortest.paths(graph, mode="in")
  if(is.null(pathway_structure_nodes)) pathway_structure_nodes <- names(V(graph))
  overlap <- intersect(target_vec, source_vec)
  target_vec <- setdiff(target_vec, overlap)
  source_vec <- setdiff(source_vec, overlap)
  test_matrix <- matrix(NA, nrow=length(target_vec[target_vec %in% pathway_structure_nodes]),
                        ncol=length(source_vec[source_vec %in% pathway_structure_nodes]))
  dim(test_matrix)
  if(nrow(test_matrix)*ncol(test_matrix)==0){
    test_matrix <- "no paths"
    warning("missing simulated states to draw paths")
  } else {
    for(ii in 1:length(target_vec[target_vec %in% pathway_structure_nodes])){
      for(jj in 1:length(source_vec[source_vec %in% pathway_structure_nodes])){
        test_matrix[ii,jj] <- test.structure(graph, target_vec[target_vec %in% pathway_structure_nodes][ii],
                                             source_vec[source_vec %in% pathway_structure_nodes][jj],
                                             shortest.paths.in = shortest.paths.in)
      }
      print(ii)
    }
    rownames(test_matrix) <- target_vec[target_vec %in% pathway_structure_nodes]
    colnames(test_matrix) <- source_vec[source_vec %in% pathway_structure_nodes]
  }
  print(table(test_matrix))
  return(test_matrix)
}

##' @rdname pathway.structure
##' @examples
##'
##' #test pathway structure between two vectors of points with permutations
##' permutation.structure(g1, letters[1:5], letters[5:7], letters)
##' permutation.structure(g2, letters[1:5], letters[5:7], letters)
##' @export
permutation.structure <- function(graph, target, source, universe, shortest.paths.in = NULL, reps = 1000, fixed_intersect = FALSE){
  if(is.null(shortest.paths.in)) shortest.paths.in <- shortest.paths(graph, mode="in")
  hits_up_minus_down <- rep(NA, reps)
  hits_up <- rep(NA, reps)
  hits_down <- rep(NA, reps)
  for(ii in 1:reps){
    if(fixed_intersect){
      #simulate intersect of given size
      overlap_sim <- sample(universe, length(intersect(target, source)))
      #sample remaining
      target_sim <- sample(setdiff(universe, overlap_sim), length(setdiff(target, overlap_sim)))
      source_sim <- sample(setdiff(universe, union(overlap_sim, target_sim)), length(setdiff(source, overlap_sim)))
    } else {
      #sample each independently
      target_sim <- sample(universe, length(target))
      source_sim <- sample(universe, length(source))
    }
    test_matrix <- matrix.structure(graph, target_sim, source_sim, shortest.paths.in = shortest.paths.in)
    hits <- as.list(table(test_matrix))
    hits_down[ii] <- ifelse(is.null(hits$down), 0, hits$down)
    hits_up[ii] <- ifelse(is.null(hits$up), 0, hits$up)
    hits_up_minus_down[ii] <- hits_up[ii] - hits_down[ii]
    print(ii)
  }
  test_matrix <- matrix.structure(graph, target, source, shortest.paths.in = shortest.paths.in)
  hits <- as.list(table(test_matrix))
  if(is.null(hits$up)) hits$up <- 0
  if(is.null(hits$down)) hits$down <- 0
  hits$up - hits$down
  print(paste(sum(hits$up - hits$down < hits_up_minus_down) / length(hits_up_minus_down), "target upstream"))
  print(paste(sum(hits_up_minus_down < hits$up - hits$down) / length(hits_up_minus_down), "target downstream"))
  table_obj <- list()
  table_obj$up <- sum(hits$up - hits$down < hits_up_minus_down) / length(hits_up_minus_down)
  table_obj$down <- sum(hits_up_minus_down < hits$up - hits$down) / length(hits_up_minus_down)
  table_obj$obs <- list(hits$down, hits$up)
  table_obj$exp <- list(hits_down, hits_up)
  names(table_obj$exp) <- names(table_obj$obs) <- c("down", "up")
  return(table_obj)
}