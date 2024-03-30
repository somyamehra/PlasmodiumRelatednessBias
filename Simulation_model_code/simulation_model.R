library(Pv3Rs)
library(paneljudge)
library(igraph)
library(dplyr)
library(tidyr)

# Simulate one meiosis, alternating geometric process governing crossovers
# M = mean segment length (M=1 for independence model)
# n_m = number of markers
# parents = 2-vector of parent lineages
simulate_one_meiosis <- function(parents, M, n_m) {
  rep(rep(sample(parents, 2), length.out=n_m),  rgeom(n_m, 1/M)+1)[1:n_m]
}

# Construct a population-level relatedness graph by pasting together
# subgraphs allowing for clone, stranger and sibling relationships
# If n_samples is not a multiple of MOIs, rounded up to nearest multiple
simulate_pop_level_graph <- function(MOIs, n_samples) {
  # sample many graphs
  clusters <- lapply(seq(0, n_samples-1, length(MOIs)), function(i) {
    g <- sample_RG(MOIs) # Sample graph
    # Override genotype names created by sample_RG()
    V(g)$name <- paste0("g", i+seq_along(MOIs)) 
    # Override group created by sample_RG s.t. coloured by cluster by plot_RG
    V(g)$group <- i+seq_along(MOIs) 
    return(g)
  })
  
  subgraph_ids <- lapply(clusters, function(g) {as_ids(V(g))})
  subgraph_ids <- setNames(rep(1:length(subgraph_ids), sapply(subgraph_ids, length)), 
                           sapply(subgraph_ids, c))
  
  # stick many graphs together following
  # https://stackoverflow.com/questions/48000074/combining-merging-two-graphs-in-igraph
  attrs <- do.call(rbind, lapply(clusters, function(g) igraph::as_data_frame(g, "vertices")))
  el <- do.call(rbind, lapply(clusters, igraph::as_data_frame))
  pop_graph <- graph_from_data_frame(el, directed = FALSE, vertices = attrs)
  
  return(list(pop_graph=pop_graph, subgraph_ids=subgraph_ids))
}

# Sample lineages many generations after founders -- each simulated parasite
# consists of a mosaic of up to n_lineages lineages
simulate_lineage_n_mosaic <- function(pop_graph, n_markers, n_lineages) {
  Lineages_per_parasite_per_marker <- sapply(1:n_markers, function(i) {
    sample_lineages(pop_graph, 1, outbred = "F", lineage_count = n_lineages)})
  
  rownames(Lineages_per_parasite_per_marker) <- rownames(sample_lineages(pop_graph, 1))
  colnames(Lineages_per_parasite_per_marker) <- paste0("m", 1:n_markers)
  
  return(Lineages_per_parasite_per_marker)
}

# Sample lineages a single generation after founders -- each simulated parasite 
# consists of a mosaic of up to 2 lineages --- with linkage
# RG = relationship graph with stranger/sibling/clonal edges
# n_m = number of markers
# lineage_freqs = vector of lineage frequencies
# lineage_groups = grouping of lineages into sub-clusters
# prob_cotrans =  probability of sampling parents from the same lineage group (without replacement)
#                 vs the population at large (with replacement)
#                 (sample without replacement from lineage group to weight sibling-sibling mating)
# M = mean no. loci spanning inherited segment
simulate_lineage_2_mosaic <- function(RG, n_m, lineage_freqs, lineage_groups, prob_cotrans, M) {
  
  # randomly assign each isolate in the previous generation a lineage code
  lineages <- sample(generate_lineages(length(lineage_freqs)), 
                     size = length(lineage_freqs), replace = FALSE)
  lineages_per_marker <- array(NA, dim = c(igraph::vcount(RG), n_m), 
                               dimnames = list(igraph::V(RG)$name, paste0("m", 1:n_m)))
  
  # iterate over components of relationship graph (sibling/clonal)
  gs_per_rcs <- igraph::groups(igraph::components(RG))
  
  # record which components are subject to enriched cotransmission model
  enriched_cotrans <- numeric(length(gs_per_rcs))
  
  for (rc_i in 1:length(gs_per_rcs)) {
    
    # sample 2 parents based on lineage frequencies
    if (runif(1)<prob_cotrans) {
      # note that component rc_i has parents sampled from the same subgraph in
      # the previous generation
      enriched_cotrans[rc_i] <- 1
      
      # sample a lineage group, weighted by lineage frequencies
      parent_group <- sample(x = lineage_groups, size = 1, prob = lineage_freqs)
      
      if (sum(lineage_groups==parent_group)<2) {
        print("Error: subgraph of size <2")
        break
      }
      
      # sample a pair of parents uniformly at random WITHOUT replacement within the sampled group
      parents <- sample(x = lineages[lineage_groups==parent_group], size=2, replace=F,
                        prob = lineage_freqs[lineage_groups==parent_group])
    } else {
      # sample parents uniformly at random WITH replacement from the population at large 
      parents <- sample(x = lineages, size = 2, replace = T, prob = lineage_freqs)
    }
    
    # extract subgraph
    rc <- igraph::induced_subgraph(RG, gs_per_rcs[[rc_i]])
    rc_gs <- igraph::V(rc)$name
    rc_size <- igraph::vcount(rc)
    rc_relationship_types <- unique(igraph::E(rc)$weight)
    
    # delete sibling edges edges from graph
    rcd <- igraph::delete_edges(rc, igraph::E(rc)[igraph::E(rc)$weight == 0.5])
    rcdC <- igraph::components(rcd)
    
    # sample genotypes for each clone with linkage
    sib_gts <- sapply(1:rcdC$no, function(x) {simulate_one_meiosis(parents, M, n_m)})
    
    lineages_per_marker[rc_gs, ] <- t(sib_gts[,rcdC$membership])
  }
  
  return(list(lineages_per_marker=lineages_per_marker, lineages=lineages,
              # record which sibling/clonal components have been subject to enriched cotransmission
              enriched_cotrans_components=Filter(function(x) length(x)>1,
                                                 gs_per_rcs[which(enriched_cotrans==1)])))
}

# Simulate 1 generation of recombination, given:
#   pop_graph: population-level relatedness graph with sibling, clone and stranger
#              relationships relative to previous generation
#   prev_gen_lineages: mosaic of original founder lineages for previous generation
simulate_next_generation <- function(pop_graph, prev_gen_lineages, prev_gen_lineage_groups,
                                     n_markers, prob_cotrans, M) {
  # Treat each isolate in previous generation as a distinct lineage
  n_lineages_prev <- nrow(prev_gen_lineages)
  
  # Each isolate in the current generation is a 2-mosaic of isolates in previous generation
  # First, simulate a 2-mosaic using arbitrary lineage codes, but n_lineages_prev lineages
  curr_gen_2_mosaic <- 
    simulate_lineage_2_mosaic(pop_graph, n_markers, 
                              rep(1/n_lineages_prev, n_lineages_prev), 
                              prev_gen_lineage_groups, prob_cotrans, M)
  
  # extract mapping between each isolate in the previous generation and a lineage code
  lineage_code_prev_gen <- setNames(1:n_lineages_prev, curr_gen_2_mosaic[["lineages"]])
  
  # Recover each isolate in current generation as a 2-mosaic of isolates in previous generation
  curr_gen_2_mosaic_prev_gen <- apply(curr_gen_2_mosaic[["lineages_per_marker"]], 2, 
                                      function(x) as.numeric(lineage_code_prev_gen[x]))
  
  # Recover each isolate in current generation as an n-mosaic of the original founder lineages
  curr_gen_founder_mosaic <- sapply(paste0("m", 1:n_markers), function(marker) {
    prev_gen_lineages[curr_gen_2_mosaic_prev_gen[, marker], marker]})
  rownames(curr_gen_founder_mosaic) <- paste0("g", 1:nrow(curr_gen_founder_mosaic))
  
  return(list(curr_gen_founder_mosaic=curr_gen_founder_mosaic,
              enriched_cotrans_components=curr_gen_2_mosaic[["enriched_cotrans_components"]]))
}

map_lineages_to_genotypes <- function(lineage_mosaic, alleles_per_lineage_per_marker) {
  for (pg in rownames(lineage_mosaic)) {
    for (marker in colnames(lineage_mosaic)) {
      lineage_mosaic[pg, marker] <- 
        alleles_per_lineage_per_marker[lineage_mosaic[pg, marker], marker]
    }
  }
  class(lineage_mosaic) <- "numeric"
  return(lineage_mosaic)
}

# ================ SIMULATE POPULATION-LEVEL RELATEDNESS GRAPH ==================

# Alleles for each founder lineage ~ Bernoulli(p)
simulate_inbreeding <- function(max_clonal_cluster, n_markers, n_lineages, 
                                n_samples, n_generations, prob_cotrans, M, p) {
  
  MOIs <- rep(1, max_clonal_cluster) # monoclonal only
  
  # Simulate n_generation relatedness graphs (clones/sibling/stranger edges) 
  relatedness_graphs <- lapply(1:n_generations, function(i) {
    simulate_pop_level_graph(MOIs, n_samples)})
  
  founder_lineage_mosaics <- enriched_cotrans_components <- list()
  
  # Given population graph, sample an n-mosaic of founder lineages for each isolate in gen 1
  generation_0_n_mosiacs <- 
    matrix(sample(Pv3Rs::generate_lineages(n_lineages), n_markers*n_lineages, replace=TRUE), 
           nrow=n_lineages, dimnames=list(paste0("g", 1:n_lineages), paste0("m", 1:n_markers)))
  
  founder_lineage_mosaics[[1]] <- 
    simulate_next_generation(relatedness_graphs[[1]][["pop_graph"]], generation_0_n_mosiacs, 
                             rep(1, n_lineages), n_markers, prob_cotrans, M)[["curr_gen_founder_mosaic"]]
  
  # Simulate successive generations as 2-mosaics of previous generations
  for (i in 2:n_generations) {
    next_generation <- 
      simulate_next_generation(relatedness_graphs[[i]][["pop_graph"]], founder_lineage_mosaics[[i-1]], 
                               relatedness_graphs[[i-1]][["subgraph_ids"]], n_markers, prob_cotrans, M)
    
    founder_lineage_mosaics[[i]] <- next_generation[["curr_gen_founder_mosaic"]]
    enriched_cotrans_components[[i]] <- next_generation[["enriched_cotrans_components"]]
  }
  
  # For each founder lineage, sample an allele per marker, assuming markers are biallelic
  founder_lineages <- unique(as.vector(founder_lineage_mosaics[[1]]))
  alleles_per_lineage_per_marker <- t(sapply(founder_lineages, function(lineage) {
    as.numeric(runif(n_markers)<p)}))
  colnames(alleles_per_lineage_per_marker) <- paste0("m", 1:n_markers)
  
  # Assign allelic states to each generation, ignoring genotyping error
  true_allelic_states <- lapply(founder_lineage_mosaics, map_lineages_to_genotypes, 
                                alleles_per_lineage_per_marker)
  
  # ================ CALCULATE ALLELE FREQUENCIES, P(IBS) ==================
  
  # Calculate P(IBS | not IBD)
  proportion_pairs_ibd <- lapply(founder_lineage_mosaics, function(mosaic) {
    apply(mosaic, 2, function(x) {sum(table(x)^2)/length(x)^2})})
  
  proportion_pairs_ibs <- lapply(true_allelic_states, function(mosaic) {
    apply(mosaic, 2, function(x) {sum(table(x)^2)/length(x)^2})})
  
  proportion_pairs_ibc <- lapply(1:n_generations, function(i) {
    (proportion_pairs_ibs[[i]]-proportion_pairs_ibd[[i]])/(1-proportion_pairs_ibd[[i]])})
  
  founder_prob_ibs <- apply(alleles_per_lineage_per_marker, 2, function(x) {
    sum(table(x)^2)/length(x)^2})
  
  locus_metrics <- lapply(1:n_generations, function(i) {
    data.frame(prop_pair_ibd=proportion_pairs_ibd[[i]],
               prop_pair_ibs=proportion_pairs_ibs[[i]],
               prop_pair_ibc=proportion_pairs_ibc[[i]],
               founder_prob_ibs=founder_prob_ibs)})
  
  ALLELE_CODES <- sort(unique(as.numeric(alleles_per_lineage_per_marker)))
  
  allele_freq_founder <- apply(alleles_per_lineage_per_marker, 2, function(x) {
    sapply(ALLELE_CODES, function(allele) {
      sum(x==allele, na.rm=TRUE)})/sum(!is.na(x))}) %>% t
  
  return(list(locus_metrics=locus_metrics,
              genotype_matrix=true_allelic_states,
              ancestry_matrix=founder_lineage_mosaics,
              relatedness_graphs=relatedness_graphs,
              enriched_cotrans_components=enriched_cotrans_components,
              allele_freq_founder=allele_freq_founder,
              max_clonal_cluster=max_clonal_cluster,
              n_lineages=n_lineages, M=M, 
              n_markers=n_markers,
              ALLELE_CODES=ALLELE_CODES, p=p))
}

# Classify pairwise relationships (siblings or clones +/- enriched cotransmission)
# Except generation 1
classify_pairs <- function(relationship_graph, enriched_cotrans_components) {
  
  RG_edges <- get.data.frame(relationship_graph[["pop_graph"]], what="edges") %>%
    transmute(p1=pmin(from, to), p2=pmax(from, to), 
              rship=ifelse(weight==1, "clone", ifelse(weight==0.5, "sibling", "stranger")))
  
  if (length(enriched_cotrans_components)==0) {
    RG_edges$enriched_cotrans <- 0
  } else {
    enriched_cotrans_edges <- 
      lapply(enriched_cotrans_components, function(x) t(combn(x, 2))) %>% 
      do.call(rbind, .) %>% as.data.frame %>% 
      transmute(p1=pmin(V1, V2), p2=pmax(V1, V2), enriched_cotrans=1)
    RG_edges <- merge(RG_edges, enriched_cotrans_edges, all.x=TRUE) %>%
      mutate(enriched_cotrans=ifelse(is.na(enriched_cotrans), 0, 1))
  }
  
  return(RG_edges)
}


# loci_to_keep must be a subset of SORTED marker names (mX)
estimate_relatedness_sim <- function(inbreeding_multigen, loci_to_keep, generations) {
  
  marker_names <- rownames(inbreeding_multigen[["locus_metrics"]][[1]])
  n_generations <- length(inbreeding_multigen[["locus_metrics"]]) 
  
  if (length(setdiff(loci_to_keep, marker_names))>0) {
    print("Chosen loci do not exist!")
    break
  }
  
  M <- inbreeding_multigen[["M"]]
  ALLELE_CODES <- inbreeding_multigen[["ALLELE_CODES"]]
  true_allelic_states <- inbreeding_multigen[["genotype_matrix"]]
  founder_lineage_mosaics <- inbreeding_multigen[["ancestry_matrix"]]
  proportion_pairs_ibs <- lapply(inbreeding_multigen[["locus_metrics"]], 
                                 function(x) {setNames(x$prop_pair_ibs, marker_names)})
  proportion_pairs_ibc <- lapply(inbreeding_multigen[["locus_metrics"]], 
                                 function(x) {setNames(x$prop_pair_ibc, marker_names)})
  
  ibd_ibs_est <- list()
  
  for (i in generations) {
    polymorphic <- names(which(apply(true_allelic_states[[i]], 2, 
                                     function(x) {length(unique(x[!is.na(x)]))})>1))
    polymorphic_loci_to_keep <- intersect(polymorphic, loci_to_keep)
    
    distances <- c(-diff(grep("m", polymorphic_loci_to_keep))*log(1-2/M), Inf)
    
    allelic_states <- true_allelic_states[[i]][, polymorphic_loci_to_keep]
    lineage_states <- founder_lineage_mosaics[[i]][, polymorphic_loci_to_keep]
    
    allele_freq_naive <- apply(allelic_states, 2, function(x) {
      sapply(ALLELE_CODES, function(allele) {
        sum(x==allele, na.rm=TRUE)})/sum(!is.na(x))}) %>% t
    
    sample_pairs <- combn(rownames(allelic_states), m=2) %>% t %>% as.data.frame
    
    summary_metrics <- sapply(1:nrow(sample_pairs), function(j) {
      ibs <- concordance(allelic_states[sample_pairs[j, 1],],
                         allelic_states[sample_pairs[j, 2],])
      
      realised_IBD <- concordance(lineage_states[sample_pairs[j, 1],],
                                  lineage_states[sample_pairs[j, 2],])
      
      naive_r <- mle_r_poisson_binomial(proportion_pairs_ibs[[i]][polymorphic_loci_to_keep], 
                                        allelic_states[sample_pairs[j, 1],],
                                        allelic_states[sample_pairs[j, 2],])
      
      corrected_r <- mle_r_poisson_binomial(proportion_pairs_ibc[[i]][polymorphic_loci_to_keep], 
                                            allelic_states[sample_pairs[j, 1],],
                                            allelic_states[sample_pairs[j, 2],])
      
      rhat_naive_indep <- estimate_ibd_paneljudge(allele_freq_naive,
                                                  allelic_states[sample_pairs[j, 1],],
                                                  allelic_states[sample_pairs[j, 2],],
                                                  rep(Inf, length(polymorphic_loci_to_keep)))
      
      rhat_naive_hmm <- estimate_ibd_paneljudge(allele_freq_naive,
                                                allelic_states[sample_pairs[j, 1],],
                                                allelic_states[sample_pairs[j, 2],],
                                                distances)
      
      N_comparable <- sum(!is.na(allelic_states[sample_pairs[j, 1],]) &
                            !is.na(allelic_states[sample_pairs[j, 2],]))
      
      
      return(c(ibs=ibs, realised_IBD=realised_IBD, 
               naive_r=naive_r, corrected_r=corrected_r,
               rhat_naive_indep=rhat_naive_indep[["rhat"]], 
               rhat_naive_hmm=rhat_naive_hmm[["rhat"]],
               khat_naive_indep=rhat_naive_indep[["khat"]], 
               khat_naive_hmm=rhat_naive_hmm[["khat"]],
               N_comparable=N_comparable))
    })
    
    ibd_ibs_est[[i]] <- cbind(sample_pairs, t(summary_metrics)) %>% as.data.frame
  }
  return(list(ibd_ibs_est=bind_rows(ibd_ibs_est, .id="generation"),
              loci_to_keep=loci_to_keep))
}