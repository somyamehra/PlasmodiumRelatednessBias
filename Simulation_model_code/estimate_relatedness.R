library(paneljudge)
library(poisbinom)

# estimate relatedness under a simple Poisson binomial model, assuming
# no linkage between loci (MLE)
mle_r_poisson_binomial <- function(prob_ibs_given_not_ibd, sample1, sample2) {
  
  comparable <- which(!is.na(sample1) & !is.na(sample2))
  
  prob_ibs_given_not_ibd <- prob_ibs_given_not_ibd[comparable]
  
  ibs_states <- as.numeric(sample1[comparable]==sample2[comparable])
  
  log_likelihood <- function(r) {
    -sum(ibs_states * log(r + prob_ibs_given_not_ibd*(1-r)) + 
           (1-ibs_states) * log((1-prob_ibs_given_not_ibd)*(1-r)))
  }
  
  MLE <- optim(0.25, log_likelihood, lower=0, upper=1, method="Brent")$par
  
  return(MLE)
}

# estimate relatedness using paneljudge, allowing for linkage between loci
estimate_ibd_paneljudge <- function(allele_freq, sample1, sample2, amplicon_distances, rho=1) {
  
  comparable <- which(!is.na(sample1) & !is.na(sample2))
  
  genotype_pair <- cbind(sample1[comparable], sample2[comparable])
  
  krhat <- estimate_r_and_k(fs=allele_freq[comparable,], 
                            ds=amplicon_distances[comparable], 
                            Ys=genotype_pair, epsilon=0, rho=rho,
                            warn_fs = FALSE)
  
  return(krhat)
}

# calculate pairwise concordance -- for calculating realised IBD + IBS
concordance <- function(sample1, sample2) {
  comparable <- which(!is.na(sample1) & !is.na(sample2))
  return(mean(sample1[comparable]==sample2[comparable]))
}

# poisbinom implements the exact DFT method of Hong (2013)
posterior_ibs <- function(ibs_given_ibd, ibs_given_not_ibd, r_vals) {
  n <- length(ibs_given_not_ibd)
  
  posterior_pdf_ibs <- sapply(r_vals, function(r) {
    p <- r*ibs_given_ibd + (1-r)*ibs_given_not_ibd
    dpoisbinom(0:n, p, log_d = FALSE)})
  colnames(posterior_pdf_ibs) <- r_vals
  posterior_pdf_ibs <- posterior_pdf_ibs %>% cbind(IBS=0:n) %>% 
    as.data.frame %>% reshape2::melt(id="IBS") 
  colnames(posterior_pdf_ibs) <- c("IBS", "r", "posterior_pdf")
  
  return(posterior_pdf_ibs)
}

generate_relatedness_summary <- function(genotype_matrix, filtered_pos, rho=1) {
  
  # basic locus metrics
  proportion_pairs_ibs <- apply(genotype_matrix, 2, function(x) {
    sum(table(x)^2)/sum(!is.na(x))^2})
  
  max_allele_code <- max(genotype_matrix, na.rm = TRUE)
  
  allele_freq <- apply(genotype_matrix, 2, function(x) {
    sapply(0:max_allele_code, function(allele) {
      sum(x==allele, na.rm=TRUE)})/sum(!is.na(x))}) %>% t
  
  distances <- c(diff(filtered_pos$V2), Inf)
  distances[distances<0] <- Inf
  
  # ================================ PAIRWISE IBD ================================
  
  sample_pairs <- combn(rownames(genotype_matrix), m=2) %>% t %>% as.data.frame
  
  summary_metrics <- sapply(1:nrow(sample_pairs), function(j) {
    #print(j/nrow(sample_pairs)*100)
    ibs <- concordance(genotype_matrix[sample_pairs[j, 1],],
                       genotype_matrix[sample_pairs[j, 2],])
    
    rhat_indep_ibs <- mle_r_poisson_binomial(proportion_pairs_ibs, 
                                             genotype_matrix[sample_pairs[j, 1],],
                                             genotype_matrix[sample_pairs[j, 2],])
    
    rhat_indep_allele <-  estimate_ibd_paneljudge(allele_freq,
                                                  genotype_matrix[sample_pairs[j, 1],],
                                                  genotype_matrix[sample_pairs[j, 2],],
                                                  rep(Inf, ncol(genotype_matrix)), rho=rho)
    
    rhat_hmm_allele <-  estimate_ibd_paneljudge(allele_freq,
                                                genotype_matrix[sample_pairs[j, 1],],
                                                genotype_matrix[sample_pairs[j, 2],],
                                                distances, rho=rho)
    
    N_comparable <- sum(!is.na(genotype_matrix[sample_pairs[j, 1],]) &
                          !is.na(genotype_matrix[sample_pairs[j, 2],]))
    
    return(c(N_comparable=N_comparable, ibs=ibs, 
             rhat_indep_ibs=rhat_indep_ibs,
             rhat_indep_allele=rhat_indep_allele[["rhat"]], 
             rhat_hmm_allele=rhat_hmm_allele[["rhat"]],
             khat_hmm_allele=rhat_hmm_allele[["khat"]]))
  })
  
  summary_metrics <- cbind(sample_pairs, t(summary_metrics)) %>% as.data.frame
  
  # ========================= POSTERIOR VS EMPIRICAL IBS =========================
  expected_ibs_unrelated <- posterior_ibs(1, proportion_pairs_ibs, 0) %>% dplyr::select(-r)
  
  collated_results <- list(genotype_matrix=genotype_matrix, 
                           filtered_pos=filtered_pos,
                           summary_metrics=summary_metrics, 
                           expected_ibs_unrelated=expected_ibs_unrelated)
  
  return(collated_results)
}

