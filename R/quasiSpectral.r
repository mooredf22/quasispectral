#' Compare mutant to wild type spectral count data, using "QuasiSeq" methodology
#'
#' @param spectralDataAll data frame of gene/protein identifiers (geneName) and their spectral counts in mutant and wildtype groups.
#' @param n.mut number of mutant animals
#' @param n.wt number of wild type animals
#' @return A dataframe with the log2 ratios, p-values, and q-values. Also eta0, the estimated proportion of null samples

quasiSpectral <- function(spectralDataAll, n.mut, n.wt) {
  # the file "spectralDataAll" must be formatted as follows:
  #have gene names in the first column, named "geneName"
  # next n.mut columns must contain spectral counts for mutant (or comparative) group
  # next n.wt columns must contain specctral counts for wild type (or control) group

  spectralDataUse <- spectralDataAll[,1 + 1:(n.mut + n.wt)]  # just use the data
  spectralDataUse[is.na(spectralDataUse)] <- 0
  trt <- c(rep(0, n.mut), rep(1, n.wt))    # assume mutants are first, followed by wild type

  # set up design matrix for QL.fit
  design.list<-vector("list",2)
  design.list[[1]]<-model.matrix(~as.factor(trt))
  design.list[[2]]<- matrix(1,length(trt))

  # define offset, which allows one to compare columns adjusted for the total counts
  log.offset <- log(apply(spectralDataUse,2,sum))   # offset is log of column totals

  #log.offset<-log(apply(spectralDataUse,2,quantile,.75))

  ### Analyze using QL, QLShrink and QLSpline methods applied to quasi-Poisson model

  fit <- QL.fit(spectralDataUse, design.list,log.offset=log.offset, Model="Poisson")
  results <- QL.results(fit)


  #######
  # recommend using Model="Poisson" instead of the following:
  #fit.nb <- QL.fit(spectralDataUse, design.list,log.offset=log.offset, Model="NegBin")
  #results.nb <- QL.results(fit.nb)

  #results <- results.nb
  #fit <- fit.nb


  ### How many significant genes at FDR=.05 from QLSpline method?
  #apply(results$Q.values[[3]]<.05,2,sum)

  ### Indexes for Top 10 most significant genes from QLSpline method
  #head(order(results$P.values[[3]]), 10)

  # convert coefficients to log base 2
  coef.main <- fit$coefficients[,2] / log(2)
  #set maximum and minimum
  max.value <- 16           # 2^4
  coef.main[coef.main > 16] <- 16
  coef.main[coef.main < -16] <- -16

  p.values <- as.numeric(results$P.values[[3]])            # QLSpline method is third component
  q.values <- results$Q.values[[3]]

  p.values.bonf <- p.adjust(p=p.values, method="bonferroni")
  p.values.holm <- p.adjust(p=p.values, method="holm")
  q.values.BH <- p.adjust(p=p.values, method="BH") # Benjamini and Hochberg

  library(fdrtool)   # must be downloaded and installed
  fdrout <- fdrtool(p.values, statistic="pvalue")
  q.values.strimmer <- fdrout$qval
  eta0 <- fdrout$param[3]  # proportion of null samples in population
  eta0

  geneName <- as.character(spectralDataAll$geneName)
  QLfit <- data.frame(geneName, spectralDataUse, coef.main, p.values, p.values.bonf, p.values.holm, q.values.BH, q.values.strimmer)

  results <- list(QLfit=QLfit, eta0=eta0)
  results
}
