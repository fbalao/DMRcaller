#' Checks whether the passed parameter is the statistical test
#'
#' @title Validate statistial test
#' @param test the statistical test used to call DMRs (\code{"betareg"} for
#' Beta regression.
#'
#' @author Alessandro Pio Greco and Nicolae Radu Zabet
.validateStatisticalTestReplicates <- function(test){
  .stopIfNotAll(c(!is.null(test), is.character(test), length(test) == 1, test %in% "betareg"),
                " test needs to be one of the following \"betareg\" for beta regression")
}

#' Checks whether the passed parameter is a valid condition for the data
#'
#' @title Validate the condition of the experiment
#' @param condition1 The vector containing the two conditions for the experiment.
#' @param condition2 The vector containing the two conditions for the experiment.
#' @author Alessandro Pio Greco, Nicolae Radu Zabet, Fran Balao
.validateConditionReplicates <- function(condition1, condition2, methylationData){
  .stopIfNotAll(c(!is.null(condition1), is.character(condition1), length(condition1) ==
                    length(grep("readsM", names(mcols(methylationData)))), 
                   !is.null(condition2), is.character(condition2), length(condition2) ==
                   length(grep("readsM", names(mcols(methylationData))))),
                " conditions need to be a vector of characters and the length needs to much the number of replicates")
  .stopIfNotAll(c(length(levels(as.factor(condition1))) == 2,length(levels(as.factor(condition2))) == 2),
                " only two levels are allowed for the condition1 and condition2")
}

#' This function computes the differentially methylated regions between
#' replicates with two conditions.
#'
#' @title Compute DMRs
#' @param methylationData the methylation data containing all the conditions
#' for all the replicates.
#' @param condition1 a vector of strings indicating the conditions for each
#' sample in \code{methylationData}. Two different values are allowed
#' (for the two conditions).
#' @param condition2 a vector of strings indicating the conditions for each
#' sample in \code{methylationData}. Two different values are allowed
#' (for the two conditions).
#' @param regions a \code{\link{GRanges}} object with the regions where to
#' compute the DMRs. If \code{NULL}, the DMRs are computed genome-wide.
#' @param context the context in which the DMRs are computed (\code{"CG"},
#' \code{"CHG"} or \code{"CHH"}).
#' @param method the method used to compute the DMRs \code{"neighbourhood"}
#' or \code{"bins"}). The \code{"neighbourhood"} method computates
#' differentially methylated cytosines. Finally, the \code{"bins"}
#' method partiones the genome into equal sized tilling bins and performs the
#' statistical test between the two conditions in each bin. For all three
#' methods, the cytosines or bins are then merged into DMRs without affecting
#' the inital parameters used when calling the differentiall methylated
#' cytosines/bins (p-value, difference in methylation levels, minimum number of
#' reads per cytosine).
#' @param binSize the size of the tiling bins in nucleotides. This parameter is
#' required only if the selected method is \code{"bins"}.
#' @param test the statistical test used to call DMRs (\code{"betareg"} for
#' Beta regression).
#' @param pseudocountM numerical value to be added to the methylated reads
#' before modelling beta regression.
#' @param pseudocountN numerical value to be added to the total reads
#' before modelling beta regression.
#' @param pValueThreshold DMRs with p-values (when performing the statistical
#' test; see \code{test}) higher or equal than \code{pValueThreshold} are
#' discarded. Note that we adjust the p-values using the Benjamini and
#' Hochberg's method to control the false discovery rate.
#' @param minCytosinesCount DMRs with less cytosines in the specified context
#' than \code{minCytosinesCount} will be discarded.
#' @param minProportionDifference DMRs where the difference in methylation
#' proportion between the two conditions is lower than
#' \code{minProportionDifference} are discarded.
#' @param minGap DMRs separated by a gap of at least \code{minGap} are not
#' merged. Note that only DMRs where the change in methylation is in the same
#' direction are joined.
#' @param minSize DMRs with a size smaller than \code{minSize} are discarded.
#' @param minReadsPerCytosine  DMRs with the average number of reads lower than
#' \code{minReadsPerCytosine} are discarded.
#' @param cores the number of cores used to compute the DMRs.
#' @return the DMRs stored as a \code{\link{GRanges}} object with the following
#' metadata columns:
#' \describe{
#'  \item{direction}{a number indicating whether the region lost (-1)  or gain
#'  (+1) methylation in condition 2 compared to condition 1.}
#'  \item{context}{the context in which the DMRs was computed (\code{"CG"},
#'  \code{"CHG"} or \code{"CHH"}).}
#'  \item{sumReadsM1}{the number of methylated reads in condition 1.}
#'  \item{sumReadsN1}{the total number of reads in condition 1.}
#'  \item{proportion1}{the proportion methylated reads in condition 1.}
#'  \item{sumReadsM2}{the number of methylated reads in condition 2.}
#'  \item{sumReadsN2}{the total number reads in condition 2.}
#'  \item{proportion2}{the proportion methylated reads in condition 2.}
#'  \item{direction2}{a number indicating whether the region lost (-1)  or gain
#'  (+1) methylation in condition 2.2 compared to condition 1.2.}
#'  \item{sumReadsM1.2}{the number of methylated reads in condition 1.2.}
#'  \item{sumReadsN1.2}{the total number of reads in condition 1.2.}
#'  \item{proportion1.2}{the proportion methylated reads in condition 1.2.}
#'  \item{sumReadsM2.2}{the number of methylated reads in condition 2.2.}
#'  \item{sumReadsN2.2}{the total number reads in condition 2.2.}
#'  \item{proportion2.2}{the proportion methylated reads in condition 2.2}
#'  \item{cytosinesCount}{the number of cytosines in the DMR.}
#'  \item{regionType}{a string indicating whether the region lost (\code{"loss"})
#'  or gained (\code{"gain"}) methylation in condition 2 compared to condition 1.}
#'  \item{pValue}{the p-value (adjusted to control the false discovery rate with
#'  the Benjamini and Hochberg's method) of the statistical test when the DMR was
#'  called.}
#' }
#'
#' @examples
#'
#' \dontrun{
#' # starting with data joined using joinReplicates
#' data("syntheticDataReplicates")
#'
#' # compute the DMRs in CG context with neighbourhood method
#'
#' # creating condition vector
#' condition1 <- c("a", "a", "b", "b")
#' condition2 <- c("a", "b", "a", "b")
#' # computing DMRs using the neighbourhood method
#' DMRsReplicatesNeighbourhood <- computeDMRsReplicates(methylationData = methylationData,
#'                                                      condition1 = condition1,
#'                                                      condition2 = condition2,
#'                                                      regions = NULL,
#'                                                      context = "CHH",
#'                                                      method = "neighbourhood",
#'                                                      test = "betareg",
#'                                                      pseudocountM = 1,
#'                                                      pseudocountN = 2,
#'                                                      pValueThreshold = 0.01,
#'                                                      minCytosinesCount = 4,
#'                                                      minProportionDifference = 0.4,
#'                                                      minGap = 200,
#'                                                      minSize = 50,
#'                                                      minReadsPerCytosine = 4,
#'                                                      cores = 1)
#' }
#' @author Alessandro Pio Greco and Nicolae Radu Zabet
#' @export
computeDMRsReplicates <- function(methylationData,
                                  condition1 = NULL,
                                  condition2 = NULL,
                                  regions = NULL,
                                  context = "CG",
                                  method= "neighbourhood",
                                  binSize = 100,
                                  test = "betareg",
                                  pseudocountM = 1,
                                  pseudocountN = 2,
                                  pValueThreshold = 0.01,
                                  minCytosinesCount = 4,
                                  minProportionDifference = 0.4,
                                  minGap = 200,
                                  minSize = 50,
                                  minReadsPerCytosine = 4,
                                  cores = 1) {
  #Parameters checking
  cat("Parameters checking ...\n")

  .validateMethylationData(methylationData, variableName="methylationData")

  .validateConditionReplicates(condition1, condition2, methylationData)

  regions <- .validateGRanges(regions, methylationData)

  .validateContext(context)

  .stopIfNotAll(c(!is.null(method),
                  all(is.character(method)),
                  length(method) == 1,
                  all(method %in% c("neighbourhood","bins"))),
                " method can be only neighbourhood or bins")


  if(method == "bins"){
    .stopIfNotAll(c(.isInteger(binSize, positive=TRUE)),
                  " the bin size used by the method is an integer higher than 0")

  }

  .validateStatisticalTestReplicates(test)

# same pValueThreshold for both conditions FB
  .stopIfNotAll(c(!is.null(pValueThreshold),
                  is.numeric(pValueThreshold),
                  pValueThreshold > 0,
                  pValueThreshold < 1),
                " the p-value threshold needs to be in the interval (0,1)")

  .stopIfNotAll(c(.isInteger(minCytosinesCount, positive=TRUE)),
                " the minCytosinesCount is an integer higher or equal to 0")

#  same minProportionDifference 
  .stopIfNotAll(c(!is.null(minProportionDifference),
                  is.numeric(minProportionDifference),
                  minProportionDifference > 0,
                  minProportionDifference < 1),
                " the minimum difference in methylation needs to be in the interval (0,1)")

  .stopIfNotAll(c(.isInteger(minGap, positive=TRUE)),
                " the minimum gap between DMRs is an integer higher or equal to 0")

  .stopIfNotAll(c(.isInteger(minSize, positive=TRUE)),
                " the minimum size of a DMR is an integer higher or equal to 0")

  .stopIfNotAll(c(.isInteger(minReadsPerCytosine, positive=TRUE)),
                " the minimum average number of reads in a DMR is an integer higher or equal to 0")

  .stopIfNotAll(c(.isInteger(cores, positive=TRUE)),
                " the number of cores to used when computing the DMRs needs to be an integer  hirger or equal to 1.")

  computedDMRs <- GRanges()

  if(method == "neighbourhood"){
    computedDMRs <- .computeDMRsReplicatesNeighbourhood(methylationData = methylationData,
                                                        condition1 = condition1,
                                                        condition2 = condition2,
                                                        regions = regions,
                                                        context = context,
                                                        pseudocountM = pseudocountM,
                                                        pseudocountN = pseudocountN,
                                                        pValueThreshold = pValueThreshold,
                                                        minCytosinesCount = minCytosinesCount,
                                                        minProportionDifference =minProportionDifference,
                                                        minGap = minGap,
                                                        minSize = minSize,
                                                        minReadsPerCytosine = minReadsPerCytosine,
                                                        cores = cores)
  } else if(method == "bins"){
    computedDMRs <- .computeDMRsReplicatesBins(methylationData = methylationData,
                                               condition1 = condition1,
                                               condition2 = condition2,
                                               regions = regions,
                                               context = context,
                                               binSize = binSize,
                                               pseudocountM = pseudocountM,
                                               pseudocountN = pseudocountN,
                                               pValueThreshold = pValueThreshold,
                                               minCytosinesCount = minCytosinesCount,
                                               minProportionDifference =minProportionDifference,
                                               minGap = minGap,
                                               minSize = minSize,
                                               minReadsPerCytosine = minReadsPerCytosine,
                                               cores = cores)

  } else{
    cat("Unknown method: ",method," \n")
  }

  return(computedDMRs)
}

#' This function computes the differentially methylated regions between replicates
#' using the neighbourhood method.
.computeDMRsReplicatesNeighbourhood <- function(methylationData,
                                                condition1,
                                                condition2,
                                                regions = NULL,
                                                context = "CG",
                                                pseudocountM = 1,
                                                pseudocountN = 2,
                                                pValueThreshold = 0.01,
                                                minCytosinesCount = 4,
                                                minProportionDifference = 0.4,
                                                minGap = 200,
                                                minSize = 50,
                                                minReadsPerCytosine = 4,
                                                cores = 1){
  condition1 <- as.factor(condition1)
  condition2 <- as.factor(condition2)
  regions <- reduce(regions)
  regionsList <- .splitGRangesEqualy(regions, cores)

  # extract the methylation data in the correct context
  cat("Extract methylation in the corresponding context \n")
  localContextMethylationData <- methylationData[methylationData$context%in%context]
  localContextMethylationData <- localContextMethylationData[queryHits(findOverlaps(
    localContextMethylationData, regions))]


  # create dataframe row wise with proportions
  m <- grep("readsM", names(mcols(localContextMethylationData)))
  n <- grep("readsN", names(mcols(localContextMethylationData)))
  

  # inner loop function for parallel::mclapply
  .computeDMPsReplicatesNeighbourhoodLoop = function(i){
    computedDMPs <- GRanges()
    for(index in 1:length(regionsList[[i]])){
      cat("Computing DMRs \n")
      currentRegion <- regionsList[[i]][index]

      overlapsCs <- findOverlaps(localContextMethylationData, currentRegion)

      if(length(overlapsCs) > 0){
        localMethylationData <- localContextMethylationData[queryHits(overlapsCs)]
        
        

        proportions <-proportions <- (as.matrix(mcols(localMethylationData)[,m]) + pseudocountM) /
          (as.matrix(mcols(localMethylationData)[,n]) + pseudocountN)

        DMPs <- GRanges()
        if(length(localMethylationData) > 0){
          DMPs <- localMethylationData
          DMPs$pValue1 <- .computeAdjuestedPValuesReplicates(proportions, condition1, condition2, cores=1)[1]
          DMPs$pValue2 <- .computeAdjuestedPValuesReplicates(proportions, condition1, condition2, cores=1)[2]
          DMPs$pValue12 <- .computeAdjuestedPValuesReplicates(proportions, condition1, condition2, cores=1)[3]
          DMPs <- DMPs[!is.na(DMPs$pValue)]
          if(length(DMPs) > 0){
            DMPs$sumReadsM1 <- apply(mcols(DMPs)[m[which(condition1 == unique(condition1)[1])]],1,sum)
            DMPs$sumReadsN1 <- apply(mcols(DMPs)[n[which(condition1 == unique(condition1)[1])]],1,sum)
            DMPs$proportion1 <- DMPs$sumReadsM1 / DMPs$sumReadsN1
            DMPs$sumReadsM2 <- apply(mcols(DMPs)[m[which(condition1 == unique(condition1)[2])]],1,sum)
            DMPs$sumReadsN2 <- apply(mcols(DMPs)[n[which(condition1 == unique(condition1)[2])]],1,sum)
            DMPs$proportion2 <- DMPs$sumReadsM2 / DMPs$sumReadsN2
            DMPs$cytosinesCount <- 1
            DMPs$direction <- sign(DMPs$proportion2 - DMPs$proportion1)
            
            DMPs$sumReadsM1.2 <- apply(mcols(DMPs)[m[which(condition2 == unique(condition2)[1])]],1,sum)
            DMPs$sumReadsN1.2 <- apply(mcols(DMPs)[n[which(condition2 == unique(condition2)[1])]],1,sum)
            DMPs$proportion1.2 <- DMPs$sumReadsM1.2 / DMPs$sumReadsN1.2
            DMPs$sumReadsM2.2 <- apply(mcols(DMPs)[m[which(condition2 == unique(condition2)[2])]],1,sum)
            DMPs$sumReadsN2.2 <- apply(mcols(DMPs)[n[which(condition2 == unique(condition2)[2])]],1,sum)
            DMPs$proportion2.2 <- DMPs$sumReadsM2.2 / DMPs$sumReadsN2.2
            DMPs$cytosinesCount <- 1
            DMPs$direction <- sign(DMPs$proportion2.2 - DMPs$proportion1.2)
          }
        }
        if(length(DMPs) > 0){
          bufferIndex <- DMPs$pValue < pValueThreshold &
            abs(DMPs$proportion2 - DMPs$proportion1) >= minProportionDifference &
            DMPs$sumReadsN1 >=minReadsPerCytosine &
            DMPs$sumReadsN2 >=minReadsPerCytosine &
            DMPs$sumReadsN1.2 >=minReadsPerCytosine &
            DMPs$sumReadsN2.2 >=minReadsPerCytosine 
          DMPs <- DMPs[bufferIndex]
          strand(DMPs) <- "*"
        }

        if(is.null(DMPs)){
          DMPs <- GRanges()
        }


        # append current DMRs to the global list of DMRs
        if(length(computedDMPs) == 0){
          computedDMPs <- DMPs
        } else{
          computedDMPs <- c(computedDMPs,DMPs)
        }
      }
    }
    return(computedDMPs)
  }


  # compute the DMRs
  if(cores > 1){
    cat("Compute the DMRs using ", cores, "cores\n")
    computedDMPs <- parallel::mclapply(1:length(regionsList), .computeDMPsReplicatesNeighbourhoodLoop, mc.cores = cores)
  } else {
    computedDMPs <- lapply(1:length(regionsList), .computeDMPsReplicatesNeighbourhoodLoop)
  }

 computedDMRs <- GRanges()
  computedDMRs <- subset(computedDMPs, !sapply(computedDMPs, is.null))
  if(length(computedDMRs) > 0){
    computedDMRs <- unlist(GRangesList(computedDMRs))
    if(length(computedDMRs) > 0){
      cat("Merge adjacent DMRs\n")
      computedDMRs <- computedDMRs[order(computedDMRs)]
      computedDMRs$pValue <- .computeaAjustedPValuesInDMRsReplicates(localContextMethylationData ,computedDMRs, condition, cores, m, n, pseudocountM, pseudocountN)
      cat("Merge DMRs iteratively\n")
      # Get rid of small gaps between DMRs.
      if(minGap > 0){
        computedDMRs <- .smartMergeDMRsReplicates(computedDMRs,
                                                  minGap = minGap,
                                                  respectSigns = TRUE,
                                                  methylationData = localContextMethylationData,
                                                  minProportionDifference=minProportionDifference,
                                                  minReadsPerCytosine = minReadsPerCytosine,
                                                  pseudocountM = pseudocountM,
                                                  pseudocountN = pseudocountN,
                                                  pValueThreshold=pValueThreshold,
                                                  condition1 = condition1,
                                                  condition2 = condition2,
                                                  m = m,
                                                  n = n,
                                                  cores = cores)
      }


      cat("Filter DMRs \n")
      if(length(computedDMRs) > 0){
        #remove small DMRs
        computedDMRs <- computedDMRs[width(computedDMRs) >= minSize]
        if(length(computedDMRs) > 0){
          #remove DMRs with few cytosines
          computedDMRs <- computedDMRs[computedDMRs$cytosinesCount >= minCytosinesCount]
          #recompute the adjusted p-values
          if(length(computedDMRs) > 0){
            computedDMRs$pValue <- .computeaAjustedPValuesInDMRsReplicates(localContextMethylationData ,computedDMRs, condition, cores, m, n, pseudocountM, pseudocountN)
            computedDMRs$regionType <- rep("loss", length(computedDMRs))
            computedDMRs$regionType[which(computedDMRs$proportion1 < computedDMRs$proportion2)] <- "gain"

          }
        }
      }
    }
  }

  return(computedDMRs)
}

#' This function computes the differentially methylated regions between replicates
#' using the bins method.
.computeDMRsReplicatesBins <- function(methylationData,
                                       condition1 = condition1,
                                       condition2 = condition2,
                                       regions = NULL,
                                       context = "CG",
                                       binSize = 100,
                                       pseudocountM = 1,
                                       pseudocountN = 2,
                                       pValueThreshold = 0.01,
                                       minCytosinesCount = 4,
                                       minProportionDifference = 0.4,
                                       minGap = 200,
                                       minSize = 50,
                                       minReadsPerCytosine = 4,
                                       cores = 1) {

  condition1 <- as.factor(condition1)
  condition2 <- as.factor(condition2)

  regions <- reduce(regions)

  # extract the methylation data in the correct context
  cat("Extract methylation in the corresponding context \n")
  localContextMethylationData <- methylationData[methylationData$context%in%context]
  localContextMethylationData <- localContextMethylationData[queryHits(findOverlaps(
    localContextMethylationData, regions))]

  regionsList <- .splitGRangesEqualy(regions, cores)


  m <- grep("readsM", names(mcols(localContextMethylationData))) # == Number of samples
  n <- grep("readsN", names(mcols(localContextMethylationData))) # == Number of samples

  # inner loop function for parallel::mclapply
  .computeDMRsReplicatesBinsLoop = function(i){
    computedDMRs <- GRanges()
    for(index in 1:length(regionsList[[i]])){
      currentRegion <- regionsList[[i]][index]

      cat("Computing DMRs at ",.printGenomicRanges(currentRegion),"\n")

      seqs <- seq(start(currentRegion), (end(currentRegion)-binSize), by = binSize);

      bins <- GRanges(seqnames(currentRegion), IRanges(seqs, (seqs+binSize-1)))

      overlapsBins <- findOverlaps(localContextMethylationData, currentRegion)

      if(length(overlapsBins) > 0){
        localMethylationData <- localContextMethylationData[queryHits(overlapsBins)]

        cat("Count inside each bin...\n")
        #bins <- .analyseReadsInsideRegions(localMethylationData, bins, context, cores)
        bins <- .analyseReadsInsideBinsReplicates(localMethylationData, bins, currentRegion, condition1, condition2, pseudocountM, pseudocountN)


        cat("Filter the bins...\n")
        # Get rid of the bins with fewer than minCytosinesCount cytosines inside.
        bins  <- bins[bins$cytosinesCount >= minCytosinesCount]

        # Get rid of the bins with fewer than minReadsPerCytosine reads per cytosine.
        bins  <- bins[(bins$sumReadsN1/bins$cytosinesCount >= minReadsPerCytosine) &
                        (bins$sumReadsN2/bins$cytosinesCount >= minReadsPerCytosine)]

        # Get rid of the bins with small difference in proportion of methylation
        bins  <- bins[(abs(bins$proportion1 - bins$proportion2) >= minProportionDifference)]

        proportions <- as.matrix(mcols(bins)[grep("proportionsR", names(mcols(bins)))])

        if(nrow(proportions) > 0){
          cat("Identifying DMRs...\n")
          pValue<- .computeAdjuestedPValuesReplicates(proportions, condition1, condition2, cores)
          mcols(bins)[grep("proportionsR", names(mcols(bins)))] <- NULL
          bins <- bins[!is.na(pValue) & pValue < pValueThreshold ]
          if(length(bins) > 0){
            bins$context <- rep(paste(context, collapse = "_"), length(bins))
            bins$direction <- rep(NA, length(bins))
            bins$direction <- sign(bins$proportion2 - bins$proportion1)
            # Select the crude list of DMRs
            DMRs <- bins[!is.na(bins$direction) & (bins$direction == 1 | bins$direction == -1)]
          }
        } else{
          DMRs <- bins
        }

        if(is.null(DMRs)){
          DMRs <- GRanges()
        }

        # append current DMRs to the global list of DMRs
        if(length(computedDMRs) == 0){
          computedDMRs <- DMRs
        } else{
          computedDMRs <- c(computedDMRs,DMRs)
        }
      }
    }
    return(computedDMRs)
  }

  # compute the DMRs
  if(cores > 1){
    cat("Compute the DMRs using ", cores, "cores\n")
    computedDMRs <- parallel::mclapply(1:length(regionsList), .computeDMRsReplicatesBinsLoop, mc.cores = cores)
  } else {
    computedDMRs <- lapply(1:length(regionsList), .computeDMRsReplicatesBinsLoop)
  }
  if(length(computedDMRs) > 0){
    computedDMRs <- subset(computedDMRs, !sapply(computedDMRs, is.null))
    if(length(computedDMRs) > 0){
      computedDMRs <- unlist(GRangesList(computedDMRs))
    } else{
      computedDMRs <- GRanges()
    }


    if(length(computedDMRs) > 0){

      cat("Merge adjacent DMRs\n")
      computedDMRs <- computedDMRs[order(computedDMRs)]
      computedDMRs$pValue <- .computeaAjustedPValuesInDMRsReplicates(localContextMethylationData ,computedDMRs, condition1, condition2, cores, m, n, pseudocountM, pseudocountN)
      cat("Merge DMRs iteratively\n")
      # Get rid of small gaps between DMRs.
      if(minGap > 0){
        computedDMRs <- .smartMergeDMRsReplicates(computedDMRs,
                                                  minGap = minGap,
                                                  respectSigns = TRUE,
                                                  methylationData = localContextMethylationData,
                                                  minProportionDifference=minProportionDifference,
                                                  minReadsPerCytosine = minReadsPerCytosine,
                                                  pseudocountM = pseudocountM,
                                                  pseudocountN = pseudocountN,
                                                  pValueThreshold=pValueThreshold,
                                                  condition1 = condition1,
                                                  consition2 =condition2,
                                                  m = m,
                                                  n = n,
                                                  cores = cores)
      }

      cat("Filter DMRs \n")
      if(length(computedDMRs) > 0){
        #remove small DMRs
        computedDMRs <- computedDMRs[width(computedDMRs) >= minSize]

        if(length(computedDMRs) > 0){
          #remove DMRs with few cytosines
          computedDMRs <- computedDMRs[computedDMRs$cytosinesCount >= minCytosinesCount]
          #recompute the adjusted p-values
          if(length(computedDMRs) > 0){
            computedDMRs$pValue <- .computeaAjustedPValuesInDMRsReplicates(localContextMethylationData ,computedDMRs, condition, cores, m, n, pseudocountM, pseudocountN)
            computedDMRs$regionType <- rep("loss", length(computedDMRs))
            computedDMRs$regionType[which(computedDMRs$proportion1 < computedDMRs$proportion2)] <- "gain"
          }
        }
      }
    }
  } else{
    computedDMRs <- GRanges()
  }
  return(computedDMRs)
}

#' This function extracts the p-value for each cytosine from the betareg object for each condition and the interaction
.convertResult <- function(x){
  if(any(is.na(x))){
    result <- c(NA, NA, NA)
  } else{
    library(car)
    result <- Anova(x, type=3)[,3] # Perform a type 3 contrasta and collect p-values from condition 1 , conditio2 and condition1:2
  }
  return(result)
}

# formula for eliminating extreme values (0,1) from the proportion matrix
.convertProportions <- function(x){
  result <- (x *(length(x) - 1) + 0.5)/length(x)
  return(result)
}

#' This function computes the p-values of the beta regression test
#'
#' @title Beta regression
#' @param y vector containing the proportions (readsM/readsN) for a single cytosine.
#' @param condition1 vector containing the condition 1 of the experiment.
#' @param condition2 vector containing the condition 2 of the experiment.
#' @return The p-values of the beta regression test.
#'
#' @author Alessandro Pio Greco and Nicolae Radu Zabet
.computeBetaRegSingle <- function(y, condition1, condition2){
  if(!.suitableForBetaReg(y)){
    result <- NA
  } else{
   # result <- betareg(.convertProportions(y) ~ condition1 + condition2 + condition1:condition2, na.action = na.pass)
    result <- betareg(.convertProportions(y) ~ condition1 + condition2 , na.action = na.pass)
    
     }
  return(result)
}


#https://stats.stackexchange.com/questions/220868/questionable-beta-regression-results
#https://stats.stackexchange.com/questions/89999/how-to-replicate-statas-robust-binomial-glm-for-proportion-data-in-r/205040#205040
.suitableForBetaReg <- function(y){
  if(any(is.na(y))){
    return(FALSE)
  } else if(length(unique(y)) < 4 | any(y==0) | any(y==1)){
    return(FALSE)
  }
  return(TRUE)
}






#' This function applies the .computeBetaRegSingle to the entire dataset
.computeBetaReg <- function(y,  condition1, condition2, cores){
  if(is.null(nrow(y))){
    result <- .computeBetaRegSingle(y,  condition1, condition2)
    result <- .convertResult(result)
  } else {
    suitableForBetaReg <- apply(y, 1, .suitableForBetaReg)
    result <- rep(NA, nrow(y))
    if(sum(suitableForBetaReg) >= 1){
      ySuitable <- matrix(y[suitableForBetaReg,], ncol=length(condition1))
      
      #buffer <- mclapply(1:nrow(ySuitable), FUN = function(x) .computeBetaRegSingle(ySuitable[x,], condition), mc.cores = cores)
      buffer <- apply(ySuitable, 1, .computeBetaRegSingle, condition1=condition1, condition2=condition2)
      buffer <- vapply(buffer, FUN=.convertResult, FUN.VALUE = 0.0)
      #buffer <- unlist(buffer)
      result[suitableForBetaReg] <- buffer
    }
  }
  return(result)
}

#' This function computes the adjusted p-values (using Benjamini & Hochberg
#' method)
#'
#'
#' @title Compute adjusted p-values
#' @param proportions matrix containing the proportions (readsM/readsN) for all the cytosines.
#' @param condition1 vector containing the conditions of the experiment.
#' @param condition2 vector containing the conditions of the experiment.
#' @param cores the number of cores to be used to fit the model.
#' @return The adjusted p-values of the statistical test.
#' @author  Nicolae Radu Zabet and Alessandro Pio Greco
.computeAdjuestedPValuesReplicates <- function(proportions, condition1, condition2, cores){

  pValue <- .computeBetaReg(proportions, condition1, condition2, cores)
  # convert p-values to FDR
  adjPValue1 <- rep(NA, times=length(pValue[,1]))

  adjPValue1[which(!is.na(pValue[,1]))] <- p.adjust(pValue[which(!is.na(pValue[,1])),1], method="fdr")

  # convert p-values to FDR
  adjPValue2 <- rep(NA, times=length(pValue[,2]))
  
  adjPValue2[which(!is.na(pValue[,2]))] <- p.adjust(pValue[which(!is.na(pValue[,2])),2], method="fdr")
  
  adjPValue3 <- rep(NA, times=length(pValue[,3]))
  
  adjPValue2[which(!is.na(pValue[,3]))] <- p.adjust(pValue[which(!is.na(pValue[,3])),3], method="fdr")
  
  adjPValue<- cbind(adjPValue1,adjPValue2,adjPValue3)
  
  return(adjPValue)
}

#' This function computes the adjusted p-values (using Benjamini & Hochberg
#' method)
#'
#'
#' @title Compute adjusted p-values
#' @param proportions matrix containing the proportions (readsM/readsN) for all the cytosines.
#' @param condition1 vector containing the conditions 1 of the experiment.
#' @param condition2 vector containing the conditions 2 of the experiment.
#' @param cores the number of cores to be used to fit the model.
#' @return The adjusted p-values of the statistical test.
#' @author  Nicolae Radu Zabet and Alessandro Pio Greco
.computeaAjustedPValuesInDMRsReplicates <- function(methylationData,
                                                    DMRs,
                                                    condition1,
                                                    condition2,
                                                    cores,
                                                    indexM,
                                                    indexN,
                                                    pseudocountM,
                                                    pseudocountN){

  M <- sapply(1:length(DMRs),
              function(x){.computeProportionsInDMRs(DMRs[x],
                                                    methylationData = methylationData,
                                                    col_indexes = indexM)})
  N <- sapply(1:length(DMRs),
              function(x){.computeProportionsInDMRs(DMRs[x],
                                                    methylationData = methylationData,
                                                    col_indexes = indexN)})
  

  proportions <- (M + pseudocountM) / (N + pseudocountN)
  proportions <- t(proportions)



  pValue <- .computeBetaReg(proportions, condition1, condition2, cores)

  # convert p-values to FDR
  adjPValue1 <- rep(NA, times=length(pValue[,1]))
  
  adjPValue1[which(!is.na(pValue[,1]))] <- p.adjust(pValue[which(!is.na(pValue[,1])),1], method="fdr")
  
  # convert p-values to FDR
  adjPValue2 <- rep(NA, times=length(pValue[,2]))
  
  adjPValue2[which(!is.na(pValue[,2]))] <- p.adjust(pValue[which(!is.na(pValue[,2])),2], method="fdr")
  
  adjPValue3 <- rep(NA, times=length(pValue[,3]))
  
  adjPValue2[which(!is.na(pValue[,3]))] <- p.adjust(pValue[which(!is.na(pValue[,3])),3], method="fdr")
  
  adjPValue<- cbind(adjPValue1,adjPValue2,adjPValue3)
  
  return(adjPValue)
}

#' This function recalculate proportions between methylated and total reads in a new region.
.computeProportionsInDMRs <- function(DMRs, methylationData,col_indexes){
  regions <- methylationData[queryHits(findOverlaps(methylationData, DMRs))]
  computedDMRsproportions <- apply(mcols(regions)[col_indexes], 2, sum)
  return(computedDMRsproportions)
}

#' Performs the analysis in equal width regions of an \code{\link{GRanges}}
#' object
#'
#' @title Analyse reads inside regions
#' @param methylationData a \code{\link{GRanges}} object with five metadata
#' columns; see \code{\link{methylationDataList}}
#' @param currentRegion a \code{\link{GRanges}} object with the identified regions
#' @param condition The vector containing the two conditions for the experiment.
#' @param pseudocountM numerical value to be added to the methylated reads
#' before modelling beta regression.
#' @param pseudocountN numerical value to be added to the total reads
#' before modelling beta regression.
#' @return a \code{\link{GRanges}} object with eaual sized tiles of the regions.
#' The object consists of the following metadata
#' \describe{
#'  \item{sumReadsM1}{the number of methylated reads in condition 1}
#'  \item{sumReadsN1}{the total number of reads in condition 1}
#'  \item{proportion1}{the proportion of methylated reads in condition 1}
#'  \item{sumReadsM2}{the number of methylated reads in condition 2}
#'  \item{sumReadsN2}{the total number of reads in condition 2}
#'  \item{proportion2}{the proportion of methylated reads in condition 2}
#'  \item{cytosinesCount}{the number of cytosines in the correct context}
#' }
#'
#' @author Nicolae Radu Zabet and Alessandro Pio Greco
.analyseReadsInsideBinsReplicates <- function(methylationData, bins, currentRegion,
                                              condition1, condition2, pseudocountM, pseudocountN){


  binSize <- min(unique(width(bins)))
  #Rcpp
  m <- grep("readsM", names(mcols(methylationData)))
  n <- grep("readsN", names(mcols(methylationData)))
  readsM <- matrix(0, ncol = length(m), nrow=length(bins))
  for(i in 1:length(m)){
    test <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), mcols(methylationData)[[m[i]]], windowSize = binSize)
    readsM[,i] <- test[seq(1,length(test)-binSize, by=binSize)]
  }

  readsN <- matrix(0, ncol = length(n), nrow=length(bins))
  for(i in 1:length(n)){
    test2 <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), mcols(methylationData)[[n[i]]], windowSize = binSize)
    readsN[,i] <- test2[seq(1,length(test)-binSize, by=binSize)]
  }

  proportions <- (readsM + pseudocountM)/ (readsN + pseudocountN)
  proportions <- as.data.frame(proportions)
  names_prop <- paste0("proportionsR", 1:ncol(proportions))
  colnames(proportions) <- names_prop

  m1 <- m[which(condition1 == unique(condition1)[1])]
  n1 <- n[which(condition1 == unique(condition1)[2])]
  m2 <- m[which(condition2 == unique(condition2)[1])]
  n2 <- n[which(condition2 == unique(condition2)[2])]

  readsM1 <- readsM[,which(condition1 == unique(condition1)[1])]
  # readsM1 <- matrix(0, ncol = length(m1), nrow=length(bins))
  # for(i in 1:length(m1)){
  #   test <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), mcols(methylationData)[[m1[i]]], windowSize = binSize)
  #   readsM1[,i] <- test[seq(1,length(test)-binSize, by=binSize)]
  # }
  sumReadsM1 <-  apply(readsM1,1,sum)



  readsN1 <- readsN[,which(condition1 == unique(condition1)[1])]
  # readsN1 <- matrix(0, ncol = length(n1), nrow=length(bins))
  # for(i in 1:length(n1)){
  #   test <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), mcols(methylationData)[[n1[i]]], windowSize = binSize)
  #   readsN1[,i] <- test[seq(1,length(test)-binSize, by=binSize)]
  # }
  sumReadsN1 <- apply(readsN1,1,sum)
  proportion1 <- (sumReadsM1 + pseudocountM)/ (sumReadsN1 + pseudocountN)


  readsM2 <- readsM[,which(condition2 == unique(condition2)[1])]
  # readsM2 <- matrix(0, ncol = length(m2), nrow=length(bins))
  # for(i in 1:length(m2)){
  #   test <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), mcols(methylationData)[[m2[i]]], windowSize = binSize)
  #   readsM2[,i] <- test[seq(1,length(test)-binSize, by=binSize)]
  # }
  sumReadsM2 <- apply(readsM2,1,sum)

  readsN2 <- readsN[,which(condition2 == unique(condition2)[2])]
  # readsN2 <- matrix(0, ncol = length(n2), nrow=length(bins))
  # for(i in 1:length(n2)){
  #   test <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), mcols(methylationData)[[n2[i]]], windowSize = binSize)
  #   readsN2[,i] <- test[seq(1,length(test)-binSize, by=binSize)]
  # }
  sumReadsN2 <- apply(readsN2,1,sum)

  #
  proportion2 <- (sumReadsM2 + pseudocountM)/ (sumReadsN2 + pseudocountN)

  cytosines <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), rep(1, length(start(methylationData))), windowSize = binSize)
  cytosinesCount <- cytosines[seq(1,length(cytosines)-binSize, by=binSize)]


  bins$sumReadsM1 <- sumReadsM1
  bins$sumReadsN1 <- sumReadsN1
  bins$proportion1 <- proportion1
  bins$sumReadsM2 <- sumReadsM2
  bins$sumReadsN2 <- sumReadsN2
  bins$proportion2 <- proportion2
  bins$cytosinesCount <- cytosinesCount
  mcols(bins) <- cbind(mcols(bins), proportions)
  return(bins)

}

#' Performs the analysis in all regions in a \code{\link{GRanges}} object
#'
#' @title Analyse reads inside regions
#' @param methylationData a \code{\link{GRanges}} object with five metadata
#' columns see \code{\link{methylationDataList}}
#' @param regions a \code{\link{GRanges}} object with the identified regions
#' @param condition The vector containing the two conditions for the experiment.
#' @param m indexes of methylated reads for creation of the proportions matrix
#' @param n indexes of total reads for creation of the proportions matrix
#' @return a \code{\link{GRanges}} object with eaual sized tiles of the regions.
#' The object consists of the following metadata
#' \describe{
#'  \item{sumReadsM1}{the number of methylated reads in condition 1}
#'  \item{sumReadsN1}{the total number of reads in condition 1}
#'  \item{proportion1}{the proportion of methylated reads in condition 1}
#'  \item{sumReadsM2}{the number of methylated reads in condition 2}
#'  \item{sumReadsN2}{the total number of reads in condition 2}
#'  \item{proportion2}{the proportion of methylated reads in condition 2}
#'  \item{cytosinesCount}{the number of cytosines in the correct context}
#' }
#'
#' @author Nicolae Radu Zabet and Alessandro Pio Greco
.analyseReadsInsideRegionsReplicates <- function(methylationData, regions, condition, m, n){
  overlaps <- findOverlaps(methylationData, regions, ignore.strand = TRUE)
  methylationDataContextList <- IRanges::splitAsList(methylationData[queryHits(overlaps)],  subjectHits(overlaps))
  regionsIndexes <- as.integer(names(methylationDataContextList))
  # methylationDataContextList <-
  # regions <- methylationData[queryHits(findOverlaps(methylationData, DMRs))]

  regions$sumReadsM1 <- rep(0, times=length(regions))
  regions$sumReadsN1 <- rep(0, times=length(regions))
  regions$proportion1 <- rep(0, times=length(regions))
  regions$sumReadsM2 <- rep(0, times=length(regions))
  regions$sumReadsN2 <- rep(0, times=length(regions))
  regions$proportion2 <- rep(0, times=length(regions))
  regions$cytosinesCount <- rep(0, times=length(regions))

  if(length(regionsIndexes) > 0){
    regions$sumReadsM1[regionsIndexes] <- sum(apply(mcols(unlist(methylationDataContextList))[m[which(condition == unique(condition)[1])]],1,sum))
    regions$sumReadsM2 <- sum(apply(mcols(unlist(methylationDataContextList))[m[which(condition == unique(condition)[2])]],1,sum))
    regions$sumReadsN1 <- sum(apply(mcols(unlist(methylationDataContextList))[n[which(condition == unique(condition)[1])]],1,sum))
    regions$sumReadsN2 <- sum(apply(mcols(unlist(methylationDataContextList))[n[which(condition == unique(condition)[2])]],1,sum))
    regions$cytosinesCount[regionsIndexes] <- sapply(methylationDataContextList,length)

    valid <- regions$cytosinesCount[regionsIndexes] > 0
    regions$proportion1[regionsIndexes[valid]] <- regions$sumReadsM1[regionsIndexes[valid]]/regions$sumReadsN1[regionsIndexes[valid]]
    regions$proportion2[regionsIndexes[valid]] <- regions$sumReadsM2[regionsIndexes[valid]]/regions$sumReadsN2[regionsIndexes[valid]]
  }
  return(regions)
}



.joinDMRsReplicates <- function(DMRs,
                                minGap = minGap,
                                respectSigns = TRUE,
                                methylationData = methylationData,
                                minProportionDifference=0.4,
                                minReadsPerCytosine = 4,
                                pseudocountM = pseudocountM,
                                pseudocountN = pseudocountN,
                                pValueThreshold=0.01,
                                condition = condition,
                                m = m,
                                n = n,
                                cores = 1
){

  #are within the requeired distance
  if(length(reduce(DMRs, drop.empty.ranges=TRUE, min.gapwidth=minGap,
                   ignore.strand=TRUE)) == 1){
    #all have the same direction
    if(!respectSigns | length(unique(DMRs$direction)) == 1){
      direction <- unique(DMRs$direction)
      if(length(unique(strand(DMRs))) == 1){
        localDMR <- DMRs[1]
        end(localDMR) <- max(end(DMRs))
        start(localDMR) <- min(start(DMRs))
        localDMR <- .analyseReadsInsideRegionsReplicates(methylationData, localDMR, condition, m, n)
        localDMR$pValue <- .computeaAjustedPValuesInDMRsReplicates(methylationData, localDMR, condition, cores, m, n, pseudocountM, pseudocountN)
        if(!is.na(localDMR$pValue)){
          if(abs(localDMR$proportion1 - localDMR$proportion2) >= minProportionDifference &
             localDMR$pValue <= pValueThreshold &
             (localDMR$sumReadsN1 / localDMR$cytosinesCount) >= minReadsPerCytosine &
             (localDMR$sumReadsN2 / localDMR$cytosinesCount) >= minReadsPerCytosine){
            DMRs <- localDMR
          }
        }
      }
    }
  }
  return(DMRs)
}

.getLongestDMRsReplicates <- function(DMRs,
                                      minGap = minGap,
                                      respectSigns = TRUE,
                                      methylationData = methylationData,
                                      minProportionDifference=0.4,
                                      minReadsPerCytosine = 4,
                                      pseudocountM = pseudocountM,
                                      pseudocountN = pseudocountN,
                                      pValueThreshold=0.01,
                                      condition = condition,
                                      m = m,
                                      n = n,
                                      cores = 1){

  newDMRs <-.joinDMRsReplicates(DMRs,
                                minGap = minGap,
                                respectSigns = respectSigns,
                                methylationData = methylationData,
                                minProportionDifference=minProportionDifference,
                                minReadsPerCytosine = minReadsPerCytosine,
                                pseudocountM = pseudocountM,
                                pseudocountN = pseudocountN,
                                pValueThreshold=pValueThreshold,
                                condition = condition,
                                m = m,
                                n = n,
                                cores = 1)
  if(length(newDMRs) == 1){
    result <- newDMRs
  } else{
    result <- .mergeDMRsIterativelyReplicates(DMRs,
                                              minGap = minGap,
                                              respectSigns = respectSigns,
                                              methylationData = methylationData,
                                              minProportionDifference=minProportionDifference,
                                              minReadsPerCytosine = minReadsPerCytosine,
                                              pseudocountM = pseudocountM,
                                              pseudocountN = pseudocountN,
                                              pValueThreshold=pValueThreshold,
                                              condition = condition,
                                              m = m,
                                              n = n,
                                              cores = 1)

  }
  return(result)
}

#' This takes a list of DMRs and attempts to merge DMRs while keeping the new
#' DMRs statistically significant
.smartMergeDMRsReplicates <- function(DMRs,
                                      minGap,
                                      respectSigns = TRUE,
                                      methylationData,
                                      minProportionDifference=0.4,
                                      minReadsPerCytosine = 4,
                                      pseudocountM = 1,
                                      pseudocountN = 2,
                                      pValueThreshold=0.01,
                                      condition = condition,
                                      m = m,
                                      n = n,
                                      cores = 1){



  overlaps <- countOverlaps(DMRs, DMRs, maxgap = minGap, ignore.strand = TRUE)
  notToJoin <- DMRs[overlaps == 1]

  DMRsToJoin <- DMRs[overlaps > 1]


  if(length(DMRsToJoin) > 0){
    overlaps <- findOverlaps(DMRsToJoin,
                             reduce(DMRsToJoin, min.gapwidth = minGap,
                                    ignore.strand=TRUE),
                             maxgap = minGap, ignore.strand = TRUE)
    DMRsList <- IRanges::splitAsList(DMRsToJoin[queryHits(overlaps)],
                                     subjectHits(overlaps))



    if(cores > 1){
      bufferDMRs <- parallel::mclapply(1:length(DMRsList), function(i){ .getLongestDMRsReplicates(DMRsList[[i]],
                                                                                                  minGap = minGap,
                                                                                                  respectSigns = respectSigns,
                                                                                                  methylationData = methylationData,
                                                                                                  minProportionDifference=minProportionDifference,
                                                                                                  minReadsPerCytosine = minReadsPerCytosine,
                                                                                                  pseudocountM = pseudocountM,
                                                                                                  pseudocountN = pseudocountN,
                                                                                                  pValueThreshold=pValueThreshold,
                                                                                                  condition = condition,
                                                                                                  m = m,
                                                                                                  n = n,
                                                                                                  cores = 1)},
                                       mc.cores = cores)

      bufferDMRs <- unlist(GRangesList(bufferDMRs))
    } else{
      bufferDMRs <- GRanges()
      for(i in 1:length(DMRsList)){
        bufferDMRs <- c(bufferDMRs, .getLongestDMRsReplicates(DMRsList[[i]],
                                                              minGap = minGap,
                                                              respectSigns = respectSigns,
                                                              methylationData = methylationData,
                                                              minProportionDifference=minProportionDifference,
                                                              minReadsPerCytosine = minReadsPerCytosine,
                                                              pseudocountM = pseudocountM,
                                                              pseudocountN = pseudocountN,
                                                              pValueThreshold=pValueThreshold,
                                                              condition = condition,
                                                              m = m,
                                                              n = n,
                                                              cores = 1))
      }
    }

    joinedDMRs <- c(bufferDMRs, notToJoin)
  } else{
    joinedDMRs <- notToJoin
  }

  joinedDMRs <- joinedDMRs[order(joinedDMRs)]

  return(joinedDMRs)
}

.mergeDMRsIterativelyReplicates <- function(DMRs,
                                            minGap,
                                            respectSigns = TRUE,
                                            methylationData,
                                            minProportionDifference=0.4,
                                            minReadsPerCytosine = 4,
                                            pseudocountM = pseudocountM,
                                            pseudocountN = pseudocountN,
                                            pValueThreshold=0.01,
                                            condition = condition,
                                            m = m,
                                            n = n,
                                            cores = 1){

  overlaps <- countOverlaps(DMRs, DMRs, maxgap = minGap, ignore.strand = TRUE)
  notToJoin <- DMRs[overlaps == 1]
  DMRs <- DMRs[overlaps > 1]

  joinedAny <- TRUE
  iteration <- 1
  while(joinedAny){
    joinedAny <- FALSE
    bufferDMRs <-GRanges()
    index <- 1
    localIndex <- index
    joinedCount <- 0
    while(localIndex < length(DMRs)){
      localDMRs <- DMRs[index]
      canJoin <- TRUE
      while(canJoin){
        localIndex <- localIndex + 1
        newDMRs <-.joinDMRsReplicates(c(localDMRs, DMRs[localIndex]),
                                      minGap = minGap,
                                      respectSigns = respectSigns,
                                      methylationData = methylationData,
                                      minProportionDifference=minProportionDifference,
                                      minReadsPerCytosine = minReadsPerCytosine,
                                      pseudocountM = pseudocountM,
                                      pseudocountN = pseudocountN,
                                      pValueThreshold=pValueThreshold,
                                      condition = condition,
                                      m = m,
                                      n = n,
                                      cores = 1)

        if(length(newDMRs) == 1){
          joinedCount <- joinedCount + 1
          localDMRs <- newDMRs
          joinedAny <- TRUE
          if(localIndex == length(DMRs)){
            canJoin <- FALSE
            bufferDMRs <- c(bufferDMRs, newDMRs)
          }
        } else{
          bufferDMRs <- c(bufferDMRs, localDMRs)
          canJoin <- FALSE
          index <- localIndex
          if(localIndex == length(DMRs)){
            bufferDMRs <- c(bufferDMRs, DMRs[localIndex])
          }
        }
      }
    }
    DMRs <- bufferDMRs
    iteration <- iteration + 1
  }

  DMRs <- c(DMRs, notToJoin)

  DMRs <- DMRs[order(DMRs)]

  return(DMRs)
}
