######
# normalized.base.measures: get delta and sigma using quantile normalied
# probe intensities for each base with bad amplicons removed
#####

normalized.base.measures <- function(chipType,dataDir,
                                     chipOrder,
                                     ampmap,badamp,exonList,
                                     chip.normalize=FALSE,
                                     amp.normalize=TRUE){

    celint.qn <- get.intensity(chipType,dataDir,chipOrder,
                        ampmap,badamp,exonList,
                        chip.normalize,
                        amp.normalize)

    delta.sigma.qn <- get.delta.sigma(chipType=chipType,
                                         dataDir=dataDir,
                                         chipOrder=chipOrder,
                                         celint=celint.qn,
                                         exonList=exonList)

    return(delta.sigma.qn)
}


