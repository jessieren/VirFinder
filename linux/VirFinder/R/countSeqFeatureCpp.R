countSeqFeatureCpp <-
function(RseqDNA, k) {
    .Call('VirFinder_countSeqFeatureCpp', PACKAGE = 'VirFinder', RseqDNA, k)
}
