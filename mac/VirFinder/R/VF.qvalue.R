VF.qvalue <-
function(pvalue, fdr.level = NULL, pfdr = FALSE, ...)
{
qvalue(pvalue)$qvalue
}
