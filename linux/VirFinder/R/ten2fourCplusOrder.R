ten2fourCplusOrder <- function(ten, WWW, ZI)
{
  getwd<-ten-1
	code<-rep(0,WWW)
	for(ijk in 1:WWW)
	{
		code[WWW-ijk+1] <- getwd%%ZI+1
		getwd <- getwd%/%ZI
	}
	code
}
