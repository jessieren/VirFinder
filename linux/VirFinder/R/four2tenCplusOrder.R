four2tenCplusOrder <- function(four, ZI)
{
  four <- rev(four)
  WWW <- length(four)
	
	index<-four[1]
	if(WWW >1)
	{
		for(ijk in 1:(WWW-1))
		{
			index<-index+(four[1+ijk]-1)*ZI^(ijk)
		}
	}
	index
}
