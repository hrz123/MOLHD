cpf2 <- function(newdes, newpfval, curdes, curpfval)
{
	pfdes=curdes
	pfval=curpfval
	temppfval=newpfval
	temppfdes=newdes
	ncdes=dim(temppfdes)[2]/dim(temppfval)[1]
	for(t in 1:dim(temppfval)[1])
    {
        temp=checkon2(matrix(temppfval[t,],nrow=1),temppfdes[,(t-1)*ncdes+1:ncdes],pfval,pfdes)
        pfval=temp[[1]]
        pfdes=temp[[2]]
    }
	return(list("pfdes"=pfdes,"pfvals"=pfval))
}


cpf3 <- function(newdes, newpfval, curdes, curpfval)
{
	pfdes=curdes
	pfval=curpfval
	temppfval=newpfval
	temppfdes=newdes
	ncdes=dim(temppfdes)[2]/dim(temppfval)[1]
	for(t in 1:dim(temppfval)[1])
    {
        temp=checkon3(matrix(temppfval[t,],nrow=1),temppfdes[,(t-1)*ncdes+1:ncdes],pfval,pfdes)
        pfval=temp[[1]]
        pfdes=temp[[2]]
    }
	return(list("pfdes"=pfdes,"pfvals"=pfval))
}