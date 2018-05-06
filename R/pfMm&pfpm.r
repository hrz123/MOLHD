pfMm <- function(n, p, 
				 crlim,
				 num,
				 nstarts=1,
				 times = 300,
				 maxiter = 1e+06, 
				 temp0=0,
				 wtset=cbind(c(1,0),c(0.8,0.2),c(0.6,0.4),c(0.4,0.6),c(0.2,0.8),c(0,1)))
{
	time0=Sys.time()
	time00=Sys.time()
	
	temp=alg4(n, p, wtset, crlim, num, temp0, times, maxiter)
	paretovals=matrix(temp[[1]],ncol=2)
	paretodes=temp[[2]]

	for(icount in 2:nstarts)
	{
		temp=alg4(n, p, wtset, crlim, num, temp0, times, maxiter)
		temp[[1]]=matrix(temp[[1]],ncol=2)
	
		for(i in 1:nrow(temp[[1]]))
		{
			newpt = matrix(temp[[1]][i,],nrow=1)
			newdes = temp[[2]][,p*(i-1)+1:p]
			temppf=checkon2(newpt,newdes,paretovals,paretodes)
			paretovals=temppf[[1]]
			paretodes=temppf[[2]]
		}
		if(icount%%100==0)
		{
			time1=sprintf('100 simulations done in %5.3f Minutes\n  ***====================================================***\n',difftime(Sys.time(),time0,units='mins'))
			cat(time1)
			time0=Sys.time()
		}
	}
	rtime=difftime(Sys.time(),time00,units='mins')
	return(list("pfdes"=paretodes, "pfvals"=paretovals, "time_rec"=rtime))
}

alg4 <- function(n, p, wtset, crlim, num, temp0=0,times = 300, maxiter=1e+06)
{
	if (temp0 == 0) 
	{
        avgdist1 = 1/(n - 1)
        avgdist2 = (1/((n - 1)^(p - 1) * (n - 2)))^(1/p)
        delta = avgdist2 - avgdist1
        temp0 = -delta/log(0.99)
    }
	c=0.95
	tempMin=1e-100
	
	desmat = LHD(n,p)[[2]]
	## generate a random LHD
	
	pfdes = desmat0 = desmat
	desval0 = Mm(desmat0)
	
	I=permutations(num,p,replace = TRUE)
	nodes=(I-0.5)/num # generate grid nodes matrix
	
	nodeDist=rdist(nodes,desmat0)
	dist=apply(nodeDist,1,min)
	fd=max(dist)
	fd=round(fd,4)
	desval0=c(desval0,fd)
	
	pfvals = matrix(desval0, nrow=1, ncol=2)
	
	for(s in 1:ncol(wtset))
	{
		desmat = desmat0
		bestvals = desval0
		
		i = 0
		count = 0
		temp=temp0
		while(i < maxiter & count < times)
		{
			newdesmat = desmat
			exchange_rows=sample(n,2)
			a=exchange_rows[1]
			b=exchange_rows[2]
			column=sample(p-1,1)+1 # The first column does not need to be exchanged
		
			newdesmat[a,column]=desmat[b,column]
			newdesmat[b,column]=desmat[a,column]
		
			nodeDist0=nodeDist
			nodeDist0[,c(a,b)]=rdist(nodes,newdesmat[c(a,b),])
		
			dist=apply(nodeDist0,1,min)
			fd0=max(dist)
			fd0=round(fd0,4)
			newbestval = c(Mm(newdesmat),fd0)
			
			if(sum(is.na(newbestval))==0)
			{
				wf <- function(x)
				{
					if(sum(is.na(x))==0)
					{
						xs=c((x[1]-crlim[1,1])/(crlim[1,1]-crlim[2,1]),(x[2]-crlim[1,2])/(crlim[1,2]-crlim[2,2]))
						## The sign needs to be changed based on the interest
						## of maximizing or minimizing the criteria
						temp=sum(xs*wtset[,s])
					}else
					{
						temp=NA
					}
					return(temp)
				}
				
				if (wf(newbestval)>wf(bestvals))
				{
					bestvals = newbestval
					desmat = newdesmat
					count = 0
				}
				else
				{
					if(temp > tempMin & runif(1)<exp((wf(newbestval)-wf(bestvals))/temp))
					{
						bestvals = newbestval
						desmat = newdesmat
						count = 0
					}
				}
				# only when the criterion is smaller than max crlim then consider adding it to parent front
				if (newbestval[1] < crlim[2,1] & newbestval[2] < crlim[2,2])
				{				
					temppf = checkon2(t(newbestval),newdesmat,pfvals,pfdes)
					pfvals = temppf[[1]]
					pfdes=temppf[[2]]
				}
			}
			i = i+1
			count = count+1
			temp=temp*c
		}
	}
	return(list(pfvals,pfdes))
}

pfpm <- function(n, p, 
				 crlim,
				 num,
				 nstarts=1,
				 times = 300,
				 maxiter = 1e+06,
				 temp0=0,				 
				 wtset=cbind(c(1,0),c(0.8,0.2),c(0.6,0.4),c(0.4,0.6),c(0.2,0.8),c(0,1)))
{
	time0=Sys.time()
	time00=Sys.time()
	
	temp=alg5(n, p, wtset, crlim, num, temp0, times, maxiter)
	paretovals=matrix(temp[[1]],ncol=2)
	paretodes=temp[[2]]

	for(icount in 2:nstarts)
	{
		temp=alg5(n, p, wtset, crlim, num, temp0, times, maxiter)
		temp[[1]]=matrix(temp[[1]],ncol=2)
	
		for(i in 1:nrow(temp[[1]]))
		{
			newpt = matrix(temp[[1]][i,],nrow=1)
			newdes = temp[[2]][,p*(i-1)+1:p]
			temppf=checkon2(newpt,newdes,paretovals,paretodes)
			paretovals=temppf[[1]]
			paretodes=temppf[[2]]
		}
		if(icount%%100==0)
		{
			time1=sprintf('100 simulations done in %5.3f Minutes\n  ***====================================================***\n',difftime(Sys.time(),time0,units='mins'))
			cat(time1)
			time0=Sys.time()
		}
	}
	rtime=difftime(Sys.time(),time00,units='mins')
	return(list("pfdes"=paretodes, "pfvals"=paretovals, "time_rec"=rtime))
}

alg5 <- function(n, p, wtset, crlim, num, temp0=0,times = 300, maxiter=1e+06)
{
	if (temp0 == 0) 
	{
        avgdist1 = 1/(n - 1)
        avgdist2 = (1/((n - 1)^(p - 1) * (n - 2)))^(1/p)
        delta = avgdist2 - avgdist1
        temp0 = -delta/log(0.99)
    }
	c=0.95
	tempMin=1e-100
	
	desmat = LHD(n,p)[[2]]
	## generate a random LHD
	
	pfdes = desmat0 = desmat
	desval0 = mp(desmat0)
	
	I=permutations(num,p,replace = TRUE)
	nodes=(I-0.5)/num # generate grid nodes matrix
	
	nodeDist=rdist(nodes,desmat0)
	dist=apply(nodeDist,1,min)
	fd=max(dist)
	fd=round(fd,4)
	desval0=c(desval0,fd)
	
	pfvals = matrix(desval0, nrow=1, ncol=2)
	
	for(s in 1:ncol(wtset))
	{
		desmat = desmat0
		bestvals = desval0
		
		i = 0
		count = 0
		temp=temp0
		while(i < maxiter & count < times)
		{
			newdesmat = desmat
			exchange_rows=sample(n,2)
			a=exchange_rows[1]
			b=exchange_rows[2]
			column=sample(p-1,1)+1 # The first column does not need to be exchanged
		
			newdesmat[a,column]=desmat[b,column]
			newdesmat[b,column]=desmat[a,column]
		
			nodeDist0=nodeDist
			nodeDist0[,c(a,b)]=rdist(nodes,newdesmat[c(a,b),])
		
			dist=apply(nodeDist0,1,min)
			fd0=max(dist)
			fd0=round(fd0,4)
			newbestval = c(mp(newdesmat),fd0)
			
			if(sum(is.na(newbestval))==0)
			{
				wf <- function(x)
				{
					if(sum(is.na(x))==0)
					{
						xs=c((x[1]-crlim[1,1])/(crlim[1,1]-crlim[2,1]),(x[2]-crlim[1,2])/(crlim[1,2]-crlim[2,2]))
						## The sign needs to be changed based on the interest
						## of maximizing or minimizing the criteria
						temp=sum(xs*wtset[,s])
					}else
					{
						temp=NA
					}
					return(temp)
				}
				
				if (wf(newbestval)>wf(bestvals))
				{
					bestvals = newbestval
					desmat = newdesmat
					count = 0
				}
				else
				{
					if(temp > tempMin & runif(1)<exp((wf(newbestval)-wf(bestvals))/temp))
					{
						bestvals = newbestval
						desmat = newdesmat
						count = 0
					}
				}
				# only when the criterion is smaller than max crlim then consider adding it to parent front
				if (newbestval[1] < crlim[2,1] & newbestval[2] < crlim[2,2])
				{				
					temppf = checkon2(t(newbestval),newdesmat,pfvals,pfdes)
					pfvals = temppf[[1]]
					pfdes=temppf[[2]]
				}
			}
			i = i+1
			count = count+1
			temp=temp*c
		}
	}
	return(list(pfvals,pfdes))
}