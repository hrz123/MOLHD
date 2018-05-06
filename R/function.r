rdist <- fields::rdist
permutations <- arrangements::permutations
runif <- stats::runif


x2D <- function(x, D)
{
    if(!is.vector(x))
	{
		stop("Error: x is not a vector")
	}
    if(!is.matrix(D))
	{
		stop("Error: D is not a matrix")
	}
    if(length(x)!=ncol(D))
	{
		stop("Error: x and D are not in the same dimension")
	}
	x=matrix(x,1) # transform x to matrix
	dist=rdist(x,D)
	d=min(dist)
	n=sum(dist==d)
	result=list("MinimumDistance"=d,"number"=n)
	return(result)
}



md <- function(D)
{
	if(!is.matrix(D))
	{
		stop("Error: D is not a matrix")
	}
	dist=rdist(D)
	dist=dist[lower.tri(dist)]
	d=min(dist)
	n=sum(dist==d)
	result=list("MinimumDistance"=d,"number"=n)
	return(result)
}


miM <- function(D, num=50)
{
	if(!is.matrix(D))
	{
		stop("Error: D is not a matrix")
	}
	p=ncol(D)
	I=permutations(num,p,replace = TRUE)
	nodes=(I-0.5)/num # generate grid nodes matrix

	nodeDist=rdist(nodes,D)
	dist=apply(nodeDist,1,min)
	d=max(dist)
	return(round(d,4)) # four decimals
}


LHD <- function(n, p)
{
	if(n%%1!=0 | p%%1!=0)
	{
		stop("Error: n and p must be integers")
	}
	if(n<=0 | p<=0)
	{
		stop("Error: n and p must be positive")
	}
	design=matrix(c(1:n,unlist(replicate((p-1),sample(n)))),n,p)
	standDesign=(design-0.5)/n
	return(list(design=design,standDesign=standDesign))
}


miMLHD <- function(n, p, num=50, temp0=0,nstarts=1, times=300, maxiter=1e+06)
{
	time0 = Sys.time()
	if(maxiter<1)
	{
		stop("Error: The maximum number of iterations must be at least 1.")
	}
	if(n<=1)
	{
		stop("Error: D has less than or equal to 1 point")
	}
	if (temp0 == 0)
	{
        avgdist1 = 1/(n - 1)
        avgdist2 = (1/((n - 1)^(p - 1) * (n - 2)))^(1/p)
        delta = avgdist2 - avgdist1
        temp0 = -delta/log(0.99)
    }
	c=0.95
	tempMin=1e-100
	des=NULL
	val=NULL
	iter=0
	for (s in 1:nstarts)
	{
		D = LHD(n,p)[[2]]

		I=permutations(num,p,replace = TRUE)
		nodes=(I-0.5)/num # generate grid nodes matrix

		nodeDist=rdist(nodes,D)
		dist=apply(nodeDist,1,min)

		fd=max(dist)
		fd=round(fd,4)

		i=0
		count=0
		temp=temp0
		while(i<maxiter & count<times)
		{
			D0=D
			exchange_rows=sample(n,2)
			a=exchange_rows[1]
			b=exchange_rows[2]
			column=sample(p-1,1)+1 # The first column does not need to be exchanged

			D0[a,column]=D[b,column]
			D0[b,column]=D[a,column]

			nodeDist0=nodeDist
			nodeDist0[,c(a,b)]=rdist(nodes,D0[c(a,b),])

			dist=apply(nodeDist0,1,min)
			fd0=max(dist)
			fd0=round(fd0,4)

			if(fd0 < fd)
			{
				D=D0;fd=fd0;nodeDist=nodeDist0;count=0
			}
			else
			{
				if(temp > tempMin & runif(1)<exp((fd-fd0)/temp))
				{
					D=D0;fd=fd0;nodeDist=nodeDist0;count=0
				}
			}
			i=i+1
			count=count+1
			temp=temp*c
		}
		des=cbind(des,D)
		val=c(val,fd)
		iter=iter+i
	}
	fd=min(val)
	D=des[,(order(val)[1]-1)*p+1:p]
	rtime=difftime(Sys.time(),time0,units='secs')
	result=list("design"=D, "criterion"=fd, "iterations"=iter, "time_rec"=rtime)
	return(result)
}


Mm <- function(D, power=100)
{
	if(!is.matrix(D))
	{
		stop("Error: D is not a matrix")
	}
	dist=rdist(D)
	dist=dist[lower.tri(dist)]
	d=sum(dist^(-power))^(1/power)
	return(round(d,4)) # four decimals
}


MmLHD <- function(n, p, power=100, temp0=0, nstarts=1, times=300, maxiter=1e+06)
{
	time0 = Sys.time()
	if(maxiter<1)
	{
		stop("Error: The maximum number of iterations must be at least 1.")
	}
	if(n<=1)
	{
		stop("Error: D has less than or equal to 1 point")
	}
	if (temp0 == 0)
	{
        avgdist1 = 1/(n - 1)
        avgdist2 = (1/((n - 1)^(p - 1) * (n - 2)))^(1/p)
        delta = avgdist2 - avgdist1
        temp0 = -delta/log(0.99)
    }
	c=0.95
	tempMin=1e-100
	des=NULL
	val=NULL
	iter=0
	for (s in 1:nstarts)
	{
		D = LHD(n,p)[[2]]
		i=0
		count=0
		temp=temp0
		while(i<maxiter & count<times)
		{
			D0=D
			exchange_rows=sample(n,2)
			a=exchange_rows[1]
			b=exchange_rows[2]
			column=sample(p-1,1)+1 # The first column does not need to be exchanged

			D0[a,column]=D[b,column]
			D0[b,column]=D[a,column]

			D1=D[-c(a,b),]
			distnew=rdist(D0[c(a,b),],D1)
			distcur=rdist(D[c(a,b),],D1)
			Cdiff=0
			if(sum(distnew[1,]==distcur[1,])!=p & sum(distnew[1,]==distcur[2,])!=p)
			{
				Cdiff=sum(distnew^(-power))-sum(distcur^(-power))
			}

			if(Cdiff < 0)
			{
				D=D0;count=0
			}
			else
			{
				if(temp > tempMin & runif(1)<exp(-Cdiff^(1/power)/temp))
				{
					D=D0;count=0
				}
			}
			i=i+1
			count=count+1
			temp=temp*c
		}
		criterion=Mm(D,power)
		des=cbind(des,D)
		val=c(val,criterion)
		iter=iter+i
	}
	criterion=min(val)
	D=des[,(order(val)[1]-1)*p+1:p]
	rtime=difftime(Sys.time(),time0,units='secs')
	result=list("design"=D, "criterion"=criterion, "iterations"=iter, "time_rec"=rtime)
	return(result)
}

mp <- function(D)
{
	if(!is.matrix(D))
	{
		stop("Error: D is not a matrix")
	}
	n=nrow(D)
	p=ncol(D)
	d=0
	for(i in 1:(n-1))
	{
		D0=D[-(1:i),]
		if(i == n-1)
		{
			D0=matrix(D0,1)
		}
		D0=D0-matrix(rep(D[i,],n-i),n-i,byrow=TRUE)
		D0=D0^(-2)
		d0=apply(D0,1,prod)
		d=d+sum(d0)
	}
	d=(d/choose(n,2))^(1/p)
	return(round(d,4)) # four decimals
}


mpLHD <- function(n, p, temp0=0,nstarts=1, times=300, maxiter=1e+06)
{
	time0 = Sys.time()
	if(maxiter<1)
	{
		stop("Error: The maximum number of iterations must be at least 1.")
	}
	if(n<=1)
	{
		stop("Error: D has less than or equal to 1 point")
	}
	if (temp0 == 0)
	{
        avgdist1 = 1/(n - 1)
        avgdist2 = (1/((n - 1)^(p - 1) * (n - 2)))^(1/p)
        delta = avgdist2 - avgdist1
        temp0 = -delta/log(0.99)
    }
	c=0.95
	tempMin=1e-100
	des=NULL
	val=NULL
	iter=0
	for (s in 1:nstarts)
	{
		D = LHD(n,p)[[2]]
		i=0
		count=0
		temp=temp0
		while(i<maxiter & count<times)
		{
			D0=D
			exchange_rows=sample(n,2)
			column=sample(p-1,1)+1 # The first column does not need to be exchanged
			a=exchange_rows[1]
			b=exchange_rows[2]

			D0[a,column]=D[b,column]
			D0[b,column]=D[a,column]

			D1=D[-c(a,b),-column]
			D2=D[-c(a,b),]
			D3=matrix(rep(D[a,-column],n-2),n-2,byrow=TRUE)-D1
			D4=matrix(rep(D[b,-column],n-2),n-2,byrow=TRUE)-D1
			D3=D3^(-2)
			D4=D4^(-2)
			v1=apply(D3,1,prod)
			v2=apply(D4,1,prod)
			v3=(D[b,column]-D2[,column])^(-2)-(D[a,column]-D2[,column])^(-2)
			Cdiff=sum((v1-v2)*v3)
			if(Cdiff < 0)
			{
				D=D0;count=0
			}
			else
			{
				if( temp > tempMin & runif(1)<exp(-Cdiff^(1/p)/temp))
				{
					D=D0;count=0
				}
			}
			i=i+1
			count=count+1
			temp=temp*c
		}
		des=cbind(des,D)
		criterion=mp(D)
		val=c(val,criterion)
		iter=iter+i
	}
	criterion=min(val)
	D=des[,(order(val)[1]-1)*p+1:p]
	rtime=difftime(Sys.time(),time0,units='secs')
	result=list("design"=D, "criterion"=criterion, "iterations"=iter, "time_rec"=rtime)
	return(result)
}



pfMp <- function(n, p,
				crlim,
				nstarts=1,
				times = 300,
				maxiter = 1e+06,
				temp0=0,
				wtset=cbind(c(1,0),c(0.8,0.2),c(0.6,0.4),c(0.4,0.6),c(0.2,0.8),c(0,1)))
{
	time0=Sys.time()
	time00=Sys.time()


	temp=alg2(n, p, wtset, crlim, temp0, times, maxiter)
	paretovals=matrix(temp[[1]],ncol=2)
	paretodes=temp[[2]]

	for(icount in 2:nstarts)
	{
		temp=alg2(n, p, wtset, crlim, temp0, times, maxiter)
		temp[[1]]=matrix(temp[[1]],ncol=2)

		for(i in 1:dim(temp[[1]])[1])
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


pfMpm <- function(n, p,
				 crlim,
				 num,
				 nstarts=1,
				 times = 300,
				 maxiter = 1e+06,
				 temp0=0,
				 wtset=cbind(c(1,0,0),c(0.5,0.5,0),c(0.5,0,0.5),c(0,0.5,0.5),c(0,1,0),c(0,0,1),c(1/3,1/3,1/3)))
{
	time0=Sys.time()
	time00=Sys.time()

	temp=alg3(n, p, wtset, crlim, num, temp0, times, maxiter)
	paretovals=matrix(temp[[1]],ncol=3)
	paretodes=temp[[2]]

	for(icount in 2:nstarts)
	{
		temp=alg3(n, p, wtset, crlim, num, temp0, times, maxiter)
		temp[[1]]=matrix(temp[[1]],ncol=3)

		for(i in 1:nrow(temp[[1]]))
		{
			newpt = matrix(temp[[1]][i,],nrow=1)
			newdes = temp[[2]][,p*(i-1)+1:p]
			temppf=checkon3(newpt,newdes,paretovals,paretodes)
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


## Update the current Pareto front and Pareto set for any newly generated design for optimization based on three criteria


checkon3 <- function(newpt, newdes, curpf, curpfdes)
{
	## this needs to be changed based on the interest
	## of maximizing or minimizing the criteria
	g1=newpt[1,1]<curpf[,1] # If you want to minimize criterion1 set it to "<"
	g2=newpt[1,2]<curpf[,2]
	g3=newpt[1,3]<curpf[,3]

	ge1=newpt[1,1]<=curpf[,1]
	ge2=newpt[1,2]<=curpf[,2]
	ge3=newpt[1,3]<=curpf[,3]

	l1=newpt[1,1]>curpf[,1]
	l2=newpt[1,2]>curpf[,2]
	l3=newpt[1,3]>curpf[,3]

	le1=newpt[1,1]>=curpf[,1]
	le2=newpt[1,2]>=curpf[,2]
	le3=newpt[1,3]>=curpf[,3]

	eq1=newpt[1,1]==curpf[,1]
	eq2=newpt[1,2]==curpf[,2]
	eq3=newpt[1,3]==curpf[,3]

	cond1=(g1*ge2*ge3+g2*ge1*ge3+g3*ge1*ge2)==0 # vector of logic values, true means the current point on pareto front is not dominated by the new point
	cond2=sum(l1*le2*le3+l2*le1*le3+l3*le1*le2+eq1*eq2*eq3) # cond2=0 means add the point
	cond3=seq(1,dim(curpf)[1])[cond1]
	ndes=dim(newdes)[2]
	cond4=rep(0,ndes*length(cond3))
	if(length(cond3)>0) # remove the points that are dominated by the new point
	{
		for(i in 1:length(cond3))
		{
			cond4[(i-1)*ndes+1:ndes]=(cond3[i]-1)*ndes+1:ndes
		}
	}

	newpf=curpf[cond1,]
	newpfdes=curpfdes[,cond4]
	if(cond2==0) # add the new point
	{
		newpf=rbind(newpf,newpt)
		newpfdes=cbind(newpfdes,newdes)
	}
	return(list(matrix(newpf,ncol=3),newpfdes))
}


## Update the current Pareto front and Pareto set for any newly generated design for optimization based on two criteria


checkon2 <- function(newpt, newdes, curpf, curpfdes)
{
	## this needs to be changed based on the interest
	## of maximizing or minimizing the criteria
	g1=newpt[1,1]<curpf[,1] # If you want to minimize criterion1 set it to "<"
	g2=newpt[1,2]<curpf[,2]

	ge1=newpt[1,1]<=curpf[,1]
	ge2=newpt[1,2]<=curpf[,2]

	l1=newpt[1,1]>curpf[,1]
	l2=newpt[1,2]>curpf[,2]

	le1=newpt[1,1]>=curpf[,1]
	le2=newpt[1,2]>=curpf[,2]

	eq1=newpt[1,1]==curpf[,1]
	eq2=newpt[1,2]==curpf[,2]

	cond1=(g1*ge2+g2*ge1)==0 # vector of true and false, true means the current point on pareto front is not dominated by the new point
	cond2=sum(l1*le2+l2*le1+eq1*eq2) # cond2=0 means add the point
	cond3=seq(1,dim(curpf)[1])[cond1]
	ndes=dim(newdes)[2]
	cond4=rep(0,ndes*length(cond3))
	if(length(cond3)>0) # remove the points that are dominated by the new point
	{
		for(i in 1:length(cond3))
		{
			cond4[(i-1)*ndes+1:ndes]=(cond3[i]-1)*ndes+1:ndes
		}
	}

	newpf=curpf[cond1,]
	newpfdes=curpfdes[,cond4]
	if(cond2==0) # add the new point
	{
		newpf=rbind(newpf,newpt)
		newpfdes=cbind(newpfdes,newdes)
	}
	return(list(matrix(newpf,ncol=2),newpfdes))
}


evalfunc2 <- function(mat)
{
	maximinLHDCriterion = Mm(mat)
	maxproCriterion = mp(mat)

	return(c(maximinLHDCriterion,maxproCriterion))
}


alg3 <- function(n, p, wtset, crlim, num, temp0=0,times = 300, maxiter=1e+06)
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
	desval0 = evalfunc2(desmat0)

	I=permutations(num,p,replace = TRUE)
	nodes=(I-0.5)/num # generate grid nodes matrix

	nodeDist=rdist(nodes,desmat0)
	dist=apply(nodeDist,1,min)
	fd=max(dist)
	fd=round(fd,4)
	desval0=c(desval0,fd)

	pfvals = matrix(desval0, nrow=1, ncol=3)

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
			newbestval = c(evalfunc2(newdesmat),fd0)

			if(sum(is.na(newbestval))==0)
			{
				wf <- function(x)
				{
					if(sum(is.na(x))==0)
					{
						xs=c((x[1]-crlim[1,1])/(crlim[1,1]-crlim[2,1]),(x[2]-crlim[1,2])/(crlim[1,2]-crlim[2,2]),(x[3]-crlim[1,3])/(crlim[1,3]-crlim[2,3]))
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
				if (newbestval[1] < crlim[2,1] & newbestval[2] < crlim[2,2] & newbestval[3] < crlim[2,3])
				{
					temppf = checkon3(t(newbestval),newdesmat,pfvals,pfdes)
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



alg2 <- function(n, p, wtset, crlim, temp0=0, times = 300, maxiter=1e+06)
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
	desval0 = evalfunc2(desmat0)
	pfvals = matrix(desval0, nrow=1, ncol=2)

	for(s in 1:dim(wtset)[2])
	{
		desmat = desmat0
		bestvals = desval0

		i=0
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

			newbestval = evalfunc2(newdesmat)

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
					count = 0
					desmat = newdesmat
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
