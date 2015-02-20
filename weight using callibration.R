######## We want to develop a meth by weighting adjustment for samples in rotation samples.
### We are interested in estimating population mean and median and Gini coeffs.
### Functions and variables explained.

# library(snowfall)

library(reldist)
library(matrixStats)


### Generate Finite Population
## N:population Size
## n:sample size
## n_g:rotation sample size
## y: responses
## t:time, 1,2,3,4,5.
### For t=1

###################################################################
#### Function for sampling
##### Use SRS to get the samples
#### In this function we input population size N and sample size n and partition size ng
### This function will return a matrix: the row is the sample index ,the column is the time t.
sampling=function(N,n,ng)
{
  A=sample(1:N,n,replace=F)
  if(n%%ng !=0)
    print("Be Careful for n_gs are not eaqual!")
  t=n/ng
  S=NULL
  for(i in 1:t)
  {
    S=cbind(S,sample(A[! A %in%S],ng,replace=F))
  }
  At=NULL
  for(i in 1:t)
  {
    At=cbind(At,A[! A%in%S[,i]])
  }
  return(At)
}

######################################################################
## Function to category response
### We want try different methods to group
### Like kmeans quantile and uniform
### delta is the difference between to closing groups

#### This function in put a variable a and the points of segment that is a vector.
### This function return the index which group the a belongs
group=function(a,point)
{
  p=length(point)
  if(a<=point[1])
    return(1)
  for(i in 2:p)
  {
    if(a>point[i-1] & a<=point[i])
      return(i)
  }
  if(a>point[p])
    return(p+1)
}
#########################
#### This function inputs the data set and the number of groups
### This function returns a list, the first of the list id the grouped response data frame. 
### The second of the list is the points used to category
category=function(y,G)
{
  delta=diff(range(y))/G
  temp=NULL
  for(i in 1:t)
    temp=c(temp,y[,i])
  ### seprate points
  point=quantile(temp,probs=seq(1/G,1-1/G,by=1/G))
  yc=y
  for(i in 1:dim(y)[1])
    for(j in 1:dim(y)[2])
      yc[i,j]=group(y[i,j],point)
  return(list(group=yc,point=point))
}



category1=function(y,G)
{
  temp=NULL
  for(i in 1:t)
    temp=rbind(temp,quantile(y[,i],probs=seq(1/G,1-1/G,by=1/G)))
  ### seprate points
  yc=y
  for(i in 1:dim(y)[1])
    for(j in 1:dim(y)[2])
      yc[i,j]=group(y[i,j],temp[j,])
  return(list(group=yc,point=temp))
}

category2=function(y,G)
{
  temp=NULL
  for(i in 1:t)
    temp=cbind(temp,kmeans(y[,i],G)$cluster)
  yc=temp
  return(list(group=yc,point=NULL))
}
####
#### Or use kmeans to group 
#### Input the all response as a vector and the number of the group
#### return a matrix that is groups we have done
group_kmeans=function(y_vector,G)
{
  y_vector=log(y_vector)
  r=kmeans(y_vector,G)$cluster
  matrix(r,N,t)
}

#####################################################################
### Function to do calibrartion
### X is always observed
## S is the sample index
## Input the auxilary information and all sample index
## Return the weights 
calibration=function(X,A)
{
  n=length(A)
  N=length(X)
  d=N/n
  X_total=sum(X)
  X_HT=sum(X[A]*d)
  Dom=sum(X[A]*X[A]*d)
  weight=NULL
  for(i in 1:n)
    weight=c(weight,d*(1+(X_total-X_HT)/Dom*X[A[i]]))
  w=data.frame(A=A,weight=weight)
  w
  
}

### Actually we don't use in this simulation beause all weights are equal by SRS
###################################################################
####  function: COunt the number of pattern
### Input the the pattern we want to search That is a vector
### return an iteger that is a count
count=function(y_prime)
{
  if(any(is.na(y_prime)))
  {
    p1=which(is.na(y_prime))
    index=1:length(y_prime)
    index=index[-p1]
    data=sample_data[apply(sample_data[,-p1],1,function(x) !any(is.na(x))),-p1]
    l=sum(apply(data,1,function(x) all(x==y_prime[index])))
    return(l)
    
  }
  else
  {
    index=1:length(y_prime)
    data=sample_data[apply(sample_data,1,function(x) !any(is.na(x))),]
    l=sum(apply(data,1,function(x) all(x==y_prime[index])))
    return(l)
  }
}

### To check the function we have another function to do the same thing.

### Input a pattern that is a vector and a matrix which is the sample in each year: row is the sample, the colnum is time
### return a lits. The first of the list is the count of the pattern and the second is the index of the same pattern in the population 
count1=function(y_prime,At)
{
  l=length(y_prime)
  A=unique(as.numeric(At))
  if(!is.na(y_prime[1]))
    temp=At[pop_group[At[,1],1]==y_prime[1],1] ### index: to store index which is possible to satisfy this pattern
  else temp=A[!A%in%At[,1]]
  for(i in 2:l)
  {
    if(!is.na(y_prime[i]))
    {
      o=At[pop_group[At[,i],i]==y_prime[i],i]
      temp=o[o%in%temp]
    }
    else
    {
      S=A[!A%in%At[,i]]
      temp=S[S%in% temp]
    }
    
  }
  return(list(number=length(temp),index=temp))
}

############################################################################################

### function for kinds of pattern
### Input G which is  the number of the group. t is the time 
### return a matrix that is all possible pattern by the group and the time 
pattern=function(G,t)
{
  if(t==1)
    return(matrix(1:G,G,1))
  else
  {
    return(cbind(rep(1:G,each=dim(pattern(G,t-1))[1]),rep(1,G)%x%pattern(G,t-1)))
  }
}

#### Use this function we can generate all possible patterns 

######################################################################################
### Function to find the poosible patterns in sample. 
### This help to reduce the heavy of computation
### Input the sample data which is the sample data we have grouped the miss place we use NA
### Return a complete matrix which is all possible patterns
possible_pat=function(sample_data)
{
  G=max(sample_data,na.rm=T)
  t=dim(sample_data)[2]
  pat=NULL
  for( i in 1:dim(sample_data)[1])
  {
    x=as.numeric(sample_data[i,])
    j=which(is.na(x))
    r=matrix(rep(x,G),ncol=t,byrow=T)
    r[,j]=1:G
    colnames(r)=paste("y",1:5,sep="")
    pat=rbind(pat,r)
  }
  pat
}

######################################################################################
### function for EM
### Input maxiter the number we wanto iteration and pop_group is the population we have grouped 
### Input At is the sample index which is a matrix . the row is the index in each time and the colnum is time
### joint_prob is the all poosible patterns and the intial probs. Which is a data frame
### Return a list which is the first is joint probs for all possibale pattern thye prob is unscaled.
## The second of the list is the trace of the function. It is a data frame. each colnum is the probs in each iteration.

EM=function(maxiter,pop_group,At,joint_prob)
{
  G=max(pop_group)
  N=dim(pop_group)[1]
  time=dim(pop_group)[2]
  p=dim(joint_prob)[1]
  weight=1
  trace=NULL
  ### We have 5^5=3125 patterns
  for(t in 1:maxiter)
  {
    joint=sapply(1:p,function(i) 
    {
      y_prime=as.numeric(joint_prob[i,1:time])
      nt=count(y_prime)*weight
      for(j in 1:time)
      {
        if(joint_prob[i,time+1]==0)
          break
        index=1:time
        y_prime=as.numeric(joint_prob[i,1:time])
        y_prime[j]=NA
        index=index[-j]
        if(count(y_prime)!=0)
        {
          nplus=count(y_prime)*weight*joint_prob[i,time+1]/ifelse(
            sum(joint_prob[apply(joint_prob[,index],1,function(x) all(x==y_prime[index])),time+1])==0,1,
            sum(joint_prob[apply(joint_prob[,index],1,function(x) all(x==y_prime[index])),time+1]))
          
        }
        else nplus=0
        nt=nt+nplus
      }
      return(nt)
    })
    #     if(sum((joint-joint_prob$prob)^2)<0.1)
    #       break
    joint_prob$prob=joint
    trace=cbind(trace,joint)
    cat("total=",sum(joint),"\n")
    cat("t=",t,"\n")
  }
  return(list(prob=joint_prob,trace=trace))
}

####################################################################################################

### Function for weights after the poststratification
### To use the poststratification to adjust weights.
## Input o is the joint probs and At is the sample index matrix and pop_group is the population 
## return a list the first of the list is weight which is a matrix and the second of the list is the the margin of the dists

final_weight=function(o,At,pop_group)
{
  time=length(At[1,])
  G=max(pop_group)
  weight=matrix(0,dim(At)[1],dim(At)[2])
  margin=matrix(0,G,dim(At)[2])
  for(i in 1:G)
    for(j in 1:time)
      margin[i,j]=sum(o[o[,j]==i,"prob"])
  for (t in 1:time)
  {
    s=At[,t]
    A=pop_group[s,t]
    for(j in 1:length(A))
    {
      weight[j,t]=margin[A[j],t]/sum(A==A[j])
    }
    
  }
  return(list(weight=weight,margin=margin))
}

#######################################################################
est=function(At,y,finalweight,G)
{
  m=NULL
  for( t in 1:5)
    m=c(m,weighted.mean(y[At[,t],t],w=finalweight[,t]))
  
  #### weighted median  
  med=NULL
  for( t in 1:5)
    med=c(med,weightedMedian(y[At[,t],t],w=finalweight[,t],ties="max"))
  
  gini_coef=NULL
  for( t in 1:5)
    gini_coef=c(gini_coef,gini(y[At[,t],t],weights=finalweight[,t])) 
  
  r=rbind(m,med,gini_coef)
  rownames(r)=paste(rownames(r),G,sep="")
  r
}

#######################################################################################

##############################################################################################
### Simmulate the data set
### Generate Finite Population
## N:population Size
## n:sample size
## n_g:rotation sample size
## y: responses
## t:time, 1,2,3,4,5.
### For t=1

N=10000
t=5
y=matrix(0,N,t)
y[,1]=exp(rnorm(N))
y_vector=y[,1]
for(t in 2:5)
{
  y[,t]=exp(rnorm(N))*y[,t-1]
  y_vector=c(y_vector,y[,t])
}

y=data.frame(y)
names(y)=c("y1","y2","y3","y4","y5")
head(y)
write.csv(y,"population.csv")
### Population mean
pop_mean=apply(y,2,mean)

#### Population median
pop_median=apply(y,2,median)

###Population Gini coeffs
pop_gini=apply(y,2,function(x) gini(x))
r=rbind(pop_mean,pop_median,pop_gini)

#### sample 
n=500
ng=100

At=sampling(N,n,ng)
### Total sample
A=unique(as.numeric(At))

######################################################
### Method 1: Direct method
direct=function(y,At)
{
  est_mean=NULL
  est_median=NULL
  est_gini=NULL
  
  for(i in 1:t)
  {
    yt=y[At[,i],i]
    est_mean=c(est_mean,mean(yt))
    est_median=c(est_median,median(yt))
    est_gini=c(est_gini,gini(yt))
  }
  
  r=rbind(est_mean,est_median,est_gini)
  r
  
}

r=rbind(r,direct(y,At))

#################################################################
## Proposed method
### G=5
### quantile seperate
y_group=category(y,5)
y_group=category1(y,5)
y_group=category2(y,5)
### So after the category, the population is:
pop_group=y_group$group
#### To look at the sample data set
sample_data=pop_group[A,]
for( i in 1:t)
{
  sample_data[!A%in%At[,i],i]=NA
}
head(sample_data)

write.csv(sample_data,"sampledata5.csv")
### Initial probs
pat=data.frame(possible_pat(sample_data))
joint_prob=data.frame(pat,prob=rep(1/dim(pat)[1],dim(pat)[1]))
result=EM(40,pop_group,At,joint_prob)
o=result$prob
o$prob=o$prob/500
write.csv(o,"prob5.csv")
write.csv(result$trace,"trace5.csv")
finalweight=final_weight(o,At,pop_group)$weight

##################
d=est(At,y,finalweight,5)

r=rbind(r,d)

#### G=10
y_group=category(y,10)
pop_group=y_group$group
#### To look at the sample data set
sample_data=pop_group[A,]
for( i in 1:t)
{
  sample_data[!A%in%At[,i],i]=NA
}
head(sample_data)
write.csv(sample_data,"sampledata10.csv")
### Initial probs
pat=data.frame(possible_pat(sample_data))
joint_prob=data.frame(pat,prob=rep(1/dim(pat)[1],dim(pat)[1]))
result=EM(40,pop_group,At,joint_prob)
o2=result$prob
o2$prob=o2$prob/500
finalweight=final_weight(o2,At,pop_group)$weight
r=rbind(r,est(At,y,finalweight,10))
write.csv(o2,"prob10.csv")
write.csv(result$trace,"trace10.csv")


#### We try log trans and kmeans group
y_group=group_kmeans(y_vector,5)
pop_group=y_group
#### To look at the sample data set
sample_data=pop_group[A,]
for( i in 1:t)
{
  sample_data[!A%in%At[,i],i]=NA
}
head(sample_data)
write.csv(sample_data,"sampledata_k5.csv")
### Initial probs
pat=data.frame(possible_pat(sample_data))
joint_prob=data.frame(pat,prob=rep(1/dim(pat)[1],dim(pat)[1]))
result=EM(40,pop_group,At,joint_prob)
o=result$prob
finalweight=final_weight(o,At,pop_group)$weight
d1=est(At,y,finalweight,10)
rownames(d1)=paste(rownames(d1),"Kmeans",sep="_")
r=rbind(r,d1)
write.csv(o,"prob_k5.csv")
write.csv(result$trace,"trace_k5.csv")
####### 
G=10
y_group=group_kmeans(y_vector,10)
pop_group=y_group
#### To look at the sample data set
sample_data=pop_group[A,]
for( i in 1:t)
{
  sample_data[!A%in%At[,i],i]=NA
}
head(sample_data)
write.csv(sample_data,"sampledata_k10.csv")
### Initial probs
pat=data.frame(possible_pat(sample_data))
joint_prob=data.frame(pat,prob=rep(1/dim(pat)[1],dim(pat)[1]))
result=EM(40,pop_group,At,joint_prob)
o=result$prob
finalweight=final_weight(o,At,pop_group)$weight
d1=est(At,y,finalweight,10)
rownames(d1)=paste(rownames(d1),"Kmeans",sep="_")
r=rbind(r,d1)
write.csv(o,"prob_k10.csv")
write.csv(result$trace,"trace_k10.csv")
write.csv(r,"finalresults.csv")