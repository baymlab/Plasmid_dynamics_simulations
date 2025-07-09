library("VGAM")

###create directories to store these files
home<-"/simulation"
prog<-"/simulation/progress"
setwd("/Users/fernandorossine/Desktop/plasmid/simulation")

home.mat<-"/simulation/mat"


daughter0<-function(N0,mu,s,phi,S_mother)  
{
  S<-mu+(S_mother-mu)*phi+rnorm(n=1,sd=s)
  N<-exp(S)
  N<-floor(N)+((N-floor(N))>runif(n=1))
  N<-max(N0,N)
  Nf<-rbinom(n=1,size=N,prob=0.5)
  return(c(Nf,S,N))
}

daughter<-function(pls,mu,s,phi,S_mother,eclipse=1,segregation="bin",r_eff=1)
{
  N0<-sum(pls)
  
  b<-pls[1]
  r<-pls[2]
  
  S<-mu+(S_mother-mu)*phi+rnorm(n=1,sd=s)
  N<-exp(S)
  N<-floor(N)+((N-floor(N))>runif(n=1))
  N<-max(N0,N)
  
  be<-0
  re<-0
  
  while (sum(b,r,be,re)<N)#####plasmids replicate up to N copies
  {
    re_f<-re*r_eff
    r_f<-r*r_eff
    id_vec<-c(b,(be*eclipse)+b,r_f+(be*eclipse)+b,((re_f+be)*eclipse)+r_f+b)/(((re_f+be)*eclipse)+r_f+b)
    id<-sum(runif(n=1)>id_vec)+1  
    if(id==1)### b replicates
    {
      b<-b-1
      be<-be+2
    }
    if(id==2)###  eclipse be replicates
    {
      be<-be+1
    }
    if(id==3)### r replicates
    {
      r<-r-1
      re<-re+2
    }
    if(id==4)### eclipsed re replicates
    {
      re<-re+1
    }
  
  }
  
  ####after replication cycle ends, returning eclipsed plasmids to replication pool
  b<-b+be
  r<-r+re
  
  
  
  #####binomial segregation
  if(segregation=="bin")
  {
    bf<-rbinom(n=1,size=b,prob=0.5)
    rf<-rbinom(n=1,size=r,prob=0.5)
  }
  
  
  #####equal copy number segregation
  Nf<-floor(N/2)+((N/2)-floor(N/2)>runif(n=1))
  if(segregation=="eq")
  {
    bf<-sum(sample(c(rep(1,b),rep(0,r)),size=Nf))
    rf<-Nf-bf
  }
 
  
  
  if(bf+rf>0)
  {
   return(list(c(bf,rf),S,N))
  }else{
    return(list(c(b,r),S,N))
  }
}


phi1<-0.2 ###autocorrelation between 0 and 1
s1<-10 ### stationary standard deviation larger or equal than 0
M<-40 ### mean plasmid copy number at cell division
r_effective<-1 ### within-cell fitness avantage of a plasmid. larger or equal tahn 1
ecl<-0.95 ###weight of eclipsed plasmids

####auxiliary mean and variance in log space
mu_ln<-log((M^2)/sqrt((M^2)+s1^2))
sig_ln<-sqrt(log(1+s1^2/(M^2))*(1-phi1^2))




###bimomial plasmid segregation or equal plasmid segregation
segr<-"bin"
###segr<-"eq"




S_0<-M
N_0<-mu_ln

####generating stationary palsmid copy number before initializing cells with equal amounts of each plasmid
for (i in 1:100000) {
  temp<-daughter0(N0=N_0[i],mu=mu_ln,phi = phi1,S_mother = S_0[i],s=sig_ln)
  N_0[i+1]<-temp[1]
  S_0[i+1]<-temp[2]
}

###Removing begining so that only stationary distribution is used
N_0<-N_0[10000:length(N_0)]
S_0<-S_0[10000:length(S_0)]



    
    n.hz<-50 ##number of time steps used to calculate heteroplasmy decay
    n.samp<-40000 ## number of cells sampled to calculate heteroplasmy decay
    mat.hz<-matrix(0,n.hz,n.samp)
    
    
    for(i in 1:n.samp)
    {
      ind<-sample(size =1,length(N_0))
      pls1<-c(max(N_0[ind],1),max(N_0[ind],1)) ####initializing a cell with equal amounts of both plasmids. the amount is sampled from the stationary distribution
      S_time<-S_0[ind]
      hz<-0.5
      for (j in 1:n.hz)
      {
        temp<-daughter(pls=pls1,mu=mu_ln,phi = phi1,S_mother = S_time[j],s=sig_ln,eclipse = ecl,segregation = segr)
        pls1<-temp[[1]]
        S_time[j+1]<-temp[[2]]
        hz[j]<-2*(pls1[1]*pls1[2]/(sum(pls1)^2))
      }
      mat.hz[,i]<-hz
    }
    
    setwd(home.mat)
    
    save(mat.hz,file = paste0("MatHZ_phi_",phi1,"sd_",s1,"M",M,"ecl_",ecl, "seg_",segr,".RData"))
  
  


setwd(home)



yy<-log(apply(FUN=mean,mat.hz,1),base=2)
xx<-1:(length(xx))

###plotting log-heteroplasmy as a functionof time (in cell divisions), for the set parameters
plot(yy,xx)





