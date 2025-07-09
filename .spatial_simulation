

library("soundgen")
  library(tripack)
  library(plotrix)
  library(png)
  library(abind)

  
  range<-1000
  

max_x<-500 #spatial dimensions in pixels. one cell per pixel
max_y<-500



t_max<-100000 # number of cell divisions in simulation



daughter2<-function(pls,mu,s,phi,S_mother,eclipse=1,segregation="bin",r_eff=1)
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
  
  
  
 
    return(list(c(bf,rf),S,N,c(b-bf,r-rf)))
  
}

daughter0<-function(N0,mu,s,phi,S_mother)
{
  S<-mu+(S_mother-mu)*phi+rnorm(n=1,sd=s)
  N<-exp(S)
  N<-floor(N)+((N-floor(N))>runif(n=1))
  N<-max(N0,N)
  Nf<-rbinom(n=1,size=N,prob=0.5)
  return(c(Nf,S,N))
}

### division timing parms
sd_time<-0.2  ###cell time untill division is lognormally distributed
mean_time<-1

mu_time<-log((1^2)/sqrt((1^2)+sd_time^2))
sig_time<-sqrt(log(1+sd_time^2/(1^2)))

n_time<-function(mean_t,sd_t) ###samples a division time for a cell
{
  mu_time<-log((mean_t^2)/sqrt((mean_t^2)+sd_t^2))
  sig_time<-sqrt(log(1+sd_t^2/(mean_t^2)))
  return(exp(rnorm(mean=mu_time,sd=sig_time,n=1)))
}
  
a<-0.05 ### dominance factor of a plasmid encoded trait

delr1000<-0.2277206 ###relative fitness advantage of a cell that is all-beneficial plasmid (blue) versus a cell that is all-nonbeneficial plasmid (red)


cell_fit<-function(blues,reds) # function generating the fitness of cells with intermediate plasmid compositions
{
  if(reds+blues!=0)
  {
  prop<-blues/(reds+blues)
  if(prop==0)
  {return(1)}
  
  del<-(prop/(a+prop))*(1+a)*delr1000
  return(log(2)/(log(2)+del))
  }else{
    return(1000)
  }
}
  

r_effective<-1.15 ### within-cell fitness of red plasmid

phi1<-0
s1<-7
M<-40
ecl<-0.05
segr<-"bin"
init_prop<-0.5 ###initial proportion of blue plasmid in each cell

mu_ln<-log((M^2)/sqrt((M^2)+s1^2))
sig_ln<-sqrt(log(1+s1^2/(M^2))*(1-phi1^2))



oc_mat<-matrix(0,max_x,max_y) ####occupation matrix 0:free, 1:occupied, 2:blocked
oc_mat[1,]<-2
oc_mat[max_x,]<-2
oc_mat[,1]<-2
oc_mat[,max_y]<-2

nb_mat<-matrix(0,max_x,max_y)   ###neighboor matrix
nb_mat[2,]<-nb_mat[2,]+3
nb_mat[max_x-1,]<-nb_mat[max_x,]+3
nb_mat[,2]<-nb_mat[,2]+3
nb_mat[,max_y-1]<-nb_mat[,max_y]+3

b_mat<-matrix(0,max_x,max_y)   ###
r_mat<-matrix(0,max_x,max_y)   ###
S_mat<-matrix(0,max_x,max_y)   ###

time_mat<--matrix(0,max_x,max_y) ### matrix with times for next replication

####generating initializations
S_0<-log(M)
N_0<-M
for (i in 1:100000) {
  temp<-daughter0(N0=N_0[i],mu=mu_ln,phi = phi1,S_mother = S_0[i],s=sig_ln)
  N_0[i+1]<-temp[1]
  S_0[i+1]<-temp[2]
}
N_0<-N_0[10000:length(N_0)]
S_0<-S_0[10000:length(S_0)]


####initialization
j<-10
i<-10

i_s<-2:(max_x-1)
j_s<-(round(max_y/2)-1):(round(max_y/2)+1)

for (i in i_s) {
  
for (j in j_s) {
  
init_NS<-sample(length(N_0),size=1)
  
oc_mat[i,j]<-1
b_mat[i,j]<-round(N_0[init_NS]*init_prop*2)
r_mat[i,j]<-round(N_0[init_NS]*(1-init_prop)*2)
S_mat[i,j]<S_0[init_NS]

time_mat[i,j]<-n_time(mean_t=mean_time,sd_t=sd_time)

nb_mat[(i-1):(i+1),(j-1):(j+1)]<-nb_mat[(i-1):(i+1),(j-1):(j+1)]+1
nb_mat[i,j]<-nb_mat[i,j]-1
}}

image(oc_mat)





for (t in 1:t_max) {
  candidates<-which(nb_mat<8&oc_mat==1,arr.ind = T)
  ###fit_mat[candidates]
  ###picking cell with shoertest wait time to division
  mother<-candidates[which(time_mat[candidates]==min(time_mat[candidates]))[1],]
  i_m<-mother[1]
  j_m<-mother[2]
  
  ### positioning new cell
  candidates<-which(oc_mat[(i_m-1):(i_m+1),(j_m-1):(j_m+1)]==0,arr.ind=T)
  candidates[,1]<-candidates[,1]+i_m-2
  candidates[,2]<-candidates[,2]+j_m-2

  new_cell<-candidates[sample(length(candidates[,1]),size=1),]
  i<-new_cell[1]
  j<-new_cell[2]
  
  oc_mat[i,j]<-1
  
  nb_mat[(i-1):(i+1),(j-1):(j+1)]<-nb_mat[(i-1):(i+1),(j-1):(j+1)]+1
  nb_mat[i,j]<-nb_mat[i,j]-1
  
  
  ###setting plasmids for old and new cell
  S_temp<-S_mat[i_m,j_m]
  pls_temp<-c(b_mat[i_m,j_m],r_mat[i_m,j_m])
  
  temp_cell<-daughter2(pls=pls_temp,mu=mu_ln,phi = phi1,S_mother = S_temp,s=sig_ln,eclipse = ecl,r_eff = r_effective,segregation = segr)
  pls_old<-temp_cell[[1]]
  pls_new<-temp_cell[[4]]
  
  S_mat[i,j]<-temp_cell[[2]]
  S_mat[i_m,j_m]<-temp_cell[[2]]
  
  
  
  if(sum(pls_new)==0) ####cells with no plasmids die
  {
    oc_mat[i,j]<-2
  }
  
  if(sum(pls_old)==0)
  {
    oc_mat[i_m,j_m]<-2
  }
  
  ####updating plasmid compositions
  r_mat[i_m,j_m]<-pls_old[2]
  b_mat[i_m,j_m]<-pls_old[1]
  
  r_mat[i,j]<-pls_new[2]
  b_mat[i,j]<-pls_new[1]
  
  ### setting up next replication times
  mean_time<-cell_fit(b_mat[i,j],r_mat[i,j])
  time_mat[i,j]<-time_mat[i_m,j_m]+n_time(mean_t=mean_time,sd_t=sd_time)
  mean_time<-cell_fit(b_mat[i_m,j_m],r_mat[i_m,j_m])
  time_mat[i_m,j_m]<-time_mat[i_m,j_m]+n_time(mean_t=mean_time,sd_t=sd_time)
  
}

image(oc_mat)
image(nb_mat)


image(b_mat/(r_mat+b_mat))
image(r_mat/(r_mat+b_mat))
