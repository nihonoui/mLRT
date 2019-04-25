# simulation at a CpG site
rep=1000
N_sub_A=2;N_sub_B=2 # number of replicates for conditions A and B

N_sub=c(N_sub_A,N_sub_B) 
N_tab=c(N_tab_A,N_tab_B)  # N_tab_A and N_tab_B are coverages from TAB_seq at a CpG site
N_wgbs=c(N_wgbs_A,N_wgbs_B)  # N_wgbs_A and N_wgbs_B are coverages from WGBS at a CpG site


# prop=c(p1A,p1B,p2A,p2B) denotes proportions of hmC and mC under two conditions
# p1A and p1B are proportions of hmC under condition A and B, respectively.
# p2A and p2B are proportions of mC under condition A and B, respectively.
# The function create.dat generates simulated data with given proportions of hmC and mC under two conditions
create.dat = function(N_sub,N_tab,N_wgbs,prop) {
  if (length(prop)==4) {
    p1A=prop[1];p1B=prop[2];p2A=prop[3];p2B=prop[4] 
  } else {
    p1A=p1B=prop[1];p2A=p2B=prop[2]
  }
  N_sub_A=N_sub[1];N_sub_B=N_sub[2]
  N_tab_A=N_tab[1:N_sub_A];N_tab_B=N_tab[(N_sub_A+1):(N_sub_A+N_sub_B)]
  N_wgbs_A=N_wgbs[1:N_sub_A];N_wgbs_B=N_wgbs[(N_sub_A+1):(N_sub_A+N_sub_B)]
  NC_tab_A=NC_wgbs_A=NC_tab_B=NC_wgbs_B=c()
  for (i in 1:N_sub_A) {
    NC_tab_A[i]=rbinom(1,N_tab_A[i],p1A)
    NC_wgbs_A[i]=rbinom(1,N_wgbs_A[i],p1A+p2A)
  }
  for (i in 1:N_sub_B) {
    NC_tab_B[i]=rbinom(1,N_tab_B[i],p1B)
    NC_wgbs_B[i]=rbinom(1,N_wgbs_B[i],p1B+p2B)
  }
  
  type=rep(c('A','B'),N_sub)
  N_tab=c(N_tab_A,N_tab_B);N_wgbs=c(N_wgbs_A,N_wgbs_B)
  NC_tab=c(NC_tab_A,NC_tab_B);NC_wgbs=c(NC_wgbs_A,NC_wgbs_B)
  dat=data.frame(N_tab,NC_tab,N_wgbs,NC_wgbs,type)
  colnames(dat)=c('N_tab','NC_tab','N_wgbs','NC_wgbs','type')
  return(dat)
}

# mLRT: our likelihood ratio test for jointly testing the proportions of 5hmC and 5mC
mLRT = function(dat) {
  N_tab_A=dat$N_tab[dat$type=='A'];NC_tab_A=dat$NC_tab[dat$type=='A']
  N_tab_B=dat$N_tab[dat$type=='B'];NC_tab_B=dat$NC_tab[dat$type=='B']
  N_wgbs_A=dat$N_wgbs[dat$type=='A'];NC_wgbs_A=dat$NC_wgbs[dat$type=='A']
  N_wgbs_B=dat$N_wgbs[dat$type=='B'];NC_wgbs_B=dat$NC_wgbs[dat$type=='B']
  N_sub_A=length(N_tab_A);N_sub_B=length(N_tab_B)
  
  p12_H0=sum(c(NC_wgbs_A,NC_wgbs_B))/sum(c(N_wgbs_A,N_wgbs_B))
  p1_H0=sum(c(NC_tab_A,NC_tab_B))/sum(c(N_tab_A,N_tab_B))
  p12A_H1=sum(NC_wgbs_A)/sum(N_wgbs_A)
  p12B_H1=sum(NC_wgbs_B)/sum(N_wgbs_B)
  p1A_H1=sum(NC_tab_A)/sum(N_tab_A)
  p1B_H1=sum(NC_tab_B)/sum(N_tab_B)
  
  if (p1_H0>p12_H0) {
    p12_H0=p1_H0=sum(c(dat$NC_tab,dat$NC_wgbs))/sum(c(dat$N_tab,dat$N_wgbs))
  }
  
  if (p1A_H1>p12A_H1) {
    p1A_H1=p12A_H1=sum(c(NC_wgbs_A,NC_tab_A))/sum(c(N_wgbs_A,N_tab_A))
  }
  
  if (p1B_H1>p12B_H1) {
    p1B_H1=p12B_H1=sum(c(NC_wgbs_B,NC_tab_B))/sum(c(N_wgbs_B,N_tab_B))
  }
  
  L0A=L1A=L0B=L1B=c()
  for (i in 1:N_sub_A) {
    L0A[i]=p12_H0^NC_wgbs_A[i]*(1-p12_H0)^(N_wgbs_A[i]-NC_wgbs_A[i])
    L0A[i]=L0A[i]*p1_H0^NC_tab_A[i]*(1-p1_H0)^(N_tab_A[i]-NC_tab_A[i])
    L1A[i]=p12A_H1^NC_wgbs_A[i]*(1-p12A_H1)^(N_wgbs_A[i]-NC_wgbs_A[i])
    L1A[i]=L1A[i]*p1A_H1^NC_tab_A[i]*(1-p1A_H1)^(N_tab_A[i]-NC_tab_A[i])
  }
  for (j in 1:N_sub_B) {
    L0B[j]=p12_H0^NC_wgbs_B[j]*(1-p12_H0)^(N_wgbs_B[j]-NC_wgbs_B[j])
    L0B[j]=L0B[j]*p1_H0^NC_tab_B[j]*(1-p1_H0)^(N_tab_B[j]-NC_tab_B[j])
    L1B[j]=p12B_H1^NC_wgbs_B[j]*(1-p12B_H1)^(N_wgbs_B[j]-NC_wgbs_B[j])
    L1B[j]=L1B[j]*p1B_H1^NC_tab_B[j]*(1-p1B_H1)^(N_tab_B[j]-NC_tab_B[j])
  }
  
  L0=c(L0A,L0B); L1=c(L1A,L1B)
  
  # test statistic
  M=2*sum(log(L1))-2*sum(log(L0))
  pv=1-pchisq(M,df=2)  #p-value
  return(pv)
}

# The function mynaive returns results of nFET and nCHI tests
mynaive=function(dat) {
  x1=aggregate(round(dat$NC_wgbs*dat$NC_tab/dat$N_tab), 
               by=list(Category=dat$type), FUN=sum) #hmC
  tmp=ifelse(dat$NC_wgbs-round(dat$NC_wgbs*dat$NC_tab/dat$N_tab)<0,0,
             dat$NC_wgbs-round(dat$NC_wgbs*dat$NC_tab/dat$N_tab))
  tmp=data.frame(tmp,dat$type)
  colnames(tmp)=c('x','type')
  x2=aggregate(tmp$x,by=list(Category=tmp$type), FUN=sum) #mC
  x3=aggregate(dat$N_wgbs-dat$NC_wgbs, 
               by=list(Category=dat$type), FUN=sum) #C
  tab=as.matrix(cbind(x1$x,x2$x,x3$x))
  #print(tab)
  pv1=fisher.test(tab)$p.value #nFET
  pv2=chisq.test(tab)$p.value  #nCHI
  return(c(pv1,pv2))
}


# The function myMLE returns results of mFET and mCHI tests
myMLE=function(dat) {
  N_tab_A=dat$N_tab[dat$type=='A'];NC_tab_A=dat$NC_tab[dat$type=='A']
  N_tab_B=dat$N_tab[dat$type=='B'];NC_tab_B=dat$NC_tab[dat$type=='B']
  N_wgbs_A=dat$N_wgbs[dat$type=='A'];NC_wgbs_A=dat$NC_wgbs[dat$type=='A']
  N_wgbs_B=dat$N_wgbs[dat$type=='B'];NC_wgbs_B=dat$NC_wgbs[dat$type=='B']
  levelA_hmC=ifelse(sum(NC_wgbs_A)/sum(N_wgbs_A)-sum(NC_tab_A)/sum(N_tab_A)<0,
                    sum(c(NC_wgbs_A,NC_tab_A))/sum(c(N_wgbs_A,N_tab_A)),
                    sum(NC_tab_A)/sum(N_tab_A))
  levelA_mC=max(sum(NC_wgbs_A)/sum(N_wgbs_A)-sum(NC_tab_A)/sum(N_tab_A),0)
  levelB_hmC=ifelse(sum(NC_wgbs_B)/sum(N_wgbs_B)-sum(NC_tab_B)/sum(N_tab_B)<0,
                    sum(c(NC_wgbs_B,NC_tab_B))/sum(c(N_wgbs_B,N_tab_B)),
                    sum(NC_tab_B)/sum(N_tab_B))
  levelB_mC=max(sum(NC_wgbs_B)/sum(N_wgbs_B)-sum(NC_tab_B)/sum(N_tab_B),0)
  N.A=round(apply(cbind(N_wgbs_A,N_tab_A),1,mean));N.B=round(apply(cbind(N_wgbs_B,N_tab_B),1,mean))
  #N.A=N_wgbs_A;N.B=N_wgbs_B
  #N.A=apply(cbind(N_wgbs_A,N_tab_A),1,max);N.B=apply(cbind(N_wgbs_B,N_tab_B),1,max)
  N.hmC=floor(c(levelA_hmC*N.A,levelB_hmC*N.B))
  N.mC=floor(c(levelA_mC*N.A,levelB_mC*N.B))
  N.C=c(N.A,N.B)-N.hmC-N.mC
  type=rep(c('A','B'),c(length(N.A),length(N.B)))
  dat1=data.frame( N.hmC, N.mC,N.C,type)
  colnames(dat1)=c('hmC','mC','C','type')
  x1=aggregate(dat1$hmC, by=list(Category=dat1$type), FUN=sum) # hmC
  x2=aggregate(dat1$mC, by=list(Category=dat1$type), FUN=sum) # mC
  x3=aggregate(dat1$C,by=list(Category=dat1$type), FUN=sum) # unmthylated C
  tab=as.matrix(cbind(x1$x,x2$x,x3$x))
  pv1=fisher.test(tab)$p.value #mFET
  pv2=chisq.test(tab)$p.value  #mCHI
  return(c(pv1,pv2))
}


p12=0.8 # p12=p1+p2 sum of proportions of hmC and mC
ncol=10
prop=matrix(NA,nrow=ncol+1,ncol=4)
pv_mLRT=pv_mFET=pv_nFET.naive=pv_FET=pv_mCHI=pv_nCHI=matrix(NA,nrow=rep,ncol=ncol+1)
for (j in 0:ncol) {
  prop[j+1,1]=0.1
  prop[j+1,2]=0.1+0.02*j
  prop[j+1,3]=p12-prop[j+1,1]
  prop[j+1,4]=p12-prop[j+1,2]
  for (i in 1:rep) {
    dat1=create.dat(N_sub,N_tab,N_wgbs,prop[j+1,])
    pv_mLRT[i,j+1]=mytest(dat1)
    re.MLE=myMLE(dat1)
    re.naive=mynaive(dat1)
    pv_nFET[i,j+1]=re.naive[1]
    pv_nCHI[i,j+1]=re.naive[2]
    pv_mFET[i,j+1]=re.MLE[1]
    pv_mCHI[i,j+1]=re.MLE[2]
    x1=aggregate(dat1$N_wgbs, by=list(Category=dat1$type), FUN=sum)
    x2=aggregate(dat1$NC_wgbs, by=list(Category=dat1$type), FUN=sum)
    tab=as.matrix(cbind(x1$x-x2$x,x2$x))
    pv_FET[i,j+1]=fisher.test(tab)$p.value
  }
}


sim.power_mLRT=apply(pv_mLRT<0.05,2,sum)/rep
sim.power_mFET=apply(pv_mFET<0.05,2,sum)/rep
sim.power_nFET=apply(pv_nFET<0.05,2,sum)/rep
sim.power_FET=apply(pv_FET<0.05,2,sum)/rep
#Pearson’s chi-squared test is not applicable when count values in a column of the contingency table are all zero. 
sim.power_mCHI=apply(pv_mCHI<0.05,2,sum)/apply(pv_mCHI,2,length) 
sim.power_nCHI=apply(pv_nCHI<0.05,2,sum)/rep