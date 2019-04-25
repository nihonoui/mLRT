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
  
  L0=c(L0A,L0B); L1=c(L1A,L1B) #likelihood function
  
  # test statistic
  M=2*sum(log(L1))-2*sum(log(L0))
  pv=1-pchisq(M,df=2)  #p-value
  return(pv)
}
