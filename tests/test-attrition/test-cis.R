
library(testthat)
library(attrition)

n1_c_s <- 100   # sum(R1==1 & Z==0)
n1_t_s <- 100   # sum(R1==1 & Z==1)
n1_c <- 400     # sum(Z==0) # clean up
n1_t <- 400     # sum(Z==1)

p1_c <- n1_c_s/n1_c
p1_t <- n1_t_s/n1_t

y1m_c <- 0 # mean(Y[R1==1 & Z==0])
y1m_t <- 0 # mean(Y[R1==1 & Z==1])

n2_c <- 5 # sum(Attempt ==1 & Z ==0)
n2_t <- 5 # sum(Attempt ==1 & Z ==1)

# This is the thing to vary
p2_c <- 0.05 # sum(R2==1 & Z==0)/n2_c
p2_t <- 0.05 # sum(R2==1 & Z==1)/n2_t

y2m_nm_c <- 0 # mean(Y[R2==1 & Z==0])
y2m_nm_t <- 0 # mean(Y[R2==1 & Z==1])

s1_c <- 1 # sd(Y[R1==1 & Z==0])
s1_t <- 1 # sd(Y[R1==1 & Z==1])
s2_nm_c <- 1 # sd(Y[R2==1 & Z==0])
s2_nm_t <- 1 # sd(Y[R2==1 & Z==1])

cis_out_m <-
manski_cis(n1_t = n1_t, n1_c = n1_c,
           n1_t_s = n1_t_s, n1_c_s = n1_c_s,
           p1_t = p1_t, p1_c = p1_c,
           y1m_t = y1m_t, y1m_c = y1m_c,
           s1_t = s1_t, s1_c = s1_c,
           minY = 0,maxY = 1,alpha = 0.05)

p2_c <- p2_t <- 1

cis_out_ds <-
  ds_manski_cis_2s(
    n1_t=n1_t,n2_t=n2_t,
    n1_c=n1_c,n2_c=n2_c,
    p1_t=p1_t,p2_t=p2_t,
    s1_t=s1_t,s1_c=s1_c,
    s2_nm_t=s2_nm_t,
    s2_nm_c=s2_nm_c,
    y1m_t=y1m_t,y1m_c=y1m_c,
    y2m_nm_t=y2m_nm_t,
    y2m_nm_c=y2m_nm_c,
    c1a_t=c1a_t,c1r_t=c1r_t,c2a_t=c2a_t,c2r_t=c2r_t,
    p1_c=p1_c,p2_c=p2_c,
    minY=0,maxY=1,alpha=0.05)


cis_out_m
cis_out_ds
