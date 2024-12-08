rm(list = ls(all = TRUE)) # Clear workspace

gg=function(l){
  
# Load necessary package
library("VGAM")

# Function to generate sample data
samp <- function(n, mu, sgm, lmd) {
  TT <- abs(rnorm(n, mu, sgm));TT # Time-to-event
  TT[TT < 0] <- 0
  
  S <- runif(n, 0, 5);S           # Monitoring time
  cc <- ifelse(TT < S, 1, 0)      # Event Status
  
  ss <- S[cc == 1]
  tt <- TT[cc == 1]
  
  p1 <- exp(-(ss - tt) * lmd)     # Recall probability
  zz <- rbinom(length(tt), 1, p1)
  
  t1 <- tt[zz == 1]
  s1 <- ss[zz == 1]
  s2 <- ss[zz == 0]
  s3 <- S[cc == 0]
  
  return(list(t1 = t1, s1 = s1, s2 = s2, s3 = s3))
}

# Set parameters and generate sample data
n <- 50; mu <- 2.2; sgm <- 0.9; lmd <- 0.50

ss <- samp(n, mu, sgm, lmd)

z1 <- unlist(ss$t1)
S1 <- unlist(ss$s1)
S2 <- unlist(ss$s2)
S3 <- unlist(ss$s3)

n1 <- length(z1);n1
n2 <- length(S2);n2
n3 <- length(S3);n3


#===================================================================#
#                                                                   #
#                     ML Estimate (Using Optim)                     #------------
#                                                                   #
#===================================================================#  

# Likelihood function
ll <- function(par) {
  if (par[1] > 0 & par[2] > 0 & par[3] > 0) {
    m.th1 <- par[1]
    m.th2 <- par[2]
    m.th3 <- par[3]
    
    # Likelihood components
    l1 <- (-n2) * log(m.th2) - (0.5 / m.th2^2) * sum((z1 - m.th1)^2) - m.th3 * sum(S1 - z1)
    
    integ_d1 <- function(u, s, m.th2, m.th3) {
      exp(-(0.5 / m.th2^2) * (s - u)^2) * (1 - exp(-m.th3 * (s - u)))
    }
    
    integd1 <- sapply(S2, function(s) {
      tryCatch({
        integrate(integ_d1, lower = 0, upper = s, s = s, m.th2 = m.th2, m.th3 = m.th3)$value
      }, error = function(e) Inf)  # Handle non-finite cases
    })
    
    l2 <- sum(log(integd1))
    l3 <- sum(1 - pnorm((S3 - m.th1) / m.th2))
    
    kk <- l1 + l2 + l3
  } else {
    kk <- -Inf
  }
  
  return(-kk)
}

# Empirical initial values
mu_init <- mean(z1)
sgm_init <- sd(z1)
lmd_init <- 1 / mean(S1-z1)

par <- c(mu_init, sgm_init, lmd_init)

# Optimization using L-BFGS-B with constraints
nl <- optim(
  par, ll, method = "L-BFGS-B",
  lower = c(0.1, 0.1, 0.1), 
  upper = c(3, 2, 1),
  control = list(maxit = 10000, trace = 1),
  hessian = TRUE
)

# Results
#ml_opt<-nl$par;ml_opt  # Estimated parameters
ml_opt<-c(nl$par[1],sqrt(nl$par[2]),nl$par[3])
#nl$value  # Minimum value of the negative log-likelihood

mse_opt=c((ml_opt[1]-mu)^2,(ml_opt[2]-sgm)^2,(ml_opt[3]-lmd)^2);mse_opt
bias_opt=c((ml_opt[1]-mu),(ml_opt[2]-sgm),(ml_opt[3]-lmd));bias_opt

#variance-covariance    
var_cov=solve(nl$hessian);var_cov

ci_l=ci_u=rep();alpha=0.05
for (i in 1:3) {
  ci_l[i]=ml_opt[i]-qnorm(1-(alpha/2))*sqrt((diag(var_cov)[i]))
  ci_u[i]=ml_opt[i]+qnorm(1-(alpha/2))*sqrt((diag(var_cov)[i]))
}

#limits
#lmt_opt<-c(ci_l[1],ci_u[1],ci_l[2],ci_u[2],ci_l[3],ci_u[3]);lmt_opt

#length    
l1=c(ci_u[1]-ci_l[1]);l1
l2=c(ci_u[2]-ci_l[2]);l2
l3=c(ci_u[3]-ci_l[3]);l3
leng_opt=c(l1,l2,l3);leng_opt

#coverage probability
cp1=ifelse(ci_l[1]<mu && mu<ci_u[1],1,0);cp1
cp2=ifelse(ci_l[2]<sgm && sgm<ci_u[2],1,0);cp2
cp3=ifelse(ci_l[3]<lmd && lmd<ci_u[3],1,0);cp3
cp_opt=c(cp1,cp2,cp3);cp_opt

NR<-c(ml_opt,mse_opt,bias_opt,leng_opt,cp_opt);NR


#===================================================================#
#                                                                   #
#                        ML Estimate (Using EM)                     #------------
#                                                                   #
#===================================================================#  

# Empirical initialization
m.th1 <- mean(z1)                   # Mean of observed events
m.th2 <- sd(z1)                     # Standard deviation of observed events
m.th3 <- 1 / mean(S1 - z1)          # Approximate initial value for lambda

it=50
# Vectors to store parameter estimates
m.thh1 <- m.thh2 <- m.thh3 <- numeric(it)  
#=========loop for alpha and beta====================#

for(j in 1:it){
  
  num1=dnorm((S2-m.th1)/m.th2)-dnorm(-m.th1/m.th2);num1 #left-censored
  den1=pnorm((S2-m.th1)/m.th2)-pnorm(-m.th1/m.th2);den1
  
  num2=dnorm((S3-m.th1)/m.th2);num2 #right
  den2=1-pnorm((S3-m.th1)/m.th2);den2
  
  #--------------------------------------------------#
  xi1=m.th1-m.th2*(num1/den1);xi1
  xi2=m.th1^2+m.th2^2-m.th2*(m.th1+S2)*(num1/den1);xi2
  xi3=m.th1+m.th2*(num2/den2);xi3
  xi4=m.th1^2+m.th2^2+m.th1*(m.th2+S3)*(num2/den2);xi4
  
  #M step of E-M 
  m.thh1[j]=m.th1;m.thh2[j]=m.th2
  
  m.thh1[j+1]=(1/n)*(sum(z1)+sum(xi1)+sum(xi3))
  m.thh2[j+1]=sqrt((1/n)*(sum(z1^2)+sum(xi2)+sum(xi4)-2*m.thh1[j+1]*(sum(z1)+sum(xi1)+sum(xi3))+n*m.thh1[j+1]^2))
  
  m.th1=m.thh1[j+1];m.th2=m.thh2[j+1] #Storing the final value of estimates
}
m.thh1;m.thh2

# par(mfrow=c(1,3))   #checking the convergence through plots
# plot(m.thh1,type='l');plot(m.thh2,type='l')

#=============================================#
#loop for lambda
for(j in 1:it){
  
  num11=dnorm((S2-m.thh1[it])/m.thh2[it])-dnorm(-m.thh1[it]/m.thh2[it]);num11 #left-censored
  den11=pnorm((S2-m.thh1[it])/m.thh2[it])-pnorm(-m.thh1[it]/m.thh2[it]);den11
  
  xi0=m.thh1[it]-m.thh2[it]*(num11/den11);xi0
  cc1=m.th3*(S2-xi0);cc1
  
  m.thh3[j]=m.th3 
  xi5=(1-((1+cc1)*exp(-cc1)))/(m.thh3[j]*(1-exp(-cc1)));xi5
  m.thh3[j+1]=n2/(sum(S1-z1)+sum(xi5))
  m.th3=m.thh3[j+1]
}
m.thh3
#plot(m.thh3,type='l')

# par(mfrow=c(2,2))   #checking the convergence through plots
# plot(m.thh1,type='l');plot(m.thh2,type='l'); plot(m.thh3,type='l')

ml_em=c(m.thh1[it],m.thh2[it],m.thh3[it]);ml_em  #store the final estimates
mse_em=c((ml_em[1]-mu)^2,(ml_em[2]-sgm)^2,(ml_em[3]-lmd)^2);mse_em #calculate the Mean Square Error
bias_em=c((ml_em[1]-mu),abs(ml_em[2]-sgm),abs(ml_em[3]-lmd));bias_em #calculate the Bias

EM<-c(ml_em,mse_em,bias_em);EM


#===================================================================#
#                                                                   #
#                        Bayesian Estimation                        #------------------
#                                                                   #
#===================================================================#
# Load necessary package
# install.packages("LaplacesDemon") # Uncomment to install if necessary
library(LaplacesDemon);library(coda)

# Prior parameters for Set 1
mu0 <- 2.00        # Prior mean for mu
sig0 <- 1.00       # Prior mean for sigma^2
nu0 <- 4.0         # Prior variance for mu
k0 <- 1.0          # Shape parameter for inverse chi-squared prior
l0 <- 1.0          # Shape parameter for lambda
m0 <- 1.0          # Shape parameter for lambda

# For Set 2, Set 3, or Set 4, you can adjust these values as needed
# Example for Set 2:
# mu0 <- 2.20
# sig0 <- 0.64
# nu0 <- 4.0
# k0 <- 1.0
# l0 <- 1.0
# m0 <- 1.0

# Number of iterations for Gibbs sampling
it1 <- 80000

# Define storage for parameter samples
mu.th <- numeric(it1)
sgm.th <- numeric(it1)
lmd.th <- numeric(it1)

# Initial values for Gibbs sampler
mu.th[1] <- ml_em[1]    # Initialize mu with ML estimate
sgm.th[1] <- ml_em[2]   # Initialize sigma with ML estimate
lmd.th[1] <- ml_em[3]   # Initialize lambda with ML estimate

## Gibbs Sampling Loop
for (j in 1:(it1 - 1)) {
  
  # Sampling t_li (left censored times) for observed censoring points (S2)
  tt1 <- qnorm(
    pnorm(-mu.th[j] / sgm.th[j]) + runif(n2, 0, 1) * (pnorm((S2 - mu.th[j]) / sgm.th[j]) - pnorm(-mu.th[j] / sgm.th[j])),
    mean = mu.th[j],
    sd = sgm.th[j]
  )
  
  # Sampling t_ri (right censored times) for censored events (S3)
  tt2 <- qnorm(
    pnorm((S3 - mu.th[j]) / sgm.th[j]) + runif(n3, 0, 1) * (1 - pnorm((S3 - mu.th[j]) / sgm.th[j])),
    mean = mu.th[j],
    sd = sgm.th[j]
  )
  
  # Sampling w_i
  if (lmd.th[j] > 0) {
    w_i <- -(1 / lmd.th[j]) * log(1 - runif(n2, 0, 1) * (1 - exp(-lmd.th[j] * (S2 - tt1))))
  } else {
    w_i <- rep(0, n2)  # Default if lambda is invalid
  }
  
  # Update for mu (posterior mean and variance)
  kn <- n1 + n2 + n3 + k0
  mu_n <- (sum(z1) + sum(tt1) + sum(tt2) + (k0 * mu0)) / kn
  mu.th[j + 1] <- rnorm(1, mean = mu_n, sd = sgm.th[j] / sqrt(kn))
  
  # Update for sigma (posterior of sigma using inverse chi-squared)
  ss <- sum((z1 - mu.th[j + 1])^2) + sum((tt1 - mu.th[j + 1])^2) + sum((tt2 - mu.th[j + 1])^2)
  scale_post <- (nu0 * sig0^2 + ss) / (nu0 + n1 + n2 + n3)
  sgm.th[j + 1] <- sqrt(rinvchisq(1, df = nu0 + n1 + n2 + n3, scale = scale_post))
  
  # Update for lambda (posterior for lambda using gamma distribution)
  lmd.th[j + 1] <- rgamma(1, shape = n2 + l0, rate = sum(S1 - z1) + sum(w_i) + m0)
}
head(mu.th);head(sgm.th);head(lmd.th)


# After the Gibbs sampling loop, you can inspect the posterior samples and perform analysis
# Trace-plots for convergence inspection
# par(mfrow = c(1, 3))
# plot(mu.th, type = 'l', main = "", ylab = expression(mu))
# plot(sgm.th, type = 'l', main = "", ylab = expression(sigma))
# plot(lmd.th, type = 'l', main = "", ylab = expression(lambda))
# 
# # Histograms for posterior distributions
# # par(mfrow = c(1, 3))
# hist(mu.th, breaks = 30, main = "", xlab = expression(mu))
# hist(sgm.th, breaks = 30, main = "", xlab = expression(sigma))
# hist(lmd.th, breaks = 30, main = "", xlab = expression(lambda))
# 
# 
# # Plot ACF for mu, sigma, and lambda
# # par(mfrow = c(1, 3))
# acf(mu.th, main = "", lag.max = 50, ylab = expression("ACF"(mu)))
# acf(sgm.th, main = "", lag.max = 50, ylab = expression("ACF"(sigma)))
# acf(lmd.th, main = "", lag.max = 50, ylab = expression("ACF"(lambda)))

# Burn-in and thinning
ch1 = mu.th[1001:it1]
ch2 = sgm.th[1001:it1]
ch3 = lmd.th[1001:it1]
zz1 = seq(15, length(ch1), 15)  # Thinning the chain
ch11 = ch1[zz1]
ch22 = ch2[zz1]
ch33 = ch3[zz1]

#===========Bayes estimates under SELF--------------
mu.tth=sgm.tth=lmd.tth=rep()
mu.tth=mean(ch11);sgm.tth=mean(ch22);lmd.tth=mean(ch33);
b.th=c(mu.tth,sgm.tth,lmd.tth);b.th 
self_estimates<-b.th
mse_self=c((b.th[1]-mu)^2,(b.th[2]-sgm)^2,(b.th[3]-lmd)^2);mse_self
bias_self=c((b.th[1]-mu),(b.th[2]-sgm),(b.th[3]-lmd));bias_self


#==================LINEX Loss functions------------
# LINEX Loss (example with 'a' = 0.5 for each parameter)
#1. If overestimation is costly (e.g., in budgeting or resource allocation scenarios where over-forecasting is risky), choose a positive 
#2. If underestimation is risky (e.g., underestimating disease risk in public health), a negative may be better.

#Case 1.
a = 1  # Modify 'a' as needed for asymmetry
mu_LINEX = -1/a * log(mean(exp(-a * ch11)))
sgm_LINEX = -1/a * log(mean(exp(-a * ch22)))
lmd_LINEX = -1/a * log(mean(exp(-a * ch33)))

linex_estimates = c(mu_LINEX, sgm_LINEX, lmd_LINEX)
linex_estimates

mse_linex=c((linex_estimates[1]-mu)^2,(linex_estimates[2]-sgm)^2,(linex_estimates[3]-lmd)^2);mse_linex
bias_linex=c((linex_estimates[1]-mu),(linex_estimates[2]-sgm),(linex_estimates[3]-lmd));bias_linex

#Case 2. 
a1 = -1  # Modify 'a' as needed for asymmetry
mu_LINEX1 = -1/a1 * log(mean(exp(-a1 * ch11)))
sgm_LINEX1 = -1/a1 * log(mean(exp(-a1 * ch22)))
lmd_LINEX1 = -1/a1 * log(mean(exp(-a1 * ch33)))

linex_estimates1 = c(mu_LINEX1, sgm_LINEX1, lmd_LINEX1)
linex_estimates1

mse_linex1=c((linex_estimates1[1]-mu)^2,(linex_estimates1[2]-sgm)^2,(linex_estimates1[3]-lmd)^2);mse_linex1
bias_linex1=c((linex_estimates1[1]-mu),(linex_estimates1[2]-sgm),(linex_estimates1[3]-lmd));bias_linex1


#===============HPD Intervals-------------------
h.th1=HPDinterval(mcmc(ch11));h.th2=HPDinterval(mcmc(ch22));h.th3=HPDinterval(mcmc(ch33))
#HPD_lim=c(h.th1[,1],h.th1[,2],h.th2[,1],h.th2[,2],h.th3[,1],h.th3[,2]);HPD_lim

h.th1_l=h.th1[,2]-h.th1[,1];h.th2_l=h.th2[,2]-h.th2[,1];h.th3_l=h.th3[,2]-h.th3[,1]
HPD_L=c(h.th1_l,h.th2_l,h.th3_l);HPD_L                                   #length of HPD

coverage_HPD1=ifelse(mu>h.th1[,1] && mu<h.th1[,2],1,0)
coverage_HPD2=ifelse(sgm>h.th2[,1] && sgm<h.th2[,2],1,0)                  #HPD interval Estimate 
coverage_HPD3=ifelse(lmd>h.th3[,1] && lmd<h.th3[,2],1,0)
coverage_HPD=c(coverage_HPD1,coverage_HPD2,coverage_HPD3);coverage_HPD    #coverage HPD

Bayes=c(self_estimates,mse_self,bias_self,
        linex_estimates,mse_linex,bias_linex,
        linex_estimates1,mse_linex1, bias_linex1,
        HPD_L,coverage_HPD)


return(c(NR,EM,Bayes))

}
#gg()

#Parallel Computation of Simulation 
library(parallel)

#Determine number of cores to use
num_cores <- detectCores() - 1

#Create a cluster
cl <- makeCluster(num_cores)

start=Sys.time()
#Run the simulation in parallel
l=50# number simulations to run (say; 500, 1000, 2000)
results <- parLapply(cl, 1:l, gg)

#Stop the cluster
stopCluster(cl)
end=Sys.time()
end-start
mm=unlist(results)
rr1=matrix(mm,nrow = l,ncol = 57,byrow = T) # here ncol will be equal to number of values in return 
rr2=round(colMeans(rr1),2);rr2 #estimates based on desired number of iterations

names(rr2)<-c(   
                 #NR  
                 "nr_mu", "nr_sgm", "nr_lmd", 
                 "nr_mu_mse", "nr_sgm_mse", "nr_lmd_mse",
                 "nr_mu_bias", "nr_sgm_bias", "nr_lmd_bias",
                 "nr_mu_al", "nr_sgm_al", "nr_lmd_al",
                 "nr_mu_cp", "nr_sgm_cp", "nr_lmd_cp",
                 
                 #EM
                 "em_mu", "em_sgm", "em_lmd", 
                 "em_mu_mse", "em_sgm_mse", "em_lmd_mse",
                 "em_mu_bias", "em_sgm_bias", "em_lmd_bias",
                 
                 #SELF
                 "self_mu", "self_sgm", "self_lmd", 
                 "self_mu_mse", "self_sgm_mse", "self_lmd_mse",
                 "self_mu_bias", "self_sgm_bias", "self_lmd_bias", 
                 
                 #linex with a=1
                 "linex_mu", "linex_sgm", "linex_lmd",
                 "linex_mu_mse", "linex_sgm_mse", "linex_lmd_mse", 
                 "linex_mu_bias", "linex_sgm_bias", "linex_lmd_bias", 
                 
                 #linex with a=-1
                 "linex_mu1", "linex_sgm1", "linex_lmd1", 
                 "linex_mu_mse1", "linex_sgm_mse1", "linex_lmd_mse1",
                 "linex_mu_bias1", "linex_sgm_bias1", "linex_lmd_bias1", 
                 
                 #HPd
                 "b_mu_al", "sgm_al", "b_lmd_al",
                 "b_mu_cp", "b_sgm_cp", "b_lmd_cp")


library(gt)

# Assuming rr2 is already defined with the names as provided
# Prepare a data frame that includes all the relevant parameters
data <- data.frame(
  Parameter = c("mu", "sgm", "lmd"),
  
  # NR
  Estimate_NR = c(rr2["nr_mu"], rr2["nr_sgm"], rr2["nr_lmd"]),
  MSE_NR = c(rr2["nr_mu_mse"], rr2["nr_sgm_mse"], rr2["nr_lmd_mse"]),
  Bias_NR = c(rr2["nr_mu_bias"], rr2["nr_sgm_bias"], rr2["nr_lmd_bias"]),
  AL_NR = c(rr2["nr_mu_al"], rr2["nr_sgm_al"], rr2["nr_lmd_al"]),
  CP_NR = c(rr2["nr_mu_cp"], rr2["nr_sgm_cp"], rr2["nr_lmd_cp"]),
  
  # EM
  Estimate_EM = c(rr2["em_mu"], rr2["em_sgm"], rr2["em_lmd"]),
  MSE_EM = c(rr2["em_mu_mse"], rr2["em_sgm_mse"], rr2["em_lmd_mse"]),
  Bias_EM = c(rr2["em_mu_bias"], rr2["em_sgm_bias"], rr2["em_lmd_bias"]),
  
  # SELF
  Estimate_SELF = c(rr2["self_mu"], rr2["self_sgm"], rr2["self_lmd"]),
  MSE_SELF = c(rr2["self_mu_mse"], rr2["self_sgm_mse"], rr2["self_lmd_mse"]),
  Bias_SELF = c(rr2["self_mu_bias"], rr2["self_sgm_bias"], rr2["self_lmd_bias"]),
  
  # Linex (a=1)
  Estimate_Linex1 = c(rr2["linex_mu"], rr2["linex_sgm"], rr2["linex_lmd"]),
  MSE_Linex1 = c(rr2["linex_mu_mse"], rr2["linex_sgm_mse"], rr2["linex_lmd_mse"]),
  Bias_Linex1 = c(rr2["linex_mu_bias"], rr2["linex_sgm_bias"], rr2["linex_lmd_bias"]),
  
  # Linex (a=-1)
  Estimate_Linex1_neg = c(rr2["linex_mu1"], rr2["linex_sgm1"], rr2["linex_lmd1"]),
  MSE_Linex1_neg = c(rr2["linex_mu_mse1"], rr2["linex_sgm_mse1"], rr2["linex_lmd_mse1"]),
  Bias_Linex1_neg = c(rr2["linex_mu_bias1"], rr2["linex_sgm_bias1"], rr2["linex_lmd_bias1"]),
  
  # HPD
  AL_Bayes = c(rr2["b_mu_al"], rr2["sgm_al"], rr2["b_lmd_al"]),
  CP_Bayes = c(rr2["b_mu_cp"], rr2["b_sgm_cp"], rr2["b_lmd_cp"])
)

# Create the gt table
gt_table <- data %>%
  gt() %>%
  tab_header(
    title = "Comparison of NR, EM, SELF, Linex, and Bayes Methods"
  ) %>%
  cols_label(
    Parameter = "Parameter",
    
    Estimate_NR = md("*Estimate*"),
    MSE_NR = md("*MSE*"),
    Bias_NR = md("*Bias*"),
    AL_NR = md("*AL*"),
    CP_NR = md("*CP*"),
    
    Estimate_EM = md("*Estimate*"),
    MSE_EM = md("*MSE*"),
    Bias_EM = md("*Bias*"),
    
    Estimate_SELF = md("*Estimate*"),
    MSE_SELF = md("*MSE*"),
    Bias_SELF = md("*Bias*"),
    
    Estimate_Linex1 = md("*Estimate*"),
    MSE_Linex1 = md("*MSE*"),
    Bias_Linex1 = md("*Bias*"),
    
    Estimate_Linex1_neg = md("*Estimate*"),
    MSE_Linex1_neg = md("*MSE*"),
    Bias_Linex1_neg = md("*Bias*"),
    
    AL_Bayes = md("*AL*"),
    CP_Bayes = md("*CP*")
  ) %>%
  tab_spanner(
    label = "NR",
    columns = vars(Estimate_NR, MSE_NR, Bias_NR, AL_NR, CP_NR)
  ) %>%
  tab_spanner(
    label = "E – M",
    columns = vars(Estimate_EM, MSE_EM, Bias_EM)
  ) %>%
  tab_spanner(
    label = "SELF",
    columns = vars(Estimate_SELF, MSE_SELF, Bias_SELF)
  ) %>%
  tab_spanner(
    label = "Linex (a=1)",
    columns = vars(Estimate_Linex1, MSE_Linex1, Bias_Linex1)
  ) %>%
  tab_spanner(
    label = "Linex (a=-1)",
    columns = vars(Estimate_Linex1_neg, MSE_Linex1_neg, Bias_Linex1_neg)
  ) %>%
  tab_spanner(
    label = "HPD",
    columns = vars(AL_Bayes, CP_Bayes)
  ) %>%
  fmt_number(
    columns = vars(MSE_NR, Bias_NR, AL_NR, CP_NR, MSE_EM, Bias_EM, MSE_SELF, Bias_SELF, MSE_Linex1, Bias_Linex1, MSE_Linex1_neg, Bias_Linex1_neg, CP_Bayes),
    decimals = 4
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  )

# Print the table
print(gt_table)


