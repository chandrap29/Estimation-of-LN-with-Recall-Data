    rm(list=ls(all=TRUE))
    #===================================================================#
    #                                                                   #
    #             Importing and Filtering Data                          #
    #                                                                   #
    #===================================================================#  
    
    
    library("VGAM");library("zipfR");library("coda") # loading necessary packages
    Data_2<-read.csv("C:/Users/Chandra Prakash/Dropbox/LN with CS/Real Data Analysis/MeNcom.csv",header=T) #importing data
    #View(Data_2)
    group<-Data_2$ToR;group #extract type of recall
    
    #different categories of recall in data
    ind1<-which(group=="Day Recall");length(ind1)        # position for Day recall   #1
    ind2<-which(group=="Month Recall");length(ind2)      # position for Month recall #2
    ind3<-which(group=="No Recall");length(ind3)         # position for No recall    #3  
    ind4<-which(group=="Not Happened");length(ind4)      # position for Not happened #4
    ind5<-which(group=="Year Recall");length(ind5)       # position for Year recall  #5  
    
    n=length(ind1)+length(ind2)+length(ind3) +length(ind4) +length(ind5);n #total sample size 
    
    #Indicator Variables for different categories
    delta=ifelse(Data_2$ToR=="Not Happened",0,1);delta        #Censoring indicator
    delta0=which(delta==0);length(delta0)                     #position of censored
    delta1=which(delta==1);length(delta1)                     #position of non-censored
    NC=length(delta1)/(length(delta0)+length(delta1));NC      #0.8784461: proportion of non-censored
    CN=length(delta0)/(length(delta0)+length(delta1));CN      #0.1215539: proportion of censored
    
    #Here, we considered "Day Recall" in recall category and merged rest as Month and Year levels recall into non-recall
    da_nr=c(ind2,ind3,ind5)       # index of non recall respect to day recall    
    
    #find position corresponding to recall 
    rem1<-which(Data_2$cause=="1st born"&Data_2$ToR=="Day Recall") ;rem1
    rem2<-which(Data_2$cause=="Later born"&Data_2$ToR=="Day Recall");rem2
    rem3=c(rem1,rem2);rem3
    
    S1<-Data_2[rem3,2];S1                     # monitoring time for recall category
    z1<-Data_2[rem3,9];z1                     # time for recall category      
    
    S2<-Data_2[da_nr,2];S2                    # monitoring time for non-recall
    S3=Data_2[ind4,2]                         # monitoring time for censored 
    
    n1<-length(z1);n1                        #sample size for recall
    n2<-length(S2);n2                        #sample size for non-recall   
    n3<-length(S3);n3                        #sample size for censored
    n<-n1+n2+n3;n                            #sample size
    
    #summary(S1);summary(S2);summary(S3)     #summary of monitoring times
    aa=15
    z1=z1/aa;S1=S1/aa;S2=S2/aa;S3=S3/aa  # Scaling of data
    
    
    #===================================================================#
    #                                                                   #
    #                        ML Estimate (Using EM)                     #
    #                                                                   #
    #===================================================================#  
    
    #Fit Weibull distribution on time to event to initialize parameters
    #install.packages("fitdistrplus")
    library(fitdistrplus)
    est=fitdist(z1,"norm")
    sc=est$estimate;sc
    a1=sc[1];b1=sc[2]
    a1;b1
    
    lower1=min(S2)/10;lower1 #set the lower limit for non-recall
    
    
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
    mu_init <- a1
    sgm_init <- b1
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
    ml_opt<-c(nl$par[1],sqrt(nl$par[2]),nl$par[3])
    
    #variance-covariance    
    var_cov=solve(nl$hessian);var_cov
    
    ci_l=ci_u=rep();alpha=0.05
    for (i in 1:3) {
      ci_l[i]=ml_opt[i]-qnorm(1-(alpha/2))*sqrt((diag(var_cov)[i]))
      ci_u[i]=ml_opt[i]+qnorm(1-(alpha/2))*sqrt((diag(var_cov)[i]))
    }
    
    #limits
    lmt_opt<-c(ci_l[1],ci_u[1],ci_l[2],ci_u[2],ci_l[3],ci_u[3]);lmt_opt
    
    #length    
    l1=c(ci_u[1]-ci_l[1]);l1
    l2=c(ci_u[2]-ci_l[2]);l2
    l3=c(ci_u[3]-ci_l[3]);l3
    leng_opt=c(l1,l2,l3);leng_opt
    
    ml=ml_opt;ml #ML estimates
    
    # Mean and Median duration of menarche in years
    mean_nr=ml[1]*15;mean_nr
    med_nr=ml[1]*15;med_nr
    
    # Confidence Intervals for Mean and Median
    mean_ci <- c(ci_l[1] * 15, ci_u[1] * 15)
    median_ci <- mean_ci  # Identical for normal distribution
    
    NR<-c(ml_opt,lmt_opt, leng_opt, mean_nr, med_nr, mean_ci, median_ci);NR
    
    
    #===================================================================#
    #                                                                   #
    #                        ML Estimate (Using EM)                     #------------
    #                                                                   #
    #===================================================================#  
    
    # Empirical initialization
    m.th1 <- mu_init                   # Mean of observed events
    m.th2 <- sgm_init                     # Standard deviation of observed events
    m.th3 <- lmd_init          # Approximate initial value for lambda
    
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
    
    ml_em=c(m.thh1[it],m.thh2[it],m.thh3[it]);ml_em  #store the final estimates
    
    # Mean and Median duration of menarche in years
    mean_em=ml_em[1]*15;mean_em
    med_em=ml_em[1]*15;med_em
    
    # Variance-Covariance Matrix
    # Fisher Information Approximation
    fisher_info <- matrix(0, 3, 3)
    fisher_info[1, 1] <- n / (ml_em[2]^2)
    fisher_info[2, 2] <- (2 * n) / (ml_em[2]^2)
    fisher_info[3, 3] <- n / (ml_em[3]^2)
    
    # Invert Fisher Information to get Variance-Covariance Matrix
    var_cov_em <- solve(fisher_info)
    
    # Confidence Intervals for Parameters
    alpha <- 0.05
    ci_l_em <- numeric(3)
    ci_u_em <- numeric(3)
    
    for (i in 1:3) {
      ci_l_em[i] <- ml_em[i] - qnorm(1 - (alpha / 2)) * sqrt(diag(var_cov_em)[i])
      ci_u_em[i] <- ml_em[i] + qnorm(1 - (alpha / 2)) * sqrt(diag(var_cov_em)[i])
    }
    
    # Scaling the Mean and Median CIs
    mean_ci_em <- c(ci_l_em[1] * 15, ci_u_em[1] * 15)
    median_ci_em <- mean_ci_em  # Identical for normal distribution
    
    EM<-c(ml_em, mean_em, med_em, mean_ci_em, median_ci_em);EM
    
    #===================================================================#
    #                                                                   #
    #                        Bayesian Estimation                        #------------------
    #                                                                   #
    #===================================================================#
    # Load necessary package
    # install.packages("LaplacesDemon") # Uncomment to install if necessary
    library(LaplacesDemon);library(coda)
    
    # Priors derived from table of medians
    mu0 <- 11.93      # Mean of medians as prior for mu
    sig0 <- 0.25      # Increased prior mean for variance
    nu0 <- 2.0        # Reduced prior variance to increase flexibility
    k0 <- 0.5         # Shape parameter for inverse chi-squared prior
    l0 <- 1.0         # Shape parameter for lambda (can remain as is)
    m0 <- 1.0         # Shape parameter for lambda (can remain as is)
    
    
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

    # Histograms for posterior distributions
    # par(mfrow = c(1, 3))
    # hist(mu.th, breaks = 30, main = "", xlab = expression(mu))
    # hist(sgm.th, breaks = 30, main = "", xlab = expression(sigma))
    # hist(lmd.th, breaks = 30, main = "", xlab = expression(lambda))


    # Plot ACF for mu, sigma, and lambda
    # par(mfrow = c(1, 3))
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
    
    #quantile plots for chains
    par(mfrow=c(1,3))
    gx1=as.vector(ch11);gx1
    ii=1:length(gx1)
    q1.25=sapply(ii,function(ii) quantile((gx1[1:ii]), probs = c(.25)));q1.25
    q1.50=sapply(ii,function(ii) quantile((gx1[1:ii]), probs = c(.5)));q1.50
    q1.75=sapply(ii,function(ii) quantile((gx1[1:ii]), probs = c(.75)));q1.75
    
    gx2=na.omit(as.vector(ch22));gx2
    ii=1:length(gx2)
    q2.25=sapply(ii,function(ii) quantile((gx2[1:ii]), probs = c(.25),na.rm=TRUE));q2.25
    q2.50=sapply(ii,function(ii) quantile((gx2[1:ii]), probs = c(.5),na.rm=TRUE));q2.50
    q2.75=sapply(ii,function(ii) quantile((gx2[1:ii]), probs = c(.75),na.rm=TRUE));q2.75
    
    gx3=as.vector(ch33)
    ii=1:length(gx3)
    q3.25=sapply(ii,function(ii) quantile((gx3[1:ii]), probs = c(.25)));q3.25
    q3.50=sapply(ii,function(ii) quantile((gx3[1:ii]), probs = c(.5)));q3.50
    q3.75=sapply(ii,function(ii) quantile((gx3[1:ii]), probs = c(.75)));q3.75
    
    #Trace Plots
    matplot(ii,cbind(q1.25,q1.50,q1.75),main="",type="l",col=1,xlab="iteration",ylab=expression(mu))
    matplot(ii,cbind(q2.25,q2.50,q2.75),main="",type="l",col=1,xlab="iteration",ylab=expression(sigma))
    matplot(ii,cbind(q3.25,q3.50,q3.75),main="",type="l",col=1,xlab="iteration",ylab=expression(lambda))
    
    #ACF Plots
    acf(ch11,main="",col="black",xlab="Lag",ylab=expression("ACF"(mu)))
    acf(ch22,main="",col="black",xlab="Lag",ylab=expression("ACF"(sigma)))
    acf(ch33,main="",col="black",xlab="Lag",ylab=expression("ACF"(lambda)))
    
    #MAT Plots
    matplot(ch11,type="l",col="black",xlab="iterations",ylab=expression(mu),main="")
    matplot(ch22,type="l",col="black",xlab="iterations",ylab=expression(sigma),main="")
    matplot(ch33,type="l",col="black",xlab="iterations",ylab=expression(lambda),main="")
    
    #density plots
    plot(density(ch11),col="black",type="l",ylab=expression("Density"(mu)),main="")
    plot(density(ch22),col="black",type="l",ylab=expression("Density"(sigma)),main="")
    plot(density(ch33),col="black",type="l",ylab=expression("Density"(lambda)),main="")
    
    
    #===========Bayes estimates under SELF--------------
    mu.tth=sgm.tth=lmd.tth=rep()
    mu.tth=mean(ch11);sgm.tth=mean(ch22);lmd.tth=mean(ch33);
    b.th=c(mu.tth,sgm.tth,lmd.tth);b.th 
    self_estimates<-b.th

    
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
    
    
    #Case 2. 
    a1 = -1  # Modify 'a' as needed for asymmetry
    mu_LINEX1 = -1/a1 * log(mean(exp(-a1 * ch11)))
    sgm_LINEX1 = -1/a1 * log(mean(exp(-a1 * ch22)))
    lmd_LINEX1 = -1/a1 * log(mean(exp(-a1 * ch33)))
    
    linex_estimates1 = c(mu_LINEX1, sgm_LINEX1, lmd_LINEX1)
    linex_estimates1
    
    #===============HPD Intervals-------------------
    h.th1=HPDinterval(mcmc(ch11));h.th2=HPDinterval(mcmc(ch22));h.th3=HPDinterval(mcmc(ch33))
    HPD_lim=c(h.th1[,1],h.th1[,2],h.th2[,1],h.th2[,2],h.th3[,1],h.th3[,2]);HPD_lim
    
    h.th1_l=h.th1[,2]-h.th1[,1];h.th2_l=h.th2[,2]-h.th2[,1];h.th3_l=h.th3[,2]-h.th3[,1]
    HPD_L=c(h.th1_l,h.th2_l,h.th3_l);HPD_L                                   #length of HPD
    
    
     Bayes=c(self_estimates,linex_estimates,linex_estimates1,HPD_lim,HPD_L);Bayes
    
     
     
     # Assuming ch11 represents the posterior chain for the mean (mu)
     
     # ========================== SELF Loss Function ============================
     mu_self = mean(15*ch11)  # Mean of mu
     mu_median_self = median(ch11*15)  # Median of mu
     
     
     # ========================== LINEX Loss Function (a = 1) ===================
     a = 1  # Modify 'a' as needed for asymmetry
     mu_linex1 = -1/a * log(mean(exp(-a * 15*ch11)))
     
     mu_median_linex1 = median(mu_linex1)
     
     
     # ========================== LINEX Loss Function (a = -1) ===================
     a1 = -1  # Modify 'a' as needed for asymmetry
     mu_linex_neg1 = -1/a1 * log(mean(exp(-a1 * 15*ch11)))
     
     mu_median_linex_neg1 = median(mu_linex_neg1)
     
     
    #===================================================================#
    #                                                                   #
    #                        End of R code                              #
    #                                                                   #
    #===================================================================# 
    