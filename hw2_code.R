############################
# HW2 Code 
############################

# problem 5.1

# using the normal model with a conjugate prior 
# a - compute posterior means and 95% CI's for the mean 
# and standard dev for each school 

# reading in the data 
setwd('/home/adbucks/Downloads')
school1 <- read.table("school1-1.dat.txt")
head(school1)
school2 <- read.table("school2.dat.txt")
school3 <- read.table("school3.dat.txt")

# read in the parameters for the posterior mean estimation 
mu_0 <- 5
var_0 <- 4
k_0 <- 1
gamma_0 <- 2

# observed data 
y1 <- school1$V1
n <- length(y1)
y1bar <- mean(y1)
s21 <- var(y1)

# posterior inference 
kn <- k_0 + n
gamma_n <- gamma_0 + n
mu_n <- (k_0 * mu_0 + n*y1bar) / kn
s2n <- (gamma_0*var_0 + (n-1) * s21 + k_0*n*(y1bar - mu_0)^2/(kn))/(gamma_n)

paste0("The Posterior mean is: " , round(mu_n,2))

# repeating this for each school 
y2 <- school2$V1
n2 <- length(y2)
y2bar <- mean(y2)
s22 <- var(y2)

kn2 <- k_0 + n2
gamma_n2 <- gamma_0 + n2
mu_n2 <- (k_0 * mu_0 + n2*y2bar) / kn2
s2n2 <- (gamma_0*var_0 + (n2-1) * s22 + k_0*n2*(y2bar - mu_0)^2/(kn2))/(gamma_n2)

paste0("The Posterior Mean for School 2 is: ",
       round(mu_n2 , 2))


#### school 3 
y3 <- school3$V1
n3 <- length(y3)
y3bar <- mean(y3)
s23 <- var(y3)

kn3 <- k_0 + n3
gamma_n3 <- gamma_0 + n3
mu_n3 <- (k_0 * mu_0 + n3*y3bar) / kn3
s2n3 <- (gamma_0*var_0 + (n3-1) * s23 + k_0*n3*(y3bar - mu_0)^2/(kn3))/(gamma_n3)

paste0("The Posterior Mean for School 3 is: ",
       round(mu_n3 , 2))

# confidence intervals for each now 
margin1 <- qt(0.975 , 
              df = n - 1)*sqrt(s21)/sqrt(n)
l1 <- mu_n - margin1
u1 <- mu_n + margin1

paste0("The 95% Confidence Interval for the mean is: ", 
       "(" , round(l1,2), ", ", round(u1,2), ")")

# now for the standard dev
library(MKinfer)
(cisd1 <- sdCI(y1))

# school 2
margin2 <- qt(0.975 , 
              df = n2 - 1)*sqrt(s22)/sqrt(n2)
l2 <- mu_n2 - margin2
u2 <- mu_n2 + margin2

paste("The 95% Confidence Interval for the mean of
      school 2 is: ", "(", round(l2,2) , ", ", 
      round(u2,2), ")")

(cisd2 <- sdCI(y2))

# school 3 
margin3 <- qt(0.975, 
              df = n3 - 1)*sqrt(s23)/sqrt(n3)
l3 <- mu_n3 - margin3
u3 <- mu_n3 + margin3

paste("The 95% Confidence Interval for the mean of 
      school 3 is: ", "(", round(l3,2), ", ",
      round(u3,2), ")")

(cisd3 <- sdCI(y3))

######
# the posterior probability that the ordered 
# parameter values appear in each permutation 
# of 1,2,3
######

# want to determine the posterior distributions 
# we have the normal inverse gamma in each case 
# school 1 params 
kn <- k_0 + n
gamma_n <- gamma_0 + n
mu_n <- (k_0 * mu_0 + n*y1bar) / kn
s2n <- (gamma_0*var_0 + (n-1) * s21 + k_0*n*(y1bar - mu_0)^2/(kn))/(gamma_n)

# school 2 params
kn2 <- k_0 + n2
gamma_n2 <- gamma_0 + n2
mu_n2 <- (k_0 * mu_0 + n2*y2bar) / kn2
s2n2 <- (gamma_0*var_0 + (n2-1) * s22 + k_0*n2*(y2bar - mu_0)^2/(kn2))/(gamma_n2)

# school 3 params 
kn3 <- k_0 + n3
gamma_n3 <- gamma_0 + n3
mu_n3 <- (k_0 * mu_0 + n3*y3bar) / kn3
s2n3 <- (gamma_0*var_0 + (n3-1) * s23 + k_0*n3*(y3bar - mu_0)^2/(kn3))/(gamma_n3)


# now can work with the distributions 
mu+n