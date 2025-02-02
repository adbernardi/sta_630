#########
# Problem 3.1 Hoff 
#########

# we can simulate the probability for each possible
# theta value now 

theta_v <- seq(0, 1, by=0.1)

# hyperparameters
x <- 57
n <- 100 
likelihoods <- vector(length = 11)

for (param in theta_v){
  i = 0
  likelihoods[i] <- dbinom(x, n, prob = param)
  i = i + 1
}

plot(x = theta_v , 
     y = likelihoods , 
     main = "Likelihood for Different Theta Values",
     xlab = "Theta Values" , 
     ylab = "Likelihood")
