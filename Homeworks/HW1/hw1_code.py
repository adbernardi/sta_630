###########################
# code for homework 1 problem 3  
########################### 

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

# setting up our initial variables 
theta_space = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
x = 57 
n = 100 # for binomial probabilities 

# calculating the binomial probabilities for each theta in the theta spapce 
import scipy.stats as stats 

# can try to do this quickly 
binom_probs = [stats.binom.pmf(x, n, theta) for theta in theta_space]

print(binom_probs)

# plotting the binomial probabilities 
plt.plot(theta_space, binom_probs) 
plt.xlabel('Theta') 
plt.ylabel('Binomial Probability') 
plt.title('Binomial Probabilities for x = 57, n = 100') 
plt.show() 

########
# part c 
######## 

# using bayes formula assuming equal probabilities for each theta 
prior = [1/11 for i in range(11)]
complement = [1 - i for i in prior] 
likelihood = binom_probs 

# can now do this for each theta and its complement 
posterior = [prior[i] * likelihood[i] / (prior[i] * likelihood[i] + complement[i] * (1 - likelihood[i])) for i in range(11)]
print("Posterior Probabilities: ", posterior)

# can now plot these and find the max 
plt.plot(theta_space , posterior) 
plt.xlabel('Theta') 
plt.ylabel('Posterior Probability') 
plt.title('Posterior Probabilities for x = 57, n = 100') 
plt.show() 

#########
# part d 
######### 

# trying the posterior again, with theta on an interval from [0,1]
theta_space_c = np.linspace(0,1, 1000) 
binom_probs_c = [stats.binom.pmf(x, n, theta) for theta in theta_space_c] 

# now we can plot 
plt.plot(theta_space_c, binom_probs_c) 
plt.xlabel('Theta') 
plt.ylabel('Binomial Probability') 
plt.title('Binomial Probabilities for x = 57, n = 100') 
plt.show() 

#############
# part e 
############# 

# plotting the beta posterior for each value of theta 

# can recycle our theta space 
alpha = 1 + 57 
beta = 1 + 100 - 57 

# calculating the beta posterior 
beta_posterior = [stats.beta.pdf(theta, alpha, beta) for theta in theta_space_c]

# now can plot with theta on the x axis 
plt.plot(theta_space_c, beta_posterior) 
plt.xlabel('Theta') 
plt.ylabel('Beta Posterior') 
plt.title('Beta Posterior for x = 57, n = 100') 
plt.show() 


