library(dplyr)
library(sp)
library(ggplot2)
library(gstat)
library(ggmap)
library(nimble)

data <- read.csv("from_station_summary_data.csv")

# histograms

par(mfrow=c(1,2))
hist(log(data$avg_duration), main = "Average Trip Duration (Seconds)")
plot(density(log(data$avg_duration)), main = "AverageTrip Duration (Seconds)", lwd = 1.5)

# plot data
ggplot(data, aes(x = longitude, y = latitude, col = avg_duration)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_viridis_c()

# idk if we really want to get rid of the outlier or not, lmk what you think
data1 <- data %>% 
  filter(avg_duration < 3500)

par(mfrow=c(1,2))
hist(log(data1$avg_duration), main = "Average Trip Duration (Seconds)")
plot(density(log(data1$avg_duration)), main = "AverageTrip Duration (Seconds)", lwd = 1.5)

ggplot(data1, aes(x = longitude, y = latitude, col = avg_duration)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_viridis_c()

# plot empirical semivariogram + various model fits
data_sp <- data
coordinates(data_sp) = ~longitude + latitude
emp_vario <- variogram(log(avg_duration) ~ 1, data_sp, width = 0.002)

emp_vario_fit <- fit.variogram(emp_vario, model = vgm('Sph'))
plot(emp_vario, emp_vario_fit)

emp_vario_fit <- fit.variogram(emp_vario, model = vgm('Exp'))
plot(emp_vario, emp_vario_fit)

emp_vario_fit <- fit.variogram(emp_vario, model = vgm('Gau'))
plot(emp_vario, emp_vario_fit)

# I think the exponential looks best, wouldn't recommend increasing the cutoff
#exponential is definitely best


# the exponential covariance
exp_cov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0), sigma2 = double(0)) {
    returnType(double(2))
    
    nrow <- dim(dists)[1]
    ncol <- dim(dists)[2]
    result <- matrix(nrow = nrow, ncol = ncol, init = FALSE)
    
    for(i in 1:nrow) {
      for(j in 1:ncol) {
        result[i, j] <- sigma2 * exp(-dists[i,j] * phi)
      }
    }
    
    return(result)
  })

set.seed(8472)

model_code <- nimbleCode({
  
  # likelihood
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], var = tau2)
    
    mu[i] <- inprod(X[i,1:p], beta[1:p]) + W[i]
    
  }
  
  # spatial process
  Sigma_W[1:n, 1:n] <- exp_cov(dists[1:n, 1:n], phi, sigma2)
  W[1:n] ~ dmnorm(zeros[1:n], cov = Sigma_W[1:n, 1:n])
  
  # priors
  sigma2 ~ dinvgamma(2, 0.05) 
  phi ~ dunif(3/0.06, 3/0.04)
  tau2 ~ dinvgamma(2, 0.01)
  
  for (j in 1:p) {
    beta[j] ~ dnorm(0, sd = 100)
  }
  
})

#subscriber count + customer count to proportions
data$subscriber_prop <- data$subscribers/(data$subscribers + data$customers)

#also analyzing gender as proportion
data$male_prop <- data$male/(data$male + data$female)

X <- cbind(1,
           data$trip_count,
           data$capacity,
           data$subscriber_prop
           data$male_prop
           data$avg_age)

dist_mat <- as.matrix(dist(data[,c('longitude','latitude')], method = "manhattan"))

data_list <- list(y = log(data$avg_duration))
constants_list <- list(n = nrow(data),
                       X = X,
                       p = ncol(X),
                       dists = dist_mat,
                       zeros = rep(0,nrow(data)))
inits_list <- list(phi = 3/0.05,
                   sigma2 = 0.05,
                   tau2 = 0.01,
                   beta = rep(0,constants_list$p))

model1 <- nimbleModel(model_code,
                      data = data_list,
                      constants = constants_list,
                      inits = inits_list)

config1 <- configureMCMC(model1, enableWAIC = TRUE)
MCMC1 <- buildMCMC(config1)
compiled1 <- compileNimble(model1, MCMC1)

model_out1 <- runMCMC(compiled1$MCMC1,
                      nburnin = 10000,
                      niter = 50000,
                      nchains = 3,
                      samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE)

plot(model_out1$samples[,c('beta[1]', 'phi', 'sigma2', 'tau2')], density = F)
plot(model_out1$samples[,c('beta[2]', 'beta[3]', 'beta[4]', 'beta[5]')], density = F)
plot(model_out1$samples[,c('beta[6]', 'beta[7]', 'beta[8]')], density = F)
summary(model_out1$samples[,c('beta[1]',
                              'phi',
                              'sigma2',
                              'tau2',
                              'beta[2]',
                              'beta[3]',
                              'beta[4]',
                              'beta[5]',
                              'beta[6]',
                              'beta[7]',
                              'beta[8]')])

samples_array <- as.array(model_out1$samples)
post_means <- colMeans(samples_array)
post_sd <- apply(samples_array, 2, sd)
post_ci <- apply(samples_array, 2, quantile, probs = c(0.025, 0.975))
out_tab <- data.frame(post_mean = post_means,
                      post_SD = post_sd,
                      post_CI = paste0('(',sprintf('%.2f', post_ci[1,]), ', ',sprintf('%.2f', post_ci[2,]), ')')
)
knitr::kable(out_tab, digits = 2)

model_out1$WAIC

##################

# model 2

set.seed(8472)

model_code2 <- nimbleCode({
  
  # likelihood
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], var = tau2)
    
    mu[i] <- beta0 + W[i]
    
  }
  
  # spatial process
  Sigma_W[1:n, 1:n] <- exp_cov(dists[1:n, 1:n], phi, sigma2)
  W[1:n] ~ dmnorm(zeros[1:n], cov = Sigma_W[1:n, 1:n])
  
  # priors
  sigma2 ~ dinvgamma(2, 0.05) 
  phi ~ dunif(3/0.06, 3/0.04)
  tau2 ~ dinvgamma(2, 0.01)
  beta0 ~ dnorm(0,100)
  
})

data_list2 <- list(y = log(data$avg_duration))
constants_list2 <- list(n = nrow(data),
                       dists = dist_mat,
                       zeros = rep(0,nrow(data)))
inits_list2 <- list(phi = 3/0.05,
                   sigma2 = 0.05,
                   tau2 = 0.01,
                   beta0 = rnorm(1, 0, 1))

model2 <- nimbleModel(model_code2,
                      data = data_list2,
                      constants = constants_list2,
                      inits = inits_list2)

config2 <- configureMCMC(model2, enableWAIC = TRUE)
MCMC2 <- buildMCMC(config2)
compiled2 <- compileNimble(model2, MCMC2)

model_out2 <- runMCMC(compiled2$MCMC2,
                      nburnin = 10000,
                      niter = 50000,
                      nchains = 3,
                      samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE)

plot(model_out2$samples[,c('beta0', 'phi', 'sigma2', 'tau2')], density = F)

samples_array2 <- as.array(model_out2$samples)
post_means2 <- colMeans(samples_array2)
post_sd2 <- apply(samples_array2, 2, sd)
post_ci2 <- apply(samples_array2, 2, quantile, probs = c(0.025, 0.975))
out_tab2 <- data.frame(post_mean = post_means2,
                      post_SD = post_sd2,
                      post_CI = paste0('(',sprintf('%.2f', post_ci2[1,]), ', ',sprintf('%.2f', post_ci2[2,]), ')')
)
knitr::kable(out_tab2, digits = 2)

model_out2$WAIC

