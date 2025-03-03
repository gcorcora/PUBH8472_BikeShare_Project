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
    beta[j] ~ dnorm(0, sd = 10) #changed sd to 10, helped with convergence
  }
  
})

#subscriber count + customer count to proportions
data$subscriber_prop <- data$subscribers/(data$subscribers + data$customers)

#also analyzing gender as proportion
data$male_prop <- data$male/(data$male + data$female)

#I think I want to add trip distance too - because slightly longer durations might happen if the bikeshare station is full, so this kind of measures that
#removed capacity because it was highly correlated with trip count

#standardized predictors because that helped convergence
data$trip_count_standard <- scale(data$trip_count)
data$capacity_standard <- scale(data$capacity)
data$subscriber_prop_standard <- scale(data$subscriber_prop)
data$male_prop_standard <- scale(data$male_prop)
data$avg_age_standard <- scale(data$avg_age)
data$avg_trip_distance_standard <- scale(data$avg_trip_distance)


X <- cbind(1,
           data$trip_count_standard,
           data$subscriber_prop_standard,
           data$male_prop_standard,
           data$avg_age_standard,
           data$avg_trip_distance_standard
          )

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
                      nburnin = 50000,
                      niter = 200000,
                      nchains = 3,
                      samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE)

#updating betas
plot(model_out1$samples[,c('beta[1]', 'phi', 'sigma2', 'tau2')], density = F)
plot(model_out1$samples[,c('beta[2]', 'beta[3]', 'beta[4]', 'beta[5]')], density = F)
plot(model_out1$samples[,c('beta[6]')], density = F)
summary(model_out1$samples[,c('beta[1]',
                              'phi',
                              'sigma2',
                              'tau2',
                              'beta[2]',
                              'beta[3]',
                              'beta[4]',
                              'beta[5]',
                              'beta[6]')])

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

#grid of points
full_grid <- expand.grid(longitude = seq(min(data$longitude), max(data$longitude), length.out = 20),
                         latitude = seq(min(data$latitude), max(data$latitude), length.out = 20))


dist_mat <- as.matrix(dist(data[,c('longitude','latitude')], method = "manhattan"))

#combining for prediciton points
combined_dists <- rbind(full_grid, data[,c('longitude', 'latitude')])

dists_obs <- as.matrix(dist(data[,c('longitude', 'latitude')]), method = "manhattan")

dists_all <- as.matrix(dist(rbind(full_grid,data[, c('longitude', 'latitude')]), method = "manhattan"))




model_code2 <- nimbleCode({

  # likelihood
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i], var = tau2)

    mu[i] <- beta0 + W[i]

  }

  # spatial process
  Sigma_W[1:n, 1:n] <- exp_cov(dists_obs[1:n, 1:n], phi, sigma2)
  W[1:n] ~ dmnorm(zeros[1:n], cov = Sigma_W[1:n, 1:n])

  # priors
  sigma2 ~ dinvgamma(2, 0.05)
  phi ~ dunif(3/0.06, 3/0.04)
  tau2 ~ dinvgamma(2, 0.01)
  beta0 ~ dnorm(0, sd = 10)

  #posterior prediction

  #full covar among unobserved and observedl ocations
  H[1:(r + n), 1:(r+n)] <- exp_cov(dists_all[1:(r+n), 1:(r+n)], phi, sigma2)

  #covar among unobserved lcoations
  Sigma_Y0[1:r, 1:r] <- H[1:r, 1:r] + tau2*diag(r)

  #covar among observed locations
  Sigma_Y[1:n, 1:n] <- H[(r+1):(r+n), (r+1):(r+n)] + tau2*diag(n)

  #covar among observed and unobserved locations
  K[1:n, 1:r] <- H[(r+1):(r+n), 1:r]

  # computing  the  elements  of the  conditional  normal distribution
  mu_Y0_pp[1:r] <- beta0 + t(K[1:n, 1:r]) %*% inverse(Sigma_Y[1:n, 1:n]) %*%  (y[1:n] - beta0)

  Sigma_Y0_pp[1:r, 1:r] <- Sigma_Y0[1:r, 1:r] -
    t(K[1:n, 1:r]) %*% inverse(Sigma_Y[1:n, 1:n]) %*% K[1:n, 1:r]

  # sampling  from  the  posterior predictive distribution
  Y0[1:r] ~ dmnorm(mu_Y0_pp[1:r], cov = Sigma_Y0_pp[1:r, 1:r])

})

data_list2 <- list(y = log(data$avg_duration))
constants_list2 <- list(n = nrow(data),
                       dists_obs = dists_obs,
                       dists_all = dists_all,
                       r = nrow(full_grid),
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
config2$addMonitors(c('Y0'))
MCMC2 <- buildMCMC(config2)
compiled2 <- compileNimble(model2, MCMC2)

model_out2 <- runMCMC(compiled2$MCMC2,
                      nburnin = 50000,
                      niter = 200000,
                      nchains = 1,
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

samples_matrix <- as.matrix(model_out2$samples)

post_mean_Y0 <- colMeans(model_out2$samples[,grep('Y0', colnames(model_out2$samples))])

to_plot <- cbind.data.frame(full_grid, Y0 = post_mean_Y0)

library(ggplot2)

ggplot(to_plot, aes(x = s1, y = s2, fill = Y0)) + 
  geom_tile() + 
  scale_fill_viridis_c() + 
  theme_minimal()

ggplot(to_plot, aes(x=longitude, y = latitude, fill = Y0)) + 
  geom_tile() +
  scale_fill_viridis_c(name = "Y0 = Log(Avg. Trip Duration)") +
  ggtitle("Bike Share Kriging - Intercept Only Model")
  theme_minimal()
  

model_out2$WAIC

### Model without outlier

data <- read.csv("/Users/gretchen/Desktop/from_station_summary_data.csv")

data1 <- data %>%
  filter(avg_duration < 3500)

data_sp1 <- data1
coordinates(data_sp1) = ~longitude + latitude
emp_vario <- variogram(log(avg_duration) ~ 1, data_sp1, width = 0.001)

emp_vario_fit <- fit.variogram(emp_vario, model = vgm('Exp'))
plot(emp_vario, emp_vario_fit)

library(sf)

#exporting so can load in QGIS map

sf_data <- st_as_sf(data1, coords = c("longitude", "latitude"), crs = 4326)
st_write(sf_data, "bikeshare_data_no_outlier.shp", delete_layer = TRUE)


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
    beta[j] ~ dnorm(0, sd = 10)
  }

})

#subscriber count + customer count to proportions
data$subscriber_prop <- data$subscribers/(data$subscribers + data$customers)

#also analyzing gender as proportion
data$male_prop <- data$male/(data$male + data$female)

#I think I want to add trip distance too - because slightly longer durations might happen if the bikeshare station is full, so this kind of measures that
#trying to remove capacity

#standardizing predictors
data$trip_count_std <- scale(data$trip_count)
data$subscriber_prop_std <- scale(data$subscriber_prop)
data$male_prop_std <- scale(data$male_prop)
data$avg_age_std <- scale(data$avg_age)
data$avg_trip_distance_std <- scale(data$avg_trip_distance)

X <- cbind(1,
           data$trip_count_std,
           data$subscriber_prop_std,
           data$male_prop_std,
           data$avg_age_std,
           data$avg_trip_distance_std
)

dist_mat <- as.matrix(dist(data1[,c('longitude','latitude')], method = "manhattan"))

data_list <- list(y = log(data1$avg_duration))
constants_list <- list(n = nrow(data1),
                       X = X,
                       p = ncol(X),
                       dists = dist_mat,
                       zeros = rep(0,nrow(data1)))
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
                      nburnin = 50000,
                      niter = 200000,
                      nchains = 3,
                      samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE)

#adding convergence diag - if more than 1 chain
coda::gelman.diag(model_out1$samples)

plot(model_out1$samples[,c('beta[1]')], density = F)
plot(model_out1$samples[,c('beta[2]', 'beta[3]', 'beta[4]', 'beta[5]')], density = F)
plot(model_out1$samples[,c('beta[6]')], density = F)
summary(model_out1$samples[,c('beta[1]',
                              'phi',
                              'sigma2',
                              'tau2',
                              'beta[2]',
                              'beta[3]',
                              'beta[4]',
                              'beta[5]',
                              'beta[6]')])

samples_array <- as.array(model_out1$samples)
post_means <- colMeans(samples_array)
post_sd <- apply(samples_array, 2, sd)
post_ci <- apply(samples_array, 2, quantile, probs = c(0.025, 0.975))
out_tab <- data.frame(post_mean = post_means,
                      post_SD = post_sd,
                      post_CI = paste0('(',sprintf('%.2f', post_ci[1,]), ', ',sprintf('%.2f', post_ci[2,]), ')')
)
knitr::kable(out_tab, digits = 2, caption = "Posterior Means and Credible Intervals - Without Outlier")
library(gridExtra)

table_plot <- tableGrob(out_tab)
ggsave("posterior_summary_table.png", plot = table_plot, width = 8, height = 4, dpi = 300)


WAIC_no_outlier = model_out1$WAIC

with_outlier_WAIC_value <- with_outlier_WAIC$WAIC
WAIC_no_outlier_value <- WAIC_no_outlier$WAIC
waic_values <- c(as.numeric(with_outlier_WAIC_value), as.numeric(WAIC_no_outlier_value))

# Create the WAIC comparison table
waic_table <- data.frame(
  Model = c("With Outlier", "Without Outlier"),
  WAIC = c(with_outlier_WAIC$WAIC, WAIC_no_outlier$WAIC),
  lppd = c(with_outlier_WAIC$lppd, WAIC_no_outlier$lppd),
  pWAIC = c(with_outlier_WAIC$pWAIC, WAIC_no_outlier$pWAIC)
)
