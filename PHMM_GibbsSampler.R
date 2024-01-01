library(readxl)
library(lubridate) 
library(dplyr) 
library(jagsUI)
library(MCMCpack)
library(coda)
library(mcmcse)
library(HiddenMarkov)
library(writexl)
library(nortest)

# MEMBACA DATA
aceh <- read_xlsx("nama_dataset.xlsx") # Ubah 'nama_dataset' menjadi dataset yang digunakan

# PHMM 2 HIDDEN STATES
x <- aceh$Count 
n <- length(x)
m <- 2
data <- list(x = x, n = n, m = m) 

# MODEL DALAM BAHASA BUGS 
modfile <- tempfile()
writeLines("
model{
  for (i in 1:m){
    delta[i] <- 1/m
    v[i] <- 1
  }
  s[1]~dcat(delta[])
  for (i in 2:n){
    s[i]~dcat(Gamma[s[i-1],])
  }
  states[1]~dcat(Gamma[s[n],1])
  x[1]~dpois(lambda[states[1]])
  for(i in 2:n){
    states[i]~dcat(Gamma[states[i-1],])
    x[i]~dpois(lambda[states[i]])
  }
  for (i in 1:m){
    tau[i]~dgamma(0.02848238, 0.04503295)
    Gamma[i,1:m]~ddirich(v[])
  }
  lambda[1] <- tau[1]
  
  for (i in 2:m){
    lambda[i] <- lambda[i-1] + tau[i]
  }
}", con=modfile)

Tau21 <- dgamma(0.001, 0.001)
Tau22 <- dgamma(0.001, 0.001)
Tau20 <- matrix(c(Tau21, Tau22), byrow = TRUE, ncol = 2)
Tau20
Lambda21 <- Tau21
Lambda22 <- Tau21 + Tau22
Lambda20 <- matrix(c(Lambda21, Lambda22), byrow = TRUE, ncol = 2)
Lambda20
Gamma20 <- rdirichlet(2, c(1,1))
Gamma20

initial20 <- function(){  
  list(Tau20 = Tau20, Lambda20 = Lambda20, Gamma20 = Gamma20)  
}

# PARAMETER YANG AKAN DI MONITOR
params <- c('tau','lambda','Gamma')

tau1 = NULL
tau2 = NULL
lambda1 = NULL
lambda2 = NULL
Gamma11 = NULL
Gamma12 = NULL
Gamma21 = NULL
Gamma22 = NULL

for(i in 1:50){
  eq.sims20 <- jags(data = data,
                    inits = initial20,
                    parameters.to.save = params,
                    model.file = modfile,
                    n.chains = 1,
                    n.iter = 100000,
                    n.burnin = 5000,
                    n.thin = 60)
  
  tau1 = c(tau1, eq.sims20$mean$tau[1])
  tau2 = c(tau2, eq.sims20$mean$tau[2])
  lambda1 = c(lambda1, eq.sims20$mean$lambda[1])
  lambda2 = c(lambda2, eq.sims20$mean$lambda[2])
  Gamma11 = c(Gamma11, eq.sims20$mean$Gamma[1,1])
  Gamma12 = c(Gamma12, eq.sims20$mean$Gamma[1,2])
  Gamma21 = c(Gamma21, eq.sims20$mean$Gamma[2,1])
  Gamma22 = c(Gamma22, eq.sims20$mean$Gamma[2,2])
}

tau1
tau2
lambda1
lambda2
Gamma11
Gamma12
Gamma21
Gamma22

mean(tau1)
mean(tau2)
mean(lambda1)
mean(lambda2)
mean(Gamma11)
mean(Gamma12)
mean(Gamma21)
mean(Gamma22)

# UJI NORMALITAS 50 SAMPEL POSTERIOR
data <- c(tau1) # Ubah sesuai parameter yang akan diuji normalitas
shapiro_result <- shapiro.test(data)
shapiro_result

# GRAFIK JEJAK RANTAI (TRACE PLOTS) DAN GRAFIK ACF 
par(mfrow=c(1,1))
plot(eq.sims20[["sims.list"]][["tau"]][,1], type = 'l') # Ubah 'tau' dan '[,1]' sesuai dengan parameter plot yang dibutuhkan
autocorr.plot(eq.sims20, col = "blue")

# HASIL ESTIMASI PARAMETER
Gamma2 <- eq.sims20[["mean"]][["Gamma"]]
Lambda2 <- eq.sims20[["mean"]][["lambda"]]
delta2 <- Gamma2
for (i in 1:n)
  delta2 <- delta2%*%Gamma2
Delta2 <- delta2[1,]

Gamma2
Lambda2
Delta2

# BANGKITAN PHMM TELAH DICOCOKAN
phmm <- dthmm(aceh$Count,
              Gamma2,
              Delta2,
              "pois",
              list(lambda = Lambda2),
              discrete = TRUE)
n = length(aceh$Count)
x <- simulate(phmm, nsim = n, seed = 10)
y <- residuals(x)
w <- aceh$Count + y
DataDuga2MCMC = ifelse(w <= 0, 0, round(w))
DataDuga2MCMC

# TES AKURASI PHMM YANG TELAH DICOCOKAN
rmse_mcmc <- sqrt(mean((aceh$Count - DataDuga2MCMC)^2))
rmse_mcmc

abs <- abs(y)
Totalresiduals <- sum(abs)
MAE2 <- Totalresiduals/n
MAE2


# PHMM 3 HIDDEN STATES
x <- aceh$Count
n <- length(x)
m <- 3
data <- list(x = x, n = n, m = m) 

# MODEL DALAM BAHASA BUGS
modfile <- tempfile()
writeLines("
model{
  for (i in 1:m){
    delta[i] <- 1/m
    v[i] <- 1
  }
  s[1]~dcat(delta[])
  for (i in 2:n){
    s[i]~dcat(Gamma[s[i-1],])
  }
  states[1]~dcat(Gamma[s[n],1])
  x[1]~dpois(lambda[states[1]])
  for(i in 2:n){
    states[i]~dcat(Gamma[states[i-1],])
    x[i]~dpois(lambda[states[i]])
  }
  for (i in 1:m){
    tau[i]~dgamma(0.02848238, 0.04503295)
    Gamma[i,1:m]~ddirich(v[])
  }
  lambda[1] <- tau[1]
  
  for (i in 2:m){
    lambda[i] <- lambda[i-1] + tau[i]
  }
}", con=modfile)

Tau31 <- dgamma(0.001, 0.001)
Tau32 <- dgamma(0.001, 0.001)
Tau33 <- dgamma(0.001, 0.001)
Tau30 <- matrix(c(Tau31, Tau32, Tau33), byrow = TRUE, ncol = 3)
Tau30
Lambda31 <- Tau31
Lambda32 <- Tau31 + Tau32 
Lambda33 <- Tau31 + Tau32 + Tau33
Lambda30 <- matrix(c(Lambda31, Lambda32, Lambda33), byrow = TRUE, ncol = 3)
Lambda30
Gamma30 <- rdirichlet(3, c(1,1,1))
Gamma30

initial30 <- function(){  
  list(Tau30=Tau30,Lambda30=Lambda30,Gamma30=Gamma30)  
}

# PARAMETER YANG AKAN DI MONITOR
params <- c('tau','lambda','Gamma')

tau1 = NULL
tau2 = NULL
tau3 = NULL
lambda1 = NULL
lambda2 = NULL
lambda3 = NULL
Gamma11 = NULL
Gamma12 = NULL
Gamma13 = NULL
Gamma21 = NULL
Gamma22 = NULL
Gamma23 = NULL
Gamma31 = NULL
Gamma32 = NULL
Gamma33 = NULL

for(i in 1:50){
  eq.sims30 <- jags(data = data,
                    inits = initial30,
                    parameters.to.save = params,
                    model.file = modfile,
                    n.chains = 1,
                    n.iter = 100000,
                    n.burnin = 5000,
                    n.thin = 60)
  
  tau1 = c(tau1, eq.sims30$mean$tau[1])
  tau2 = c(tau2, eq.sims30$mean$tau[2])
  tau3 = c(tau3, eq.sims30$mean$tau[3])
  lambda1 = c(lambda1, eq.sims30$mean$lambda[1])
  lambda2 = c(lambda2, eq.sims30$mean$lambda[2])
  lambda3 = c(lambda3, eq.sims30$mean$lambda[3])
  Gamma11 = c(Gamma11, eq.sims30$mean$Gamma[1,1])
  Gamma12 = c(Gamma12, eq.sims30$mean$Gamma[1,2])
  Gamma13 = c(Gamma13, eq.sims30$mean$Gamma[1,3])
  Gamma21 = c(Gamma21, eq.sims30$mean$Gamma[2,1])
  Gamma22 = c(Gamma22, eq.sims30$mean$Gamma[2,2])
  Gamma23 = c(Gamma23, eq.sims30$mean$Gamma[2,3])
  Gamma31 = c(Gamma31, eq.sims30$mean$Gamma[3,1])
  Gamma32 = c(Gamma32, eq.sims30$mean$Gamma[3,2])
  Gamma33 = c(Gamma33, eq.sims30$mean$Gamma[3,3])
}

tau1
tau2
tau3
lambda1
lambda2
lambda3
Gamma11
Gamma12
Gamma13
Gamma21
Gamma22
Gamma23
Gamma31
Gamma32
Gamma33

mean(tau1)
mean(tau2)
mean(tau3)
mean(lambda1)
mean(lambda2)
mean(lambda3)
mean(Gamma11)
mean(Gamma12)
mean(Gamma13)
mean(Gamma21)
mean(Gamma22)
mean(Gamma23)
mean(Gamma31)
mean(Gamma32)
mean(Gamma33)

# UJI NORMALITAS 50 SAMPEL POSTERIOR
data <- c(tau1) # Ubah sesuai parameter yang akan diuji normalitas
shapiro_result <- shapiro.test(data)
shapiro_result

# GRAFIK JEJAK RANTAI (TRACE PLOTS) DAN GRAFIK ACF 
par(mfrow=c(1,1))
plot(eq.sims30[["sims.list"]][["tau"]][,1], type = 'l') # Ubah 'tau' dan '[,1]' sesuai dengan parameter plot yang dibutuhkan
autocorr.plot(eq.sims30, col = "blue")

# HASIL ESTIMASI PARAMETER
Gamma3 <- eq.sims30[["mean"]][["Gamma"]]
Lambda3 <- eq.sims30[["mean"]][["lambda"]]
delta3 <- Gamma3
for (i in 1:n)
  delta3 <- delta3%*%Gamma3
Delta3 <- delta3[1,]

Gamma3
Lambda3
Delta3

# BANGKITAN PHMM TELAH DICOCOKAN
phmm <- dthmm(aceh$Count,
              Gamma3,
              Delta3,
              "pois",
              list(lambda = Lambda3),
              discrete = TRUE)
n = length(aceh$Count)
x <- simulate(phmm, nsim = n, seed = 10)
y <- residuals(x)
w <- aceh$Count + y
DataDuga3MCMC = ifelse(w <= 0, 0, round(w))
DataDuga3MCMC

# TES AKURASI PHMM YANG TELAH DICOCOKAN
rmse_mcmc <- sqrt(mean((aceh$Count - DataDuga3MCMC)^2))
rmse_mcmc

abs <- abs(y)
Totalresiduals <- sum(abs)
MAE3 <- Totalresiduals/n
MAE3


# PHMM 4 HIDDEN STATES
x <- aceh$Count
n <- length(x)
m <- 4
data <- list(x = x, n = n, m = m) 

# MODEL DALAM BAHASA BUGS 
modfile <- tempfile()
writeLines("
model{
  for (i in 1:m){
    delta[i] <- 1/m
    v[i] <- 1
  }
  s[1]~dcat(delta[])
  for (i in 2:n){
    s[i]~dcat(Gamma[s[i-1],])
  }
  states[1]~dcat(Gamma[s[n],1])
  x[1]~dpois(lambda[states[1]])
  for(i in 2:n){
    states[i]~dcat(Gamma[states[i-1],])
    x[i]~dpois(lambda[states[i]])
  }
  for (i in 1:m){
    tau[i]~dgamma(0.02848238, 0.04503295)
    Gamma[i,1:m]~ddirich(v[])
  }
  lambda[1] <- tau[1]
  
  for (i in 2:m){
    lambda[i] <- lambda[i-1] + tau[i]
  }
}", con=modfile)

Tau41 <- dgamma(0.001, 0.001)
Tau42 <- dgamma(0.001, 0.001)
Tau43 <- dgamma(0.001, 0.001)
Tau44 <- dgamma(0.001, 0.001)
Tau40 <- matrix(c(Tau41, Tau42, Tau43, Tau44), byrow = TRUE, ncol = 4)
Tau40
Lambda41 <- Tau41
Lambda42 <- Tau41 + Tau42 
Lambda43 <- Tau41 + Tau42 + Tau43
Lambda44 <- Tau41 + Tau42 + Tau43 + Tau44
Lambda40 <- matrix(c(Lambda41, Lambda42, Lambda43, Lambda44), byrow = TRUE, ncol = 4)
Lambda40
Gamma40 <- rdirichlet(4, c(1,1,1,1))
Gamma40

initial40 <- function(){  
  list(Tau40=Tau40,Lambda40=Lambda40,Gamma40=Gamma40)  
}

#PARAMETER YANG AKAN DI MONITOR
params <- c('tau','lambda','Gamma')

tau1 = NULL
tau2 = NULL
tau3 = NULL
tau4 = NULL
lambda1 = NULL
lambda2 = NULL
lambda3 = NULL
lambda4 = NULL
Gamma11 = NULL
Gamma12 = NULL
Gamma13 = NULL
Gamma14 = NULL
Gamma21 = NULL
Gamma22 = NULL
Gamma23 = NULL
Gamma24 = NULL
Gamma31 = NULL
Gamma32 = NULL
Gamma33 = NULL
Gamma34 = NULL
Gamma41 = NULL
Gamma42 = NULL
Gamma43 = NULL
Gamma44 = NULL

for(i in 1:50){
  eq.sims40 <- jags(data = data,
                    inits = initial40,
                    parameters.to.save = params,
                    model.file = modfile,
                    n.chains = 1,
                    n.iter = 800000,
                    n.burnin = 250000,
                    n.thin = 400)
  
  tau1 = c(tau1, eq.sims40$mean$tau[1])
  tau2 = c(tau2, eq.sims40$mean$tau[2])
  tau3 = c(tau3, eq.sims40$mean$tau[3])
  tau4 = c(tau4, eq.sims40$mean$tau[4])
  lambda1 = c(lambda1, eq.sims40$mean$lambda[1])
  lambda2 = c(lambda2, eq.sims40$mean$lambda[2])
  lambda3 = c(lambda3, eq.sims40$mean$lambda[3])
  lambda4 = c(lambda4, eq.sims40$mean$lambda[4])
  Gamma11 = c(Gamma11, eq.sims40$mean$Gamma[1,1])
  Gamma12 = c(Gamma12, eq.sims40$mean$Gamma[1,2])
  Gamma13 = c(Gamma13, eq.sims40$mean$Gamma[1,3])
  Gamma14 = c(Gamma14, eq.sims40$mean$Gamma[1,4])
  Gamma21 = c(Gamma21, eq.sims40$mean$Gamma[2,1])
  Gamma22 = c(Gamma22, eq.sims40$mean$Gamma[2,2])
  Gamma23 = c(Gamma23, eq.sims40$mean$Gamma[2,3])
  Gamma24 = c(Gamma24, eq.sims40$mean$Gamma[2,4])
  Gamma31 = c(Gamma31, eq.sims40$mean$Gamma[3,1])
  Gamma32 = c(Gamma32, eq.sims40$mean$Gamma[3,2])
  Gamma33 = c(Gamma33, eq.sims40$mean$Gamma[3,3])
  Gamma34 = c(Gamma34, eq.sims40$mean$Gamma[3,4])
  Gamma41 = c(Gamma41, eq.sims40$mean$Gamma[4,1])
  Gamma42 = c(Gamma42, eq.sims40$mean$Gamma[4,2])
  Gamma43 = c(Gamma43, eq.sims40$mean$Gamma[4,3])
  Gamma44 = c(Gamma44, eq.sims40$mean$Gamma[4,4])
}

tau1
tau2
tau3
tau4
lambda1
lambda2
lambda3
lambda4
Gamma11
Gamma12
Gamma13
Gamma14
Gamma21
Gamma22
Gamma23
Gamma24
Gamma31
Gamma32
Gamma33
Gamma34
Gamma41
Gamma42
Gamma43
Gamma44

mean(tau1)
mean(tau2)
mean(tau3)
mean(tau4)
mean(lambda1)
mean(lambda2)
mean(lambda3)
mean(lambda4)
mean(Gamma11)
mean(Gamma12)
mean(Gamma13)
mean(Gamma14)
mean(Gamma21)
mean(Gamma22)
mean(Gamma23)
mean(Gamma24)
mean(Gamma31)
mean(Gamma32)
mean(Gamma33)
mean(Gamma34)
mean(Gamma41)
mean(Gamma42)
mean(Gamma43)
mean(Gamma44)

# UJI NORMALITAS 50 SAMPEL POSTERIOR
data <- c(tau1) # Ubah sesuai parameter yang akan diuji normalitas
shapiro_result <- shapiro.test(data)
shapiro_result

# GRAFIK JEJAK RANTAI (TRACE PLOTS) DAN GRAFIK ACF 
par(mfrow=c(1,1))
plot(eq.sims40[["sims.list"]][["tau"]][,1], type = 'l') # Ubah 'tau' dan '[,1]' sesuai dengan parameter plot yang dibutuhkan
autocorr.plot(eq.sims40, col = "blue")

# HASIL ESTIMASI PARAMETER
Gamma4 <- eq.sims40[["mean"]][["Gamma"]]
Lambda4 <- eq.sims40[["mean"]][["lambda"]]
delta4 <- Gamma4
for (i in 1:n)
  delta4 <- delta4%*%Gamma4
Delta4 <- delta4[1,]

Gamma4
Lambda4
Delta4

# BANGKITAN PHMM TELAH DICOCOKAN
phmm <- dthmm(aceh$Count,
              Gamma4,
              Delta4,
              "pois",
              list(lambda = Lambda4),
              discrete = TRUE)
n = length(aceh$Count)
x <- simulate(phmm, nsim = n, seed = 10)
y <- residuals(x)
w <- aceh$Count + y
DataDuga4MCMC = ifelse(w <= 0, 0, round(w))
DataDuga4MCMC

# TES AKURASI PHMM YANG TELAH DICOCOKAN
rmse_mcmc <- sqrt(mean((aceh$Count - DataDuga4MCMC)^2))
rmse_mcmc

abs <- abs(y)
Totalresiduals <- sum(abs)
MAE4 <- Totalresiduals/n
MAE4