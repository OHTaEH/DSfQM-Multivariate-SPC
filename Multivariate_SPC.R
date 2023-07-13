library(MASS)
rm(list=ls())

############## First Multivariate SPC ##############
# mean vector
mean_v <- c(0, 0, 0)
mean_v2 <- c(1, 2, 1)
# covariance matrix
cov_v <- c(1, 0.9, 0.9, 
           0.9, 1, 0.9,
           0.9, 0.9, 1)
cov_m <- matrix(cov_v, nrow=3, byrow=TRUE)

# phase1 dataset
set.seed(221013)
data_phase1 <- mvrnorm(n = 50, mean_v, cov_m)
alpha_CL <- 0.005

# phase2 dataset
data_phase2_1 <- mvrnorm(n = 20, mean_v, cov_m) #phase2 first mean vector (0,0,0)
data_phase2_2 <- mvrnorm(n = 20, mean_v2, cov_m) #phase2 first mean vector (1,2,1)

# math
T2_all <- rep(0, nrow(data_phase1))

for(i in 1:nrow(data_phase1)){
    T2_all[i] = t(data_phase1[i, ]-mean_v) %*% solve(cov_m) %*% (data_phase1[i, ]-mean_v)
}

# UCL
alpha_CL <- 0.005
m <- nrow(data_phase1)
p <- ncol(data_phase1)
UCL_beta <- ((m-1)^2/m) * qbeta((1-alpha_CL), p/2, (m-p-1)/2)
UCL_beta <- round(UCL_beta, 3)

# phase1 plot
par=(mai= c(1,1,1,1))
plot <- plot(T2_all, ylab=expression(T^2), pch=19, type="o", ylim = c(0, 13))
abline(h=UCL_beta, lty=3)

# phase1 + phase2 dataset
data_all <- rbind(data_phase1, data_phase2_1, data_phase2_2) # rbind
T2_all <- rep(0, nrow(data_all))

## math
for(i in 1:nrow(data_all)){
    T2_all[i] = t(data_all[i, ]-mean_v) %*% solve(cov_m) %*% (data_all[i, ]-mean_v)
}
T2_all

## phase 1+2 plot
par=(mai= c(1,1,1,1))
plot <- plot(T2_all, ylab=expression(T^2), pch=19, type="o", ylim = c(0, 20))
abline(h=UCL_beta, lty=3)

## phase 1+2 plot : length of alarm
length(which(T2_all>UCL_beta))

#########################################
#Repeat 200

data_phase2_1 <- vector(mode = "list", length = 200)

data_phase2_2 <- vector(mode = "list", length = 200)

data_all <- vector(mode = "list", length = 200)

T2_all <- vector(mode = "list", length = 200)

RL <-vector(mode = "list", length = 200)
RL
for(i in 1:200){
    data_phase2_1[[i]] <- mvrnorm(n = 20, mean_v, cov_m)
    data_phase2_2[[i]] <- mvrnorm(n = 20, mean_v2, cov_m)
    data_all[[i]] <- rbind(data_phase1, data_phase2_1[[i]], data_phase2_2[[i]])# rbind
    T2_all[[i]] <- rep(0, nrow(data_all[[i]]))
    for(j in 1:nrow(data_all[[i]])){
        T2_all[[i]][j] = t(data_all[[i]][j, ]-mean_v) %*% solve(cov_m) %*% (data_all[[i]][j, ]-mean_v)
    RL[[i]] <- which(T2_all[[i]]>UCL_beta)
    }
}
which(T2_all[[3]] > UCL_beta)
RL <- unlist(RL)
RL 
ARL <- mean(RL)
print(ARL)


