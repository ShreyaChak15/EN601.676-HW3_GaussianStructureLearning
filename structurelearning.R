library(glasso)
library(pcalg)
library(network)
library(ggplot2)
library(igraph)
library(network)

gm <- read.table("data.txt", header = TRUE)
V <- colnames(gm)

#####PC ALGORITHM
a <-c(0.005, 0.001, 0.05, 0.01, 0.03, 0.5, 0.1)
for(j in a) {
  pc.fit <- pc(suffStat = list(C = cor(gm), n = nrow(gm)),
               indepTest = gaussCItest,
               alpha=j, labels = V, verbose = FALSE)
  #plot(pc.fit, main = paste("Estimated CPDAG, alpha= ", j))
}
best_alpha=0.01

###GLASSO ALGORITHM
#Output True Network
trueA <- ifelse(gm!=0 & row(gm)!=col(gm),1,0)
g <- network(trueA, directed=FALSE)
plot(gm,label=V,main="True network")

S <- cov(gm)
n <-nrow(gm)
nr  <- 100
#Change tuning parameter over a range of values
rho <- seq(0.1,1,length=nr)
bic <- rho
#Calculating BIC score for a range of values
for(j in 1:nr){
  a       <- glasso(S,rho[j])
  p_off_d <- sum(a$wi!=0 & col(S)<row(S))
  bic[j]  <- -2*(a$loglik) + p_off_d*log(n)
}

best <- which.min(bic)
plot(rho,bic)

#Plot estimated graph with best BIC score 
a <- glasso(S,rho[best])
best_rho<-rho[best]
#Covariance matrx
#round(a$wi,2)
P <- a$wi
A <- ifelse(P!=0 & row(P)!=col(P),1,0)
g <- network(A, directed=FALSE)
plot(g,label=V,main="Estimated network")

library(mice)
library(VIM)

gm_mis <- read.table("data_missing.txt", header=TRUE)
V_mis <- colnames(gm_mis)

#histogram for missing values
#mice_plot <- aggr(gm_mis, col=c('navyblue','yellow'),
#                  numbers=TRUE, sortVars=TRUE,
#                  labels=V, cex.axis=.7,
#                  gap=3, ylab=c("Missing data","Pattern"))
n <- ncol(gm_mis)
n_est<-10

#multiple imputation using PMM
adj_mat_pc <- array(data=0, dim=c(n_est, n, n))
adj_mat_gl <- array(data=0, dim=c(n_est, n, n))
imputed_Data <- mice(gm_mis, maxit= 10, m=n_est, method = 'pmm', seed = 500, verbose = FALSE)
for (i in 1:n_est){
  completeData <- complete(imputed_Data, action = i)
  pc.fit <- pc(suffStat = list(C = cor(completeData), n = nrow(completeData)),
               indepTest = gaussCItest,
               alpha=best_alpha, labels = V_mis, verbose = FALSE)
  gl <- glasso(cov(completeData), rho = best_rho)
  A <- ifelse(gl$wi!=0 & row(gl$wi)!=col(gl$wi),1,0)
  adj_mat_gl[i,,] <- as(A, "matrix")
  adj_mat_pc[i,,] <- as(pc.fit@graph, "matrix")
}

check_common <-function(n, n_est, adj_mat){
  M <-matrix(0L, nrow=n, ncol=n)
  #print(M)
  #print(adj_mat[1,,])
  for (j in 1:n){
    s <-0
    for (k in 1:n){
      for (i in 1:n_est){
        s <- s + adj_mat[i, j, k]
      }
      
      if (s > n_est/2){
        M[j, k] = adj_mat[i, j, k]
      }
    }
  }
  return(M)
}

final_pc<-check_common(n, n_est, adj_mat_gl)
final_gl<-check_common(n, n_est, adj_mat_pc)

g <- network(final_gl, directed=FALSE)
plot(g,label=V,main="Estimated Undirected Graph after Multiple Imputation")
g1 <- network(final_pc, interactive ='T')
plot(g1,label=V,main="Estimated CPDAG Graph after Multiple Imputation")
