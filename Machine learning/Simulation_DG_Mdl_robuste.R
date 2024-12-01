library(tidyverse)
library(extraDistr)
library(quantreg)
library(SparseM)
library(boot)

set.seed(123)
# Simulation n couples de v.a (Xi,Yi) 
# Yi = B1 + B2*Xi + Ei

#  Loi normale
n=1000
b1<-1
b2<-2
Xi<- runif(n,min=-1,max=1)
ein<- rnorm(n, mean=0, sd=0.8)
yi<-b1 +b2*Xi +ein

plot (Xi,yi)
abline(b1,b2,col = "red", lwd = 2)
# Ajuster le modèle de régression quantile
model <- rq(yi ~ Xi, tau = 0.5, method = "br")

# Afficher les résultats
summary(model)


# loi laplace 

n=1000
b1<-1
b2<-2
Xi_pl<- runif(n,min=-1,max=1)
ei_pl<-rlaplace(n, 0, 0.8)
yi_pl<-b1 +b2*Xi_pl +ei_pl

plot(Xi_pl,yi_pl)
abline(b1,b2,col = "blue", lwd = 2)

# Ajuster le modèle de régression quantile
model_pl <- rq(yi_pl ~ Xi_pl, tau = 0.5, method = "br")

# Afficher les résultats
summary(model_pl)





# Tracer les données
plot(Xi_pl, yi_pl, pch = 19, col = "black", main = "Régression quantile avec IC", xlab = "Xi", ylab = "yi")
abline(model_pl, col = "red", lwd = 2)  # Droite ajustée

# Tracer les droites pour les bornes des IC
abline(a = 0.91766, b = 1.89954, col = "green", lwd = 2, lty = 2)  # Borne inférieure
abline(a = 1.01071, b = 2.04745, col = "purple", lwd = 2, lty = 2)  # Borne supérieure

# Légende
legend("topleft", legend = c("Solution centrale", "IC inférieur", "IC supérieur"),
       col = c("red", "green", "purple"), lty = c(1, 2, 2), lwd = 2)


# Tracer les données
plot(Xi_pl, yi_pl, pch = 19, col = "black", main = "Régression quantile avec IC", xlab = "Xi", ylab = "yi")
abline(model_pl, col = "red")  # Droite ajustée centrale

# Tracer les droites pour les bornes des IC
abline(a = 0.91766, b = 1.89954, col = rgb(0, 1, 0, 0.5), lwd = 2,lty=1)  # Borne inférieure en vert avec transparence
abline(a = 1.01071, b = 2.04745, col = rgb(0.5, 0, 1, 0.5), lwd = 2,lty=1)  # Borne supérieure en violet avec transparence

# Ajouter une légende
legend("topleft", 
       legend = c("Solution centrale", "IC inférieur", "IC supérieur"),
       col = c("red", rgb(0, 1, 0, 0.5), rgb(0.5, 0, 1, 0.5)), 
       lty = c(1, 2, 2), 
       lwd = 2)






#Determiner les estimateurs


