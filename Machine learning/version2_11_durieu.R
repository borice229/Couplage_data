library(doParallel)
library(tidyverse)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

regressiondata <- function(n) 
{ 
  X <- sort(runif(n)) 
  m <- (sin(2*pi*X^3))^3 
  e <- rnorm(n,0,sqrt(0.02)) 
  y <- m+e 
  result <- list(x=X,y=y,m=m,e=e) 
} 

vraifonction <- function(x) 
{ 
  res <- (sin(2*pi*x^3))^3 
  return(res) 
} 

n <- 1000 
data <- regressiondata(n) 
plot(data$x,data$y) 
curve(vraifonction,add=T,col="red") 

# Noyaux utilisés 
rectangle <- function(x) 
{ 
  0.5*(abs(x)<=1) 
} 

triangle <- function(x) 
{ 
  (1-abs(x))*(abs(x)<=1) 
} 

gauss <- function(x) 
{ 
  (1/(2*pi)^0.5)*exp(-x^2/2) 
} 

epanechnikov <- function(x)
{
  ((3/4)*(1-x^2))*(abs(x)<=1)
}

# Estimateur de Nadaraya-Watson 
NW <- function(x,X,Y,h,noy="G") 
{ 
  if (noy=="G") 
    noyau=gauss 
  if (noy=="R") 
    noyau=rectangle 
  if (noy=="T") 
    noyau=triangle 
  if (noy=="E") 
    noyau=epanechnikov 
  dXh <- (x-X)/h 
  num <- sum( noyau(dXh)*Y ) 
  denom <- sum( noyau(dXh) ) 
  num/denom 
  return(num/denom) 
} 
plot(data$x,data$y) 
curve(vraifonction,add=T,col="red") 

# Calcul pour différents h (h=0.5, h=0.1 et h=0.01)
h <- 0.5
resu1 <- NULL 
for(i in 1: n) 
{ 
  print(i) 
  resu1[i] = NW(data$x[i],data$x,data$y,h) 
} 
lines(data$x,resu1,col="blue",lty=1) 

h <- 0.1
resu2 <- NULL
for(i in 1: n) 
{ 
  print(i) 
  resu2[i] = NW(data$x[i],data$x,data$y,h) 
} 
lines(data$x,resu2,col="green",lty=2) 

h <- 0.01
resu3 <- NULL
for(i in 1: n) 
{ 
  print(i) 
  resu3[i] = NW(data$x[i],data$x,data$y,h) 
} 
lines(data$x,resu3,col="pink") 

Mat <- cbind(data$x,data$y,resu1,resu2,resu3)

# Représentation graphique
# Choix de la taille de la fenêtre h par validation croisée 

CV <- function(X,Y,h,noy="G") 
{ 
  E <- NULL 
  for (i in 1:n) 
  { 
    E[i] <- (Y[i]- NW(X[i],X[-i],Y[-i],h,noy))^2 
  } 
  return(mean(E)) 
} 

hCV <- function(X,Y,noy="G") 
{ 
  CV1 <- function(h) 
  {
    CV(X,Y,h,noy="G") 
    #print(CV(X,Y,h,noy="G"))
  } 
  #grilleh=1:200/1000 
  E <- sapply(grilleh,CV1) 
  #E <- mclapply(grilleh,CV1,mc.cores=no_cores)
  return(list(h=grilleh[which.min(E)],erreur=E)) 
}

n <- length(data$x)
grilleh <- 1:200/1000 
T1<-Sys.time() 
hopt <- hCV(data$x,data$y) 
T2<-Sys.time() 
Tdiff <- difftime(T2, T1) 

plot(grilleh,as.numeric(hopt$erreur))
y <- hopt$erreur[as.numeric(hopt$erreur)==min(as.numeric(hopt$erreur))] 
abline(v=hopt$h,h=y,col="blue") 
# Choix de h par la méthode de la validation croisée 
hopt$h 

# Calcul de l'estimateur de Nadaraya-Watson pour chaque x et choix de h 
# par la méthode de la validation croisée
h <- hopt$h 
resu <-  rep(0,n) 
for(i in 1: n) 
{ 
  print(i) 
  resu[i] <- NW(data$x[i],data$x,data$y,h) 
} 
plot(data$x,data$y) 
curve(vraifonction,col="blue",add=T) 
lines(data$x,resu,col="red")



# Noyaux utilisés 
rectangle2 <- function(x) 
{ 
  (0.5*(abs(x)<=1))^2
} 

triangle2 <- function(x) 
{ 
  ((1-abs(x))*(abs(x)<=1))^2 
} 

gauss2 <- function(x) 
{ 
  ((1/(2*pi)^0.5)*exp(-x^2/2))^2 
} 

epanechnikov2 <- function(x)
{
  (((3/4)*(1-x^2))*(abs(x)<=1))^2
}

#tau

tau_g<-as.numeric(integrate(gauss2,-Inf,Inf)[1])
tau_g
tau_t<-as.numeric(integrate(triangle2,-Inf,Inf)[1])
tau_t
tau_r<-as.numeric(integrate(rectangle2,-Inf,Inf)[1])
tau_r
tau_e<-as.numeric(integrate(epanechnikov2,-Inf,Inf)[1])
tau_e

# Estimateur de la densité de x 
dens <- function(x,X,h,noy="G") 
{
  n <- length(data$x)
  if (noy=="G")
    noyau=gauss
  if (noy=="R")
    noyau=rectangle
  if (noy=="T")
    noyau=triangle
  dXh <- (x-X)/h
  num <- sum( noyau(dXh))
  res <- num/(n*h)
  return(res) 
} 


# Estimateur de la variance conditionnelle var(Y/X=x) 
varcond <- function(x,X,Y,h,noy="G") 
{
  n <- length(data$x)
  if (noy=="G")
    noyau=gauss
  if (noy=="R")
    noyau=rectangle
  if (noy=="T")
    noyau=triangle
  dXh <- (x-X)/h
  num <- sum( noyau(dXh)*(Y-NW(x,X,Y,h,noy))^2)
  res <- num/dens(x,X,h,noy)
  return(res) 
} 


# Intervalle de confiance assymptoptique

# Niveau de confiance

h <- 0.01
resu3 <- NULL
for(i in 1: n) 
{ 
  print(i) 
  resu3[i] = NW(data$x[i],data$x,data$y,h) 
} 
lines(data$x,resu3,col="pink") 

Mat <- cbind(data$x,data$y,resu1,resu2,resu3)

h<-0.01
n<-1000
noy<-"G"
alpha <- 0.05
z_alpha <- qnorm(1 - alpha / 2)
IC_lower<-NULL
IC_upper<- NULL
estNW<-NULL
sigma2<-NULL
erreur_standard<-NULL
for(i in 1: n) 
{ 
sigma2[i]<- (varcond(data$x[i],data$x,data$y,h,noy=noy)*tau_g^2)/dens(data$x[i],data$x,h,noy=noy)
erreur_standard[i]<- sqrt(sigma2[i]/n*h)
estNW[i]<- NW(data$x[i],data$x,data$y,h,noy)
IC_lower[i] <- estNW[i] - z_alpha * erreur_standard[i]
IC_upper[i] <- estNW[i] + z_alpha * erreur_standard[i]
}
MatInt <- cbind(data$x,data$y,estNW,IC_lower,IC_upper)
colnames(MatInt) <- c("x", "y", "estimationNW", "IC_lower", "IC_upper")

library(ggplot2)

# Supposons que MatInt est déjà défini avec les bonnes colonnes
colnames(MatInt) <- c("x", "y", "estimation", "IC_lower", "IC_upper")

# Convertir la matrice MatInt en DataFrame pour ggplot2
MatInt_df <- as.data.frame(MatInt)

# Tracer
ggplot(MatInt_df, aes(x = x)) +
  # Nuages de points (x, y)
  geom_point(aes(y = y, color = "Observations"), size = 2) +
  # Courbe pour l'estimateur
  geom_line(aes(y = estimation, color = "Estimateur"), size = 1) +
  # Ruban pour les intervalles de confiance
  geom_ribbon(aes(ymin = IC_lower, ymax = IC_upper, fill = "Intervalle de confiance"), alpha = 0.3) +
  # Courbes pour IC_lower et IC_upper
  geom_line(aes(y = IC_lower, color = "IC_lower"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = IC_upper, color = "IC_upper"), linetype = "dashed", size = 0.8) +
  # Labels et thème
  labs(
    title = "Estimateur avec Intervalle de Confiance et Observations",
    x = "x",
    y = "Valeur"
  ) +
  # Personnalisation des couleurs et de la légende
  scale_color_manual(
    name = "Légende",
    values = c(
      "Observations" = "blue",
      "Estimateur" = "red",
      "IC_lower" = "green",
      "IC_upper" = "purple"
    )
  ) +
  scale_fill_manual(
    name = "Intervalle",
    values = c("Intervalle de confiance" = "gray")
  ) +
  theme_minimal()






