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

tau_g<-as.numeric(integrate(gauss2,-Inf,Inf)[1])
tau_g
tau_t<-as.numeric(integrate(triangle2,-Inf,Inf)[1])
tau_t
tau_r<-as.numeric(integrate(rectangle2,-Inf,Inf)[1])
tau_r
tau_e<-as.numeric(integrate(epanechnikov2,-Inf,Inf)[1])
tau_e

# Intervalle de confiance assymptoptique

NW()






