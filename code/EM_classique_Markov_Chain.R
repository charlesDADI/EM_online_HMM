################################################
#Markov Chain 
###############################################
MarkovChain=function(N,q,mu1,mu2,sigma){
  X=NULL
  Y=NULL
  #V=NULL
  X<-matrix(nrow=N,ncol=1)
  Y<-matrix(nrow=N,ncol=1)
  #V<-matrix(nrow=N,ncol=1)
  X[1,1]=mu2
  #V=rnorm(N,0,sigma)
  u=runif(1)
  for(i in 2:N){  
    if(X[i-1,1]==mu1){
                if(u<q[1,1]){X[i,1]=mu1}else X[i,1]=mu2}
    if(X[i-1,1]==mu2){
                if(u<q[2,2]){X[i,1]=mu2}else X[i,1]=mu1}
  }
  Y=X+rnorm(N,0,sigma)
  Y
}

Q<-matrix(data=c(0.7,0.5,0.3,0.5),nrow=2,ncol=2)# to define matrix Q
mu=matrix(data=c(2,0),nrow=1,ncol=2) # to define vector mu 

Y=MarkovChain(N=1000,q=Q,mu1=mu[1],mu2=mu[2],sigma=.5)
par(mfrow=c(1,2))
ts.plot(Y)
hist(Y,col="lightblue",breaks=100)

##################################
# Algo EM  for HMM
#################################
#Définition et initialisation des parametres :
t1=proc.time()# horloge au début de l'algo.

K = 100 # on fera K itérations de l’algorithme
Q_M_STEP<-matrix(data=0,nrow=K,ncol=4)#to stock value of q_ij in (M)step
SIGMA_M_STEP<-matrix(data=0,nrow=K,ncol=2)#To stock sigma in (M) step
MU_M_STEP<-matrix(data=0,nrow=K,ncol=2)#to stock value of mu_k in (M)step
#Initialisation
MU_M_STEP[1,1] = -0.5 #à modfier pour voir les changements
MU_M_STEP[1,2] = 0.5
SIGMA_M_STEP[1,1] = 2
SIGMA_M_STEP[1,2]= 2  #à modfier pour voir les changements
#Q_M_STEP[1,]=c(0.7,0.5,0.3,0.5)
Q_M_STEP[1,]=c(0.7,0.3,0.5,0.5)

i=2
for (i in 2:K) {

  vrais11_ =  Q_M_STEP[i-1,1]*dnorm(Y,mean=MU_M_STEP[i-1,1],sd=SIGMA_M_STEP[i-1,1])
  vrais12_ =  Q_M_STEP[i-1,2]*dnorm(Y,mean=MU_M_STEP[i-1,1],sd=SIGMA_M_STEP[i-1,1])
  vrais22_ = Q_M_STEP[i-1,3]*dnorm(Y,mean=MU_M_STEP[i-1,2],sd=SIGMA_M_STEP[i-1,2])
  vrais21_ = Q_M_STEP[i-1,4]*dnorm(Y,mean=MU_M_STEP[i-1,2],sd=SIGMA_M_STEP[i-1,2])

  denominateur_=vrais11_+vrais12_+vrais21_+vrais22_ 
  
  vrais12 = vrais11_ / denominateur_ # probas a posteriori p_{i,1,1}
  vrais22 = vrais22_ / denominateur_ # probas a posteriori p_{i,2}
  
  ## Mise a jour de lambda1 = P(Z=1 | X,Thteta) :
  #Q_M_STEP[i,1] = mean(vrais12)         #q_11
  #Q_M_STEP[i,2] = 1-Q_M_STEP[i-1,1]     #q_12
  #Q_M_STEP[i,3] = mean(vrais22)         #q_21
  #Q_M_STEP[i,4] = 1-Q_M_STEP[i-1,4]     #q_22
  
  Q_M_STEP[i,1] = mean(vrais12)         #q_11
  Q_M_STEP[i,2] = 1-Q_M_STEP[i,1]       #q_12
  Q_M_STEP[i,3] = mean(vrais22)         #q_21
  Q_M_STEP[i,4] = 1-Q_M_STEP[i,3]       #q_22  
  
  ## Mise  à jour de mu1 et mu2 :
  MU_M_STEP[i,1]= sum(vrais12*Y)/sum(vrais12) #mu_1
  MU_M_STEP[i,2]= sum(vrais22*Y)/sum(vrais22) #mu_2
  
  ## Mise a jour de sigma1 et sigma2 :
  SIGMA_M_STEP[i,1] = sqrt(sum(vrais12*(Y-MU_M_STEP[i,1])^2)/(sum(vrais12)))
  SIGMA_M_STEP[i,2] = sqrt(sum(vrais22*(Y-MU_M_STEP[i,2])^2)/(sum(vrais22)))
}
par(mfrow=c(2,2))
plot(1:K, SIGMA_M_STEP[,1],  pch = 16, type='o', col = "red", cex = 0.6)
plot(1:K, Q_M_STEP[,1],pch = 16, type='o', col = "red",cex = 0.6)
plot(1:K, MU_M_STEP[,2],pch = 16, type='o', col = "red",cex = 0.6)

t2=proc.time()-t1 # calcul le temps pour effecuter algo EM 
t2

resultats = data.frame(SIGMA_M_STEP,MU_M_STEP, Q_M_STEP)
names(resultats) = c("sigma1", "sigma2", "mu1", "mu2", "q11", "q12", "q21", "q22")
View(resultats)