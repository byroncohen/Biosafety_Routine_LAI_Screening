##PaperOne_4_2_2024_freq7##
####Install/load Packages####

library(tidyverse)
library(purrr)
library(statnet)
library(reshape2)


#Create the Contact Matrix Plotting Function using geometry
plot.cm= function(CM){
  N= dim(CM)[1] #obtain number of rows in CM
  theta= seq(0,2*pi, length= N+1)
  x= cos(theta[1:N])
  y= sin(theta[1:N])
  symbols(x,y, fg=0, circles = rep(1,N),
          inches=0.1, bg=1, xlab="", ylab="")
  segx1= as.vector(matrix(x, ncol=length(x),
                          nrow= length(x),byrow=TRUE))
  segx2= as.vector(matrix(x, ncol=length(x),
                          nrow= length(x),byrow=FALSE))
  segy1= as.vector(matrix(y, ncol=length(x),
                          nrow= length(x),byrow=TRUE))
  segy2= as.vector(matrix(y, ncol=length(x),
                          nrow= length(x),byrow=FALSE))
  segments(segx1, segy1, segx2, segy2, lty = as.vector(CM))
}


####Create a summary table function that also provides a plotting option####
summary.cm = function(x, plot=FALSE){ #x is the contact matrix. The default plot argument is that plot=FALSE
  x= table(apply(x,2,sum)) #apply sum function over the columns of the x matrix, present in table
  res=data.frame(n=x) #store result as dataframe
  names(res)= c("degree","freq") #label column names of res
  if (plot) #if plot argument is specified as plot=TRUE...
    barplot(x,xlab="degree") #create a barplot
  return(res) #return the result
}

##BANetFunction

#The Watts-Strogatz model can have at most Poisson-like variance in degree distribution (mean degree= degree variance). 
#As a result, it cannot mimic heavy-tailed distributions that are seen in "scale-free networks". 
#Scale-free networks are those whose degree distribution follows a power law. 
#Scale-free networks are characterized by preferential attachment ("rich get richer") dynamics, where certain nodes acquire 
#vastly more connections than most other nodes, acting as "hubs". Examples of scale-free networks include 
#certain social networks (such as academic citation networks) as well as computer networks (such as the hyperlink
#connections between different pages on the internet).Barabasi-Albert Networks can account for 
#preferential attachment in scale-free networks.

BarabasiAlbert= function(N,K){
  CM= matrix(0, ncol = N, nrow=N) #create empty CM
  CM[1,2]=1 #assign row 1, column
  CM[2,1]=1  
  for(i in 3:N){
    probs=apply(CM,1,sum) #sum probabilities across the rows of the contract matrix
    link= unique(sample(c(1:N)[-i], #extract unique elements from vector 1:N[-i]
                        size = min(c(K, i-1)), prob=probs[-i]))
    CM[i,link]=CM[link,i]=1
  }
  class(CM)="cm"
  return(CM)
}


#To visualize networks more clearly, we can use the plotting functions of the Statnet package.
#The network function of the Statnet package converts the contact matrix (currently in class CM) to a network class object.
#install.packages("statnet")

BAnetwork = BarabasiAlbert(100,1)
plot(network(BAnetwork, directed=FALSE)) #create statnet network plot
plot(BAnetwork)#create regular network plot


##Create Network
cm6=BarabasiAlbert(N=1000,K=4) #mean degree=8


##Summary func TripleSEIIIRiso 

#Let's define a summary function for the TripleSEIIIRO class.
summary.TripleSEIIIRiso= function(x){
  t=dim(x$S)[2]# get number of columns in x$S
  SL=apply(x$SL,2,sum) # S= sum of x$SL over all columns
  SG=apply(x$SG,2,sum) # S= sum of x$SG over all columns
  ENL=apply(x$ENL,2,sum) # ENL= sum of x$ENL over all columns
  ENG=apply(x$ENG,2,sum) # ENG= sum of x$ENG over all columns
  ENISO=apply(x$ENISO,2,sum) #ENISO= sum of x$ENISO over all columns
  IPL=apply(x$IPL,2,sum) # IPL= sum of x$IPL over all columns
  IPG=apply(x$IPG,2,sum) # IPG= sum of x$IPG over all columns
  IPISO=apply(x$IPISO, 2, sum) # IPISO= sum of x$IPISO over all columns
  IAL= apply(x$IAL, 2, sum) #IAL= sum of x#IAL over all columns
  IAG= apply(x$IAG, 2, sum) #IAG= sum of x#IAG over all columns
  ISL=apply(x$ISL,2,sum) # ISL= sum of x$ISL over all columns
  ISG=apply(x$ISG,2,sum) # ISG= sum of x$ISG over all columns
  ISISO=apply(x$ISISO, 2, sum) #ISISO= sum of x$ISISO over all columns
  EVISO= apply(x$EVISO,2, sum) #EVISO= sum of x$EVISO over all columns
  LWEVINF= apply(x$LWEVINF,2,sum) #LWEVINF= sum of x$LWEVINF over all columns
  GPEVINF= apply(x$GPEVINF,2,sum) #GPEVINF= sum of x$GPEVINF over all columns
  Incidence= apply(x$Incidence,2,sum) #Incidence= sum of x$Incidence over all columns
  R=apply(x$R,2,sum) # RL= sum of x$RL over all columns
  freq= mean(x$freq)
  sigma= mean(x$sigma)
  isolationdelay= mean(x$isolationdelay)
  res= data.frame(SL=SL, SG=SG, ENL=ENL, ENG=ENG, ENISO=ENISO, IPL=IPL, IPG=IPG, IPISO=IPISO, IAL=IAL, IAG=IAG, ISL=ISL, ISG=ISG, ISISO=ISISO, EVISO=EVISO, LWEVINF= LWEVINF,GPEVINF=GPEVINF, Incidence=Incidence, R=R, freq=freq, sigma=sigma, isolationdelay=isolationdelay) #save results in a dataframe
  res$Sum = SL + SG + ENL+ ENG+ ENISO+ IPL+IPG+IPISO+IAL+IAG+ISL+ISG+ISISO+R
  res$Time = rownames(res)
  res$Outbreak100 <- ifelse(max(res$R)>=100,"yes","no")
  res$Outbreak50 <- ifelse(max(res$R)>=50,"yes","no")
  res$Outbreak10<- ifelse(max(res$R)>=10,"yes","no")
  res$everinfected<- max(res$R)
  res$peakNinfected<- max(res$ENL+ res$ENG +res$ENISO + res$IPL + res$IPG + res$IPISO+ res$IAL +res$IAG+ res$ISL + res$ISG + res$ISISO)
  res$peakNinfectiousnotiso<- max(res$IPL + res$IPG + res$IAL +res$IAG+ res$ISL + res$ISG)
  res$shareinfectedisolated<- max(res$EVISO)/max(res$R)
  res$Ninfectiousnotiso<- res$IPL + res$IPG + res$IAL +res$IAG+ res$ISL + res$ISG
  res$PTinfectiousnotiso<- sum(res$Ninfectiousnotiso)
  res$Ninfectiousiso<- res$IPISO+res$ISISO
  res$PTinfectiousiso<- sum(res$Ninfectiousiso)
  res$PTallinfectious <- res$PTinfectiousiso+ res$PTinfectiousnotiso 
  res$shareinfectiousPTiso<-  res$PTinfectiousiso/ res$PTallinfectious 
  return(res)
}

##Plot funct TripleSEIIIRiso

#Let's define a plot function for the tripleSEIIIRiso class
plot.TripleSEIIIRiso= function(x){
  y=summary(x)
  plot(y$SL, type="l",xlab="time", ylab="", ylim = c(0,400))
  lines(y$SG, type="l", col="brown")
  lines(y$ENL, type="l", col="green")
  lines(y$ENG, type="l", col="blue")
  lines(y$ENISO, type="l",col="grey")
  lines(y$IPL, type="l", col="purple")
  lines(y$IPG, type="l", col="orange")
  lines(y$IPISO, type="l", col="yellow")
  lines(y$IAL, type="l", col="#00CC99")
  lines(y$IAG, type="l", col="#99CC00")
  lines(y$ISL, type="l", col="red")
  lines(y$ISG, type="l", col="pink")
  lines(y$ISISO,type="l", col="gold")
  lines(y$EVISO, type= "l", col= "turquoise")
  lines(y$R, type = "l", col= "black")
  legend("right", legend= c("SL","SG","ENL", "ENG","ENISO", "IPL", "IPG","IPISO","ISL","ISG","IAL","IAG","ISISO","EVISO","R"),
         lty = c(1,1,1), pch= c(1,1,1),
         col= c("black", "brown","green","blue", "grey", "purple", "orange", "yellow","#00CC99","#99CC00", "red", "pink", "gold", "turquoise", "black"))
}




##CreateModelTripleSEIIIRiso: Creating a new Model with a Completely Asymptomatic 
# Infectious (IA) state, faulty isolation, partial post-symptom isolation, 
# dynamic infectiousness, dynamic test sensitivity, a dynamic testing regime, 
#and isolation delays: TripleSEIIIRiso 

########################################################################################
#The TripleSEIIIRiso function we will define below will simulate an SEIIIR epidemic with isolation on an arbitrary contact matrices with two linked networks and return an object with the class TripleSEIIIRiso
#CM= contact matrix
#tau= probability that an infection is transmitted across a given S-I edge per time step t
#gamma= prob of removal of infected node per time step
TripleSEIIIRiso= function(CM, tauL, epsilonL, thetaL, delta, gammaL, piL, tauG,thetaG,omegaG, chiG, gammaG, piG, epsilonIS, gammaIS, epsilonIAL, gammaIAL, gammaIAG, iota, psi, zeta, sigma,mu, nu, isolationdelay) { {
  #set.seed(1)
  N= dim(CM)[1] # get no. of indivs by taking # of rows in CM
  NL= 20 #100 lab workers
  NG= N-NL #400 general population
  SL= as.data.frame(matrix(c(rep(1,NL),rep(0,NG)), nrow=N, ncol = 1)) #First susceptible lab worker: create single column matrix where first 100 people are susceptible lab workers
  ENL= as.data.frame(matrix(rep(0,N), nrow=N, ncol=1)) #First exposed lab worker: create single column matrix where no one is in the exposed but not yet infectious lab worker category.
  IPL= as.data.frame(matrix(rep(0,N), nrow=N, ncol=1)) #First exposed lab worker: create single column matrix where no one is in the infectious, pre-symptomatic lab worker category.
  ISL= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) #First infectious lab worker: create single column matrix where no one is in the infectious, symptomatic lab worker category.
  SG= as.data.frame(matrix(c(rep(0,NL),rep(1,NG)),nrow=N, ncol = 1)) #First susceptible general population: create single column matrix where the last 400 people are susceptibles in the general population
  ENG= as.data.frame(matrix(rep(0,N), nrow=N, ncol=1)) #First exposed general population: create single column matrix where no one is in the exposed but not yet infectious general population category.
  IPG= as.data.frame(matrix(rep(0,N), nrow=N, ncol=1)) #First exposed general population: create single column matrix where no one is in the infectious, pre-symptomatic general population category.
  ISG= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) #First infectious general population: create single column matrix where no one is in the infectious, symptomatic general population category.
  ENISO= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) #First exposed, noninfectious isolated person:single column matrix where no one is in the exposed, non-infectious category.
  IPISO= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) #First infectious, presymptomatic isolated person:single column matrix where no one is in the infectious, pre-symptomatic isolated person
  ISISO= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) #First infectious, symptomatic isolated person:single column matrix where no one is in the infectious, symptomatic isolated person
  IAL= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) # First infectious, permanently asymptomatic lab worker: single column matrix where no lab worker is in the infectious, permanently asymptomatic category
  IAG= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) # First infectious, permanently asymptomatic general population person: single column matrix where no general population person is in the infectious, permanently asymptomatic category
  EVISO= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) # Single column matrix for people who have ever been isolated (starts at 0)
  LWEVINF= as.data.frame(matrix(rep(1,N), nrow=N, ncol = 1)) # Single column matrix for lab workers who have ever been infected (starts at 1)
  GPEVINF= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) # Single column matrix for general population who have ever been infected (starts at 0)
  Incidence= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1))  #Single column matrix for new cases at each time step. Start it at 0 for now.
  R= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) #First removed: create single column matrix where no people are recovered/dead
  ENL1= sample(1:NL, size=1) #Pick row number containing the first random exposed lab worker, name it ENL1
  ENL[ENL1,1]=1 #make ENL=1 in row ENL1 (this is assigning the first exposed lab worker)
  SL[ENL1,1]=0 #make S=0 in row ENL1 (this is accounting for the fact that the first infected person is no longer susceptible)
  Incidence[ENL1,1]=1 #make Incidence=1 in row ENL1 (this is the very first incident case)
  freq=thetaL
  t=1
  while(t<100){ #as long as t<100, and...
    print(t)
    #browser() #:only use when debugging
    #print(EVISO[1,])
    if(t==1 | sum(EVISO[,t-1])==0) { #as long as no one has been isolated, OR if you're at the initial time step (t==1)...
      tauL = tauL
      epsilonL = epsilonL 
      thetaL = thetaL 
      delta = delta
      gammaL = gammaL
      piL = piL 
      tauG = tauG
      thetaG= thetaG
      omegaG= omegaG
      chiG =chiG
      gammaG=gammaG
      piG= piG
      epsilonIS = epsilonIS
      gammaIS =gammaIS
      epsilonIAL= epsilonIAL
      gammaIAL= gammaIAL
      gammaIAG =gammaIAG
      iota= iota
      psi= psi 
      zeta= zeta
      sigma =sigma
      mu= mu
      nu= nu
      isolationdelay = isolationdelay
    }else if(sum(EVISO[,t-1])>0 & sum(EVISO[,t-1])<10){
      thetaL=1
    }else {
      chiG=1/7
      thetaL=1
    }
    #######################ACTION CODE#################################
    t=t+1 #move to the next time step
    infneighISL=CM%*%ISL[,t-1] #number of infectious, symptomatic lab workers is product of multiplying the CM matrix by ISL matrix, using all rows, and data from previous time-step
    infneighIPL= CM%*%IPL[,t-1] #number of exposed, infectious lab workers is product of multiplying the CM matrix by ISL matrix, using all rows, and data from previous time-step
    infneighISG=CM%*%ISG[,t-1] #number of infectious, symptomatic general population is product of multiplying the CM matrix by ISG matrix, using all rows, and data from previous time-step
    infneighIPG= CM%*%IPG[,t-1] #number of exposed, infectious, general population is product of multiplying the CM matrix by ISG matrix, using all rows, and data from previous time-step
    infneighIAL= CM%*%IAL[,t-1] # number of infectious, permanently asymptomatic lab workers is product of multiplying the CM matrix by the IAL matrix, using all rows, and data from previous time-step
    infneighIPISO= CM%*%IPISO[,t-1] # number of infectious, pre-symptomatic isolated people is the product of multiplying the CM matrix by the IPISO matrix, using all rows, and data from previous time-step
    infneighISISO= CM%*%ISISO[,t-1] # number of infectious, symptomatic isolated people is the product of multiplying the CM matrix by the ISISO matrix, using all rows, and data from the previous time-step
    infneighIAG= CM%*%IAG[,t-1] # number of infectious, permanently asymptomatic GPs is product of multiplying the CM matrix by the IAG matrix, using all rows, and data from previous time-step
    pinfL= 1-(1-tauL)^(infneighISL+infneighIPL*psi+ infneighISG+infneighIPG*psi + infneighIAL*zeta +infneighIAG*zeta + infneighIPISO*psi*iota + infneighISISO*iota) #define the probability of any given susceptible lab worker becoming infected
    pinfG= 1-(1-tauL)^(infneighISL+infneighIPL*psi+ infneighISG+infneighIPG*psi + infneighIAL*zeta +infneighIAG*zeta + infneighIPISO*psi*iota + infneighISISO*iota) #define the probability of any given general population person becoming infected
    
    ########################Define the Flows##############
    
    #Create the infection flows
    newENL= rbinom(N, SL[,t-1],pinfL) #random draw from previous time-step's susceptible lab workers to select the new exposed lab workers. Run SL[,t-1] trials with pinfL probability of infection for each observation in 1:N.
    newENG= rbinom(N, SG[,t-1],pinfG) #random draw from previous time-step's susceptible general population to select the new exposed general population. Run SG[,t-1] trials with pinfG probability of infection for each observation in 1:N.
    #Create empty vectors for new non-infection flows
    newENISOfromENL=rep(0,N)
    newIPL= rep(0,N)
    newIALfromENL = rep(0,N) #new
    newENISOfromENG= rep(0,N)
    newIPG= rep(0,N)
    newIAGfromENG = rep(0,N) #new
    newIPISOfromENISO= rep(0,N)
    newIPISOfromIPL= rep(0,N)
    newISL= rep(0,N)
    newIPISOfromIPG= rep(0,N)
    newISG=rep(0,N)
    newISISOfromIPISO= rep(0,N)
    newISISOfromISL= rep(0,N)
    newISISOfromIAL= rep(0,N) #new
    newRecoveredfromISL= rep(0,N)
    newISISOfromISG= rep(0,N)
    newRecoveredfromISG= rep(0,N)
    newRecoveredfromISISO=rep(0,N)
    newRecoveredfromIAL= rep(0,N) #new
    newRecoveredfromIAG= rep(0,N) #new
    
    #Transition rates out of ENL with delay
    ENLtoENISO<- 1/(1/(thetaL*sigma*mu)+isolationdelay)
    ENLtoIPL<- (1-ENLtoENISO)*epsilonL
    ENLtoIAL<- (1-ENLtoENISO)*(1-epsilonL)*epsilonIAL
    ENLstay<- (1-ENLtoENISO)*(1-epsilonL)*(1-epsilonIAL)
    
    #Code for flows out of ENL: 
    leavingENL<- sapply(1:N, function(x) rmultinom(1, ENL[x,t-1], prob = c(ENLstay,#Stay in ENL
                                                                           ENLtoIPL,#to IPL
                                                                           ENLtoIAL,#to IAL
                                                                           ENLtoENISO))) #to ENISO
    
    #extract each row of leavingENL, save as a vector
    #skip the first row of leavingENL, because those people are actually staying in ENL  
    newIPL<- leavingENL[2,]
    newIPL<- as.vector(newIPL)
    
    newIALfromENL<- leavingENL[3,]
    newIALfromENL<- as.vector(newIALfromENL)
    
    newENISOfromENL<- leavingENL[4,]
    newENISOfromENL<- as.vector(newENISOfromENL)
    
    #Transition rates out of ENG with delay
    ENGtoENISO<- 1/(1/(thetaG*sigma*mu)+isolationdelay)
    ENGtoIPG<- (1-ENGtoENISO)*epsilonL
    ENGtoIAG<- (1-ENGtoENISO)*(1-epsilonL)*epsilonIAL
    ENGstay<- (1-ENGtoENISO)*(1-epsilonL)*(1-epsilonIAL)
    
    #Code for flows out of ENG: 
    leavingENG<- sapply(1:N, function(x) rmultinom(1, ENG[x,t-1], prob = c(ENGstay,#Stay in ENG
                                                                           ENGtoIPG,#to IPG
                                                                           ENGtoIAG,#to IAG
                                                                           ENGtoENISO))) #to ENISO
    
    #extract each row of leavingENG, save as a vector
    #skip the first row of leavingENG, because those people are actually staying in ENG  
    newIPG<- leavingENG[2,]
    newIPG<- as.vector(newIPG)
    
    newIAGfromENG<- leavingENG[3,]
    newIAGfromENG<- as.vector(newIAGfromENG)
    
    newENISOfromENG<- leavingENG[4,]
    newENISOfromENG<- as.vector(newENISOfromENG)
    
    #Code for flow out of ENISO:
    newIPISOfromENISO= rbinom(N, ENISO[,t-1], epsilonIS) 
    
    #Code for flows out of IPL: 
    newIPISOfromIPL= rbinom(N, IPL[,t-1], 1/(1/(thetaL*sigma*nu)+isolationdelay))
    newISL= rbinom(N,IPL[,t-1]*(1-newIPISOfromIPL), delta)
    
    #Code for flows out of IPG: 
    newIPISOfromIPG= rbinom(N, IPG[,t-1], 1/(1/(omegaG*sigma*nu) +isolationdelay)) #select isolated people
    newISG= rbinom(N, IPG[,t-1]*(1-newIPISOfromIPG), delta) #from non-isolated people, select progressors
    
    #Code for flow out of IPISO:
    newISISOfromIPISO= rbinom(N, IPISO[,t-1], delta) #select those isolated who will progress to symptoms
    
    #Code for flows out of ISL: 
    newISISOfromISL= rbinom(N, ISL[,t-1], 1/(1/(thetaL*sigma*(1-piL) +piL)+isolationdelay)) #select the isolated
    newRecoveredfromISL= rbinom(N, ISL[,t-1]*(1-newISISOfromISL), gammaL) #select the recovered LWs who didn't get isolated
    
    #Code for flows out of ISG: 
    newISISOfromISG= rbinom(N, ISG[,t-1], 1/(1/(chiG*sigma*(1-piG)+piG) +isolationdelay)) #select the isolated
    newRecoveredfromISG= rbinom(N, ISG[,t-1]*(1-newISISOfromISG), gammaG) #select the recovered GPs who didn't get isolated    
    
    #Code for flows out of IAL:
    newISISOfromIAL= rbinom(N, IAL[,t-1], thetaL *sigma*0.8)
    newRecoveredfromIAL= rbinom(N, IAL[,t-1]*(1-newISISOfromIAL), gammaIAL)
    
    #Code for flows out of IAG:
    newRecoveredfromIAG= rbinom(N, IAG[,t-1], gammaIAG)
    
    #Code for flows out of ISISO:
    newRecoveredfromISISO= rbinom(N, ISISO[,t-1], gammaIS)
    
    #Code for counting up all the people who have been newly isolated in a given time-step
    newIsolatedAll= newENISOfromENL+ newENISOfromENG+ newIPISOfromIPL+ newIPISOfromIPG+ newISISOfromISL+ newISISOfromIAL+ newISISOfromISG
    
    #######Define the states at the next time-step#######
    
    nextSL= SL[,t-1]-newENL #next step's susceptible lab workers is the previous step's susceptible lab workers, minus the newly exposed but not infectious lab workers.
    nextSG= SG[,t-1]-newENG #next step's susceptible general population is the previous step's susceptible general population, minus the newly exposed but not infectious general population.
    nextENL= ENL[,t-1]+newENL -newIPL -newENISOfromENL- newIALfromENL #next step's exposed but not infectious lab workers is the previous step's exposed but not infectious lab workers, plus the newly exposed but not infectious lab workers, minus the newly exposed and infectious lab workers, minus the newly isolated exposed but non infectious lab workers, minus the new infected perm-asymptomatic LWs
    nextENG= ENG[,t-1]+newENG -newIPG -newENISOfromENG -newIAGfromENG #next step's exposed but not infectious general population is the previous step's exposed but not infectious general population, plus the newly exposed but not infectious general population, minus the newly exposed and infectious general population, minus the newly isolated exposed but non infectious general population, minus the new infected perm-asymptomatic GPs.
    
    nextENISO= ENISO[,t-1] + newENISOfromENL + newENISOfromENG  - newIPISOfromENISO #next step's exposed, noninfectious isolateds includes the previous step's exposed, noninfectious isolateds, plus the new exposed, noninfectious isolateds from ENL and ENG, minus the newly infectious, pre-symptomatic isolateds
    
    nextIPL= IPL[,t-1]+newIPL -newISL -newIPISOfromIPL #next step's exposed and infectious lab workers is the previous step's exposed and infectious lab workers, plus the newly exposed and infectious lab workers, minus the newly infectious and symptomatic lab workers, minus the newly isolated exposed, infectious lab workers.
    nextIPG= IPG[,t-1]+newIPG -newISG -newIPISOfromIPG #next step's exposed and infectious general population is the previous step's exposed and infectious general population, plus the newly exposed and infectious general population, minus the newly infectious and symptomatic general population, minus the newly isolated exposed, infectious general population
    nextIPISO= IPISO[,t-1]  + newIPISOfromIPL +newIPISOfromIPG +newIPISOfromENISO- newISISOfromIPISO
    #next step's infectious,pre-symptomatic isolateds (IPISO) is previous step's IPISO, plus newly isolated from infectious pre-symptomatic LWs, plus newly ioslated from infectious presymptomatic GPs, plus those who progressed from exposed, noninfectious isolated status, minus those who progress to infectious symptomatic isolated status
    nextISL= ISL[,t-1] +newISL -newRecoveredfromISL - newISISOfromISL # next step's infectious and symptomatic is previous step's infectious and symptomatic LWs, plus newly infectious and symptomatic LWs, minus the newly recovered LWs, minus the newly isolated LWs.
    nextISG= ISG[,t-1] +newISG -newRecoveredfromISG - newISISOfromISG # next step's infectious and symptomatic is previous step's infectious and symptomatic GPs, plus newly infectious and symptomatic GPs, minus the newly recovered GPs, minus the newly isolated GPs
    nextIAL= IAL[,t-1] + newIALfromENL -newRecoveredfromIAL -newISISOfromIAL # next step's infectious perm-asymptomatic LWs (IALs) is previous step's IALs, plus new IALs who were previously exposed non-infectious LWs, minues newly isolated IALs, minus newly recovered IALs
    nextIAG= IAG[,t-1] +newIAGfromENG -newRecoveredfromIAG  # next step's infectious perm-asymptomatic GPs (IAGs) is previous step's IAGs, plus new IAGs who were previously exposed non-infectious GPs, minus newly recovered IAGs
    nextISISO= ISISO[,t-1] + newISISOfromISL + newISISOfromISG +newISISOfromIAL+ newISISOfromIPISO -newRecoveredfromISISO
    # next step's infectious symptomatic isolateds (ISISO) is the previous step's ISISO, plus the infected symptomatic LWs who are isolated, plus the infected symptomatic GPs who are isolated, plus the permanently asymptomatic infected LWs who are isolated, plus the isolated who progress from the infectious presymptomatic stage, minus the newly recovered ISISO
    nextEVISO= EVISO[,t-1]+newIsolatedAll #Next step's total number of ever isolated is the previous step's number of ever isolated plus all of the newly isolated people
    nextLWEVINF= LWEVINF[,t-1] + newENL # Next steps's LWs who were ever infected is the previous step's LWs who were ever infected plus the newly infected LWs
    nextGPEVINF= GPEVINF[,t-1] + newENG # Next steps's GPs who were ever infected is the previous step's GPs who were ever infected plus the newly infected GPs
    nextstepcases= newENL + newENG #Next step's newly infected cases
    nextR= R[,t-1]+ newRecoveredfromISL + newRecoveredfromISG+ newRecoveredfromISISO+ newRecoveredfromIAL +newRecoveredfromIAG   #next step's removed people is previous step's removed, plus the newly recovered infectious symptomatic LWs, plus the newly recovered infectious symptomatic GPs, plus the newly recovered infectious symptomatic LWs, plus the newly recovered isolated infectious symptomatic GPs, plus the newly recovered infectious but asymptomatic LWs, plus the newly recovered infectious but asymptomatic GPs
    SL=cbind(SL, nextSL) #bind SL, nextSL matrices together
    SG=cbind(SG, nextSG) #bind SG, nextSG matrices together
    ENL=cbind(ENL, nextENL) #bind ENL, nextENL matrices together
    ENG=cbind(ENG, nextENG) #bind ENG, nextENG matrices together
    ENISO= cbind(ENISO, nextENISO) #bind ENISO, nextENISO matrices together
    IPL=cbind(IPL, nextIPL) #bind IPL, nextIPL matrices together
    IPG=cbind(IPG, nextIPG) #bind IPG, nextIPG matrices together
    IPISO= cbind(IPISO, nextIPISO) #bind IPISO, nextIPISO together
    IAL= cbind(IAL, nextIAL) #bind IAL, nextIAL together
    IAG= cbind(IAG, nextIAG) #bind IAG, nextIAG together
    ISL=cbind(ISL, nextISL) #bind ISL, nextISL matrices together
    ISG=cbind(ISG, nextISG) #bind ISG, nextISG matrices together
    ISISO= cbind(ISISO, nextISISO) #bind ISISO, nextISISO together
    EVISO= cbind(EVISO, nextEVISO) #bind EVISO, nextEVISO together
    LWEVINF= cbind(LWEVINF, nextLWEVINF) #bind LWEVINF, nextLWEVINF together
    GPEVINF= cbind(GPEVINF, nextGPEVINF) #bind GPEVINF, nextGPEVINF together
    Incidence= cbind(Incidence, nextstepcases)
    R= cbind(R, nextR) #bind R, nextR matrices together
  }
}
  freq= rep(freq, t)
  sigma= rep(sigma, t)
  isolationdelay= rep(sigma, t)
  res=list(SL=SL, SG=SG, ENL=ENL, ENG=ENG, ENISO=ENISO, IPL=IPL, IPG=IPG, IPISO=IPISO, IAL=IAL, IAG=IAG, ISL=ISL, ISG=ISG, ISISO=ISISO, EVISO=EVISO, LWEVINF= LWEVINF,GPEVINF=GPEVINF,Incidence=Incidence, R=R,freq= freq, sigma=sigma, isolationdelay=isolationdelay) #save results in a list of matrices SL,SG,ENL,ENG,ENISO, IPL,IPG, IPISO, IAL, IAG, ISL,ISG, ISISO, R
  class(res)="TripleSEIIIRiso" #set class of results as "TripleSEIIIRiso"
  return(res) #return results
}


##CreateFuncMinimalSaving

########################################################################################
#The TripleSEIIIRiso function we will define below will simulate an SEIIIR epidemic with isolation on an arbitrary contact matrices with two linked networks and return an object with the class TripleSEIIIRiso
#CM= contact matrix
#tau= probability that an infection is transmitted across a given S-I edge per time step t
#gamma= prob of removal of infected node per time step
TripleSEIIIRiso2= function(CM, tauL, epsilonL, thetaL, delta, gammaL, piL, tauG,thetaG,omegaG, chiG, gammaG, piG, epsilonIS, gammaIS, epsilonIAL, gammaIAL, gammaIAG, iota, psi, zeta, sigma,mu, nu, isolationdelay) { {
  #set.seed(1)
  N= dim(CM)[1] # get no. of indivs by taking # of rows in CM
  NL= 20 #100 lab workers
  NG= N-NL #400 general population
  SL= as.data.frame(matrix(c(rep(1,NL),rep(0,NG)), nrow=N, ncol = 1)) #First susceptible lab worker: create single column matrix where first 100 people are susceptible lab workers
  ENL= as.data.frame(matrix(rep(0,N), nrow=N, ncol=1)) #First exposed lab worker: create single column matrix where no one is in the exposed but not yet infectious lab worker category.
  IPL= as.data.frame(matrix(rep(0,N), nrow=N, ncol=1)) #First exposed lab worker: create single column matrix where no one is in the infectious, pre-symptomatic lab worker category.
  ISL= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) #First infectious lab worker: create single column matrix where no one is in the infectious, symptomatic lab worker category.
  SG= as.data.frame(matrix(c(rep(0,NL),rep(1,NG)),nrow=N, ncol = 1)) #First susceptible general population: create single column matrix where the last 400 people are susceptibles in the general population
  ENG= as.data.frame(matrix(rep(0,N), nrow=N, ncol=1)) #First exposed general population: create single column matrix where no one is in the exposed but not yet infectious general population category.
  IPG= as.data.frame(matrix(rep(0,N), nrow=N, ncol=1)) #First exposed general population: create single column matrix where no one is in the infectious, pre-symptomatic general population category.
  ISG= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) #First infectious general population: create single column matrix where no one is in the infectious, symptomatic general population category.
  ENISO= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) #First exposed, noninfectious isolated person:single column matrix where no one is in the exposed, non-infectious category.
  IPISO= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) #First infectious, presymptomatic isolated person:single column matrix where no one is in the infectious, pre-symptomatic isolated person
  ISISO= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) #First infectious, symptomatic isolated person:single column matrix where no one is in the infectious, symptomatic isolated person
  IAL= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) # First infectious, permanently asymptomatic lab worker: single column matrix where no lab worker is in the infectious, permanently asymptomatic category
  IAG= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) # First infectious, permanently asymptomatic general population person: single column matrix where no general population person is in the infectious, permanently asymptomatic category
  EVISO= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) # Single column matrix for people who have ever been isolated (starts at 0)
  LWEVINF= as.data.frame(matrix(rep(1,N), nrow=N, ncol = 1)) # Single column matrix for lab workers who have ever been infected (starts at 1)
  GPEVINF= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) # Single column matrix for general population who have ever been infected (starts at 0)
  Incidence= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1))  #Single column matrix for new cases at each time step. Start it at 0 for now.
  R= as.data.frame(matrix(rep(0,N), nrow=N, ncol = 1)) #First removed: create single column matrix where no people are recovered/dead
  ENL1= sample(1:NL, size=1) #Pick row number containing the first random exposed lab worker, name it ENL1
  ENL[ENL1,1]=1 #make ENL=1 in row ENL1 (this is assigning the first exposed lab worker)
  SL[ENL1,1]=0 #make S=0 in row ENL1 (this is accounting for the fact that the first infected person is no longer susceptible)
  Incidence[ENL1,1]=1 #make Incidence=1 in row ENL1 (this is the very first incident case)
  freq=thetaL
  t=1
  while(t<100){ #as long as t<100, and...
    print(t)
    #browser() #:only use when debugging
    #print(EVISO[1,])
    if(t==1 | sum(EVISO[,t-1])==0) { #as long as no one has been isolated, OR if you're at the initial time step (t==1)...
      tauL = tauL
      epsilonL = epsilonL 
      thetaL = thetaL 
      delta = delta
      gammaL = gammaL
      piL = piL 
      tauG = tauG
      thetaG= thetaG
      omegaG= omegaG
      chiG =chiG
      gammaG=gammaG
      piG= piG
      epsilonIS = epsilonIS
      gammaIS =gammaIS
      epsilonIAL= epsilonIAL
      gammaIAL= gammaIAL
      gammaIAG =gammaIAG
      iota= iota
      psi= psi 
      zeta= zeta
      sigma =sigma
      mu= mu
      nu=nu
      isolationdelay = isolationdelay
    }else if(sum(EVISO[,t-1])>0 & sum(EVISO[,t-1])<10){
      thetaL=1
    }else {
      chiG=1/7
      thetaL=1
    }
    #######################ACTION CODE#################################
    t=t+1 #move to the next time step

    #Factorize the matrix multiplication for efficiency
    Ninfneighbors<-(CM%*%(ISL[,t-1] + IPL[,t-1]*psi +ISG[,t-1] +IPG[,t-1]*psi +IAL[,t-1]*zeta +IAG[,t-1]*zeta +IPISO[,t-1]*psi*iota + ISISO[,t-1]*iota))
    pinfL<- 1-(1-tauL)^Ninfneighbors
    pinfG<- 1-(1-tauL)^Ninfneighbors
    
    
    ########################Define the Flows##############
    
    #Create the infection flows
    newENL= rbinom(N, SL[,t-1],pinfL) #random draw from previous time-step's susceptible lab workers to select the new exposed lab workers. Run SL[,t-1] trials with pinfL probability of infection for each observation in 1:N.
    newENG= rbinom(N, SG[,t-1],pinfG) #random draw from previous time-step's susceptible general population to select the new exposed general population. Run SG[,t-1] trials with pinfG probability of infection for each observation in 1:N.
    #Create empty vectors for new non-infection flows
    newENISOfromENL=rep(0,N)
    newIPL= rep(0,N)
    newIALfromENL = rep(0,N) #new
    newENISOfromENG= rep(0,N)
    newIPG= rep(0,N)
    newIAGfromENG = rep(0,N) #new
    newIPISOfromENISO= rep(0,N)
    newIPISOfromIPL= rep(0,N)
    newISL= rep(0,N)
    newIPISOfromIPG= rep(0,N)
    newISG=rep(0,N)
    newISISOfromIPISO= rep(0,N)
    newISISOfromISL= rep(0,N)
    newISISOfromIAL= rep(0,N) #new
    newRecoveredfromISL= rep(0,N)
    newISISOfromISG= rep(0,N)
    newRecoveredfromISG= rep(0,N)
    newRecoveredfromISISO=rep(0,N)
    newRecoveredfromIAL= rep(0,N) #new
    newRecoveredfromIAG= rep(0,N) #new
    
    #Transition rates out of ENL with delay
    ENLtoENISO<- 1/(1/(thetaL*sigma*mu)+isolationdelay)
    ENLtoIPL<- (1-ENLtoENISO)*epsilonL
    ENLtoIAL<- (1-ENLtoENISO)*(1-epsilonL)*epsilonIAL
    ENLstay<- (1-ENLtoENISO)*(1-epsilonL)*(1-epsilonIAL)
    
    #Code for flows out of ENL: 
    leavingENL<- sapply(1:N, function(x) rmultinom(1, ENL[x,t-1], prob = c(ENLstay,#Stay in ENL
                                                                           ENLtoIPL,#to IPL
                                                                           ENLtoIAL,#to IAL
                                                                           ENLtoENISO))) #to ENISO
    
    #extract each row of leavingENL, save as a vector
    #skip the first row of leavingENL, because those people are actually staying in ENL  
    newIPL<- leavingENL[2,]
    newIPL<- as.vector(newIPL)
    
    newIALfromENL<- leavingENL[3,]
    newIALfromENL<- as.vector(newIALfromENL)
    
    newENISOfromENL<- leavingENL[4,]
    newENISOfromENL<- as.vector(newENISOfromENL)
    
    #Transition rates out of ENG with delay
    ENGtoENISO<- 1/(1/(thetaG*sigma*mu)+isolationdelay)
    ENGtoIPG<- (1-ENGtoENISO)*epsilonL
    ENGtoIAG<- (1-ENGtoENISO)*(1-epsilonL)*epsilonIAL
    ENGstay<- (1-ENGtoENISO)*(1-epsilonL)*(1-epsilonIAL)
    
    #Code for flows out of ENG: 
    leavingENG<- sapply(1:N, function(x) rmultinom(1, ENG[x,t-1], prob = c(ENGstay,#Stay in ENG
                                                                           ENGtoIPG,#to IPG
                                                                           ENGtoIAG,#to IAG
                                                                           ENGtoENISO))) #to ENISO
    
    #extract each row of leavingENG, save as a vector
    #skip the first row of leavingENG, because those people are actually staying in ENG  
    newIPG<- leavingENG[2,]
    newIPG<- as.vector(newIPG)
    
    newIAGfromENG<- leavingENG[3,]
    newIAGfromENG<- as.vector(newIAGfromENG)
    
    newENISOfromENG<- leavingENG[4,]
    newENISOfromENG<- as.vector(newENISOfromENG)
    
    #Code for flow out of ENISO:
    newIPISOfromENISO= rbinom(N, ENISO[,t-1], epsilonIS) 
    
    #Code for flows out of IPL: 
    newIPISOfromIPL= rbinom(N, IPL[,t-1], 1/(1/(thetaL*sigma*nu)+isolationdelay))
    newISL= rbinom(N,IPL[,t-1]*(1-newIPISOfromIPL), delta)
    
    #Code for flows out of IPG: 
    newIPISOfromIPG= rbinom(N, IPG[,t-1], 1/(1/(omegaG*sigma*nu) +isolationdelay)) #select isolated people
    newISG= rbinom(N, IPG[,t-1]*(1-newIPISOfromIPG), delta) #from non-isolated people, select progressors
    
    #Code for flow out of IPISO:
    newISISOfromIPISO= rbinom(N, IPISO[,t-1], delta) #select those isolated who will progress to symptoms
    
    #Code for flows out of ISL: 
    newISISOfromISL= rbinom(N, ISL[,t-1], 1/(1/(thetaL*sigma*(1-piL) +piL)+isolationdelay)) #select the isolated
    newRecoveredfromISL= rbinom(N, ISL[,t-1]*(1-newISISOfromISL), gammaL) #select the recovered LWs who didn't get isolated
    
    #Code for flows out of ISG: 
    newISISOfromISG= rbinom(N, ISG[,t-1], 1/(1/(chiG*sigma*(1-piG)+piG) +isolationdelay)) #select the isolated
    newRecoveredfromISG= rbinom(N, ISG[,t-1]*(1-newISISOfromISG), gammaG) #select the recovered GPs who didn't get isolated    
    
    #Code for flows out of IAL:
    newISISOfromIAL= rbinom(N, IAL[,t-1], thetaL *sigma*0.8)
    newRecoveredfromIAL= rbinom(N, IAL[,t-1]*(1-newISISOfromIAL), gammaIAL)
    
    #Code for flows out of IAG:
    newRecoveredfromIAG= rbinom(N, IAG[,t-1], gammaIAG)
    
    #Code for flows out of ISISO:
    newRecoveredfromISISO= rbinom(N, ISISO[,t-1], gammaIS)
    
    #Code for counting up all the people who have been newly isolated in a given time-step
    newIsolatedAll= newENISOfromENL+ newENISOfromENG+ newIPISOfromIPL+ newIPISOfromIPG+ newISISOfromISL+ newISISOfromIAL+ newISISOfromISG
    
    #######Define the states at the next time-step#######
    
    nextSL= SL[,t-1]-newENL #next step's susceptible lab workers is the previous step's susceptible lab workers, minus the newly exposed but not infectious lab workers.
    nextSG= SG[,t-1]-newENG #next step's susceptible general population is the previous step's susceptible general population, minus the newly exposed but not infectious general population.
    nextENL= ENL[,t-1]+newENL -newIPL -newENISOfromENL- newIALfromENL #next step's exposed but not infectious lab workers is the previous step's exposed but not infectious lab workers, plus the newly exposed but not infectious lab workers, minus the newly exposed and infectious lab workers, minus the newly isolated exposed but non infectious lab workers, minus the new infected perm-asymptomatic LWs
    nextENG= ENG[,t-1]+newENG -newIPG -newENISOfromENG -newIAGfromENG #next step's exposed but not infectious general population is the previous step's exposed but not infectious general population, plus the newly exposed but not infectious general population, minus the newly exposed and infectious general population, minus the newly isolated exposed but non infectious general population, minus the new infected perm-asymptomatic GPs.
    
    nextENISO= ENISO[,t-1] + newENISOfromENL + newENISOfromENG  - newIPISOfromENISO #next step's exposed, noninfectious isolateds includes the previous step's exposed, noninfectious isolateds, plus the new exposed, noninfectious isolateds from ENL and ENG, minus the newly infectious, pre-symptomatic isolateds
    
    nextIPL= IPL[,t-1]+newIPL -newISL -newIPISOfromIPL #next step's exposed and infectious lab workers is the previous step's exposed and infectious lab workers, plus the newly exposed and infectious lab workers, minus the newly infectious and symptomatic lab workers, minus the newly isolated exposed, infectious lab workers.
    nextIPG= IPG[,t-1]+newIPG -newISG -newIPISOfromIPG #next step's exposed and infectious general population is the previous step's exposed and infectious general population, plus the newly exposed and infectious general population, minus the newly infectious and symptomatic general population, minus the newly isolated exposed, infectious general population
    nextIPISO= IPISO[,t-1]  + newIPISOfromIPL +newIPISOfromIPG +newIPISOfromENISO- newISISOfromIPISO
    #next step's infectious,pre-symptomatic isolateds (IPISO) is previous step's IPISO, plus newly isolated from infectious pre-symptomatic LWs, plus newly ioslated from infectious presymptomatic GPs, plus those who progressed from exposed, noninfectious isolated status, minus those who progress to infectious symptomatic isolated status
    nextISL= ISL[,t-1] +newISL -newRecoveredfromISL - newISISOfromISL # next step's infectious and symptomatic is previous step's infectious and symptomatic LWs, plus newly infectious and symptomatic LWs, minus the newly recovered LWs, minus the newly isolated LWs.
    nextISG= ISG[,t-1] +newISG -newRecoveredfromISG - newISISOfromISG # next step's infectious and symptomatic is previous step's infectious and symptomatic GPs, plus newly infectious and symptomatic GPs, minus the newly recovered GPs, minus the newly isolated GPs
    nextIAL= IAL[,t-1] + newIALfromENL -newRecoveredfromIAL -newISISOfromIAL # next step's infectious perm-asymptomatic LWs (IALs) is previous step's IALs, plus new IALs who were previously exposed non-infectious LWs, minues newly isolated IALs, minus newly recovered IALs
    nextIAG= IAG[,t-1] +newIAGfromENG -newRecoveredfromIAG  # next step's infectious perm-asymptomatic GPs (IAGs) is previous step's IAGs, plus new IAGs who were previously exposed non-infectious GPs, minus newly recovered IAGs
    nextISISO= ISISO[,t-1] + newISISOfromISL + newISISOfromISG +newISISOfromIAL+ newISISOfromIPISO -newRecoveredfromISISO
    # next step's infectious symptomatic isolateds (ISISO) is the previous step's ISISO, plus the infected symptomatic LWs who are isolated, plus the infected symptomatic GPs who are isolated, plus the permanently asymptomatic infected LWs who are isolated, plus the isolated who progress from the infectious presymptomatic stage, minus the newly recovered ISISO
    nextEVISO= EVISO[,t-1]+newIsolatedAll #Next step's total number of ever isolated is the previous step's number of ever isolated plus all of the newly isolated people
    nextLWEVINF= LWEVINF[,t-1] + newENL # Next steps's LWs who were ever infected is the previous step's LWs who were ever infected plus the newly infected LWs
    nextGPEVINF= GPEVINF[,t-1] + newENG # Next steps's GPs who were ever infected is the previous step's GPs who were ever infected plus the newly infected GPs
    nextstepcases= newENL + newENG #Next step's newly infected cases
    nextR= R[,t-1]+ newRecoveredfromISL + newRecoveredfromISG+ newRecoveredfromISISO+ newRecoveredfromIAL +newRecoveredfromIAG   #next step's removed people is previous step's removed, plus the newly recovered infectious symptomatic LWs, plus the newly recovered infectious symptomatic GPs, plus the newly recovered infectious symptomatic LWs, plus the newly recovered isolated infectious symptomatic GPs, plus the newly recovered infectious but asymptomatic LWs, plus the newly recovered infectious but asymptomatic GPs
    SL=cbind(SL, nextSL) #bind SL, nextSL matrices together
    SG=cbind(SG, nextSG) #bind SG, nextSG matrices together
    ENL=cbind(ENL, nextENL) #bind ENL, nextENL matrices together
    ENG=cbind(ENG, nextENG) #bind ENG, nextENG matrices together
    ENISO= cbind(ENISO, nextENISO) #bind ENISO, nextENISO matrices together
    IPL=cbind(IPL, nextIPL) #bind IPL, nextIPL matrices together
    IPG=cbind(IPG, nextIPG) #bind IPG, nextIPG matrices together
    IPISO= cbind(IPISO, nextIPISO) #bind IPISO, nextIPISO together
    IAL= cbind(IAL, nextIAL) #bind IAL, nextIAL together
    IAG= cbind(IAG, nextIAG) #bind IAG, nextIAG together
    ISL=cbind(ISL, nextISL) #bind ISL, nextISL matrices together
    ISG=cbind(ISG, nextISG) #bind ISG, nextISG matrices together
    ISISO= cbind(ISISO, nextISISO) #bind ISISO, nextISISO together
    EVISO= cbind(EVISO, nextEVISO) #bind EVISO, nextEVISO together
    LWEVINF= cbind(LWEVINF, nextLWEVINF) #bind LWEVINF, nextLWEVINF together
    GPEVINF= cbind(GPEVINF, nextGPEVINF) #bind GPEVINF, nextGPEVINF together
    Incidence= cbind(Incidence, nextstepcases)
    R= cbind(R, nextR) #bind R, nextR matrices together
  }
}
  freq= rep(freq, t)
  sigma= rep(sigma, t)
  isolationdelay= rep(isolationdelay, t)
  res=list(SL=SL, SG=SG, ENL=ENL, ENG=ENG, ENISO=ENISO, IPL=IPL, IPG=IPG, IPISO=IPISO, IAL=IAL, IAG=IAG, ISL=ISL, ISG=ISG, ISISO=ISISO, EVISO=EVISO, LWEVINF= LWEVINF,GPEVINF=GPEVINF,Incidence=Incidence, R=R,freq= freq, sigma=sigma, isolationdelay=isolationdelay) #save results in a list of matrices SL,SG,ENL,ENG,ENISO, IPL,IPG, IPISO, IAL, IAG, ISL,ISG, ISISO, R
  res=summary.TripleSEIIIRiso(res) #turn the complete simulation into a summary dataframe
  res= as.data.frame(res)
  #res$outbreak100<-ifelse(max(R)>=100,"yes","no")
  #class(res)="TripleSEIIIRiso" #set class of results as "TripleSEIIIRiso"
  return(res) #return results
}



##Sim Test 

sim7=TripleSEIIIRiso(cm6, tauL=0.5, epsilonL=0.5, thetaL=1, delta=0.5, gammaL=0.25, piL=0.1, tauG=0.05, thetaG=0,omegaG=0, chiG=0.01, gammaG=0.25,piG=0.1, epsilonIS=0.5, gammaIS=0.25, epsilonIAL = 0.5, gammaIAL=0.25, gammaIAG=0.25, iota=0.2, psi=0.8, zeta=0.9, sigma=0.95,mu=0.15, nu=0.5, isolationdelay = 0)
#sim7
#plot(apply(c(sim5$IL),2,sum),type="l", xlab = "Time", ylab= "Infected Lab Workers")
sim7summary<-summary.TripleSEIIIRiso(sim7) #save summary as a dataframe
sim7summary
#sim5summary<- sim5summary %>% mutate(time= row_number())
#extended.summary.TripleSEEIR(sim5)
plot.TripleSEIIIRiso(sim7)



##Simtestforminimalsavingfunc

start_time <- Sys.time()

sim8=TripleSEIIIRiso2(cm6, tauL=0.5, epsilonL=0.5, thetaL=1, delta=0.5, gammaL=0.25, piL=0.1, tauG=0.05, thetaG=0, omegaG=0, chiG=0.01, gammaG=0.25,piG=0.1, epsilonIS=0.5, gammaIS=0.25, epsilonIAL = 0.5, gammaIAL=0.25, gammaIAG=0.25, iota=0.2, psi=0.8, zeta=0.9, sigma=0.95,mu=0.15, nu=0.5, isolationdelay = 0)
sim8
end_time <- Sys.time()
end_time - start_time
# avg Time difference of 1.70 seconds when we didn't minimize matrix multiplication



##R0 Estimator

#Create Approximate R0 estimator based on network, tau, and recovery rate (Bjornstad 2018)
r0fun= function(CM, tau, recoveryrate){
  x= apply(CM,2,sum)
  (tau/(tau+recoveryrate))*(mean(x^2)-mean(x))/mean(x)
}
#Here, we estimate the recovery rate as the inverse of the sum of the inverses of epsilon, delta, and gamma... e.g. if epsilon, delta, and gamma are 0.19, 0.77, and 0.1, the sum of their inverses (which is the total recovery time) is 5.26+1.3+10= 16.56, and thus the recovery rate is 1/16.56
# 
#r0fun(cm6,0.01,1/16.56)

####Set Params####

#Use expand.grid to create par_tab, which is a df table with all of the possible  combinations of params of interest in separate rows (holding some params constant)
par_tab <- expand.grid(thetaL= 7/7, sigma= seq(0.5,1,0.1), 
                       isolationdelay= 0:3, tauL=0.0125, epsilonL=0.19, 
                       delta=0.77, gammaL=0.1, piL=0.02, tauG=0.0125, 
                       epsilonG=0.19,thetaG=0, omegaG=0, 
                       chiG=0.01, gammaG=0.1,piG=0, epsilonIS=0.19, 
                       gammaIS=0.1, epsilonIAL = 0.0396, 
                       gammaIAL=0.07, epsilonIAG=0.0396, gammaIAG=0.07, 
                       iota=0.2, psi=0.75, zeta=0.5, mu=0.15, nu=0.48)
n_sims = 1000


####Run Sims####

start_time <- Sys.time()
#Create a list to store the results
output <- list()

#Run loop where...
for(i in 1:nrow(par_tab)){ #for every row in parameter table representing a unique parameter combination [outer loop]
  output[[i]] <- list() #create a result list for that row
  for(j in 1:n_sims){ #for every simulation... [inner loop]
    set.seed(j) #set a unique simulation seed
    
    #browser()
    output[[i]][[j]] <- TripleSEIIIRiso2(CM= cm6, thetaL =par_tab$thetaL[i], sigma= par_tab$sigma[i],isolationdelay = par_tab$isolationdelay[i], par_tab$tauL[i], par_tab$epsilonL[i], par_tab$delta[i], par_tab$gammaL[i], par_tab$piL[i], par_tab$tauG[i], par_tab$thetaG[i], par_tab$omegaG[i], par_tab$chiG[i], par_tab$gammaG[i], par_tab$piG[i], par_tab$epsilonIS[i], par_tab$gammaIS[i], par_tab$epsilonIAL[i], par_tab$gammaIAL[i], par_tab$gammaIAG[i], par_tab$iota[i], par_tab$psi[i], par_tab$zeta[i], par_tab$mu[i], par_tab$nu[i])
    # Calculate unique simulation number
    unique_sim_number =7*n_sims*nrow(par_tab)+(i - 1) * n_sims + j
    # Add unique simulation number to the result
    output[[i]][[j]]$unique_sim_number <- unique_sim_number
    
    end_time <- Sys.time()
    timeelapsed<-end_time - start_time    
    
  }
}

####Unpack lists of outputs####
output<- flatten(output)

####Rename outputs####
#Rename simulation sets in output list according to the parameter values for each simulation set

#create a vector of unique names for each row in the param table (for each unique param combination)

outputname= rep(0, nrow(par_tab))
for (i in 1:nrow(par_tab)){ #for every row in param table
  outputname[i]= paste(paste("sigma=", par_tab[i,2]), paste("isolationdelay=", par_tab[i,3])) 
}


#Create new vector of output names where we repeat each unique name n_sims times
outputname<- rep(outputname, each=n_sims)
#Since we unpacked output in order, we can rename the simulations in the same order as outputname
names(output)<- outputname
#save(output, file= "output_flat.RData")


####Create Analytic DF####
analyticdf<- output[[1]]
for (i in 2:length(output)){
  analyticdf<- rbind(analyticdf, output[[i]])
  print(i)
}

####Create Mini Analytic DF####
minianalyticdf<- analyticdf[analyticdf$Time==100,]
minianalyticdf$Outbreak100<- ifelse(minianalyticdf$Outbreak100=="yes",1,0)
minianalyticdf$Outbreak50<- ifelse(minianalyticdf$Outbreak50=="yes",1,0)
minianalyticdf$Outbreak10<- ifelse(minianalyticdf$Outbreak10=="yes",1,0)
minianalyticdf<- minianalyticdf %>% mutate(Outbreak20= ifelse(R>=20,1,0))

####Rename vars####
minianalyticdf$frequency<-minianalyticdf$freq
minianalyticdf$sensitivity<-minianalyticdf$sigma
minianalyticdf$sensitivity<- minianalyticdf$sensitivity*100 #convert sensitivity to percentage point scale
####Save minianalyticdf####
save(minianalyticdf, file= "Manalyticdf_paper1_freq7.RData")
####Save analyticdf####
save(analyticdf, file= "Analyticdf_paper1_freq7.RData")

