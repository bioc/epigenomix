.bayesMixModel <- function(z, normNull=c(), expNeg=c(), expPos=c(), gamNeg=c(), gamPos=c(),
                           sdNormNullInit=c(), rateExpNegInit=c(), rateExpPosInit=c(),
                           shapeGamNegInit=c(), scaleGamNegInit=c(), shapeGamPosInit=c(), scaleGamPosInit=c(),
                           piInit, classificationsInit, dirichletParInit=1, shapeDir=1, scaleDir=1, weightsPrior="FDD", sdAlpha,
                           shapeNorm0=c(), scaleNorm0=c(), shapeExpNeg0=c(), scaleExpNeg0=c(), shapeExpPos0=c(), scaleExpPos0=c(),
                           shapeGamNegAlpha0=c(), shapeGamNegBeta0=c(), scaleGamNegAlpha0=c(), scaleGamNegBeta0=c(),
                           shapeGamPosAlpha0=c(), shapeGamPosBeta0=c(), scaleGamPosAlpha0=c(), scaleGamPosBeta0=c(),
                           itb, nmc, thin, average="mean",sdShape) {

  
k <- length(c(normNull,expNeg,expPos,gamNeg,gamPos))#length(c(normNull,normNeg,normPos,expNeg,expPos))
nnorm <- length(normNull)#length(c(normNull,normNeg,normPos))
nexp <- length(c(expNeg,expPos))
ngamma <- length(c(gamNeg,gamPos))
n <- length(z)

if (missing(piInit)) {
  piInit <- rep(1/k, k)
}

if (missing(classificationsInit)) {
  classificationsInit <- rep(floor(k/2), length(z))
}

# Umordnung der Indizes
indizes <- c(normNull,expNeg,expPos,gamNeg,gamPos)

normNull <- which(indizes %in% normNull)
expNeg <- which(indizes %in% expNeg)
expPos <- which(indizes %in% expPos)
gamNeg <- which(indizes %in% gamNeg)
gamPos <- which(indizes %in% gamPos)

if(length(shapeGamNegInit) < length(gamNeg)){stop("shapeGamNegInit is insufficiently specified.")}
if(length(scaleGamNegInit) < length(gamNeg)){stop("scaleGamNegInit is insufficiently specified.")}
if(length(shapeGamPosInit) < length(gamPos)){stop("shapeGamPosInit is insufficiently specified.")}
if(length(scaleGamPosInit) < length(gamPos)){stop("scaleGamPosInit is insufficiently specified.")}
if(length(rateExpNegInit) < length(expNeg)){stop("rateExpNegInit is insufficiently specified.")}
if(length(rateExpPosInit) < length(expPos)){stop("rateExpPosInit is insufficiently specified.")}
if(length(sdNormNullInit) != nnorm){stop("Length of vector sdNormNullInit must correspond to number of normal components.")}
if(length(piInit) != k){stop("Length of vector piInit must correspond to total number of components.")}
if(length(classificationsInit) != n){stop("Length of vector classificationsInit must correspond to length of data vector.")}

if((length(dirichletParInit) != 1) | (dirichletParInit <= 0)){stop("dirichletParInit must be a single positive real number.")}
if(length(shapeNorm0) != nnorm){stop("Length of shapeNorm0 must correspond to number of normal components.")}
if(length(scaleNorm0) != nnorm){stop("Length of scaleNorm0 must correspond to number of normal components.")}
if(length(shapeExpNeg0) < length(expNeg)){stop("shapeExpNeg0 is insufficiently specified.")}
if(length(scaleExpNeg0) < length(expNeg)){stop("scaleExpNeg0 is insufficiently specified.")}
if(length(shapeExpPos0) < length(expPos)){stop("shapeExpPos0 is insufficiently specified.")}
if(length(scaleExpPos0) < length(expPos)){stop("scaleExpPos0 is insufficiently specified.")}
if(length(shapeGamNegAlpha0) < length(gamNeg)){stop("shapeGamNegAlpha0 is insufficiently specified.")}
if(length(shapeGamNegBeta0) < length(gamNeg)){stop("shapeGamNegBeta0 is insufficiently specified.")}
if(length(scaleGamNegAlpha0) < length(gamNeg)){stop("shapeGamNegAlpha0 is insufficiently specified.")}
if(length(scaleGamNegBeta0) < length(gamNeg)){stop("shapeGamNegBeta0 is insufficiently specified.")}
if(length(shapeGamPosAlpha0) < length(gamPos)){stop("shapeGamPosAlpha0 is insufficiently specified.")}
if(length(shapeGamPosBeta0) < length(gamPos)){stop("shapeGamPosBeta0 is insufficiently specified.")}
if(length(scaleGamPosAlpha0) < length(gamPos)){stop("shapeGamPosAlpha0 is insufficiently specified.")}
if(length(scaleGamPosBeta0) < length(gamPos)){stop("shapeGamPosBeta0 is insufficiently specified.")}

##### Initialisierung der Komponenten #####
mu <- c(rep(NA,nnorm),1/rateExpNegInit,1/rateExpPosInit)
Sigma <- sdNormNullInit
shapeGam <- c(shapeGamNegInit,shapeGamPosInit)
scaleGam <- c(scaleGamNegInit,scaleGamPosInit)
Pi <- piInit
Pi <- Pi/sum(Pi); zdp <- classificationsInit

##### Initialisierung der Ketten #####
S.pi <- matrix(rep(NA,ceiling(nmc/thin)*k),ncol=k); S.mu <- matrix(rep(NA,ceiling(nmc/thin)*(nnorm+nexp)),ncol=nnorm+nexp); S.sigma <- matrix(rep(NA,ceiling(nmc/thin)*nnorm),ncol=nnorm); S.zdp <- matrix(numeric(n*k),ncol=k)
S.shapeNorm <- matrix(rep(NA,ceiling(nmc/thin)*nnorm),ncol=nnorm); S.scaleNorm <- matrix(rep(NA,ceiling(nmc/thin)*nnorm),ncol=nnorm);
S.alphadir <- rep(NA,ceiling(nmc/thin)); S.alphadirAccept <- rep(NA,ceiling(nmc/thin)); S.allocs <- matrix(rep(NA,ceiling(nmc/thin)*k),ncol=k);
if(nexp > 0){
     S.shapeExp <- matrix(rep(NA,ceiling(nmc/thin)*nexp),ncol=nexp); S.scaleExp <- matrix(rep(NA,ceiling(nmc/thin)*nexp),ncol=nexp);
}
if(ngamma > 0){
    S.shapeGam <- matrix(rep(NA,ceiling(nmc/thin)*ngamma),ncol=ngamma); S.scaleGam <- matrix(rep(NA,ceiling(nmc/thin)*ngamma),ncol=ngamma); S.acceptanceProb <- matrix(rep(NA,ceiling(nmc/thin)*ngamma),ncol=ngamma)
    S.scaleGamAlpha <- matrix(rep(NA,ceiling(nmc/thin)*ngamma),ncol=ngamma); S.scaleGamBeta <- matrix(rep(NA,ceiling(nmc/thin)*ngamma),ncol=ngamma);
}
#Initial sampling fuer DP-Mischung
alphadir <- dirichletParInit
if(weightsPrior == "FDD"){
    PiV <- .sampleMixture(k,zdp,alphadir)
    Pi <- PiV$Pi
}
if(weightsPrior == "TDP"){
       PiV <- .sampleMixtureTDP(k,zdp,alphadir)
    Pi <- matrix(PiV$Pi, nrow=1)
       V <- PiV$V
}
ind <- 1
############################## Start Gibbs sampling ##############################
    for(it in 1:(itb+nmc)){
        if(it %% 50 == 0){print(paste("Iteration:", it))}
        ################ update Theta ################
        #Allokationsvariable zdp sampeln
        zdp <- .sampleAllocations(z,n,k,mu,Sigma,shapeGam,scaleGam,Pi,normNull,expNeg,expPos,gamNeg,gamPos)
        allocs <- numeric(k)
        allocs_sparse <- table(zdp)
        dabei <- as.numeric(dimnames(allocs_sparse)$zdp)
        for(l in 1:length(dabei)){
           allocs[dabei[l]] <- allocs_sparse[l]    
        }
        #v <- unique(zdp); c <- hist(zdp,v);
    
        #Mixture sampeln;
        if(weightsPrior == "FDD"){ 
            PiV <- .sampleMixture(k,zdp,alphadir)
        Pi <- PiV$Pi
        }
        if(weightsPrior == "TDP"){  
            PiV <- .sampleMixtureTDP(k,zdp,alphadir)
        Pi <- matrix(PiV$Pi, nrow=1)
              V <- PiV$V
        }
      
        #Parameter der Mischungskomponenten sampeln
        muSigma <- .sampleComponentParameters(z,k,nnorm,nexp,ngamma,zdp,normNull,expNeg,expPos,gamNeg,gamPos,shapeNorm0,scaleNorm0,shapeExpNeg0,scaleExpNeg0,shapeExpPos0,scaleExpPos0,shapeGamNegAlpha0,shapeGamNegBeta0,scaleGamNegAlpha0,scaleGamNegBeta0,shapeGamPosAlpha0,shapeGamPosBeta0,scaleGamPosAlpha0,scaleGamPosBeta0,shapeGam,scaleGam,sdShape)
        mu <- muSigma$mu
        Sigma <- muSigma$Sigma
        shapeGam <- muSigma$shapeGam
        scaleGam <- muSigma$scaleGam
        acceptanceProb <- muSigma$acceptanceProb
        shapeNorm <- muSigma$shapeNorm
        scaleNorm <- muSigma$scaleNorm
        shapeExpNeg <- muSigma$shapeExpNeg
        scaleExpNeg <- muSigma$scaleExpNeg
        shapeExpPos <- muSigma$shapeExpPos
        scaleExpPos <- muSigma$scaleExpPos
        scaleGamNegAlpha <- muSigma$scaleGamNegAlpha
        scaleGamNegBeta <- muSigma$scaleGamNegBeta
        scaleGamPosAlpha <- muSigma$scaleGamPosAlpha
        scaleGamPosBeta <- muSigma$scaleGamPosBeta
               
        #Parameter des Dirichlezprozesses sampeln
        if(weightsPrior == "FDD"){       
         saAlpha <- .sampleAlpha(alphadir,a=shapeDir,b=scaleDir,werte=Pi,sdAlpha=sdAlpha)
            alphadir <- saAlpha$alpha_new
            alphadirAccept <- saAlpha$alpha_accept 
        }
        if(weightsPrior == "TDP"){
            alphadir <- rgamma(n=1,shape=shapeDir+k-1,scale=1/(scaleDir-sum(log(1-V[1:(k-1)]))))
        }

        #Ergebnisse nach Burnin speichern
        if (it>itb){
            if( length(intersect(it,seq(itb+1,itb+nmc,thin))) == 1){
                S.pi[ind,] <- Pi
                S.mu[ind,] <- mu
                S.sigma[ind,] <- Sigma
                S.shapeNorm[ind,] <- shapeNorm
                S.scaleNorm[ind,] <- scaleNorm
                S.allocs[ind,] <- allocs
                if(nexp > 0){
                    S.shapeExp[ind,] <- c(shapeExpNeg, shapeExpPos)
                    S.scaleExp[ind,] <- c(scaleExpNeg, scaleExpPos)
                }
                if(ngamma > 0){
                    S.shapeGam[ind,] <- shapeGam
                    S.scaleGam[ind,] <- scaleGam
                    S.acceptanceProb[ind,] <- acceptanceProb
                    S.scaleGamAlpha[ind,] <- c(scaleGamNegAlpha,scaleGamPosAlpha)
                    S.scaleGamBeta[ind,] <- c(scaleGamNegBeta,scaleGamPosBeta)
                }
                S.alphadir[ind] <- alphadir
                if(weightsPrior == "FDD"){  
                S.alphadirAccept[ind] <- alphadirAccept
             }
                for(j in 1:k){
                    S.zdp[zdp==j,j] <- S.zdp[zdp==j,j]+1
                }
                ind <- ind + 1
            }
        }    
    }
#Umordnung der Komponenten fuer Medianberechnung
reihenfolge_alles <- c(gamNeg,expNeg,normNull,expPos,gamPos)
reihenfolge_mu <- c(expNeg,normNull,expPos)
S.pi <- S.pi[,reihenfolge_alles]
S.zdp <- S.zdp[,reihenfolge_alles]
S.mu <- S.mu[,reihenfolge_mu]

gamNeg <- which(reihenfolge_alles %in% gamNeg)
expNeg <- which(reihenfolge_alles %in% expNeg)
normNull <- which(reihenfolge_alles %in% normNull)
expPos <- which(reihenfolge_alles %in% expPos)
gamPos <- which(reihenfolge_alles %in% gamPos)

if(average=="mean"){
mu_mean <- apply(S.mu,2,function(x){mean(x[!is.na(x)])})
sigma_mean <- apply(S.sigma,2,function(x){mean(x[!is.na(x)])})
if(ngamma > 0){
    shapeGam_mean <- apply(S.shapeGam,2,function(x){mean(x[!is.na(x)])})
    scaleGam_mean <- apply(S.scaleGam,2,function(x){mean(x[!is.na(x)])})
    acceptanceProb_mean <- apply(S.acceptanceProb,2,function(x){mean(x[!is.na(x)])})
}
pi_mean <- apply(S.pi,2,function(x){mean(x[!is.na(x)])})
alphadir_mean <- mean(S.alphadir[!is.na(S.alphadir)])
}
if(average=="median"){
mu_mean <- apply(S.mu,2,function(x){median(x[!is.na(x)])})
sigma_mean <- apply(S.sigma,2,function(x){median(x[!is.na(x)])})
if(ngamma > 0){
    shapeGam_mean <- apply(S.shapeGam,2,function(x){median(x[!is.na(x)])})
    scaleGam_mean <- apply(S.scaleGam,2,function(x){median(x[!is.na(x)])})
    acceptanceProb_mean <- apply(S.acceptanceProb,2,function(x){median(x[!is.na(x)])})
}
pi_mean <- apply(S.pi,2,function(x){median(x[!is.na(x)])})
alphadir_mean <- median(S.alphadir[!is.na(S.alphadir)])
}

calculateMode <- function(x){
  tabelle <- table(x)
  if(length(tabelle) > 0){
    wo <- which(tabelle == max(tabelle))
  } else {
    wo <- integer(0)
  }
  if(length(wo) == 1){erg <- as.numeric(labels(tabelle)$x[wo])}
  if(length(wo) == 0){erg <- floor(median(as.numeric(labels(tabelle)$x)))}
  if(length(wo) > 1){
    erg <- as.numeric(labels(tabelle)$x[wo])
    erg <- floor(median(erg))
  }
  erg
}

# Allokationsvariable auswerten
allocations_median <- allocations_mode <- numeric(n)
for(m in 1:n){
    alles = c()
    for(l in 1:k){
        alles <- c(alles, rep(l,S.zdp[m,l]))
    } 
    allocations_median[m] <- floor(median(alles))
    allocations_mode[m] <- calculateMode(alles)
}

values <- matrix(numeric(length(z) * k), nrow = k)
    if (nnorm > 0) {
        for (i in normNull) {
            values[i, ] <- pi_mean[i] * dnorm(x = z, mean = 0,
                sd = sigma_mean[reihenfolge_alles[i]])
        }
    }
    if (length(expNeg) > 0) {
        i <- expNeg
        values[i, ] <- pi_mean[i] * dexp(x = -z, rate = 1/mu_mean[reihenfolge_mu == reihenfolge_alles[i]])
    }
    if (length(expPos) > 0) {
        i <- expPos
        values[i, ] <- pi_mean[i] * dexp(x = z, rate = 1/mu_mean[reihenfolge_mu == reihenfolge_alles[i]])
    }
    if (length(gamNeg) > 0) {
        i <- gamNeg
        values[i, ] <- pi_mean[i] * dgamma(x = -z, shape = shapeGam_mean[i -
            max(c(0, expPos, expNeg, normNull))], scale = scaleGam_mean[i -
            max(c(0, expPos, expNeg, normNull))])
    }
    if (length(gamPos) > 0) {
        i <- gamPos
        values[i, ] <- pi_mean[i] * dgamma(x = z, shape = shapeGam_mean[i -
            max(c(0, expPos, expNeg, normNull, gamNeg))], scale = scaleGam_mean[i -
            max(c(0, expPos, expNeg, normNull, gamNeg))])
    }
    allocations_max_dens <- unlist(apply(values, 2, function(x) {
        which(x == max(x))
    }))

# construct lists with components parameters
compsInit = list()
for (ind in normNull) {
  compsInit[[ind]] = new("MixtureComponent",
             name="NormNull",
             parameters=list(mean=0, sd=sdNormNullInit[ind==normNull]),
             pdf=dnorm,
             color="blue")
}
for (ind in expNeg) { 
  compsInit[[ind]] = new("MixtureComponent",
             name="ExpNeg",
             parameters=list(rate=rateExpNegInit[ind==expNeg]),
             pdf=function(x, rate) { dexp(-x, rate) },
             color="red")
}
for (ind in expPos) {
  compsInit[[ind]] = new("MixtureComponent",
             name="ExpPos",
             parameters=list(rate=rateExpPosInit[ind==expPos]),
             pdf=dexp,
             color="green")
}
for (ind in gamNeg) {
  compsInit[[ind]] = new("MixtureComponent",
             name="GamNeg",
             parameters=list(shape=shapeGamNegInit[ind==gamNeg], scale=scaleGamNegInit[ind==gamNeg]),
             pdf=function(x, shape, scale) { dgamma(-x, shape=shape, scale=scale) },
             color="purple")
}
for (ind in gamPos) {
    compsInit[[ind]] = new("MixtureComponent",
             name="GamPos",
             parameters=list(shape=shapeGamPosInit[ind==gamPos], scale=scaleGamPosInit[ind==gamPos]),
             pdf=dgamma,
             color="cyan")
}

components_priors = list()
for (ind in normNull) {
  components_priors[[ind]] =  list(name="NormNull", shape=shapeNorm0[ind==normNull], scale=scaleNorm0[ind==normNull])
}
for (ind in expNeg) {
    components_priors[[ind]] = list(name="ExpNeg", shape=shapeExpNeg0[ind==expNeg], scale=shapeExpNeg0[ind==expNeg])
}
for (ind in expPos) {
    components_priors[[ind]] = list(name="ExpPos", shape=shapeExpPos0[ind==expPos], scale=shapeExpPos0[ind==expPos])
}
for (ind in gamNeg) {
    components_priors[[ind]] = list(name="GamNeg", shapeAlpha=shapeGamNegAlpha0[ind==gamNeg], shapeBeta=shapeGamNegBeta0[ind==gamNeg], scaleAlpha=scaleGamNegAlpha0[ind==gamNeg], scaleBeta=scaleGamNegBeta0[ind==gamNeg])
}
for (ind in gamPos) {
    components_priors[[ind]] = list(name="GamPos", shapeAlpha=shapeGamPosAlpha0[ind==gamPos], shapeBeta=shapeGamPosBeta0[ind==gamPos], scaleAlpha=scaleGamPosAlpha0[ind==gamPos], scaleBeta=scaleGamPosBeta0[ind==gamPos])
}

compsResults = list()
for (ind in normNull) {
  compsResults[[ind]] = new("MixtureComponent",
                name="NormNull",
                parameters=list(mean=0, sd=sigma_mean[ind==normNull]),
                pdf=dnorm,
                color="blue")
}
for (ind in expNeg) {
  if(average=="mean"){
  compsResults[[ind]] = new("MixtureComponent",
             name="ExpNeg",
             parameters=list(rate=mean(1/S.mu[,reihenfolge_mu == reihenfolge_alles[ind]])),
             pdf=function(x, rate) { dexp(-x, rate) },
             color="red")
  }
  if(average=="median"){
  compsResults[[ind]] = new("MixtureComponent",
             name="ExpNeg",
             parameters=list(rate=median(1/S.mu[,reihenfolge_mu == reihenfolge_alles[ind]])),
             pdf=function(x, rate) { dexp(-x, rate) },
             color="red")
  }
}
for (ind in expPos) {
  if(average=="mean"){
  compsResults[[ind]] = new("MixtureComponent",
             name="ExpPos",
             parameters=list(rate=mean(1/S.mu[,reihenfolge_mu == reihenfolge_alles[ind]])),
             pdf=dexp,
             color="green")
  }
  if(average=="median"){
  compsResults[[ind]] = new("MixtureComponent",
             name="ExpPos",
             parameters=list(rate=median(1/S.mu[,reihenfolge_mu == reihenfolge_alles[ind]])),
             pdf=dexp,
             color="green")
  }
}
for (ind in gamNeg) {
  compsResults[[ind]] = new("MixtureComponent",
                name="GamNeg",
                parameters=list(shape=shapeGam_mean[ind==c(gamNeg, gamPos)], scale=scaleGam_mean[ind==c(gamNeg, gamPos)]), #shapeAcceptanceProbability=acceptanceProb_mean[ind==c(gamNeg, gamPos)]
                pdf=function(x, shape, scale) { dgamma(-x, shape=shape, scale=scale) },
                color="purple")
}
for (ind in gamPos) {
  compsResults[[ind]] = new("MixtureComponent",
                name="GamPos",
                parameters=list(shape=shapeGam_mean[ind==c(gamNeg, gamPos)], scale=scaleGam_mean[ind==c(gamNeg, gamPos)]), #shapeAcceptanceProbability=acceptanceProb_mean[ind==c(gamNeg, gamPos)]
                pdf=dgamma,
                color="orange")
}

compsChains = list()
for (ind in normNull) {
    compsChains[[ind]] = list(name="NormNull", sd=S.sigma[,ind==normNull], precision.shape=S.shapeNorm[,ind==normNull], precision.scale=S.shapeNorm[,ind==normNull])
}
for (ind in expNeg) {
    compsChains[[ind]] = list(name="ExpNeg", rate=1/S.mu[,reihenfolge_mu == reihenfolge_alles[ind]], rate.shape=S.shapeExp[,ind==c(expNeg, expPos)], rate.scale=S.scaleExp[,ind==c(expNeg, expPos)])
}
for (ind in expPos) {
    compsChains[[ind]] = list(name="ExpPos", rate=1/S.mu[,reihenfolge_mu == reihenfolge_alles[ind]], rate.shape=S.shapeExp[,ind==c(expNeg, expPos)], rate.scale=S.scaleExp[,ind==c(expNeg, expPos)])
}
for (ind in gamNeg) {
    compsChains[[ind]] = list(name="GamNeg", shape=S.shapeGam[,ind==c(gamNeg, gamPos)], scale=S.scaleGam[,ind==c(gamNeg, gamPos)], scale.shape=S.scaleGamAlpha[,ind==c(gamNeg, gamPos)], scale.scale=S.scaleGamBeta[,ind==c(gamNeg, gamPos)], shape.acceptance=S.acceptanceProb[,ind==c(gamNeg, gamPos)])
}
for (ind in gamPos) {
    compsChains[[ind]] = list(name="GamPos", shape=S.shapeGam[,ind==c(gamNeg, gamPos)], scale=S.scaleGam[,ind==c(gamNeg, gamPos)], scale.shape=S.scaleGamAlpha[,ind==c(gamNeg, gamPos)], scale.scale=S.scaleGamBeta[,ind==c(gamNeg, gamPos)], shape.acceptance=S.acceptanceProb[,ind==c(gamNeg, gamPos)])
}
#configuration
configuration_inits <- list(components=compsInit,pi=piInit, classification=classificationsInit, dirichletParameter=dirichletParInit)
configuration_priors <- list(components=components_priors, dirichlet=list(shapealphadir=shapeDir, scalealphadir=scaleDir, weightsPrior=weightsPrior))
if(weightsPrior == "FDD"){
    if(ngamma > 0){
            configuration_chains <- list(itb=itb, thin=thin, nmc=nmc, sdAlpha=sdAlpha, sdShape=sdShape)
    } else {
            configuration_chains <- list(itb=itb, thin=thin, nmc=nmc, sdAlpha=sdAlpha)
    }
}
if(weightsPrior == "TDP"){
    if(ngamma > 0){
            configuration_chains <- list(itb=itb, thin=thin, nmc=nmc, sdShape=sdShape)
    } else {
            configuration_chains <- list(itb=itb, thin=thin, nmc=nmc)
    }
}

#results
results_classification <- list(mode=allocations_mode, median=allocations_median, maxDens=allocations_max_dens)

if(weightsPrior == "FDD"){
    mm <- new("MixModelBayes",
          mmData=z,
          configuration=list(inits=configuration_inits, priors=configuration_priors, chains=configuration_chains),
          results=list(components=compsResults, pi=pi_mean, classification=results_classification),
          chains=list(components=compsChains,pi=S.pi,dirichletParameter=S.alphadir,dirichletParameterAcceptance=S.alphadirAccept,classification=S.zdp,allocations=S.allocs))
}
if(weightsPrior == "TDP"){
    mm <- new("MixModelBayes",
          mmData=z,
          configuration=list(inits=configuration_inits, priors=configuration_priors, chains=configuration_chains),
          results=list(components=compsResults, pi=pi_mean, classification=results_classification),
          chains=list(components=compsChains,pi=S.pi,dirichletParameter=S.alphadir,classification=S.zdp,allocations=S.allocs))
}

return(mm)
}


#Berechnung der Gesamtdichte fuer alle n Punkte
.sampleAllocations <- function(x,n,k,mu,Sigma,shapeGam,scaleGam,Pi,normNull,expNeg,expPos,gamNeg, gamPos){

  mu[normNull] <- 0
  P1 <- Pi; P2 <- matrix(numeric(k*n), ncol=n)
  for(j in 1:max(normNull)){
    P2[j,] <- dnorm(x,mu[j],Sigma[j])
  }
  if(length(expNeg) > 0){
    for(j in expNeg){
    P2[j,] <- dexp(-x,1/mu[j])
  }
  }
  if(length(expPos) > 0){
    for(j in expPos){
      P2[j,] <- dexp(x,1/mu[j])
    }
  }  
  if(length(gamNeg) > 0){
    for(j in gamNeg){
      P2[j,] <- dgamma(-x,shape=shapeGam[j==c(gamNeg, gamPos)], scale=scaleGam[j==c(gamNeg, gamPos)])
    }
  }
  if(length(gamPos) > 0){
    for(j in gamPos){
      P2[j,] <- dgamma(x,shape=shapeGam[j==c(gamNeg, gamPos)], scale=scaleGam[j==c(gamNeg, gamPos)])
    }
  }    
  P <-  matrix(rep(P1,n),ncol=n);P <- P*P2
  #Normierung
  P <- P/t(matrix(rep(apply(P,2,sum),k),ncol=k))
  
  #Resampling Index z for i=1,...,n
  Prop_n <- apply(P,2,cumsum)
  rn <- runif(n);  z <- numeric(n)
  ir <- which(rn <= Prop_n[1,])
  if(is.null(ir)==FALSE){
    z[ir] <- 1
  }
  for(j in 2:k){
    ir <- which(rn>=Prop_n[j-1,] & rn < Prop_n[j,])
    if(is.null(ir)==FALSE){
      z[ir] <- j
    }
  }
  
  return(z)
}


.sampleComponentParameters <- function(x,k,nnorm,nexp,ngamma,zdp,normNull,expNeg,expPos,gamNeg,gamPos,shapeNorm0,scaleNorm0,shapeExpNeg0,scaleExpNeg0,shapeExpPos0,scaleExpPos0,shapeGamNegAlpha0,shapeGamNegBeta0,scaleGamNegAlpha0,scaleGamNegBeta0,shapeGamPosAlpha0,shapeGamPosBeta0,scaleGamPosAlpha0,scaleGamPosBeta0,shapeGam,scaleGam,sdShape){

  # Speicherobjekte initialisieren
  mu    <- rep(NA,nnorm+nexp)
  Sigma <- rep(NA,nnorm)
  shapeNorm <- rep(NA,nnorm)
  scaleNorm <- rep(NA,nnorm)

  shapeExpNeg <- rep(NA,length(expNeg))
  scaleExpNeg <- rep(NA,length(expNeg))
  shapeExpPos <- rep(NA,length(expPos))
  scaleExpPos <- rep(NA,length(expPos))

  acceptanceProb <- rep(NA,ngamma)
  shapeGamNegAlpha <- rep(NA,length(gamNeg))
  shapeGamNegBeta  <- rep(NA,length(gamNeg))
  scaleGamNegAlpha <- rep(NA,length(gamNeg))
  scaleGamNegBeta  <- rep(NA,length(gamNeg))
  shapeGamPosAlpha <- rep(NA,length(gamPos))
  shapeGamPosBeta  <- rep(NA,length(gamPos))
  scaleGamPosAlpha <- rep(NA,length(gamPos))
  scaleGamPosBeta  <- rep(NA,length(gamPos))

  for(j in 1:k){
    i <- which(zdp==j); nj <- length(i)
    # wenn die Komponente nicht "leer" ist 
    if(nj>0){
      # Update der Verteilungen
      if(length(intersect(j, normNull)) == 1){
       # Normalverteilungen mit festem Erwartungswert 0 
        werte <- x[i];
        xquer <- mean(werte); SS <- sum((werte-xquer)^2)
        # Updates
        shapeNorm[j] <- shapeNorm0[j] + nj/2;
        scaleNorm[j] <- 1/(1/scaleNorm0[j] + sum(werte^2)/2 )
      } else if(length(intersect(j, expNeg)) == 1){
      # Exponentialverteilungen fuer negative Werte
        # Einschraenkung auf negative Werte
        werte <- x[i]
        unter_null_ind <- which(werte < 0)
        xsum <- sum(abs(werte[unter_null_ind]))
        # Updates fuer Exponentialverteilung fuer negative Werte
        shapeExpNeg[j-max(normNull)] <- shapeExpNeg0[j-max(normNull)] + nj
        scaleExpNeg[j-max(normNull)] <- scaleExpNeg0[j-max(normNull)]/(1 + scaleExpNeg0[j-max(normNull)]*xsum)
      } else if(length(intersect(j, expPos)) == 1){
        # Einschraenkung auf positive Werte
        werte <- x[i]
        ueber_null_ind <- which(werte > 0)
        xsum <- sum(werte[ueber_null_ind])
        # Updates
        shapeExpPos[j-max(expNeg)] <- shapeExpPos0[j-max(expNeg)] + nj
        scaleExpPos[j-max(expNeg)] <- scaleExpPos0[j-max(expNeg)]/(1 + scaleExpPos0[j-max(expNeg)]*xsum)
      } else if(length(intersect(j, gamNeg)) == 1){
        # Einschraenkung auf negative Werte
        werte <- x[i]
        unter_null_ind <- which(werte < 0)
        S_neg <- round(sum(abs(werte[unter_null_ind])),digits=22)
        P_neg <- round(prod(abs(werte[unter_null_ind])),digits=22)
        # Updates
        scaleGamNegAlpha[j-max(c(0,expPos,expNeg,normNull))] <- shapeGam[j-max(c(0,expPos,expNeg,normNull))]*nj + scaleGamNegAlpha0[j-max(c(0,expPos,expNeg,normNull))]
        scaleGamNegBeta[j-max(c(0,expPos,expNeg,normNull))] <- scaleGamNegBeta0[j-max(c(0,expPos,expNeg,normNull))]/(1 + scaleGamNegBeta0[j-max(c(0,expPos,expNeg,normNull))]*S_neg)
      } else if(length(intersect(j, gamPos)) == 1){
        # Einschraenkung auf positive Werte
        werte <- x[i]
        ueber_null_ind <- which(werte > 0)
        S_pos <- round(sum(werte[ueber_null_ind]),digits=22)
        P_pos <- round(prod(werte[ueber_null_ind]),digits=22)
        # Updates
        scaleGamPosAlpha[j-max(c(0,expPos,expNeg,normNull,gamNeg))] <- shapeGam[j-max(c(0,expPos,expNeg,normNull))]*nj + scaleGamPosAlpha0[j-max(c(0,expPos,expNeg,normNull,gamNeg))]
        scaleGamPosBeta[j-max(c(0,expPos,expNeg,normNull,gamNeg))] <- scaleGamPosBeta0[j-max(c(0,expPos,expNeg,normNull,gamNeg))]/(1 + scaleGamPosBeta0[j-max(c(0,expPos,expNeg,normNull,gamNeg))]*S_pos)
      }
        
      # anschliessend Ziehen aus den upgedateten Verteilungen
      if(length(intersect(j, normNull)) == 1){
        # Ziehen von Werten aus Normalverteilung mit festem Erwartungswert
        mu[j] <- 0
        Sigma[j] <- sqrt(1/rgamma(n=1,shape=shapeNorm[j],scale=scaleNorm[j]))
      } else if( length(intersect(j, expNeg)) == 1){
        # Ziehen von Werten aus Exponentialverteilungen fuer negative Werte
        mu[j] <- 1/rgamma(n=1,shape=shapeExpNeg[j-max(normNull)],scale=scaleExpNeg[j-max(normNull)])
      } else if( length(intersect(j, expPos)) == 1){
        # Ziehen von Werten aus Exponentialverteilungen fuer positive Werte
        mu[j] <- 1/rgamma(n=1,shape=shapeExpPos[j-max(expNeg)],scale=scaleExpPos[j-max(expNeg)])
      } else if( length(intersect(j, gamNeg)) == 1){
        # Ziehen von Werten aus Gammaverteilungen fuer positive Werte
        scaleGam[j-max(c(0,expPos,expNeg,normNull))] <- 1/rgamma(n=1,shape=scaleGamNegAlpha[j-max(c(0,expPos,expNeg,normNull))],scale=scaleGamNegBeta[j-max(c(0,expPos,expNeg,normNull))])
        saSh <- .sampleShape(alpha_alt=shapeGam[j-max(c(0,expPos,expNeg,normNull))],beta_r=scaleGam[j-max(c(0,expPos,expNeg,normNull))],a=shapeGamNegAlpha0[j-max(c(0,expPos,expNeg,normNull))],b=shapeGamNegBeta0[j-max(c(0,expPos,expNeg,normNull))],werte=abs(werte[unter_null_ind]),sdShape=sdShape)
        shapeGam[j-max(c(0,expPos,expNeg,normNull))] <- saSh$alpha_new
        acceptanceProb[j-max(c(0,expPos,expNeg,normNull))] <- saSh$alpha_accept
      } else if(length(intersect(j, gamPos)) == 1){
        # Ziehen von Werten aus Gammaverteilungen fuer positive Werte
        scaleGam[j-max(c(0,expPos,expNeg,normNull))] <- 1/rgamma(n=1,shape=scaleGamPosAlpha[j-max(c(0,expPos,expNeg,normNull,gamNeg))],scale=scaleGamPosBeta[j-max(c(0,expPos,expNeg,normNull,gamNeg))])
        saSh <- .sampleShape(alpha_alt=shapeGam[j-max(c(0,expPos,expNeg,normNull))],beta_r=scaleGam[j-max(c(0,expPos,expNeg,normNull))],a=shapeGamPosAlpha0[j-max(c(0,expPos,expNeg,normNull,gamNeg))],b=shapeGamPosBeta0[j-max(c(0,expPos,expNeg,normNull,gamNeg))],werte=werte[ueber_null_ind],sdShape=sdShape)
        shapeGam[j-max(c(0,expPos,expNeg,normNull))] <- saSh$alpha_new
        acceptanceProb[j-max(c(0,expPos,expNeg,normNull))] <- saSh$alpha_accept
      }
        
    } else {
   
   # wenn Komponente "leer" ist: ziehe aus Prior
      if(length(intersect(j, normNull)) == 1){
        # Normalverteilungen mit festgelegtem Erwartungswert
        mu[j] <- 0
        Sigma[j] <- sqrt(1/rgamma(n=1,shape=shapeNorm0[j],scale=scaleNorm0[j])) 
      } else if( length(intersect(j, expNeg)) == 1){
        # Exponentialverteilungen fuer negative Werte
        mu[j] <- 1/rgamma(n=1,shape=shapeExpNeg0[j-max(normNull)],scale=scaleExpNeg0[j-max(normNull)])    
      } else if( length(intersect(j, expPos)) == 1){
        # Exponentialverteilungen fuer positive Werte
        mu[j] <- 1/rgamma(n=1,shape=shapeExpPos0[j-max(expNeg)],scale=scaleExpPos0[j-max(expNeg)])  
      } else if( length(intersect(j, gamNeg)) == 1){
        # Gammaverteilungen fuer positive Werte
        scaleGam[j-max(c(0,expPos,expNeg,normNull))] <- 1/rgamma(n=1,shape=scaleGamNegAlpha0[j-max(c(0,expPos,expNeg,normNull))],scale=scaleGamNegBeta0[j-max(c(0,expPos,expNeg,normNull))])
        shapeGam[j-max(c(0,expPos,expNeg,normNull))] <- rgamma(n=1,shape=shapeGamNegAlpha0[j-max(c(0,expPos,expNeg,normNull))],scale=shapeGamNegBeta0[j-max(c(0,expPos,expNeg,normNull))])
      } else if( length(intersect(j, gamPos)) == 1){
        # Gammaverteilungen fuer positive Werte
        scaleGam[j-max(c(0,expPos,expNeg,normNull))] <- 1/rgamma(n=1,shape=scaleGamPosAlpha0[j-max(c(0,expPos,expNeg,normNull,gamNeg))],scale=scaleGamPosBeta0[j-max(c(0,expPos,expNeg,normNull,gamNeg))])
        shapeGam[j-max(c(0,expPos,expNeg,normNull))] <- rgamma(n=1,shape=shapeGamPosAlpha0[j-max(c(0,expPos,expNeg,normNull,gamNeg))],scale=shapeGamPosBeta0[j-max(c(0,expPos,expNeg,normNull,gamNeg))])
      }    
    }    
  }
  RVAL <- list(mu=mu,Sigma=Sigma,scaleGam=scaleGam,shapeGam=shapeGam,acceptanceProb=acceptanceProb,shapeNorm=shapeNorm,scaleNorm=scaleNorm,shapeExpNeg=shapeExpNeg,scaleExpNeg=scaleExpNeg,shapeExpPos=shapeExpPos,scaleExpPos=scaleExpPos,scaleGamNegAlpha=scaleGamNegAlpha,scaleGamNegBeta=scaleGamNegBeta,scaleGamPosAlpha=scaleGamPosAlpha,scaleGamPosBeta=scaleGamPosBeta)

  return(RVAL)
}


.sampleMixture <- function(k,zdp,alphadir){
 
  n1 <- rep(1,k)
  for(j in 1:k){
    i=which(zdp==j); n1[j]=length(i)
  }
  
  Pi <- as.numeric(rdirichlet(n=1,alpha=rep(alphadir/k,k)+n1))
  #Pi <- rdirichlet(n=1,alpha=rep(alphadir,k)+n1)

  return(list(Pi=Pi,V=NULL))
}

.sampleMixtureTDP <- function(k,zdp,alphadir){
  n1 <- rep(1,k)
  for(j in 1:k){
    i=which(zdp==j); n1[j]=length(i)
  }

  V <- Pi <- numeric(k)
  for(j in 1:(k-1)){
    gamma1 <- 1+n1[j]
    gamma2 <- alphadir+sum(n1[(j+1):k])
    V[j]   <- rbeta(n=1,shape1=gamma1,shape2=gamma2)
    Pi[j]  <- V[j]*prod(1-V[1:(j-1)])
  }
  Pi[k] <- 1-sum(Pi[1:(length(Pi)-1)])
  
  return(list(Pi=Pi,V=V))
}

.sampleShape <- function(alpha_alt,beta_r,a,b,werte,sdShape){
  ############### pi_Y #################
  likelihood_alt <- sum(log(dgamma(x=werte,shape=alpha_alt,scale=beta_r)))
  #2) Faktor alpha0
  alpha_alt_prior <- log(dgamma(x=alpha_alt,shape=a,scale=b))
  pi_X_log <- likelihood_alt+alpha_alt_prior
  
  ############################ Proposal q(Y|X) ############################
  alpha_neu_prop <- abs(rnorm(n=1,mean=alpha_alt, sd=sdShape))
  q_yx_log <- log(dnorm(x=alpha_neu_prop,mean=alpha_alt, sd=sdShape))
  
  ############### pi_Y #################
  #%1) Faktor Distanzen
  likelihood_neu <- sum(log(dgamma(x=werte,shape=alpha_neu_prop,scale=beta_r)))
  #%2) Faktor alpha0
  alpha_neu_prior <- log(dgamma(x=alpha_neu_prop,shape=a,scale=b))
  pi_Y_log <- likelihood_neu+alpha_neu_prior
  
  ######### Proposal q(X|Y) ############
  alpha_alt_prop <- abs(rnorm(n=1,mean=alpha_neu_prop,sd=sdShape))
  q_xy_log <- log(dnorm(x=alpha_alt_prop,mean=alpha_neu_prop, sd=sdShape))
  
  #% Akzeptanzschritt
  shapeAlpha_akzep <- min(0, pi_Y_log+q_xy_log-pi_X_log-q_yx_log)
  u <- log(runif(1))

  if(!is.na(shapeAlpha_akzep) && u < shapeAlpha_akzep){
    alpha_neu <- alpha_neu_prop
    alpha_accept <- 1
  } else {
    alpha_neu <- alpha_alt
    alpha_accept <- 0
  }
  RVAL <- list(alpha_new=alpha_neu,alpha_accept=alpha_accept)
}

.sampleAlpha <- function(alpha_alt,a,b,werte,sdAlpha){ # werte = p's (Anteile)
  
  ##############################################################
  
  ############ Proposal q(Y|X) ###############
  alpha_neu_prop <- abs(rnorm(n=1,mean=alpha_alt, sd=sdAlpha))
  q_yx_log <- log(dnorm(x=alpha_neu_prop,mean=alpha_alt,sd=sdAlpha))
  
  ############# pi_X ###############
  likelihood_alt <- sum(log(ddirichlet(x=werte,alpha=rep(alpha_alt,length(werte)))))
  
  alpha_alt_prior <- log(dgamma(x=alpha_alt,shape=a,scale=b))
  pi_X_log <- likelihood_alt+alpha_alt_prior
  
  #############################################################

  ############ Proposal q(X|Y) ###############
  alpha_alt_prop <- abs(rnorm(n=1,mean=alpha_neu_prop,sd=sdAlpha))
  q_xy_log <- log(dnorm(x=alpha_alt_prop,mean=alpha_neu_prop, sd=sdAlpha))
  
  ############# pi_Y ###############
  likelihood_neu <- sum(log(ddirichlet(x=werte,alpha=rep(alpha_neu_prop,length(werte)))))
  
  alpha_neu_prior <- log(dgamma(x=alpha_neu_prop,shape=a,scale=b))
  pi_Y_log <- likelihood_neu+alpha_neu_prior
  
  ############################################# Akzeptanzschritt #############################################
  DirAlpha_akzep <- min(0, pi_Y_log+q_xy_log-pi_X_log-q_yx_log)
  u <- log(runif(1))

  if(!is.na(DirAlpha_akzep) && u < DirAlpha_akzep){
    alpha_neu <- alpha_neu_prop
    alpha_accept <- 1
  } else {
    alpha_neu <- alpha_alt
    alpha_accept <- 0
  }
  RVAL <- list(alpha_new=alpha_neu,alpha_accept=alpha_accept)
}


setMethod("bayesMixModel", signature=c(z="numeric"), .bayesMixModel)
