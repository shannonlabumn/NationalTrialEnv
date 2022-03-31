getAF <- function(X, x){ # X nxm matrix, x ploidy
	return(apply(X, 2, function(y){ mean(y, na.rm=T)/x}))
}

whichMA <- function(X){ # X nxm marker matrix, # copies of the minor allele by indiv
		ma <- as.numeric(apply(X, 2, function(x) {y <- mean(x, na.rm=T); 
		if(y >2){ paste(0)
			}else{ 
				if(y<2){ paste(4)
					}else{ paste(NA)}}}))
	#return(ma)
	
	
	MAct <- matrix(NA, nrow=dim(X)[1], ncol=dim(X)[2]); rownames(MAct) <- rownames(X); colnames(MAct)<- colnames(X)
	for(i in 1:length(ma)){
		if(!is.na(ma[i])){
				if (ma[i]==4) { MAct[,i] <-X[,i]
					}else{MAct[, i] <-4-X[,i]
						} 
				}
			}
			
	
 return(list(ma, MAct))
}

# position is vector same length of marker trait, size is width in Mbp of the bin
binpos <- function(pos, trait, size=1){
	pos2 <- pos*1e-6
	right.break <- seq(size, max(pos2), by=size)
	bins <- c(0, right.break)
	
	median <- rep(NA, length(bins))
	quant25 <- rep(NA, length(bins))
	quant75 <- rep(NA, length(bins))
	

	for(n in 2:length(bins)){
		median[n] <- median(trait[intersect(which(pos2 > bins[n-1]), which(pos2 < bins[n]))] )
		quant25[n]<- quantile(trait[intersect(which(pos2 > bins[n-1]), which(pos2 < bins[n]))])[2]
		quant75[n] <- quantile(trait[intersect(which(pos2 > bins[n-1]), which(pos2 < bins[n]))])[4]
	}
	
	ANS <- list(bins, median, quant25, quant75)
	return(ANS)
	
}

# sliding window
winpos <- function(pos, trait, size=1, step=0.5){
	pos2 <- pos*1e-6
	breaks <- seq(from=min(pos2), to=max(pos2)-size, by=step)
	median <- rep(NA, length(breaks))
	sd <- rep(NA, length(breaks))
	for(i in 1:length(breaks)){
		median[i] <- median(trait[intersect(which(pos2 > breaks[i]), which(pos2 < breaks[i]+size))] )
		sd[i]<- sd(trait[intersect(which(pos2 > breaks[i]), which(pos2 < breaks[i]+size))], na.rm=T)
	}

	
	ANS <- list(breaks, median, sd)
	return(ANS)
	
}



# tetraploid rho 
# Meirmans PG, Liu S (2018) Analysis of Molecular Variance (AMOVA) for Autopolyploids. Frontiers in Ecology and Evolution 6:66
# X n x m genotype matrix, pop list of population index, x ploidy for all if length(x)=1 else vector of length n of ploidy for each indiv
get.rho <- function(X, pop, x=4){
	X2 <- X/x # scale genotypes 0-1, X2[n, m]= frequency of allele m in individual n
	X3 <- 1-X2
	M <- dim(X)[2]


	rho <- rep(NA, M)
	Fst <- rep(NA, M)
	Ho <- rep(NA, M)
	He <- vector("list", M)
	Ht <- rep(NA, M)
	Hs <- rep(NA, M)
		
	for(m in 1:M){
		ix <- which(!is.na(X2[ , m])) # non-missing genotypes only
		N <- length(ix) # N is the number of individuals without missing data
		
		HO <- 1/N* sum( 1-(X2[ ix, m]^2 +  X3[ix, m]^2) ) # (8) from Meirmans and Liu 2018
    
		# each pop gets an HE, where length(n) is the number of indiv without missing data in the subpop (9) from Meirmans and Liu
		HE <- lapply(pop, function(p){ n <- intersect(which(rownames(X2) %in% p), ix);
			1- ( ((sum(X2[n , m])/ length(n))^2) + ((sum(X3[n , m])/ length(n))^2)) }) 
		pop.size <- unlist(lapply(pop, function(foo){length(intersect(which(rownames(X2) %in% foo), ix))}))
		# HT is the HE for all individuals at once 
		HT <- 1- ( (sum(X2[ix , m])/ N)^2 + (sum(X3[ix , m])/ N)^2) 

		# HS is the average HE for all subpops
		HS <- sum(unlist(HE)*pop.size)/sum(pop.size)
		

		rho[m] <- 	(HT-HS)/(HT -HO) # (14)
		Fst[m] <- (HT-HS)/HT
		Ho[m] <- HO
		He[[m]] <- unlist(HE)
		Ht[m] <- HT
		Hs[m] <- HS
	}
	names(rho) <- colnames(X)
	return(list(rho=rho, Fst=Fst, Ho=Ho, He=He, Ht=Ht, Hs=Hs))
}
	
# a function allowing for comparisons across ploidy, eq 11 and 12 for Meirmans and Liu
# X is nxm, x is a nx1 vector of ploidy for each individual
get.divers <- function(X, pop, x){ 
	X2 <- X/x # scale genotypes 0-1, X2[n, m]= frequency of allele m in individual n
	X3 <- 1-X2
	M <- dim(X)[2]


	rho <- rep(NA, M)
	Fst <- rep(NA, M)
	Ho <- rep(NA, M)
	He <- vector("list", M)
	Ht <- rep(NA, M)
	Hs <- rep(NA, M)
		
	for(m in 1:M){
		ix <- which(!is.na(X2[ , m])) # non-missing genotypes only
		N <- length(ix) # N is the number of individuals without missing data
		
		HO <- 1/N* sum( (1-(X2[ ix, m]^2 +  X3[ix, m]^2)* x[ix]/(x[ix]-1) )) # (11) from Meirmans and Liu 2018
	
		# each pop gets an HE, where length(n) is the number of indiv without missing data in the subpop (12) from Meirmans and Liu
		HE <- lapply(pop, function(p){ n <- intersect(which(rownames(X2) %in% p), ix);
			1- ( (sum(x[n]*X2[n , m])/ sum(x[n]))^2 + (sum(x[n]*X3[n , m])/ sum(x[n]))^2) }) 
		
		# HT is the HE for all individuals at once 
		HT <- 1- ( (sum(x[ix]*X2[ix , m])/ sum(x[ix]))^2 + (sum(x[ix]*X3[ix , m])/ sum(x[ix]))^2) 

		# HS is the average HE for all subpops
		HS <- mean(unlist(HE))
		

		rho[m] <- 	(HT-HS)/(HT -HO) # (14)
		Fst[m] <- (HT-HS)/HT
		Ho[m] <- HO
		He[[m]] <- unlist(HE)
		Ht[m] <- HT
		Hs[m] <- HS
	}
	names(rho) <- colnames(X)
	return(list(rho=rho, Fst=Fst, Ho=Ho, He=He, Ht=Ht, Hs=Hs))

}





