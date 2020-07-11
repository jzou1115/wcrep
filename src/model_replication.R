library(mvtnorm)

args = commandArgs(trailingOnly=TRUE)
filename<-args[1]
n1 <- as.numeric(args[2])
n2 <- as.numeric(args[3])
outfile <- args[4]
data<-read.table(filename, header=T)



#estimating variance of s1, accounting for missing data
MLE<-function(var, s1){
	estimate<-0
	for(i in s1){
		estimate<-estimate+log(dnorm(i, mean=0, sd=sqrt(var)))
	}

	return(estimate)
}

#variance explained by genetics
expected_mean_ratio<-function(sigma_g_2, maxVar, n1, n2){
	sqrt(n1)*sqrt(n2)*sigma_g_2/maxVar
}

#estimate variance components
estimate_parameters <- function(s1, s2, n1, n2){
	#sigma_g, sigma_c1, sigma_c2 initial values
	parameters <- c(1e-28, 1e-28, 1e-28)

	######### estimate total variance in s1  #############
	max<- -Inf
	maxVar<- -Inf
	for(i in seq(from=0.1, to=100, by=.01)){
		temp<-MLE(i, s1)
		if (temp>max) {
			max<-temp
 			maxVar<-i
 		}
	}



	######### estimate sigma_g  #############
	min_rms <- Inf
	var_g_est2 <- 1e-28


	for(i in seq(from=1e-28,to=1, by=0.0001)){
		ratio<-expected_mean_ratio(i, maxVar, n1, n2)
		expected_s2 <- s1*ratio
  		cur_rms <- sqrt(sum((expected_s2-s2)^2))
  		if (cur_rms < min_rms){
    			min_rms<- cur_rms
    			var_g_est2 <- i
  		}
	}	
	#set sigma_g parameter
	parameters[1] <- var_g_est2

	######### estimate sigma_c1 #############
	c1_est2<-(maxVar-1-n1*var_g_est2)/n1


	if(c1_est2<0){
 		c1_est2<-1e-28
 		parameters[1] <- maxVar-1 #c1_est2 sampling error, so set to zero
	}
	#set sigma_c1 parameter
	parameters[2] <- c1_est2

	######### estimate sigma_c2 #############
	max<- - Inf
	c2_est<- 1e-28
	stats <- cbind(s1, s2)
	lik <- c()
	for(i in seq(from=1e-28,to=1, by=0.0001)){
		temp <- MLE_joint_probability(var_g_est2, c1_est2, i, n1, n2, stats)
		if(temp>max){
			max<-temp
			c2_est<-i
		}
		lik <- c(lik, temp)
	}

	#set sigma_c2 parameter
	parameters[3] <- c2_est

	return(parameters)
}


MLE_joint_probability<-function(var_g, var_c1, var_c2, n1, n2, stats){
	cov_matrix=matrix(data=NA, nrow=2, ncol=2)
	cov_matrix[1,1]=n1*var_g+n1*var_c1+1
	cov_matrix[1,2]=sqrt(n1)*var_g*sqrt(n2)
	cov_matrix[2,1]=sqrt(n1)*var_g*sqrt(n2)
	cov_matrix[2,2]=n2*var_g+n2*var_c2+1

	mean_matrix=matrix(data=NA, nrow=2, ncol=1)
	mean_matrix[1,1]=0
	mean_matrix[2,1]=0

	estimate<-0

	prob<-dmvnorm(x=stats, mean = mean_matrix, sigma = cov_matrix, log = TRUE)

	for(i in prob){
 		estimate<-estimate + i
  	}
  	return(estimate)
}



predict_replication <- function(mean, sd, z){
	#calculate predicted replication rate
	lower <-  pnorm(z, mean, sd)
	upper <- 1- pnorm(-z, mean, sd)
	return(lower + upper)
}


#calculate conditional probability of s2 > t | s1 = x under model with no confounding
calc_conditional_no_confounding <- function(s1, n1, n2, var_g, z){


	mean <- (sqrt(n1)*sqrt(n2)*var_g)/(1+n1*var_g) *s1
	sd <- sqrt(n2*var_g+1-((n1*var_g*var_g*n2)/(n1*var_g+1)))
	lower <-  pnorm(z, mean, sd)
	upper <- 1- pnorm(-z, mean, sd)
	return(lower + upper)
}

#calculate conditional probability of s2 > t | s1 = x under model with confounding
calc_conditional_with_confounding <- function(s1, n1, n2, var_g, var_c1, var_c2, z){

	mean <- (sqrt(n1)*sqrt(n2)*var_g)/(1+n1*var_g+n1*var_c1) *s1
	sd <- sqrt(n2*var_g+n2*var_c2+1-((n1*var_g^2*n2)/(n1*var_g+n1*var_c1+1)))
	lower <-  pnorm(z, mean, sd)
	upper <- 1-pnorm(-z, mean, sd)
	return(lower + upper)

}

#calculate predicted replication rate under model with no confounding
predict_no_confounding <- function(data, z, sigma_g){
	
	predicted_replication <- 0
	for(i in 1:nrow(data)){
		predicted_replication <- predicted_replication +  calc_conditional_no_confounding(data[i,1], n1, n2, sigma_g, z)

	}

	return(predicted_replication/nrow(data))
}

#calculate predicted replication rate under model with confounding
predict_with_confounding <- function(data, z, sigma_g, sigma_c1, sigma_c2){

	predicted_replication <- 0
	for(i in 1:nrow(data)){
		 predicted_replication <- predicted_replication +  calc_conditional_with_confounding(data[i,1], n1, n2, sigma_g, sigma_c1, sigma_c2, z)

	}

	return(c(predicted_replication/nrow(data), sigma_g, sigma_c1, sigma_c2))
}

calcReplication <- function(data, z){
	count <- 0
	for(i in 1:nrow(data)){
		if(abs(data[i,2])>abs(z)){
			count  <- count + 1
		}
	}
	return(count/nrow(data))
}




run_mouse_analysis <- function(sig){

#	stopifnot(nrow(sig)>1) #cannot apply method if there is not at least one significant variant

	#z-score significance threshold for replication study
	threshold2 <- 0.05/nrow(sig)
	z2 <- qnorm(threshold2/2)
	#true replication rate
	r_true <- calcReplication(sig, z2)

	
	#estimate parameters
	s1 <- sig[,1]
	s2 <- sig[,2]

	parameters <- estimate_parameters(s1, s2, n1, n2)
	sigma_g <- parameters[1]
	sigma_c1 <- parameters[2]
	sigma_c2 <- parameters[3]

	#predicted replication rate with no confounding
	pr_no_confounding <- predict_no_confounding(sig, z2, sigma_g)

	#predicted replication rate with confounding, sigma_g1, sigma_c1, and sigma_c2
	pr_with_confounding <- predict_with_confounding(sig, z2, sigma_g, sigma_c1, sigma_c2)


	#output predicted replication rates and parameters
	write.table(c(r_true, pr_no_confounding, pr_with_confounding), file=outfile, row.names=F, col.names=F, quote=F)


}

run_mouse_analysis(data)
	
