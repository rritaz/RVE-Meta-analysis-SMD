###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
### Functions to compute the degrees of freedom for HC2 
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###

var_hc2_term1 <- function(c_hc2, betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + (c_hc2[i]^2)*((1-2*betas[i])^2)*(sigmas[i]+tauHat_value)^2
    }
    return(as.numeric(2*numerator))
} # term 1

var_hc2_term2_1 <- function(c_hc2){
    numerator <- 0 
    for(i in 1:length(c_hc2)){
        numerator <- numerator + c_hc2[i]^2
    }
    return(as.numeric(2*numerator))
} # term 2_1

var_hc2_term2_2 <- function(c_hc2, betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
    }
    return(as.numeric(numerator^2))
} # term 2_2

var_hc2_term3_1 <- function(c_hc2, betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + c_hc2[i]^2*(1-2*betas[i])*(sigmas[i]+tauHat_value)
    }
    return(as.numeric(4*numerator))
} # term 3_1

var_hc2_term3_2 <- function(betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
    }
    return(as.numeric(numerator))
} # term 3_2

var_hc2_term4_1_1 <- function(c_hc2){
    numerator <- 0 
    for(i in 1:length(c_hc2)){
        numerator <- numerator + c_hc2[i]
    }
    return(as.numeric(numerator^2))
} # term 4_1_1

var_hc2_term4_1_2 <- function(c_hc2, kVal){
    numerator <- 0 
    for(i in 1:length(c_hc2)){
        numerator <- numerator + c_hc2[i]^2
    }
    return(as.numeric(numerator))
} # term 4_1_2

var_hc2_term4_2 <- function(betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
    }
    return(as.numeric(2*numerator^2))
} # term 4_2

var_hc2_term5_1 <- function(c_hc2, betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + c_hc2[i]*betas[i]^2*(sigmas[i]+tauHat_value)^2
    }
    return(as.numeric(4*numerator))
} # term 5_1

var_hc2_term5_2 <- function(c_hc2, sigmas){
    numerator <- 0
    for (i in 1:length(sigmas)){
        numerator <- numerator + c_hc2[i]
    }
    return(as.numeric(numerator))
} # term 5_2

var_hc2_term6 <- function(c_hc2, betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + c_hc2[i]*betas[i]*(sigmas[i]+tauHat_value)
    }
    return(as.numeric(4*numerator^2))
} # term 6

var_hc2_term7 <- function(c_hc2, betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + c_hc2[i]^2*betas[i]^2*(sigmas[i]+tauHat_value)^2
    }
    return(as.numeric(8*numerator))
} # term 7

var_hc2_term8_1 <- function(betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
    }
    return(as.numeric(8*numerator))
} # term 8_1

var_hc2_term8_2_1 <- function(c_hc2, betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + c_hc2[i]*betas[i]*(sigmas[i]+tauHat_value)
    }
    return(as.numeric(numerator))
} # term 8_2_1

var_hc2_term8_2_2 <- function(c_hc2){
    numerator <- 0
    for (i in 1:length(c_hc2)){
        numerator <- numerator + c_hc2[i]
    }
    return(as.numeric(numerator))
} # term 8_2_2

var_hc2_term8_3 <- function(c_hc2, betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + c_hc2[i]^2*betas[i]*(sigmas[i]+tauHat_value)
    }
    return(as.numeric(numerator))
} # term 8_3

mean_hc2_term1 <- function(c_hc2, betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + c_hc2[i]*(1-2*betas[i])*(sigmas[i]+tauHat_value)
    }
    return(as.numeric(numerator))
} # term 1

mean_hc2_term2_1 <- function(c_hc2, kVal){
    numerator <- 0 
    for(i in 1:length(c_hc2)){
        numerator <- numerator + c_hc2[i]
    }
    return(as.numeric(numerator))
} # term 2_1

mean_hc2_term2_2 <- function(betas, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(betas)){
        numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
    }
    return(as.numeric(numerator))
} # term 2_2


