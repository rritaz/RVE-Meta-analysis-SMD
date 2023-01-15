# load libraries
library(tidyverse)
library(readxl)

# Read in excel file containing n values
n_values <- read_excel("~/Desktop/n_values.xlsx", col_names = FALSE) 
delta <- 0.5
A <- 0.95
B <- 1.05

var_cohenD <- function(n) {
    output <- as.numeric((2 / n) + ((delta^2)/(4*as.numeric(n))))
    return (output)
} # variance of estimate of cohen's d

get_delta <- function(tau_sq) {
    return(rnorm(1, delta, sqrt(tau_sq)))
} # generate delta

get_d <- function(n, tau_sq, delta) {
    return(rnorm(1, delta, sqrt(var_cohenD(n))))
} # generate d

#----------------------------------------------------------------------------#

# Generating data

generate_deltas_forCol <- function(kVal, currentTau){
    deltas <- vector()
    for(i in 1:kVal) {
        deltas[i] <- get_delta(currentTau)
    }
    return(deltas)
}

generate_d_forCol <- function(colNum, kVal, currentTau, deltas){
    ds <- vector()
    for(i in 1:kVal) {
        ds[i] <- get_d(n_values[i, colNum], currentTau, deltas[i])
    }
    return(ds)
}

generate_sigma_forCol <- function(colNum, kVal, ds){
    sigmas <- vector()
    for(i in 1:kVal) {
        sigmas[i] <- (2/as.numeric(n_values[i, colNum]))+(ds[i]^2/(4*as.numeric(n_values[i, colNum])))
    }
    return(sigmas)
}

#----------------------------------------------------------------------------#

# Meta-analysis 

weight <- function(sigma) {
    output <- 1 / sigma
    return(output)
} # fixed effect weights

weightTrue <- function(tau_sqVal, sigma) {
    output <- 1 / (sigma + tau_sqVal)
    return(output)
} # random effect weights

sumW <- function(kVal, powToRaise, sigmas) {
    sum_w <- 0
    for (i in 1:kVal) {
        sum_w <- sum_w + weight(sigmas[i])^powToRaise
    }
    return(as.numeric(sum_w))
} # sum of fixed effect weights

sumTrueW <- function(kVal, tauHat_value, sigmas) {
    sum_w <- 0
    for (i in 1:kVal) {
        sum_w <- sum_w + weightTrue(tauHat_value, sigmas[i])
    }
    return(as.numeric(sum_w))
} # sum of random effect weights

tau_sq_fun <- function(colNum, kVal) {
    total <- 0
    for(i in 1:kVal) {
        total <- total + var_cohenD(n_values[i, colNum])
    }
    average <- total/kVal
    return(as.numeric(average))
} # true heterogeneity = arithmetic average of the within study variances

delta_hat <- function(kVal, tauHat_value, ds, sigmas) {
    sum_Numerator <- 0
    for (i in 1:kVal) {
        sum_Numerator <- sum_Numerator + weightTrue(tauHat_value, sigmas[i]) * ds[i]
    }
    return (as.numeric(sum_Numerator/sumTrueW(kVal, tauHat_value, sigmas)))
} # calculating meta-analytic estimate in RE procedures

#----------------------------------------------------------------------------#

Q <- function(ds, deltaHat_FE_val, kVal, sigmas){
    numerator <- 0
    for (i in 1:length(ds)){
        numerator <- numerator + weight(sigmas[i])*(ds[i]-deltaHat_FE_val)^2
    }
    return(numerator)
} # weighted sum of squares Q using fixed effect weights

delta_hatFE <- function(kVal, ds, sigmas) {
    sum_Numerator <- 0
    for (i in 1:kVal) {
        sum_Numerator <- sum_Numerator + weight(sigmas[i]) * ds[i]
    }
    return (as.numeric(sum_Numerator/sumW(kVal, 1, sigmas)))
} # delta hat with fixed effect weights

denomConst <- function(kVal, sigmas) {
    toReturn <- (sumW(kVal, 1, sigmas) - sumW(kVal, 2, sigmas)/sumW(kVal, 1, sigmas))
    return(as.numeric(toReturn))
} # constant function of weights in DL denominator of tau hat squared

tauHat_DL <- function(ds, deltaHat_FE_val, kVal, sigmas){
    toReturn <- ((Q(ds, deltaHat_FE_val, kVal, sigmas)-(kVal-1))/denomConst(kVal, sigmas))
    if (toReturn<0){
        return(0)
    }else{
        return(as.numeric(toReturn))
    }
} # DerSimonian and Laird tau hat

betas <- function(kVal, tauHat_value, sigmas){
    beta <- vector()
    for(i in 1:kVal) {
        beta[i] <- weightTrue(tauHat_value, sigmas[i])/sumTrueW(kVal, tauHat_value, sigmas)
    }
    return(as.numeric(beta))
} # normalized weights

sumOfBetasSQ <- function(betas, kVal){
    sum_betas <- 0
    for(i in 1:kVal){
        sum_betas <- sum_betas + betas[i]^2
    }
    return(as.numeric(sum_betas))
} # sum of normalized weights

ksi_i_Beta <- function(betas, kVal){
    ksi_i_Beta <- vector()
    for(i in 1:kVal){
        ksi_i_Beta[i] <- betas[i]+((-betas[i]+(betas[i])^2)/(1-sumOfBetasSQ(betas, kVal)))
    }
    return(as.numeric(ksi_i_Beta))
}

QofBeta_term2 <- function(kVal, sigmas, ksi_i_Beta){
    term2 <- 0
    for(i in 1:kVal){
        term2 <- term2 + ksi_i_Beta[i]*sigmas[i]
    }
    return(as.numeric(term2))
} 

Q_RE <- function(ds, deltaHat_val, kVal, betas){
    numerator <- 0
    for (i in 1:length(ds)){
        numerator <- numerator + betas[i]*(ds[i]-deltaHat_val)^2
    }
    return(numerator)
} # weighted sum of squares - using normalized weights

RofBeta <- function(betas, sigmas){
    numerator <- 0
    for (i in 1:length(sigmas)){
        numerator <- numerator + betas[i]^2*sigmas[i]
    }
    return(as.numeric(numerator))
} 

hc1 <- function(ds, deltaHat_val, kVal, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(ds)){
        numerator <- numerator + weightTrue(tauHat_value, sigmas[i])^2*(kVal/(kVal-1))*
            (ds[i]-deltaHat_val)^2
    }
    return(numerator/sumTrueW(kVal, tauHat_value, sigmas)^2)
} 

hc2 <- function(ds, deltaHat_val, kVal, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(ds)){
        numerator <- numerator + weightTrue(tauHat_value, sigmas[i])^2*(1/(1-(weightTrue(tauHat_value, sigmas[i])/
                                                                                  sumTrueW(kVal, tauHat_value, sigmas))))*
            (ds[i]-deltaHat_val)^2
    }
    return(numerator/sumTrueW(kVal, tauHat_value, sigmas)^2)
} 

hc3 <- function(ds, deltaHat_val, kVal, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(ds)){
        numerator <- numerator + weightTrue(tauHat_value, sigmas[i])^2*(1/(1-(weightTrue(tauHat_value, sigmas[i])/
                                                                                  sumTrueW(kVal, tauHat_value, sigmas))))^2*
            (ds[i]-deltaHat_val)^2
    }
    return(numerator/sumTrueW(kVal, tauHat_value, sigmas)^2)
} 

c_hc2 <- function(kVal, betas){
    c <- vector()
    for (i in 1:kVal){
        c[i] <- betas[i]^2/(1-betas[i])
    }
    return(c)
} # c for HC2

c_hc3 <- function(kVal, betas){
    c <- vector()
    for (i in 1:kVal){
        c[i] <- betas[i]^2/(1-betas[i])^2
    }
    return(c)
} # c for HC3

main <- function(){
    kVal <- 2 # number of studies
    tau_sq <- c(0, 0.5, 1, 2) # heterogeneity
    totalNumCol <- length(n_values[1,]) 
    
    mainResults <- matrix(NA, nrow=length(tau_sq)*totalNumCol, ncol=22)
    currentRow <- 1
            for(j in 1:totalNumCol) {
                for(k in 1:length(tau_sq)){
                    currentTau <- tau_sq[k]*tau_sq_fun(j, kVal)
                    deltas <- generate_deltas_forCol(kVal, currentTau)
                    ds <- generate_d_forCol(j, kVal, currentTau, deltas)
                    sigmas <- generate_sigma_forCol(j, kVal, ds)

                    deltaHat_FE_val <- delta_hatFE(kVal, ds, sigmas) # delta hat value
                    tauHat_value <- tauHat_DL(ds, deltaHat_FE_val, kVal, sigmas)

                    deltaHat_val <- delta_hat(kVal, tauHat_value, ds, sigmas) # random effect meta-analytic estimate
                    
                    betas <- betas(kVal, tauHat_value, sigmas)
                    c_hc2 <- c_hc2(kVal, betas)
                    c_hc3 <- c_hc3(kVal, betas)
                    sumOfBetasSQ <- sumOfBetasSQ(betas, kVal)
                    lambdaOfBeta <- sumOfBetasSQ(betas, kVal)/(1-sumOfBetasSQ(betas, kVal))
                    ksi_i_Beta <- ksi_i_Beta(betas, kVal)
                    
                    QofBeta <- lambdaOfBeta*Q_RE(ds, deltaHat_val, kVal, betas) + QofBeta_term2(kVal, sigmas, ksi_i_Beta)
                    RofBeta <- RofBeta(betas, sigmas)
                    LofBeta <- min(1, max(0, ((QofBeta/RofBeta) - A)/(B-A)))
                    qofBeta <- LofBeta*QofBeta + (1-LofBeta)*RofBeta # truncated Hartung variance estimator
                    
                    Var_QofBeta <- (lambdaOfBeta^2)*(1/sumTrueW(kVal, tauHat_value, sigmas))^2*(2*(kVal-1))
                    Var_qofBeta <- LofBeta^2*Var_QofBeta # ignoring variance of within study variances
                    df <- 2*((qofBeta^2)/Var_qofBeta)
                    
                    varDL_deltaHat <- 1/sumTrueW(kVal, tauHat_value, sigmas) # standard variance estimator
                    varHK_deltaHat <- (1/((kVal-1)))*Q_RE(ds, deltaHat_val, kVal, betas)
                    
                    hc1 <-hc1(ds, deltaHat_val, kVal, sigmas, tauHat_value)
                    hc2 <-hc2(ds, deltaHat_val, kVal, sigmas, tauHat_value)
                    hc3 <-hc3(ds, deltaHat_val, kVal, sigmas, tauHat_value)
                    
                    #-----------------------------------------------------------------------------------#
                    
                    # Degrees of freedom adjustment for HC2
                    # Variance of HC2
                    
                    var_hc2_term1 <- function(c_hc2, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc2[i]^2*(1-2*betas[i])^2*(sigmas[i]+tauHat_value)^2
                        }
                        return(2*numerator)
                    } # term 1
                    
                    var_hc2_term2_1 <- function(kVal, c_hc2){
                        numerator <- 0 
                        for(i in 1:length(kVal)){
                            numerator <- numerator + c_hc2[i]^2
                        }
                        return(2*numerator)
                    } # term 2_1
                    
                    var_hc2_term2_2 <- function(c_hc2, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
                        }
                        return(numerator^2)
                    } # term 2_2
                    
                    var_hc2_term2 <- var_hc2_term2_1(kVal, c_hc2)*var_hc2_term2_2(c_hc2, betas, sigmas, tauHat_value)
                    
                    var_hc2_term3_1 <- function(c_hc2, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc2[i]^2*(1-betas[i])*(sigmas[i]+tauHat_value)
                        }
                        return(4*numerator)
                    } # term 3_1
                    
                    var_hc2_term3_2 <- function(betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
                        }
                        return(numerator)
                    } # term 3_2
                    
                    var_hc2_term3 <- var_hc2_term3_1(c_hc2, betas, sigmas, tauHat_value)*
                        var_hc2_term3_2(betas, sigmas, tauHat_value)
                    
                    
                    var_hc2_term4_1_1 <- function(c_hc2, kVal){
                        numerator <- 0 
                        for(i in 1:length(kVal)){
                            numerator <- numerator + c_hc2[i]
                        }
                        return(numerator^2)
                    } # term 4_1_1
                    
                    var_hc2_term4_1_2 <- function(c_hc2, kVal){
                        numerator <- 0 
                        for(i in 1:length(kVal)){
                            numerator <- numerator + c_hc2[i]^2
                        }
                        return(numerator)
                    } # term 4_1_2
                    
                    var_hc2_term4_1 <- 2*(var_hc2_term4_1_1(c_hc2, kVal)-var_hc2_term4_1_2(c_hc2, kVal))
                    
                    var_hc2_term4_2 <- function(betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
                        }
                        return(numerator^2)
                    } # term 4_2
                    
                    var_hc2_term4<- var_hc2_term4_1*var_hc2_term4_2(betas, sigmas, tauHat_value)
                    
                    var_hc2_term5_1 <- function(c_hc2, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc2[i]*betas[i]^2*(sigmas[i]+tauHat_value)^2
                        }
                        return(4*numerator)
                    } # term 5_1
                    
                    var_hc2_term5_2 <- function(c_hc2){
                        numerator <- 0
                        for (i in 1:length(sigmas)){
                            numerator <- numerator + c_hc2[i]
                        }
                        return(numerator)
                    } # term 5_2
                    
                    var_hc2_term5 <- var_hc2_term5_1(c_hc2, betas, sigmas, tauHat_value)*
                        var_hc2_term5_2(c_hc2)
                    
                    var_hc2_term6 <- function(c_hc2, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc2[i]*betas[i]*(sigmas[i]+tauHat_value)
                        }
                        return(4*numerator^2)
                    } # term 6
                    
                    var_hc2_term7 <- function(c_hc2, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc2[i]^2*betas[i]^2*(sigmas[i]+tauHat_value)^2
                        }
                        return(8*numerator)
                    } # term 7
                    
                    var_hc2_term8_1 <- function(betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
                        }
                        return(8*numerator)
                    } # term 8_1
                    
                    var_hc2_term8_2_1 <- function(c_hc2, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc2[i]*betas[i]*(sigmas[i]+tauHat_value)
                        }
                        return(numerator)
                    } # term 8_2_1
                    
                    var_hc2_term8_2_2 <- function(c_hc2){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc2[i]
                        }
                        return(numerator)
                    } # term 8_2_2
                    
                    var_hc2_term8_2 <- var_hc2_term8_2_1(c_hc2, betas, sigmas, tauHat_value)*var_hc2_term8_2_2(c_hc2)
                    
                    var_hc2_term8_3 <- function(c_hc2, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc2[i]^2*betas[i]*(sigmas[i]+tauHat_value)
                        }
                        return(numerator)
                    } # term 8_3
                    
                    var_hc2_term8 <- var_hc2_term8_1(betas, sigmas, tauHat_value)*
                        (var_hc2_term8_2-var_hc2_term8_3(c_hc2, betas, sigmas, tauHat_value))
                    
                    var_hc2 <- var_hc2_term1(c_hc2, betas, sigmas, tauHat_value) + var_hc2_term2 + var_hc2_term3 + var_hc2_term4 +
                        var_hc2_term5 + var_hc2_term6(c_hc2, betas, sigmas, tauHat_value) - var_hc2_term7(c_hc2, betas, sigmas, tauHat_value)-
                        var_hc2_term8
                    
                    # Expected value of HC2
                    
                    mean_hc2_term1 <- function(c_hc2, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc2[i]*(1-2*betas[i])*(sigmas[i]+tauHat_value)
                        }
                        return(numerator)
                    } # term 1
                    
                    mean_hc2_term2_1 <- function(c_hc2, kVal){
                        numerator <- 0 
                        for(i in 1:length(kVal)){
                            numerator <- numerator + c_hc2[i]
                        }
                        return(numerator)
                    } # term 2_1
                    
                    mean_hc2_term2_2 <- function(betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
                        }
                        return(numerator)
                    } # term 2_2
                    
                    mean_hc2 <- mean_hc2_term1(c_hc2, betas, sigmas, tauHat_value)+
                        mean_hc2_term2_1(c_hc2, kVal)*mean_hc2_term2_2(betas, sigmas, tauHat_value) # Exp value of HC2
                    
                    df_hc2 <- 2*mean_hc2^2/var_hc2
                    
                    #-----------------------------------------------------------------------------------#
                    
                    # Degrees of freedom adjustment for HC3
                    # Variance
              
                    var_hc3_term1 <- function(c_hc3, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc3[i]^2*(1-2*betas[i])^2*(sigmas[i]+tauHat_value)^2
                        }
                        return(2*numerator)
                    } # term 1
                    
                    var_hc3_term2_1 <- function(c_hc3, kVal){
                        numerator <- 0 
                        for(i in 1:length(kVal)){
                            numerator <- numerator + c_hc3[i]^2
                        }
                        return(2*numerator)
                    } # term 2_1
                    
                    var_hc3_term2_2 <- function(c_hc3, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
                        }
                        return(numerator^2)
                    } # term 2_2
                    
                    var_hc3_term2 <- var_hc3_term2_1(c_hc3, kVal)*var_hc3_term2_2(c_hc3, betas, sigmas, tauHat_value) # term 2
                    
                    var_hc3_term3_1 <- function(c_hc3, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc3[i]^2*(1-betas[i])*(sigmas[i]+tauHat_value)
                        }
                        return(4*numerator)
                    } # term 3_1
                    
                    var_hc3_term3_2 <- function(betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
                        }
                        return(numerator)
                    } # term 3_2
                    
                    var_hc3_term3 <- var_hc3_term3_1(c_hc3, betas, sigmas, tauHat_value)*
                        var_hc3_term3_2(betas, sigmas, tauHat_value) # term 3
                    
                    
                    var_hc3_term4_1_1 <- function(c_hc3, kVal){
                        numerator <- 0 
                        for(i in 1:length(kVal)){
                            numerator <- numerator + c_hc3[i]
                        }
                        return(numerator^2)
                    } # term 4_1_1
                    
                    var_hc3_term4_1_2 <- function(c_hc3, kVal){
                        numerator <- 0 
                        for(i in 1:length(kVal)){
                            numerator <- numerator + c_hc3[i]^2
                        }
                        return(numerator)
                    } # term 4_1_2
                    
                    var_hc3_term4_1 <- 2*(var_hc3_term4_1_1(c_hc3, kVal)-var_hc3_term4_1_2(c_hc3, kVal)) 
                    
                    var_hc3_term4_2 <- function(betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
                        }
                        return(numerator^2)
                    } # term 4_2
                    
                    var_hc3_term4<- var_hc3_term4_1*var_hc3_term4_2(betas, sigmas, tauHat_value) # term 4
                    
                    var_hc3_term5_1 <- function(c_hc3, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc3[i]*betas[i]^2*(sigmas[i]+tauHat_value)^2
                        }
                        return(4*numerator)
                    } # term 5_1
                    
                    var_hc3_term5_2 <- function(c_hc3){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc3[i]
                        }
                        return(numerator)
                    } # term 5_2
                    
                    var_hc3_term5 <- var_hc3_term5_1(c_hc3, betas, sigmas, tauHat_value)*
                        var_hc3_term5_2(c_hc3) # term 5
                    
                    var_hc3_term6 <- function(c_hc3, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc3[i]*betas[i]*(sigmas[i]+tauHat_value)
                        }
                        return(4*numerator^2)
                    } # term 6
                    
                    var_hc3_term7 <- function(c_hc3, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc3[i]^2*betas[i]^2*(sigmas[i]+tauHat_value)^2
                        }
                        return(8*numerator)
                    } # term 7
                    
                    var_hc3_term8_1 <- function(betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
                        }
                        return(8*numerator)
                    } # term 8_1
                    
                    var_hc3_term8_2_1 <- function(c_hc3, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc3[i]*betas[i]*(sigmas[i]+tauHat_value)
                        }
                        return(numerator)
                    } # term 8_2_1
                    
                    var_hc3_term8_2_2 <- function(c_hc3){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc3[i]
                        }
                        return(numerator)
                    } # term 8_2_2
                    
                    var_hc3_term8_2 <- var_hc3_term8_2_1(c_hc3, betas, sigmas, tauHat_value)*var_hc3_term8_2_2(c_hc3)
                    
                    var_hc3_term8_3 <- function(c_hc3, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc3[i]^2*betas[i]*(sigmas[i]+tauHat_value)
                        }
                        return(numerator)
                    } # term 8_3
                    
                    var_hc3_term8 <- var_hc3_term8_1(betas, sigmas, tauHat_value)*
                        (var_hc3_term8_2-var_hc3_term8_3(c_hc3, betas, sigmas, tauHat_value)) # term 8 
                    
                    var_hc3 <- var_hc3_term1(c_hc3, betas, sigmas, tauHat_value) + var_hc3_term2 + var_hc3_term3 + var_hc3_term4 +
                        var_hc3_term5 + var_hc3_term6(c_hc3, betas, sigmas, tauHat_value) - var_hc3_term7(c_hc3, betas, sigmas, tauHat_value)-
                        var_hc3_term8 # variance of HC3
                    
                    #-----------------------------------------------------------------------------------#
                    
                    ### Expected value of HC3
                    
                    mean_hc3_term1 <- function(c_hc3, betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + c_hc3[i]*(1-2*betas[i])*(sigmas[i]+tauHat_value)
                        }
                        return(numerator)
                    } # term 1
                    
                    mean_hc3_term2_1 <- function(c_hc3, kVal){
                        numerator <- 0 
                        for(i in 1:length(kVal)){
                            numerator <- numerator + c_hc3[i]
                        }
                        return(numerator)
                    } # term 2_1
                    
                    mean_hc3_term2_2 <- function(betas, sigmas, tauHat_value){
                        numerator <- 0
                        for (i in 1:length(betas)){
                            numerator <- numerator + betas[i]^2*(sigmas[i]+tauHat_value)
                        }
                        return(numerator)
                    } # term 2_2
                    
                    mean_hc3 <- mean_hc3_term1(c_hc3, betas, sigmas, tauHat_value)+
                        mean_hc3_term2_1(c_hc3, kVal)*mean_hc3_term2_2(betas, sigmas, tauHat_value) # Exp value of HC3
                    
                    df_hc3 <- 2*mean_hc3^2/var_hc3
                    
                    #-----------------------------------------------------------------------------------#
                    
                    LB_hc1 <- deltaHat_val - qt(.975, kVal - 1) * sqrt(hc1)
                    UB_hc1 <- deltaHat_val + qt(.975, kVal - 1) * sqrt(hc1)
                    is_in_CI_hc1 <- ifelse((delta >= LB_hc1 & delta <= UB_hc1), 1, 0)
                    
                    LB_hc2 <- deltaHat_val - qt(.975, kVal - 1) * sqrt(hc2)
                    UB_hc2 <- deltaHat_val + qt(.975, kVal - 1) * sqrt(hc2)
                    is_in_CI_hc2 <- ifelse((delta >= LB_hc2 & delta <= UB_hc2), 1, 0)
                    
                    LB_hc3 <- deltaHat_val - qt(.975, kVal - 1) * sqrt(hc3)
                    UB_hc3 <- deltaHat_val + qt(.975, kVal - 1) * sqrt(hc3)
                    is_in_CI_hc3 <- ifelse((delta >= LB_hc3 & delta <= UB_hc3), 1, 0)
                    
                    
                    LB_hc2_newDOF_est <- deltaHat_val - qt(.975, df_hc2) * sqrt(hc2)
                    UB_hc2_newDOF_est <- deltaHat_val + qt(.975, df_hc2) * sqrt(hc2)
                    is_in_CI_hc2_newDOF_est <- ifelse((delta >= LB_hc2_newDOF_est & delta <= UB_hc2_newDOF_est), 1, 0)
                    
                    LB_hc3_newDOF_est <- deltaHat_val - qt(.975, df_hc3) * sqrt(hc3)
                    UB_hc3_newDOF_est <- deltaHat_val + qt(.975, df_hc3) * sqrt(hc3)
                    is_in_CI_hc3_newDOF_est <- ifelse((delta >= LB_hc3_newDOF_est & delta <= UB_hc3_newDOF_est), 1, 0)
                    
                    LB_DL <- deltaHat_val-qnorm(.975)*sqrt(varDL_deltaHat)
                    UB_DL <- deltaHat_val+ qnorm(.975)*sqrt(varDL_deltaHat)
                    is_in_CI_DL <- ifelse((delta >= LB_DL & delta <= UB_DL), 1, 0)
                    
                    LB_HK <- deltaHat_val- qt(.975, kVal - 1)*sqrt(varHK_deltaHat) 
                    UB_HK <- deltaHat_val+ qt(.975, kVal - 1)*sqrt(varHK_deltaHat) 
                    is_in_CI_HK <- ifelse((delta >= LB_HK & delta <= UB_HK), 1, 0)
                    
                    LB_Hartung <- deltaHat_val- qt(.975, df)*sqrt(qofBeta) 
                    UB_Hartung <- deltaHat_val+ qt(.975, df)*sqrt(qofBeta) 
                    is_in_CI_Hartung <- ifelse((delta >= LB_Hartung & delta <= UB_Hartung), 1, 0)
                    
                    mainResults[currentRow,] <- c(kVal, currentTau, j, df, df_hc2, df_hc3, deltaHat_val, tauHat_value,
                                                  # variance estimators
                                                  varDL_deltaHat, varHK_deltaHat, qofBeta, hc1, hc2, hc3,
                                                  # confidence intervals
                                                  is_in_CI_DL, is_in_CI_HK, is_in_CI_Hartung,
                                                  is_in_CI_hc1, is_in_CI_hc2, is_in_CI_hc3,
                                                  is_in_CI_hc2_newDOF_est, is_in_CI_hc3_newDOF_est)
                    currentRow <- currentRow + 1
                }
                
            }
    return(mainResults)
}

df_sims <- data.frame()
for (l in 1:10000){
    currentResults <- main()
    df_sims <- rbind(df_sims, currentResults) 
}








