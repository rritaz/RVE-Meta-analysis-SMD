var_cohenD <- function(n) {
    output <- as.numeric((2 / n) + ((delta^2)/(4*as.numeric(n))))
    return (output)
} # variance of estimator of cohen's d

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
} # weighted mean with fixed effect weights

denomConst <- function(kVal, sigmas) {
    toReturn <- (sumW(kVal, 1, sigmas) - sumW(kVal, 2, sigmas)/sumW(kVal, 1, sigmas))
    return(as.numeric(toReturn))
} # constant function of weights in DL denominator of the heterogeneity variance

tauHat_DL <- function(ds, deltaHat_FE_val, kVal, sigmas){
    toReturn <- ((Q(ds, deltaHat_FE_val, kVal, sigmas)-(kVal-1))/denomConst(kVal, sigmas))
    if (toReturn<0){
        return(0)
    }else{
        return(as.numeric(toReturn))
    }
} # DerSimonian and Laird estimate of the heterogeneity variance

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
} # variance estimator HC1

hc2 <- function(ds, deltaHat_val, kVal, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(ds)){
        numerator <- numerator + weightTrue(tauHat_value, sigmas[i])^2*(1/(1-(weightTrue(tauHat_value, sigmas[i])/
                                                                                  sumTrueW(kVal, tauHat_value, sigmas))))*
            (ds[i]-deltaHat_val)^2
    }
    return(numerator/sumTrueW(kVal, tauHat_value, sigmas)^2)
} # variance estimator HC2

hc3 <- function(ds, deltaHat_val, kVal, sigmas, tauHat_value){
    numerator <- 0
    for (i in 1:length(ds)){
        numerator <- numerator + weightTrue(tauHat_value, sigmas[i])^2*(1/(1-(weightTrue(tauHat_value, sigmas[i])/
                                                                                  sumTrueW(kVal, tauHat_value, sigmas))))^2*
            (ds[i]-deltaHat_val)^2
    }
    return(numerator/sumTrueW(kVal, tauHat_value, sigmas)^2)
} # variance estimator HC3

c_hc2 <- function(kVal, betas){
    c <- vector()
    for (i in 1:kVal){
        c[i] <- (betas[i]^2)/(1-betas[i])
    }
    return(c)
} # c for HC2

c_hc3 <- function(kVal, betas){
    c <- vector()
    for (i in 1:kVal){
        c[i] <- (betas[i]^2)/((1-betas[i])^2)
    }
    return(c)
} # c for HC3


main <- function(){
    totalNumCol <- length(n_values[1,]) 

    mainResults <- matrix(NA, nrow=length(tau_sq)*totalNumCol, ncol=30)
    ##
    currentRow <- 1
    for(j in 1:totalNumCol) {
        for(k in 1:length(tau_sq)){
            currentTau <- tau_sq[k]*tau_sq_fun(j, kValue) 
            deltas <- generate_deltas_forCol(kValue, currentTau)
            ds <- generate_d_forCol(j, kValue, currentTau, deltas)
            sigmas <- generate_sigma_forCol(j, kValue, ds)
            
            deltaHat_FE_val <- delta_hatFE(kValue, ds, sigmas) 
            tauHat_value <- tauHat_DL(ds, deltaHat_FE_val, kValue, sigmas)
            
            deltaHat_val <- delta_hat(kValue, tauHat_value, ds, sigmas)
            
            betas <- betas(kValue, tauHat_value, sigmas)
            c_hc2 <- c_hc2(kValue, betas)
            c_hc3 <- c_hc3(kValue, betas)
            sumOfBetasSQ <- sumOfBetasSQ(betas, kValue)
            lambdaOfBeta <- sumOfBetasSQ(betas, kValue)/(1-sumOfBetasSQ(betas, kValue))
            ksi_i_Beta <- ksi_i_Beta(betas, kValue)
            
            QofBeta <- lambdaOfBeta*Q_RE(ds, deltaHat_val, kValue, betas) + QofBeta_term2(kValue, sigmas, ksi_i_Beta)
            RofBeta <- RofBeta(betas, sigmas)
            LofBeta <- min(1, max(0, ((QofBeta/RofBeta) - A)/(B-A)))
            qofBeta <- LofBeta*QofBeta + (1-LofBeta)*RofBeta # truncated Hartung variance estimator
            
            Var_QofBeta <- (lambdaOfBeta^2)*(1/sumTrueW(kValue, tauHat_value, sigmas))^2*(2*(kValue-1))
            Var_qofBeta <- LofBeta^2*Var_QofBeta # ignoring variance of within study variances
            df <- 2*((qofBeta^2)/Var_qofBeta)
            
            varDL_deltaHat <- 1/sumTrueW(kValue, tauHat_value, sigmas) # DL variance estimator
            varHK_deltaHat <- (1/((kValue-1)))*Q_RE(ds, deltaHat_val, kValue, betas) # HKSJ variance estimator
            
            hc1 <-hc1(ds, deltaHat_val, kValue, sigmas, tauHat_value)
            hc2 <-hc2(ds, deltaHat_val, kValue, sigmas, tauHat_value)
            hc3 <-hc3(ds, deltaHat_val, kValue, sigmas, tauHat_value)
            
            #----------------------------------------------------------------------------#
            
            # Compute degrees of freedom adjustment for HC2 
            
            var_hc2_term1 <- var_hc2_term1(c_hc2, betas, sigmas, tauHat_value)
            var_hc2_term2 <- var_hc2_term2_1(c_hc2)*var_hc2_term2_2(c_hc2, betas, sigmas, tauHat_value)
            var_hc2_term3 <- var_hc2_term3_1(c_hc2, betas, sigmas, tauHat_value)*
                var_hc2_term3_2(betas, sigmas, tauHat_value)
            var_hc2_term4<- (var_hc2_term4_1_1(c_hc2)-var_hc2_term4_1_2(c_hc2))*
                var_hc2_term4_2(betas, sigmas, tauHat_value)
            var_hc2_term5 <- var_hc2_term5_1(c_hc2, betas, sigmas, tauHat_value)*
                var_hc2_term5_2(c_hc2, sigmas)
            var_hc2_term6 <- var_hc2_term6(c_hc2, betas, sigmas, tauHat_value)
            var_hc2_term7 <- var_hc2_term7(c_hc2, betas, sigmas, tauHat_value)
            var_hc2_term8 <- var_hc2_term8_1(betas, sigmas, tauHat_value)*((
                var_hc2_term8_2_1(c_hc2, betas, sigmas, tauHat_value)*var_hc2_term8_2_2(c_hc2))-
                    var_hc2_term8_3(c_hc2, betas, sigmas, tauHat_value))
            
            # Compute variance of HC2
            var_hc2 <- var_hc2_term1 + var_hc2_term2 + var_hc2_term3 + var_hc2_term4 +
                var_hc2_term5 + var_hc2_term6 - var_hc2_term7 - var_hc2_term8
            
            # Compute expectetation of HC2
            mean_hc2 <- mean_hc2_term1(c_hc2, betas, sigmas, tauHat_value)+
                mean_hc2_term2_1(c_hc2, kValue)*mean_hc2_term2_2(betas, sigmas, tauHat_value) # Exp value of HC2
            
            # Compute degrees of freedom of HC2
            df_hc2_est <- 2*(mean_hc2^2)/var_hc2
            
            #-----------------------------------------------------------------------------------#
            
            # Compute degrees of freedom adjustment for HC3
            
            var_hc3_term1 <- var_hc3_term1(c_hc3, betas, sigmas, tauHat_value)
            var_hc3_term2 <- var_hc3_term2_1(c_hc3)*var_hc3_term2_2(c_hc3, betas, sigmas, tauHat_value)
            var_hc3_term3 <- var_hc3_term3_1(c_hc3, betas, sigmas, tauHat_value)*
                var_hc3_term3_2(betas, sigmas, tauHat_value)
            var_hc3_term4<- (var_hc3_term4_1_1(c_hc3)-var_hc3_term4_1_2(c_hc3))*
                var_hc3_term4_2(betas, sigmas, tauHat_value)
            var_hc3_term5 <- var_hc3_term5_1(c_hc3, betas, sigmas, tauHat_value)*
                var_hc3_term5_2(c_hc3, sigmas)
            var_hc3_term6 <- var_hc3_term6(c_hc3, betas, sigmas, tauHat_value)
            var_hc3_term7 <- var_hc3_term7(c_hc3, betas, sigmas, tauHat_value)
            var_hc3_term8 <- var_hc3_term8_1(betas, sigmas, tauHat_value)*((
                var_hc3_term8_2_1(c_hc3, betas, sigmas, tauHat_value)*var_hc3_term8_2_2(c_hc3))-
                    var_hc3_term8_3(c_hc3, betas, sigmas, tauHat_value))
            
            # Compute variance of HC3
            var_hc3 <- var_hc3_term1 + var_hc3_term2 + var_hc3_term3 + var_hc3_term4 +
                var_hc3_term5 + var_hc3_term6 - var_hc3_term7 - var_hc3_term8
            
            # Compute expectetation of HC3
            mean_hc3 <- mean_hc3_term1(c_hc3, betas, sigmas, tauHat_value)+
                mean_hc3_term2_1(c_hc3, kValue)*mean_hc3_term2_2(betas, sigmas, tauHat_value) # Exp value of HC3
            
            # Compute degrees of freedom of HC3
            df_hc3_est <- 2*(mean_hc3^2)/var_hc3
            
            #-----------------------------------------------------------------------------------#
            
            # Compute confidence intervals and CI length for all methods
            
            LB_DL <- deltaHat_val-qnorm(.975)*sqrt(varDL_deltaHat)
            UB_DL <- deltaHat_val+ qnorm(.975)*sqrt(varDL_deltaHat)
            is_in_CI_DL <- ifelse((delta >= LB_DL & delta <= UB_DL), 1, 0)
            DL_CI_length <- UB_DL-LB_DL
            
            LB_HK <- deltaHat_val- qt(.975, kValue - 1)*sqrt(varHK_deltaHat) 
            UB_HK <- deltaHat_val+ qt(.975, kValue - 1)*sqrt(varHK_deltaHat) 
            is_in_CI_HK <- ifelse((delta >= LB_HK & delta <= UB_HK), 1, 0)
            HK_CI_length <- UB_HK-LB_HK
            
            LB_Hartung <- deltaHat_val- qt(.975, df)*sqrt(qofBeta) 
            UB_Hartung <- deltaHat_val+ qt(.975, df)*sqrt(qofBeta) 
            is_in_CI_Hartung <- ifelse((delta >= LB_Hartung & delta <= UB_Hartung), 1, 0)
            Hartung_CI_length <- UB_Hartung-LB_Hartung
            
            LB_hc1 <- deltaHat_val - qt(.975, kValue - 1) * sqrt(hc1)
            UB_hc1 <- deltaHat_val + qt(.975, kValue - 1) * sqrt(hc1)
            is_in_CI_hc1 <- ifelse((delta >= LB_hc1 & delta <= UB_hc1), 1, 0)
            hc1_CI_length <- UB_hc1-LB_hc1
            
            LB_hc2 <- deltaHat_val - qt(.975, kValue - 1) * sqrt(hc2)
            UB_hc2 <- deltaHat_val + qt(.975, kValue - 1) * sqrt(hc2)
            is_in_CI_hc2 <- ifelse((delta >= LB_hc2 & delta <= UB_hc2), 1, 0)
            hc2_CI_length <- UB_hc2-LB_hc2
            
            LB_hc3 <- deltaHat_val - qt(.975, kValue - 1) * sqrt(hc3)
            UB_hc3 <- deltaHat_val + qt(.975, kValue - 1) * sqrt(hc3)
            is_in_CI_hc3 <- ifelse((delta >= LB_hc3 & delta <= UB_hc3), 1, 0)
            hc3_CI_length <- UB_hc3-LB_hc3
            
            LB_hc2_newDOF_est <- deltaHat_val - qt(.975, df_hc2_est) * sqrt(hc2)
            UB_hc2_newDOF_est <- deltaHat_val + qt(.975, df_hc2_est) * sqrt(hc2)
            is_in_CI_hc2_newDOF_est <- ifelse((delta >= LB_hc2_newDOF_est & delta <= UB_hc2_newDOF_est), 1, 0)
            hc2_newDOF_est_CI_length <- UB_hc2_newDOF_est-LB_hc2_newDOF_est
            
            LB_hc3_newDOF_est <- deltaHat_val - qt(.975, df_hc3_est) * sqrt(hc3)
            UB_hc3_newDOF_est <- deltaHat_val + qt(.975, df_hc3_est) * sqrt(hc3)
            is_in_CI_hc3_newDOF_est <- ifelse((delta >= LB_hc3_newDOF_est & delta <= UB_hc3_newDOF_est), 1, 0)
            hc3_newDOF_est_CI_length <- UB_hc3_newDOF_est-LB_hc3_newDOF_est
            
            mainResults[currentRow,] <- c(kValue, currentTau, j, 
                                          df, df_hc2_est, df_hc3_est, 
                                          deltaHat_val, tauHat_value,
                                          varDL_deltaHat, varHK_deltaHat, qofBeta, hc1, hc2, hc3,
                                          is_in_CI_DL, is_in_CI_HK, is_in_CI_Hartung,
                                          is_in_CI_hc1, is_in_CI_hc2, is_in_CI_hc3,
                                          is_in_CI_hc2_newDOF_est, is_in_CI_hc3_newDOF_est,
                                          DL_CI_length, HK_CI_length, Hartung_CI_length,
                                          hc1_CI_length, hc2_CI_length, hc3_CI_length,
                                          hc2_newDOF_est_CI_length, hc3_newDOF_est_CI_length)
            currentRow <- currentRow + 1
        }
        
    }
    return(mainResults)
}
