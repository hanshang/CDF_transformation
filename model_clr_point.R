#######
# clr
#######

select_K <- function(tau, eigenvalue)
{
    k_max = length(eigenvalue)
    k_all = rep(0, k_max-1)
    for(k in 1:(k_max-1))
    {
        k_all[k] = (eigenvalue[k+1]/eigenvalue[k])*ifelse(eigenvalue[k]/eigenvalue[1] > tau, 1, 0) + ifelse(eigenvalue[k]/eigenvalue[1] < tau, 1, 0)
    }
    K_hat = which.min(k_all)
    return(K_hat)
}

# fdata: n by p data matrix
# ncomp_selection: method for selecting the number of retained components
# fh: forecast horizon
# fore_method: forecasting method
# object_interest: point or interval forecasts
# B = 399: number of bootstrap replications
# alpha: level of significance

clr_fun <- function(fdata, ncomp_selection, fh, fore_method, object_interest, B = 399, alpha)
{
    n_age = ncol(fdata)
    n_year = nrow(fdata)
    h_x_t = CLR(fdata)$LR

    SVD_decomp = svd(h_x_t)
    if(ncomp_selection == "EVR")
    {
        ncomp = select_K(tau = 0.001, eigenvalue = SVD_decomp$d^2)
    }
    else if(ncomp_selection == "fixed")
    {
        ncomp = 6
    }
    else
    {
        warning("The number of retained component must be chosen by EVR or fixed at 6.")
    }
    basis = SVD_decomp$v[,1:ncomp]
    score = t(basis) %*% t(h_x_t)
    recon = basis %*% score
    resi = t(h_x_t) - recon

    # reconstruction (model in-sample fitting)

    recon = invCLR(t(recon))

    if(object_interest == "point")
    {
        # forecasts of principal component scores

        score_fore = matrix(NA, ncomp, 1)
        for(ik in 1:ncomp)
        {
            if(fore_method == "RWF_no_drift")
            {
                score_fore[ik,] = rwf(as.numeric(score[ik,]), h = fh, drift = FALSE)$mean[fh]
            }
            else if(fore_method == "RWF_drift")
            {
                score_fore[ik,] = rwf(as.numeric(score[ik,]), h = fh, drift = TRUE)$mean[fh]
            }
            else if(fore_method == "ETS")
            {
                score_fore[ik,] = forecast(ets(as.numeric(score[ik,])), h = fh)$mean[fh]
            }
            else if(fore_method == "ARIMA")
            {
                score_fore[ik,] = forecast(auto.arima(as.numeric(score[ik,])), h = fh)$mean[fh]
            }
            else
            {
                warning("Univariate time series forecasting method is not on the list.")
            }
        }

        # obtain forecasts in real-valued space

        fore_val = basis %*% score_fore
        fore_count = invCLR(t(fore_val))
        return(list(ncomp = ncomp, fore_count = fore_count))
    }
    else if(object_interest == "interval")
    {
        # determine in-sample forecast error for principal component scores

        olivia = matrix(NA, ncomp, fh)
        if(fore_method == "ets")
        {
            for(ij in 1:ncomp)
            {
                olivia[ij,] = forecast(ets(score[ij,]), h = fh)$mean
            }
        }
        else if(fore_method == "arima")
        {
            for(ij in 1:ncomp)
            {
                olivia[ij,] = forecast(auto.arima(score[ij,]), h = fh)$mean
            }
        }
        else if(fore_method == "rwf")
        {
            for(ij in 1:ncomp)
            {
                olivia[ij,] = rwf(score[ij,], h = fh, drift = TRUE)$mean
            }
        }
        else if(fore_method == "rw")
        {
            for(ij in 1:ncomp)
            {
                olivia[ij,] = rwf(score[ij,], h = fh, drift = FALSE)$mean
            }
        }
        else
        {
            warning("Forecasting method must be ets, arima, rwf or rw.")
        }
        forerr = matrix(NA, (n_year - ncomp - fh + 1), ncomp)
        for(i in fh:(n_year - ncomp))
        {
            k = i + (ncomp - fh)
            fore = matrix(NA, 1, ncomp)
            if(fore_method == "ets")
            {
                for(j in 1:ncomp)
                {
                    fore[,j] = forecast(ets(score[j,1:k]), h = fh)$mean[fh]
                }
            }
            else if(fore_method == "arima")
            {
                for(j in 1:ncomp)
                {
                    fore[,j] = forecast(auto.arima(score[j,1:k]), h = fh)$mean[fh]
                }
            }
            else if(fore_method == "rwf")
            {
                if(k <= 2)
                {
                    for(j in 1:ncomp)
                    {
                        fore[,j] = score[j,k]
                    }
                }
                if(k > 2)
                {
                    for(j in 1:ncomp)
                    {
                        fore[,j] = rwf(score[j,1:k], h = fh, drift = TRUE)$mean[fh]
                    }
                }
            }
            else if(fore_method == "rw")
            {
                if(k == 1)
                {
                    for(j in 1:ncomp)
                    {
                        fore[,j] = score[j,1]
                    }
                }
                if(k > 1)
                {
                    for(j in 1:ncomp)
                    {
                        fore[,j] = rwf(score[j,1:k], h = fh, drift = FALSE)$mean[fh]
                    }
                }
            }
            forerr[i - fh + 1,] = score[, k + fh] - fore
        }
        # bootstrapping residuals
        K = 1

        q = array(NA, dim = c(n_age, B, K, fh))
        for(j in 1:fh)
        {
            for(i in 1:n_age)
            {
                for(k in 1:K)
                {
                    q[i,,k,j] = sample(resi[i,], size = B, replace = TRUE)
                }
            }
        }
        rm(i); rm(j); rm(k)
        # bootstrapping PC score errors
        ny = array(NA, dim = c(ncomp, B, fh))
        for(j in 1:fh)
        {
            for(i in 1:ncomp)
            {
                ny[i,,j] = sample(forerr[,i], size = B, replace = TRUE)
            }
        }
        rm(i); rm(j)
        # adding the PC score error to the predicted score
        oli = array(rep(olivia, B * fh), dim = c(ncomp, B, fh))
        fo = array(NA, dim = c(ncomp, B, fh))
        for(j in 1:fh)
        {
            for(i in 1:B)
            {
                fo[,i,j] = oli[,i,j] + ny[,i,j]
            }
        }
        rm(i); rm(j)
        # construct bootstrapped samples
        pred = array(NA, dim = c(n_age, B, K, fh))
        for(j in 1:fh)
        {
            for(i in 1:B)
            {
                for(k in 1:K)
                {
                    pred[,i,k,j] = basis %*% fo[,i,j] + q[,i,k,j]
                }
            }
        }
        rm(i); rm(j); rm(k)

        pred_resize = array(NA, dim = c(n_age, B * K, fh))
        for(j in 1:fh)
        {
            for(i in 1:B)
            {
                pred_resize[, (((i-1)*K+1):(i*K)), ] = pred[,i,,j]
            }
        }
        rm(i); rm(j)

        # transform back

        d_x_t_star_fore = array(NA, dim = c(n_age, B * K, fh))
        for(iw in 1:fh)
        {
            for(ij in 1:(B * K))
            {
                d_x_t_star_fore[,ij,iw] = invCLR(t(pred_resize[,ij,iw]))
            }
        }
        rm(iw); rm(ij)
        return(list(PI = apply(d_x_t_star_fore, c(1, 3), quantile, c(alpha/2, 1-alpha/2)), ncomp = ncomp))
    }
    else
    {
        warning("object_interest must either be point or interval.")
    }
}

#############################################################
# point forecast error for h = 1,...,16 using the clr method
#############################################################

point_fore_national_clr_EVR_err_female = point_fore_national_clr_EVR_err_male =
point_fore_national_clr_EVR_err_female_lt = point_fore_national_clr_EVR_err_male_lt =
point_fore_national_clr_fixed_err_female = point_fore_national_clr_fixed_err_male =
point_fore_national_clr_fixed_err_female_lt = point_fore_national_clr_fixed_err_male_lt = matrix(NA, 16, 2)
for(iw in 1:16)
{
    # ncomp = "EVR"

    dum = point_fore_national_cdf(fdata = female_JPN_pop, sex = "female", horizon = iw,
                                  fore_method = "CLR", variable_interest = "point",
                                  CLR_ncomp_selection = "EVR")
    point_fore_national_clr_EVR_err_female[iw,] = dum$err
    point_fore_national_clr_EVR_err_female_lt[iw,] = dum$err_lt
    rm(dum)

    dum = point_fore_national_cdf(fdata = male_JPN_pop, sex = "male", horizon = iw,
                                  fore_method = "CLR", variable_interest = "point",
                                  CLR_ncomp_selection = "EVR")
    point_fore_national_clr_EVR_err_male[iw,] = dum$err
    point_fore_national_clr_EVR_err_male_lt[iw,] = dum$err_lt
    rm(dum)

    # ncomp = 6

    dum = point_fore_national_cdf(fdata = female_JPN_pop, sex = "female", horizon = iw,
                                  fore_method = "CLR", variable_interest = "point",
                                  CLR_ncomp_selection = "fixed")
    point_fore_national_clr_fixed_err_female[iw,] = dum$err
    point_fore_national_clr_fixed_err_female_lt[iw,] = dum$err_lt
    rm(dum)

    dum = point_fore_national_cdf(fdata = male_JPN_pop, sex = "male", horizon = iw,
                                  fore_method = "CLR", variable_interest = "point",
                                  CLR_ncomp_selection = "fixed")
    point_fore_national_clr_fixed_err_male[iw,] = dum$err
    point_fore_national_clr_fixed_err_male_lt[iw,] = dum$err_lt
    print(iw); rm(iw); rm(dum)
}
colnames(point_fore_national_clr_EVR_err_female) = colnames(point_fore_national_clr_EVR_err_male) =
colnames(point_fore_national_clr_fixed_err_female) = colnames(point_fore_national_clr_fixed_err_male) = c("KLD", "JSD")

colnames(point_fore_national_clr_EVR_err_female_lt) = colnames(point_fore_national_clr_EVR_err_male_lt) =
colnames(point_fore_national_clr_fixed_err_female_lt) = colnames(point_fore_national_clr_fixed_err_male_lt) = c("RMSE", "MAE")

rownames(point_fore_national_clr_EVR_err_female) = rownames(point_fore_national_clr_EVR_err_male) =
rownames(point_fore_national_clr_fixed_err_female) = rownames(point_fore_national_clr_fixed_err_male) =
rownames(point_fore_national_clr_EVR_err_female_lt) = rownames(point_fore_national_clr_EVR_err_male_lt) =
rownames(point_fore_national_clr_fixed_err_female_lt) = rownames(point_fore_national_clr_fixed_err_male_lt) = 1:16

round(colMeans(point_fore_national_clr_EVR_err_female), 4) # 0.0128 0.0034
round(colMeans(point_fore_national_clr_EVR_err_male), 4)   # 0.0060 0.0015

round(colMeans(point_fore_national_clr_fixed_err_female), 4) # 0.0109 0.0029
round(colMeans(point_fore_national_clr_fixed_err_male), 4)   # 0.0043 0.0011
