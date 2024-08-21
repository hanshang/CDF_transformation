interval_score <- function(holdout, lb, ub, alpha)
{
    lb_ind = ifelse(holdout < lb, 1, 0)
    ub_ind = ifelse(holdout > ub, 1, 0)
    score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
    cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout)
    cpd = abs(cover - (1 - alpha))
    return(c(mean(score), cpd, cover))
}

############################################
# Multivariate functional time series model
############################################

# data_set: a list of p by n data matrix
# ncomp_method: method for selecting the number of components
# fh: forecast horizon
# fore_method: univariate time-series forecasting method, "ARIMA" or "ETS"
# object_interest: point or interval forecasts
# boot_number: number of bootstrap samples
# PI_level: prediction interval level customarily 80% or 95%

MFTS_model <- function(data_input, ncomp_method, fh, fore_method, object_interest, boot_number, PI_level)
{
    n_age  = dim(data_input)[1]
    n_year = dim(data_input)[2]
    n_pop  = dim(data_input)[3]

    # clean the data

    data_set_array = array(NA, dim = c(n_age, n_year, n_pop))
    if(any(!is.finite(data_input)))
    {
        for(iw in 1:n_pop)
        {
            for(ij in 1:n_age)
            {
                data_set_array[ij,,iw] = na.interp(data_input[ij,,iw])
            }
        }
    }
    else
    {
        data_set_array = data_input
    }

    #comb_object_raw = list()
    #for(ik in 1:n_pop)
    #{
    #    comb_object_raw[[ik]] = data_set_array[,,ik]
    #}

    rowmeans_object = sd_object = decenter_object = list()
    for(ik in 1:n_pop)
    {
        # compute mean and sd function
        rowmeans_object[[ik]] = rowMeans(data_set_array[,,ik], na.rm = TRUE)
        sd_object[[ik]] = apply(data_set_array[,,ik], 1, sd, na.rm = TRUE)

        # de-center functional data
        decenter_object[[ik]] = t(scale(t(data_set_array[,,ik]), center = TRUE, scale = TRUE))
    }

    comb_object = do.call(rbind, decenter_object)
    # comb_object = do.call(rbind, comb_object_raw)
    colnames(comb_object) = 1:ncol(comb_object)

    eigen_value = eigen(cov(t(comb_object)))$values
    if(ncomp_method == "EVR")
    {
        ncomp = select_K(tau = 10^-3, eigenvalue = eigen_value)
    }
    else if(ncomp_method == "provide")
    {
        ncomp = 6
    }
    else
    {
        warning("The number of components is required.")
    }
    fore_ftsm = forecast(ftsm(fts(1:nrow(comb_object), comb_object), order = ncomp), h = fh, method = fore_method,
                         level = PI_level, pimethod = "nonparametric", B = 399)
    if(object_interest == "point")
    {
        res_fore = as.matrix(fore_ftsm$mean$y[,fh] * do.call(c, sd_object) + do.call(c, rowmeans_object))
        #res_fore = as.matrix(fore_ftsm$mean$y[,fh])
    }
    else if(object_interest == "interval")
    {
        res_fore = matrix(NA, n_age * 2, 399)
        for(iw in 1:399)
        {
            res_fore[,iw] = as.matrix(fore_ftsm$bootsamp[,iw,fh] * do.call(c, sd_object) + do.call(c, rowmeans_object))
            rm(iw)
        }
    }
    else
    {
        warning("forecasts must either be point or interval.")
    }
    return(res_fore)
}

# data_F: female data
# data_M: male data
# fh: forecast horizon
# fmethod: forecasting method
# object_interest: point or interval forecasts
# no_boot: number of bootstrap samples
# alpha: level of significance
# method_ncomp: K = 6 or fixed

fore_national_cdf_MFTS <- function(data_F, data_M, fh, fmethod, object_interest, no_boot, alpha, method_ncomp)
{
    data_cumsum_dum_F = data_cumsum_dum_M = matrix(NA, nrow(data_F), ncol(data_F))
    for(ij in 1:nrow(data_F))
    {
        data_cumsum_dum_F[ij,] = cumsum(data_F[ij,])
        data_cumsum_dum_M[ij,] = cumsum(data_M[ij,])
        rm(ij)
    }

    if(any(data_cumsum_dum_F == 0))
    {
        data_cumsum_F = replace(data_cumsum_dum_F, which(data_cumsum_dum_F == 0), 10^-5)
    }
    else
    {
        data_cumsum_F = data_cumsum_dum_F
    }
    if(any(data_cumsum_dum_M == 0))
    {
        data_cumsum_M = replace(data_cumsum_dum_M, which(data_cumsum_dum_M == 0), 10^-5)
    }
    else
    {
        data_cumsum_M = data_cumsum_dum_M
    }
    rm(data_cumsum_dum_F); rm(data_cumsum_dum_M)

    data_cumsum_logit_F = data_cumsum_logit_M = matrix(NA, nrow(data_F), (ncol(data_F) - 1))
    for(ij in 1:nrow(data_F))
    {
        data_cumsum_logit_F[ij,] = logit(data_cumsum_F[ij, 1:(ncol(data_F) - 1)])
        data_cumsum_logit_M[ij,] = logit(data_cumsum_M[ij, 1:(ncol(data_M) - 1)])
        rm(ij)
    }
    rownames(data_cumsum_logit_F) = rownames(data_cumsum_logit_M) = years[1:nrow(data_F)]

    data_comb = array(NA, dim = c((ncol(data_F) - 1), nrow(data_F), 2))
    data_comb[,,1] = t(data_cumsum_logit_F)
    data_comb[,,2] = t(data_cumsum_logit_M)

    dum = MFTS_model(data_input = data_comb, ncomp_method = method_ncomp, fh = fh, fore_method = fmethod,
                     object_interest = object_interest, boot_number = no_boot, PI_level = (1 - alpha) * 100)

    if(object_interest == "point")
    {
        data_cumsum_logit_F_fore = dum[1:(ncol(data_F) - 1),]
        data_cumsum_logit_M_fore = dum[ncol(data_F):(2 * (ncol(data_F) - 1)),]
        rm(dum)

        data_cumsum_logit_F_fore_add = c(invlogit(data_cumsum_logit_F_fore), 1)
        data_cumsum_logit_M_fore_add = c(invlogit(data_cumsum_logit_M_fore), 1)

        data_cumsum_logit_F_fore_add_diff = c(data_cumsum_logit_F_fore_add[1], diff(data_cumsum_logit_F_fore_add))
        data_cumsum_logit_M_fore_add_diff = c(data_cumsum_logit_M_fore_add[1], diff(data_cumsum_logit_M_fore_add))
        return(list(mfts_fore_F = data_cumsum_logit_F_fore_add_diff, mfts_fore_M = data_cumsum_logit_M_fore_add_diff))
    }
    else if(object_interest == "interval")
    {
        data_cumsum_logit_F_fore = dum[1:(ncol(data_F) - 1),]
        data_cumsum_logit_M_fore = dum[ncol(data_F):(2 * (ncol(data_F) - 1)),]
        rm(dum)

        data_cumsum_logit_F_fore_add = rbind(invlogit(data_cumsum_logit_F_fore), rep(1, 399))
        data_cumsum_logit_M_fore_add = rbind(invlogit(data_cumsum_logit_M_fore), rep(1, 399))

        data_cumsum_logit_F_fore_add_diff = data_cumsum_logit_M_fore_add_diff = matrix(NA, ncol(data_F), 399)
        for(iw in 1:399)
        {
            data_cumsum_logit_F_fore_add_diff[,iw] = c(data_cumsum_logit_F_fore_add[1,iw], diff(data_cumsum_logit_F_fore_add[,iw]))
            data_cumsum_logit_M_fore_add_diff[,iw] = c(data_cumsum_logit_M_fore_add[1,iw], diff(data_cumsum_logit_M_fore_add[,iw]))
            rm(iw)
        }
        rownames(data_cumsum_logit_F_fore_add_diff) = rownames(data_cumsum_logit_M_fore_add_diff) = ages
        colnames(data_cumsum_logit_F_fore_add_diff) = colnames(data_cumsum_logit_M_fore_add_diff) = 1:399
        return(list(mfts_fore_F_PI = t(apply(data_cumsum_logit_F_fore_add_diff, 1, quantile, c(alpha/2, (1 - alpha/2)))),
                    mfts_fore_M_PI = t(apply(data_cumsum_logit_M_fore_add_diff, 1, quantile, c(alpha/2, (1 - alpha/2))))))
    }
    else
    {
        warning("Forecasts must either be point or interval.")
    }
}

# fdata_F: a data matrix of dimension (n by p)
# fdata_M: a data matrix of dimension (n by p)
# horizon: forecast horizon 1 to 16
# way_ncomp: K = 6 or EVR

point_fore_national_cdf_MFTS <- function(fdata_F, fdata_M, horizon, way_ncomp)
{
    forecast_val_F = forecast_val_M = matrix(NA, ncol(fdata_F), (17 - horizon))
    for(ij in 1:(17 - horizon))
    {
        dum <- fore_national_cdf_MFTS(data_F = fdata_F[1:(31+ij),], data_M = fdata_M[1:(31+ij),],
                                      fh = horizon, fmethod = "ets", object_interest = "point",
                                      no_boot = 399, alpha = 0.2, method_ncomp = way_ncomp)
        forecast_val_F[,ij] = dum$mfts_fore_F
        forecast_val_M[,ij] = dum$mfts_fore_M
        rm(ij)
    }

    lt_mat_F = lt_mat_M = matrix(NA, (17 - horizon), ncol(fdata))
    for(ik in 1:(17 - horizon))
    {
        lt_mat_F[ik,] = LifeTable(0:110, dx = forecast_val_F[,ik] * 10^5)$lt$ex
        lt_mat_M[ik,] = LifeTable(0:110, dx = forecast_val_M[,ik] * 10^5)$lt$ex
        rm(ik)
    }
    err_lt_F = c(ftsa:::rmse(forecast = lt_mat_F, true = female_JPN_ex[(32 + horizon):48,]),
                 ftsa:::mae(forecast = lt_mat_F,  true = female_JPN_ex[(32 + horizon):48,]))

    err_lt_M = c(ftsa:::rmse(forecast = lt_mat_M, true = male_JPN_ex[(32 + horizon):48,]),
                 ftsa:::mae(forecast = lt_mat_M,  true = male_JPN_ex[(32 + horizon):48,]))


    holdout_val_dum_F = t(matrix(fdata_F[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
    holdout_val_dum_M = t(matrix(fdata_M[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
    if(any(holdout_val_dum_F == 0))
    {
        holdout_val_F = replace(holdout_val_dum_F, which(holdout_val_dum_F == 0), 10^-5)
    }
    else
    {
        holdout_val_F = holdout_val_dum_F
    }
    if(any(holdout_val_dum_M == 0))
    {
        holdout_val_M = replace(holdout_val_dum_M, which(holdout_val_dum_M == 0), 10^-5)
    }
    else
    {
        holdout_val_M = holdout_val_dum_M
    }
    rm(holdout_val_dum_F); rm(holdout_val_dum_M)

    KL_div_val_F = JS_div_val_F = KL_div_val_M = JS_div_val_M = vector("numeric", (17 - horizon))
    for(ij in 1:(17 - horizon))
    {
        KL_div_val_F[ij] = mean(KLdiv(cbind(forecast_val_F[,ij], holdout_val_F[,ij]))[2:3])
        JS_div_val_F[ij] = mean(KLdiv(cbind(forecast_val_F[,ij], apply(cbind(forecast_val_F[,ij], holdout_val_F[,ij]), 1, geometric.mean)))[2:3])

        KL_div_val_M[ij] = mean(KLdiv(cbind(forecast_val_M[,ij], holdout_val_M[,ij]))[2:3])
        JS_div_val_M[ij] = mean(KLdiv(cbind(forecast_val_M[,ij], apply(cbind(forecast_val_M[,ij], holdout_val_M[,ij]), 1, geometric.mean)))[2:3])
    }

    err_F = c(mean(KL_div_val_F), mean(JS_div_val_F))
    err_M = c(mean(KL_div_val_M), mean(JS_div_val_M))

    return(list(forecast_pdf_F = forecast_val_F, forecast_pdf_M = forecast_val_M,
                holdout_pdf_F = holdout_val_F, holdout_pdf_M = holdout_val_M,
                err_F = err_F, err_M = err_M, err_lt_F = err_lt_F,
                err_lt_M = err_lt_M))
}

#################################################
## point forecast error for h = 1,...,16 (K = 6)
#################################################

point_fore_national_err_MFTS_F = point_fore_national_err_MFTS_M =
point_fore_national_err_MFTS_F_lt = point_fore_national_err_MFTS_M_lt = matrix(NA, 16, 2)
for(iw in 1:16)
{
    dum = point_fore_national_cdf_MFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop, horizon = iw,
                                       way_ncomp = "provide")
    point_fore_national_err_MFTS_F[iw,] = dum$err_F
    point_fore_national_err_MFTS_M[iw,] = dum$err_M
    point_fore_national_err_MFTS_F_lt[iw,] = dum$err_lt_F
    point_fore_national_err_MFTS_M_lt[iw,] = dum$err_lt_M
    print(iw); rm(iw); rm(dum)
}
colnames(point_fore_national_err_MFTS_F) = colnames(point_fore_national_err_MFTS_M) = c("KLD", "JSD (geo)")
colnames(point_fore_national_err_MFTS_F_lt) = colnames(point_fore_national_err_MFTS_M_lt) = c("RMSE", "MAE")
rownames(point_fore_national_err_MFTS_F) = rownames(point_fore_national_err_MFTS_M) =
rownames(point_fore_national_err_MFTS_F_lt) = rownames(point_fore_national_err_MFTS_M_lt) = 1:16

######
# EVR
######

point_fore_national_err_MFTS_F_EVR = point_fore_national_err_MFTS_M_EVR =
point_fore_national_err_MFTS_F_EVR_lt = point_fore_national_err_MFTS_M_EVR_lt = matrix(NA, 16, 2)
for(iw in 1:16)
{
    dum = point_fore_national_cdf_MFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop, horizon = iw,
                                       way_ncomp = "EVR")
    point_fore_national_err_MFTS_F_EVR[iw,] = dum$err_F
    point_fore_national_err_MFTS_M_EVR[iw,] = dum$err_M
    point_fore_national_err_MFTS_F_EVR_lt[iw,] = dum$err_lt_F
    point_fore_national_err_MFTS_M_EVR_lt[iw,] = dum$err_lt_M
    print(iw); rm(iw); rm(dum)
}
colnames(point_fore_national_err_MFTS_F_EVR) = colnames(point_fore_national_err_MFTS_M_EVR) = c("KLD", "JSD (geo)")
colnames(point_fore_national_err_MFTS_F_EVR_lt) = colnames(point_fore_national_err_MFTS_M_EVR_lt) = c("RMSE", "MAE")

rownames(point_fore_national_err_MFTS_F_EVR) = rownames(point_fore_national_err_MFTS_M_EVR) =
rownames(point_fore_national_err_MFTS_F_EVR_lt) = rownames(point_fore_national_err_MFTS_M_EVR_lt) = 1:16

#####################
# interval forecasts
#####################

interval_fore_national_cdf_MFTS <- function(fdata_F, fdata_M, horizon, alpha_level, way_ncomp)
{
    forecast_val_PI_F = forecast_val_PI_M = array(NA, dim = c(ncol(fdata_F), 2, (17 - horizon)))
    for(ij in 1:(17 - horizon))
    {
        dum <- fore_national_cdf_MFTS(data_F = fdata_F[1:(31+ij),], data_M = fdata_M[1:(31+ij),],
                                      fh = horizon, fmethod = "ets", object_interest = "interval",
                                      no_boot = 1000, alpha = alpha_level, method_ncomp = way_ncomp)
        forecast_val_PI_F[,,ij] = dum$mfts_fore_F_PI
        forecast_val_PI_M[,,ij] = dum$mfts_fore_M_PI
        rm(ij)
    }

    holdout_val_dum_F = t(matrix(fdata_F[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
    holdout_val_dum_M = t(matrix(fdata_M[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
    if(any(holdout_val_dum_F == 0))
    {
        holdout_val_F = replace(holdout_val_dum_F, which(holdout_val_dum_F == 0), 10^-5)
    }
    else
    {
        holdout_val_F = holdout_val_dum_F
    }
    if(any(holdout_val_dum_M == 0))
    {
        holdout_val_M = replace(holdout_val_dum_M, which(holdout_val_dum_M == 0), 10^-5)
    }
    else
    {
        holdout_val_M = holdout_val_dum_M
    }
    rm(holdout_val_dum_F); rm(holdout_val_dum_M)

    return(list(int_err_F = interval_score(holdout = holdout_val_F, lb = forecast_val_PI_F[,1,],
                                           ub = forecast_val_PI_F[,2,], alpha = alpha_level),
                int_err_M = interval_score(holdout = holdout_val_M, lb = forecast_val_PI_M[,1,],
                                           ub = forecast_val_PI_M[,2,], alpha = alpha_level)))
}

#########
## K = 6
#########

# interval forecast errors for h = 1,...,16 (alpha = 0.2)

interval_fore_national_err_MFTS_F = interval_fore_national_err_MFTS_M = matrix(NA, 16, 3)
for(iw in 1:16)
{
    dum = interval_fore_national_cdf_MFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop, horizon = iw,
                                          alpha_level = 0.2, way_ncomp = "provide")
    interval_fore_national_err_MFTS_F[iw,] = dum$int_err_F
    interval_fore_national_err_MFTS_M[iw,] = dum$int_err_M
    print(iw); rm(iw); rm(dum)
}
colnames(interval_fore_national_err_MFTS_F) = colnames(interval_fore_national_err_MFTS_M) = c("score", "CPD", "ECP")
rownames(interval_fore_national_err_MFTS_F) = rownames(interval_fore_national_err_MFTS_M) = 1:16

# interval forecast errors for h = 1,...,16 (alpha = 0.05)

interval_fore_national_err_MFTS_F_PI_95 = interval_fore_national_err_MFTS_M_PI_95 = matrix(NA, 16, 3)
for(iw in 1:16)
{
    dum = interval_fore_national_cdf_MFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop, horizon = iw,
                                          alpha_level = 0.05, way_ncomp = "provide")
    interval_fore_national_err_MFTS_F_PI_95[iw,] = dum$int_err_F
    interval_fore_national_err_MFTS_M_PI_95[iw,] = dum$int_err_M
    print(iw); rm(iw); rm(dum)
}

colnames(interval_fore_national_err_MFTS_F_PI_95) = colnames(interval_fore_national_err_MFTS_M_PI_95) = c("score", "CPD", "ECP")
rownames(interval_fore_national_err_MFTS_F_PI_95) = rownames(interval_fore_national_err_MFTS_F_PI_95) = 1:16

#######
## EVR
#######

# interval forecast errors for h = 1,...,16 (alpha = 0.2)

interval_fore_national_err_MFTS_F_EVR = interval_fore_national_err_MFTS_M_EVR = matrix(NA, 16, 3)
for(iw in 1:16)
{
    dum = interval_fore_national_cdf_MFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop, horizon = iw,
                                          alpha_level = 0.2, way_ncomp = "provide")
    interval_fore_national_err_MFTS_F_EVR[iw,] = dum$int_err_F
    interval_fore_national_err_MFTS_M_EVR[iw,] = dum$int_err_M
    print(iw); rm(iw); rm(dum)
}
colnames(interval_fore_national_err_MFTS_F_EVR) = colnames(interval_fore_national_err_MFTS_M_EVR) = c("score", "CPD", "ECP")
rownames(interval_fore_national_err_MFTS_F_EVR) = rownames(interval_fore_national_err_MFTS_M_EVR) = 1:16

# interval forecast errors for h = 1,...,16 (alpha = 0.05)

interval_fore_national_err_MFTS_F_PI_95_EVR = interval_fore_national_err_MFTS_M_PI_95_EVR = matrix(NA, 16, 3)
for(iw in 1:16)
{
    dum = interval_fore_national_cdf_MFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop, horizon = iw,
                                          alpha_level = 0.05, way_ncomp = "provide")
    interval_fore_national_err_MFTS_F_PI_95_EVR[iw,] = dum$int_err_F
    interval_fore_national_err_MFTS_M_PI_95_EVR[iw,] = dum$int_err_M
    print(iw); rm(iw); rm(dum)
}
colnames(interval_fore_national_err_MFTS_F_PI_95_EVR) = colnames(interval_fore_national_err_MFTS_M_PI_95_EVR) = c("score", "CPD", "ECP")
rownames(interval_fore_national_err_MFTS_F_PI_95_EVR) = rownames(interval_fore_national_err_MFTS_M_PI_95_EVR) = 1:16
