# data_raw: a data matrix of dimension (n by p) where n denotes the number of years and p denotes the number of ages
# ncomp_method: method for selecting the number of components
# fh: forecast horizon
# fmethod: forecasting method
# object_interest: point or interval forecasts
# no_boot: number of bootstrap samples
# alpha: level of significance between 0 and 1

fore_national_cdf <- function(data, ncomp_method, fh, fmethod, object_interest, no_boot = 1000, alpha = 0.2)
{
    data_cumsum_dum = matrix(NA, nrow(data), ncol(data))
    for(ij in 1:nrow(data))
    {
        data_cumsum_dum[ij,] = cumsum(data[ij,])
        rm(ij)
    }

    # check if any cumsum values equal to 0

    if(any(data_cumsum_dum == 0))
    {
        data_cumsum = replace(data_cumsum_dum, which(data_cumsum_dum == 0), 10^-5)
    }
    else
    {
        data_cumsum = data_cumsum_dum
    }
    rm(data_cumsum_dum)

    # logit transformation

    data_cumsum_logit = matrix(NA, nrow(data), (ncol(data) - 1))
    for(ij in 1:nrow(data))
    {
        data_cumsum_logit[ij,] = logit(data_cumsum[ij, 1:(ncol(data) - 1)])
        rm(ij)
    }
    rm(data_cumsum)
    rownames(data_cumsum_logit) = years[1:nrow(data)]

    # fitting a functional time series forecasting method

    if(ncomp_method == "EVR")
    {
        ncomp = select_K(tau = 10^-3, eigenvalue = (svd(data_cumsum_logit)$d)^2)
    }
    else if(ncomp_method == "provide")
    {
        ncomp = 6
    }
    else
    {
        warning("The number of components is required.")
    }
    data_cumsum_logit_fore = forecast(ftsm(fts(ages[1:110], t(data_cumsum_logit)), order = ncomp), h = fh,
                                      method = fmethod, pimethod = "nonparametric",
                                      level = (1 - alpha) * 100, B = no_boot)

    if(object_interest == "point")
    {
        # h-step-ahead forecast

        data_cumsum_logit_fore_add = c(invlogit(data_cumsum_logit_fore$mean$y[,fh]), 1)
        data_cumsum_logit_fore_add_diff = c(data_cumsum_logit_fore_add[1], diff(data_cumsum_logit_fore_add))
        return(data_cumsum_logit_fore_add_diff)
    }
    else if(object_interest == "interval")
    {
        # bootstrap forecasts

        data_cumsum_logit_fore_add_boot_diff = matrix(NA, ncol(data), no_boot)
        for(ik in 1:no_boot)
        {
            data_cumsum_logit_fore_add_boot = c(invlogit(data_cumsum_logit_fore$bootsamp[,ik,fh]), 1)
            data_cumsum_logit_fore_add_boot_diff[,ik] = c(data_cumsum_logit_fore_add_boot[1],
                                                          diff(data_cumsum_logit_fore_add_boot))
            rm(ik); rm(data_cumsum_logit_fore_add_boot)
        }
        return(t(apply(data_cumsum_logit_fore_add_boot_diff, 1, quantile, c(alpha/2, 1 - alpha/2))))
    }
    else
    {
        warning("object_interest must either be point or interval.")
    }
}

# fdata: a data matrix of dimension (n by p)
# sex: female or male series
# method_ncomp: way of selecting the number of components
# horizon: forecast horizon 1 to 16
# fore_method: forecasting method
# variable_interest: point or interval forecasts
# CLR_ncomp_selection: when the fore_method = "CLR", it requires a way for selecting number of components

point_fore_national_cdf <- function(fdata, sex, method_ncomp, horizon, fore_method, variable_interest,
                                    CLR_ncomp_selection)
{
    fore_val = matrix(NA, ncol(fdata), (17 - horizon))
    if(fore_method == "CDF")
    {
        for(ij in 1:(17 - horizon))
        {
            fore_val[,ij] <- fore_national_cdf(data = fdata[1:(31 + ij),], ncomp_method = method_ncomp,
                                               fh = horizon, fmethod = "ets",
                                               object_interest = variable_interest)
            rm(ij)
        }
    }
    else if(fore_method == "CLR")
    {
        for(ij in 1:(17 - horizon))
        {
            fore_val[,ij] <- as.numeric(clr_fun(fdata = fdata[1:(31 + ij),], ncomp_selection = CLR_ncomp_selection,
                                                fh = horizon, fore_method = "ETS", object_interest = "point")$fore_count)
            rm(ij)
        }
    }
    else
    {
        warning("forecasting method must either be CDF or CLR.")
    }

    lt_mat = matrix(NA, (17 - horizon), ncol(fdata))
    for(ik in 1:(17 - horizon))
    {
        lt_mat[ik,] = LifeTable(0:110, dx = fore_val[,ik] * 10^5)$lt$ex
        rm(ik)
    }
    if(sex == "female")
    {
        err_lt = c(ftsa:::rmse(forecast = lt_mat, true = female_JPN_ex[(32 + horizon):48,]),
                   ftsa:::mae(forecast = lt_mat,  true = female_JPN_ex[(32 + horizon):48,]))
    }
    else if(sex == "male")
    {
        err_lt = c(ftsa:::rmse(forecast = lt_mat, true = male_JPN_ex[(32 + horizon):48,]),
                   ftsa:::mae(forecast = lt_mat,  true = male_JPN_ex[(32 + horizon):48,]))
    }

    holdout_val_dum = t(matrix(fdata[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata)))
    if(any(holdout_val_dum == 0))
    {
        holdout_val = replace(x = holdout_val_dum, list = which(holdout_val_dum == 0), values = 10^-5)
    }
    else
    {
        holdout_val = holdout_val_dum
    }
    rm(holdout_val_dum)

    # compute the KL divergence and JS divergence

    KL_div_val = JS_div_val = vector("numeric", (17 - horizon))
    for(ij in 1:(17 - horizon))
    {
        KL_div_val[ij] = mean(KLdiv(cbind(fore_val[,ij], holdout_val[,ij]))[2:3])
        JS_div_val[ij] = mean(KLdiv(cbind(fore_val[,ij],
                              apply(cbind(fore_val[,ij], holdout_val[,ij]), 1, geometric.mean)))[2:3])
    }
    err = c(mean(KL_div_val), mean(JS_div_val))
    return(list(err = err, err_lt = err_lt, forecast_pdf = fore_val, lt_mat = lt_mat, holdout_pdf = holdout_val))
}

#####################################################################
# point forecast error for h = 1,...,16 using the CDF method (K = 6)
#####################################################################

point_fore_national_err_female = point_fore_national_err_male =
point_fore_national_err_female_lt = point_fore_national_err_male_lt = matrix(NA, 16, 2)
for(iw in 1:16)
{
    dum = point_fore_national_cdf(fdata = female_JPN_pop, sex = "female", method_ncomp = "provide",
                                  horizon = iw, fore_method = "CDF", variable_interest = "point")
    point_fore_national_err_female[iw,] = dum$err
    point_fore_national_err_female_lt[iw,] = dum$err_lt
    rm(dum)

    dum = point_fore_national_cdf(fdata = male_JPN_pop, sex = "male", method_ncomp = "provide",
                                  horizon = iw, fore_method = "CDF", variable_interest = "point")
    point_fore_national_err_male[iw,] = dum$err
    point_fore_national_err_male_lt[iw,] = dum$err_lt
    print(iw); rm(iw); rm(dum)
}
colnames(point_fore_national_err_female) = colnames(point_fore_national_err_male) = c("KLD", "JSD (geo)")
colnames(point_fore_national_err_female_lt) = colnames(point_fore_national_err_male_lt) = c("RMSE", "MAE")

rownames(point_fore_national_err_female) = rownames(point_fore_national_err_male) =
rownames(point_fore_national_err_female_lt) = rownames(point_fore_national_err_male_lt) = 1:16

######
# EVR
######

point_fore_national_err_female_EVR = point_fore_national_err_male_EVR =
point_fore_national_err_female_EVR_lt = point_fore_national_err_male_EVR_lt = matrix(NA, 16, 2)
for(iw in 1:16)
{
    dum = point_fore_national_cdf(fdata = female_JPN_pop, sex = "female", method_ncomp = "EVR",
                                  horizon = iw, fore_method = "CDF", variable_interest = "point")
    point_fore_national_err_female_EVR[iw,] = dum$err
    point_fore_national_err_female_EVR_lt[iw,] = dum$err_lt
    rm(dum)

    dum = point_fore_national_cdf(fdata = male_JPN_pop, sex = "male", method_ncomp = "EVR",
                                  horizon = iw, fore_method = "CDF", variable_interest = "point")
    point_fore_national_err_male_EVR[iw,] = dum$err
    point_fore_national_err_male_EVR_lt[iw,] = dum$err_lt
    print(iw); rm(iw); rm(dum)
}
colnames(point_fore_national_err_female_EVR) = colnames(point_fore_national_err_male_EVR) = c("KLD", "JSD (geo)")
colnames(point_fore_national_err_female_EVR_lt) = colnames(point_fore_national_err_male_EVR_lt) = c("RMSE", "MAE")

rownames(point_fore_national_err_female_EVR) = rownames(point_fore_national_err_male_EVR) =
rownames(point_fore_national_err_female_EVR_lt) = rownames(point_fore_national_err_male_EVR_lt) = 1:16

#############################
# interval forecast accuracy
#############################

interval_score <- function(holdout, lb, ub, alpha)
{
    lb_ind = ifelse(holdout < lb, 1, 0)
    ub_ind = ifelse(holdout > ub, 1, 0)
    score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
    cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout)
    cpd = abs(cover - (1 - alpha))
    return(c(mean(score), cpd, cover))
}

# fdata: a data matrix of dimension (n by p)
# sex: female or male series
# method_ncomp: K = 6 or EVR
# horizon: forecast horizon 1 to 16
# alpha_val: level of significance between 0 and 1

interval_fore_national_cdf <- function(fdata, sex, method_ncomp, horizon, alpha_val)
{
    forecast_val_lb = forecast_val_ub = matrix(NA, ncol(fdata), (17 - horizon))
    for(ij in 1:(17 - horizon))
    {
        dum = fore_national_cdf(data = fdata[1:(31+ij),], ncomp_method = method_ncomp,
                                fh = horizon, fmethod = "ets",
                                object_interest = "interval", alpha = alpha_val)
        forecast_val_lb[,ij] = dum[,1]
        forecast_val_ub[,ij] = dum[,2]
        rm(ij)
    }

    lt_mat_lb = lt_mat_ub = matrix(NA, (17 - horizon), ncol(fdata))
    for(ik in 1:(17 - horizon))
    {
        lt_mat_lb[ik,] = LifeTable(0:110, dx = forecast_val_lb[,ik] * 10^5)$lt$ex
        lt_mat_ub[ik,] = LifeTable(0:110, dx = forecast_val_ub[,ik] * 10^5)$lt$ex
        rm(ik)
    }
    if(sex == "female")
    {
        int_err_lt = interval_score(holdout = female_JPN_ex[(32 + horizon):48,],
                                    lb = lt_mat_lb, ub = lt_mat_ub, alpha = alpha_val)
    }
    else if(sex == "male")
    {
        int_err_lt = interval_score(holdout = male_JPN_ex[(32 + horizon):48,],
                                    lb = lt_mat_lb, ub = lt_mat_ub, alpha = alpha_val)
    }
    holdout_val = t(matrix(fdata[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata)))
    dx_int_err = interval_score(holdout = holdout_val, lb = forecast_val_lb, ub = forecast_val_ub,
                                alpha = alpha_val)
    return(list(dx_int_err = dx_int_err, lt_int_err = int_err_lt))
}

###################################################
# interval forecast error for h = 1,...,16 (K = 6)
###################################################

interval_fore_national_err_female_PI_80 = interval_fore_national_err_male_PI_80 =
interval_fore_national_err_female_PI_80_lt = interval_fore_national_err_male_PI_80_lt =
interval_fore_national_err_female_PI_95 = interval_fore_national_err_male_PI_95 =
interval_fore_national_err_female_PI_95_lt = interval_fore_national_err_male_PI_95_lt = matrix(NA, 16, 3)
for(iw in 1:16)
{
    # alpha = 0.2

    dum = interval_fore_national_cdf(fdata = female_JPN_pop, sex = "female", method_ncomp = "provide",
                                     horizon = iw, alpha_val = 0.2)
    interval_fore_national_err_female_PI_80[iw,] = dum$dx_int_err
    interval_fore_national_err_female_PI_80_lt[iw,] = dum$lt_int_err
    rm(dum)

    dum = interval_fore_national_cdf(fdata = male_JPN_pop, sex = "female", method_ncomp = "provide",
                                     horizon = iw, alpha_val = 0.2)
    interval_fore_national_err_male_PI_80[iw,] = dum$dx_int_err
    interval_fore_national_err_male_PI_80_lt[iw,] = dum$lt_int_err
    rm(dum)

    # alpha = 0.05

    dum = interval_fore_national_cdf(fdata = female_JPN_pop, sex = "male", method_ncomp = "provide",
                                     horizon = iw, alpha_val = 0.05)
    interval_fore_national_err_female_PI_95[iw,] = dum$dx_int_err
    interval_fore_national_err_female_PI_95_lt[iw,] = dum$lt_int_err
    rm(dum)

    dum = interval_fore_national_cdf(fdata = male_JPN_pop, sex = "male", method_ncomp = "provide",
                                     horizon = iw, alpha_val = 0.05)
    interval_fore_national_err_male_PI_95[iw,] = dum$dx_int_err
    interval_fore_national_err_male_PI_95_lt[iw,] = dum$lt_int_err
    print(iw); rm(iw); rm(dum)
}
colnames(interval_fore_national_err_female_PI_80) = colnames(interval_fore_national_err_male_PI_80) =
colnames(interval_fore_national_err_female_PI_95) = colnames(interval_fore_national_err_male_PI_95) = c("score", "CPD", "ECP")

colnames(interval_fore_national_err_female_PI_80_lt) = colnames(interval_fore_national_err_male_PI_80_lt)
colnames(interval_fore_national_err_female_PI_95_lt) = colnames(interval_fore_national_err_male_PI_95_lt)

rownames(interval_fore_national_err_female_PI_80) = rownames(interval_fore_national_err_male_PI_80) =
rownames(interval_fore_national_err_female_PI_95) = rownames(interval_fore_national_err_male_PI_95) = 1:16

######
# EVR
######

interval_fore_national_err_female_PI_80_EVR = interval_fore_national_err_male_PI_80_EVR =
interval_fore_national_err_female_PI_95_EVR = interval_fore_national_err_male_PI_95_EVR = matrix(NA, 16, 3)
for(iw in 1:16)
{
    # alpha = 0.2

    interval_fore_national_err_female_PI_80_EVR[iw,] = interval_fore_national_cdf(fdata = female_JPN_pop,
                                                                              method_ncomp = "EVR",
                                                                              horizon = iw, alpha_val = 0.2)

    interval_fore_national_err_male_PI_80_EVR[iw,] = interval_fore_national_cdf(fdata = male_JPN_pop,
                                                                            method_ncomp = "EVR",
                                                                            horizon = iw, alpha_val = 0.2)

    # alpha = 0.05

    interval_fore_national_err_female_PI_95_EVR[iw,] = interval_fore_national_cdf(fdata = female_JPN_pop,
                                                                              method_ncomp = "EVR",
                                                                              horizon = iw, alpha_val = 0.05)

    interval_fore_national_err_male_PI_95_EVR[iw,] = interval_fore_national_cdf(fdata = male_JPN_pop,
                                                                            method_ncomp = "EVR",
                                                                            horizon = iw, alpha_val = 0.05)
    print(iw); rm(iw)
}
colnames(interval_fore_national_err_female_PI_80_EVR) = colnames(interval_fore_national_err_male_PI_80_EVR) =
colnames(interval_fore_national_err_female_PI_95_EVR) = colnames(interval_fore_national_err_male_PI_95_EVR) = c("score", "CPD", "ECP")
rownames(interval_fore_national_err_female_PI_80_EVR) = rownames(interval_fore_national_err_male_PI_80_EVR) =
rownames(interval_fore_national_err_female_PI_95_EVR) = rownames(interval_fore_national_err_male_PI_95_EVR) = 1:16
