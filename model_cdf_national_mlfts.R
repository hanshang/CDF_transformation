##########################################
# Multilevel functional time series model
##########################################

# data_set: a list of p by n data matrix
# aux_var: an aggregated p by n data matrix
# ncomp_method: method for selecting the number of components
# fh: forecast horizon
# fore_method: univariate time-series forecasting method, "ARIMA" or "ETS"

MLFTS_model <- function(data_input, aux_var, ncomp_method, fh, fore_method)
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

    # compute the mean function

    mean_function_list = list()
    for(ik in 1:n_pop)
    {
        mean_function_list[[ik]] = rowMeans(data_set_array[,,ik], na.rm = TRUE)
        rm(ik)
    }

    data_set = array(NA, dim = c(n_age, n_year, n_pop))
    for(ik in 1:n_pop)
    {
        data_set[,,ik] = t(scale(t(data_set_array[,,ik]), center = TRUE, scale = FALSE))
        rm(ik)
    }

    if(missing(aux_var)|is.null(aux_var))
    {
        aggregate_data = apply(data_set, c(1, 2), mean)
    }
    else
    {
        aggregate_data = t(aux_var)
    }
    colnames(aggregate_data) = 1:n_year
    rownames(aggregate_data) = 1:n_age

    # 1st FPCA

    eigen_value_aggregate = eigen(cov(t(aggregate_data)))$values
    if(ncomp_method == "EVR")
    {
        ncomp_aggregate = select_K(tau = 10^-3, eigenvalue = eigen_value_aggregate)
    }
    else if(ncomp_method == "provide")
    {
        ncomp_aggregate = 6
    }
    ftsm_aggregate = ftsm(fts(1:n_age, aggregate_data), order = ncomp_aggregate)

    # calculate sum of lambda_k
    sum_lambda_k = sum(eigen_value_aggregate[1:ncomp_aggregate])

    # compute the residual trend
    data_residual = array(NA, dim = c(n_age, n_year, n_pop))
    for(iw in 1:n_pop)
    {
        data_residual[,,iw] = data_set[,,iw] - ftsm_aggregate$fitted$y
        colnames(data_residual[,,iw]) = 1:n_year
        rownames(data_residual[,,iw]) = 1:n_age
        rm(iw)
    }

    # 2nd FPCA

    if(ncomp_method == "EVR")
    {
        ncomp_resi = vector("numeric", n_pop)
        for(iw in 1:n_pop)
        {
            eigen_value_resi = eigen(cov(t(data_residual[,,iw])))$values
            ncomp_resi[iw] = select_K(tau = 10^-3, eigenvalue = eigen_value_resi)
        }
    }
    else if(ncomp_method == "provide")
    {
        ncomp_resi = rep(6, n_pop)
    }

    sum_lambda_l = vector("numeric", n_pop)
    for(iw in 1:n_pop)
    {
        eigen_value_resi = eigen(cov(t(data_residual[,,iw])))$values

        # calculate sum of lambda_l
        sum_lambda_l[iw] = sum(eigen_value_resi[1:(ncomp_resi[iw])])
    }

    ftsm_resi = list()
    for(iw in 1:n_pop)
    {
        ftsm_resi[[iw]] = ftsm(fts(1:n_age, data_residual[,,iw]), order = ncomp_resi[iw])
        rm(iw)
    }

    # within-cluster variability

    within_cluster_variability = vector("numeric", n_pop)
    for(iw in 1:n_pop)
    {
        within_cluster_variability[iw] = sum_lambda_k/(sum_lambda_k + sum_lambda_l[iw])
    }

    # reconstruction

    coef_fore = matrix(NA, ncomp_aggregate, fh)
    if(fore_method == "arima")
    {
        for(ik in 1:ncomp_aggregate)
        {
            coef_fore[ik,] = forecast(auto.arima(ftsm_aggregate$coeff[,ik+1]), h = fh)$mean
        }
    }
    else if(fore_method == "ets")
    {
        for(ik in 1:ncomp_aggregate)
        {
            coef_fore[ik,] = forecast(ets(ftsm_aggregate$coeff[,ik+1]), h = fh)$mean
        }
    }
    else
    {
        warning("Forecasting method can either be ARIMA or ETS.")
    }
    rownames(coef_fore) = 1:ncomp_aggregate
    colnames(coef_fore) = 1:fh

    if(ncomp_aggregate == 1)
    {
        aggregate_fore = as.matrix(ftsm_aggregate$basis[,2]) %*% matrix(coef_fore, nrow = 1)
    }
    else
    {
        aggregate_fore = ftsm_aggregate$basis[,2:(ncomp_aggregate+1)] %*% coef_fore
    }

    # residual forecasts

    coef_fore_resi_list = list()
    for(iw in 1:n_pop)
    {
        coef_fore_resi = matrix(NA, ncomp_resi[iw], fh)
        if(fore_method == "arima")
        {
            for(ik in 1:ncomp_resi[iw])
            {
                coef_fore_resi[ik,] = forecast(auto.arima(ftsm_resi[[iw]]$coeff[,ik+1]), h = fh)$mean
            }
        }
        else if(fore_method == "ets")
        {
            for(ik in 1:ncomp_resi[iw])
            {
                coef_fore_resi[ik,] = forecast(ets(ftsm_resi[[iw]]$coeff[,ik+1]), h = fh)$mean
            }
        }
        else
        {
            warning("Forecasting method can either be ARIMA or ETS.")
        }
        coef_fore_resi_list[[iw]] = coef_fore_resi
        rm(iw)
    }

    resi_fore = list()
    for(iw in 1:n_pop)
    {
        resi_fore[[iw]] = ftsm_resi[[iw]]$basis[,2:(ncomp_resi[iw] + 1)] %*% coef_fore_resi_list[[iw]]
        rm(iw)
    }

    final_fore = list()
    for(iw in 1:n_pop)
    {
        final_fore[[iw]] = mean_function_list[[iw]] + (aggregate_fore + resi_fore[[iw]])[,fh]
        rm(iw)
    }
    return(final_fore)
}

# point forecast evaluation

fore_national_cdf_MLFTS <- function(data_F, data_M, aux_variable, fh, fmethod, method_ncomp)
{
    if(missing(aux_variable)|is.null(aux_variable))
    {
        data_cumsum_F = data_cumsum_M = matrix(NA, nrow(data_F), ncol(data_F))
        for(ij in 1:nrow(data_F))
        {
            data_cumsum_F[ij,] = cumsum(data_F[ij,])
            data_cumsum_M[ij,] = cumsum(data_M[ij,])
            rm(ij)
        }

        data_cumsum_logit_F = data_cumsum_logit_M = matrix(NA, nrow(data_F), (ncol(data_F) - 1))
        for(ij in 1:nrow(data_F))
        {
            data_cumsum_logit_F[ij,] = logit(data_cumsum_F[ij, 1:(ncol(data_F) - 1)])
            data_cumsum_logit_M[ij,] = logit(data_cumsum_M[ij, 1:(ncol(data_M) - 1)])
            rm(ij)
        }
        rownames(data_cumsum_logit_F) = rownames(data_cumsum_logit_M) = years[1:nrow(data_F)]
        colnames(data_cumsum_logit_F) = colnames(data_cumsum_logit_M) = 1:(ncol(data_F) - 1)
        data_common = NULL
    }
    else
    {
        data_cumsum_F = data_cumsum_M = data_cumsum_T = matrix(NA, nrow(data_F), ncol(data_F))
        for(ij in 1:nrow(data_F))
        {
            data_cumsum_F[ij,] = cumsum(data_F[ij,])
            data_cumsum_M[ij,] = cumsum(data_M[ij,])
            data_cumsum_T[ij,] = cumsum(aux_variable[ij,])
            rm(ij)
        }

        data_cumsum_logit_F = data_cumsum_logit_M = data_cumsum_logit_T = matrix(NA, nrow(data_F), (ncol(data_F) - 1))
        for(ij in 1:nrow(data_F))
        {
            data_cumsum_logit_F[ij,] = logit(data_cumsum_F[ij, 1:(ncol(data_F) - 1)])
            data_cumsum_logit_M[ij,] = logit(data_cumsum_M[ij, 1:(ncol(data_M) - 1)])
            data_cumsum_logit_T[ij,] = logit(data_cumsum_T[ij, 1:(ncol(data_M) - 1)])
            rm(ij)
        }
        rownames(data_cumsum_logit_F) = rownames(data_cumsum_logit_M) = rownames(data_cumsum_logit_T) = years[1:nrow(data_F)]
        colnames(data_cumsum_logit_F) = colnames(data_cumsum_logit_M) = colnames(data_cumsum_logit_T) = 1:(ncol(data_F) - 1)
        data_common = data_cumsum_logit_T
    }
    data_comb = array(NA, dim = c((ncol(data_F) - 1), nrow(data_F), 2))
    data_comb[,,1] = t(data_cumsum_logit_F)
    data_comb[,,2] = t(data_cumsum_logit_M)

    # implementing the multilevel functional data model

    dum = MLFTS_model(data_input = data_comb, aux_var = data_common, ncomp_method = method_ncomp,
                      fh = fh, fore_method = fmethod)
    data_cumsum_logit_F_fore = dum[[1]]
    data_cumsum_logit_M_fore = dum[[2]]
    rm(dum)

    data_cumsum_logit_F_fore_add = c(invlogit(data_cumsum_logit_F_fore), 1)
    data_cumsum_logit_M_fore_add = c(invlogit(data_cumsum_logit_M_fore), 1)

    data_cumsum_logit_F_fore_add_diff = c(data_cumsum_logit_F_fore_add[1], diff(data_cumsum_logit_F_fore_add))
    data_cumsum_logit_M_fore_add_diff = c(data_cumsum_logit_M_fore_add[1], diff(data_cumsum_logit_M_fore_add))
    return(list(mlfts_fore_F = data_cumsum_logit_F_fore_add_diff,
                mlfts_fore_M = data_cumsum_logit_M_fore_add_diff))
}

# fdata_F: a data matrix of dimension (n by p)
# fdata_M: a data matrix of dimension (n by p)
# fdata_common: a common data matrix
# horizon: forecast horizon 1 to 16
# way_ncomp: way of selecting the number of components

point_fore_national_cdf_MLFTS <- function(fdata_F, fdata_M, fdata_common, horizon, way_ncomp)
{
    forecast_val_F = forecast_val_M = matrix(NA, ncol(fdata_F), (17 - horizon))
    for(ij in 1:(17 - horizon))
    {
        dum <- fore_national_cdf_MLFTS(data_F = fdata_F[1:(31+ij),], data_M = fdata_M[1:(31+ij),],
                                       aux_variable = fdata_common[1:(31+ij),], fh = horizon, fmethod = "ets",
                                       method_ncomp = way_ncomp)
        forecast_val_F[,ij] = dum$mlfts_fore_F
        forecast_val_M[,ij] = dum$mlfts_fore_M
        rm(ij); rm(dum)
    }
    rownames(forecast_val_F) = rownames(forecast_val_M) = 1:ncol(fdata_F)
    colnames(forecast_val_F) = colnames(forecast_val_M) = 1:(17 - horizon)

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

    holdout_val_F = t(matrix(fdata_F[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
    holdout_val_M = t(matrix(fdata_M[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))

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
                err_F = err_F, err_M = err_M, err_lt_F = err_lt_F, err_lt_M = err_lt_M))
}

################################################
# point forecast error for h = 1,...,16 (K = 6)
################################################

point_fore_national_err_MLFTS_F = point_fore_national_err_MLFTS_M =
point_fore_national_err_MLFTS_F_lt = point_fore_national_err_MLFTS_M_lt = matrix(NA, 16, 2)
for(iw in 1:16)
{
    dum = point_fore_national_cdf_MLFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop,
                                        fdata_common = NULL, horizon = iw, way_ncomp = "provide")
    point_fore_national_err_MLFTS_F[iw,] = dum$err_F
    point_fore_national_err_MLFTS_M[iw,] = dum$err_M
    point_fore_national_err_MLFTS_F_lt[iw,] = dum$err_lt_F
    point_fore_national_err_MLFTS_M_lt[iw,] = dum$err_lt_M
    print(iw); rm(iw); rm(dum)
}
colnames(point_fore_national_err_MLFTS_F) = colnames(point_fore_national_err_MLFTS_M) = c("KLD", "JSD")
colnames(point_fore_national_err_MLFTS_F_lt) = colnames(point_fore_national_err_MLFTS_M_lt) = c("RMSE", "MAE")
rownames(point_fore_national_err_MLFTS_F) = rownames(point_fore_national_err_MLFTS_M) =
rownames(point_fore_national_err_MLFTS_F_lt) = rownames(point_fore_national_err_MLFTS_M_lt) = 1:16

######
# EVR
######

point_fore_national_err_MLFTS_F_EVR = point_fore_national_err_MLFTS_M_EVR =
point_fore_national_err_MLFTS_F_EVR_lt = point_fore_national_err_MLFTS_M_EVR_lt = matrix(NA, 16, 2)
for(iw in 1:16)
{
    dum = point_fore_national_cdf_MLFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop,
                                        fdata_common = NULL, horizon = iw, way_ncomp = "EVR")
    point_fore_national_err_MLFTS_F_EVR[iw,] = dum$err_F
    point_fore_national_err_MLFTS_M_EVR[iw,] = dum$err_M
    point_fore_national_err_MLFTS_F_EVR_lt[iw,] = dum$err_lt_F
    point_fore_national_err_MLFTS_M_EVR_lt[iw,] = dum$err_lt_M
    print(iw); rm(iw); rm(dum)
}
colnames(point_fore_national_err_MLFTS_F_EVR) = colnames(point_fore_national_err_MLFTS_M_EVR) = c("KLD", "JSD")
colnames(point_fore_national_err_MLFTS_F_EVR_lt) = colnames(point_fore_national_err_MLFTS_M_EVR_lt) = c("RMSE", "MAE")
rownames(point_fore_national_err_MLFTS_F_EVR) = rownames(point_fore_national_err_MLFTS_M_EVR) =
rownames(point_fore_national_err_MLFTS_F_EVR_lt) = rownames(point_fore_national_err_MLFTS_M_EVR_lt) = 1:16

# aux_var = "total"

# point_fore_national_err_total_MLFTS_F = point_fore_national_err_total_MLFTS_M = matrix(NA, 16, 2)
# for(iw in 1:16)
# {
#    dum = point_fore_national_cdf_MLFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop,
#                                        fdata_common = total_JPN_pop, horizon = iw)
#    point_fore_national_err_total_MLFTS_F[iw,] = dum$err_F
#    point_fore_national_err_total_MLFTS_M[iw,]   = dum$err_M
#    print(iw); rm(iw); rm(dum)
# }
# colnames(point_fore_national_err_total_MLFTS_F) = colnames(point_fore_national_err_total_MLFTS_M) = c("KLD", "JSD (geo)")
# rownames(point_fore_national_err_total_MLFTS_F) = rownames(point_fore_national_err_total_MLFTS_M) = 1:16

# round(colMeans(point_fore_national_err_total_MLFTS_F), 4) # 0.0097 0.0027
# round(colMeans(point_fore_national_err_total_MLFTS_M), 4) # 0.0040 0.0011
