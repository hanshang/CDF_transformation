##########################################
# Multilevel functional time series model
##########################################

# data_set: a list of p by n data matrix
# aux_var: an aggregated p by n data matrix
# ncomp_method: method for selecting the number of components
# fh: forecast horizon
# fore_method: univariate time-series forecasting method, "ARIMA" or "ETS"
# B: number of bootstrap samples

MLFTS_model_int <- function(data_input, aux_var, ncomp_method, fh, fore_method, B)
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
    ftsm_aggregate = ftsm(fts(1:n_age, aggregate_data), order = ncomp_aggregate, mean = FALSE)

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
        data_residual_mat = data_residual[,,iw]
        colnames(data_residual_mat) = 1:n_year
        ftsm_resi[[iw]] = ftsm(fts(1:n_age, data_residual_mat), order = ncomp_resi[iw], mean = FALSE)
        rm(iw); rm(data_residual_mat)
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
            coef_fore[ik,] = forecast(auto.arima(ftsm_aggregate$coeff[,ik]), h = fh)$mean
        }
    }
    else if(fore_method == "ets")
    {
        for(ik in 1:ncomp_aggregate)
        {
            coef_fore[ik,] = forecast(ets(ftsm_aggregate$coeff[,ik]), h = fh)$mean
        }
    }
    else
    {
        warning("Forecasting method can either be ARIMA or ETS.")
    }
    rownames(coef_fore) = 1:ncomp_aggregate
    colnames(coef_fore) = 1:fh

    forerr = matrix(NA, (n_year - ncomp_aggregate - fh + 1), ncomp_aggregate)
    for(i in fh:(n_year - ncomp_aggregate))
    {
        k = i + (ncomp_aggregate - fh)
        fore = matrix(NA, 1, ncomp_aggregate)
        if(fore_method == "ets")
        {
            for(j in 1:ncomp_aggregate)
            {
                fore[,j] = forecast(ets(ftsm_aggregate$coeff[1:k,j]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "arima")
        {
            for(j in 1:ncomp_aggregate)
            {
                fore[,j] = forecast(auto.arima(ftsm_aggregate$coeff[1:k,j]), h = fh)$mean[fh]
            }
        }
        forerr[i - fh + 1, ] = ftsm_aggregate$coeff[k + fh,] - fore
    }

    coef_fore_boot = t(matrix(rep(coef_fore[,fh], B), B,,byrow = TRUE) + forerr[sample(1:nrow(forerr), B, replace = TRUE),])
    ftsm_aggregate_boot = ftsm_aggregate$basis %*% coef_fore_boot
    rm(coef_fore_boot)

    aggregate_residuals_boot = ftsm_aggregate$residuals$y[,sample(1:n_year, size = B, replace = TRUE)]

    # residual forecasts

    coef_fore_resi_list = list()
    for(iw in 1:n_pop)
    {
        coef_fore_resi = matrix(NA, ncomp_resi[iw], fh)
        if(fore_method == "arima")
        {
            for(ik in 1:ncomp_resi[iw])
            {
                coef_fore_resi[ik,] = forecast(auto.arima(ftsm_resi[[iw]]$coeff[,ik]), h = fh)$mean
            }
        }
        else if(fore_method == "ets")
        {
            for(ik in 1:ncomp_resi[iw])
            {
                coef_fore_resi[ik,] = forecast(ets(ftsm_resi[[iw]]$coeff[,ik]), h = fh)$mean
            }
        }
        else
        {
            warning("Forecasting method can either be ARIMA or ETS.")
        }
        coef_fore_resi_list[[iw]] = coef_fore_resi
        rm(iw)
    }

    # add noise

    coef_fore_resi_list_boot = list()
    for(iw in 1:n_pop)
    {
        forerr = matrix(NA, (n_year - ncomp_resi[[iw]] - fh + 1), ncomp_resi[[iw]])
        for(i in fh:(n_year - ncomp_resi[[iw]]))
        {
            k = i + (ncomp_resi[[iw]] - fh)
            fore = matrix(NA, 1, ncomp_resi[[iw]])
            if(fore_method == "ets")
            {
                for(j in 1:ncomp_resi[[iw]])
                {
                    fore[,j] = forecast(ets(ftsm_resi[[iw]]$coeff[1:k,j]), h = fh)$mean[fh]
                }
            }
            else if(fore_method == "arima")
            {
                for(j in 1:ncomp_resi[[iw]])
                {
                    fore[,j] = forecast(auto.arima(ftsm_resi[[iw]]$coeff[1:k,j]), h = fh)$mean[fh]
                }
            }
            else
            {
                warning("Forecasting method can either be ARIMA or ETS.")
            }
            forerr[i - fh + 1, ] = ftsm_resi[[iw]]$coeff[k + fh,] - fore
        }
        coef_fore_resi_list_boot[[iw]] = t(matrix(rep((coef_fore_resi_list[[iw]])[,fh], B), B,,byrow = TRUE) + forerr[sample(1:nrow(forerr), B, replace = TRUE),])
    }

    # residual bootstrap samples

    resi_fore_boot = list()
    for(iw in 1:n_pop)
    {
        resi_fore_boot[[iw]] = ftsm_resi[[iw]]$basis[,1:ncomp_resi[[iw]]] %*% coef_fore_resi_list_boot[[iw]]
        rm(iw)
    }

    # individual residual bootstrap samples

    resi_residual_boot = list()
    for(iw in 1:n_pop)
    {
        resi_residual_boot[[iw]] = ftsm_resi[[iw]]$residuals$y[,sample(1:n_year, size = B, replace = TRUE)]
        rm(iw)
    }

    # add algoether

    final_fore_boot = list()
    for(iw in 1:n_pop)
    {
        final_fore_boot[[iw]] = mean_function_list[[iw]] + ftsm_aggregate_boot + aggregate_residuals_boot +
                                              resi_fore_boot[[iw]] + resi_residual_boot[[iw]]
        #final_fore_boot[[iw]] = mean_function_list[[iw]] + ftsm_aggregate_boot + resi_fore_boot[[iw]] + resi_residual_boot[[iw]]
        rm(iw)
    }
    return(final_fore_boot)
}

# point forecast errors

fore_national_cdf_MLFTS_int <- function(data_F, data_M, aux_variable, fh, fmethod, no_boot,
                                        method_ncomp, alpha_level)
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

    dum = MLFTS_model_int(data_input = data_comb, aux_var = data_common, ncomp_method = method_ncomp,
                          fh = fh, fore_method = fmethod, B = no_boot)
    data_cumsum_logit_F_fore_boot = dum[[1]]
    data_cumsum_logit_M_fore_boot = dum[[2]]
    rm(dum)

    data_cumsum_logit_F_fore_add_boot = rbind(invlogit(data_cumsum_logit_F_fore_boot), rep(1, no_boot))
    data_cumsum_logit_M_fore_add_boot = rbind(invlogit(data_cumsum_logit_M_fore_boot), rep(1, no_boot))

    data_cumsum_logit_F_fore_add_diff_boot = data_cumsum_logit_M_fore_add_diff_boot = matrix(NA, ncol(data_F), no_boot)
    for(iw in 1:no_boot)
    {
        data_cumsum_logit_F_fore_add_diff_boot[,iw] = c(data_cumsum_logit_F_fore_add_boot[1,iw], diff(data_cumsum_logit_F_fore_add_boot[,iw]))
        data_cumsum_logit_M_fore_add_diff_boot[,iw] = c(data_cumsum_logit_M_fore_add_boot[1,iw], diff(data_cumsum_logit_M_fore_add_boot[,iw]))
        rm(iw)
    }
    data_cumsum_logit_F_fore_add_diff_PI = t(apply(data_cumsum_logit_F_fore_add_diff_boot, 1, quantile, c(alpha_level/2, (1-alpha_level/2))))
    data_cumsum_logit_M_fore_add_diff_PI = t(apply(data_cumsum_logit_M_fore_add_diff_boot, 1, quantile, c(alpha_level/2, (1-alpha_level/2))))

    return(list(mlfts_fore_F_boot = data_cumsum_logit_F_fore_add_diff_PI,
                mlfts_fore_M_boot = data_cumsum_logit_M_fore_add_diff_PI))
}

# interval forecast errors

int_fore_national_cdf_MLFTS <- function(fdata_F, fdata_M, fdata_common, horizon, way_ncomp, boot_number, level_alpha)
{
    forecast_val_F = forecast_val_M = array(NA, dim = c(ncol(fdata_F), 2, (17 - horizon)))
    for(ij in 1:(17 - horizon))
    {
        dum <- fore_national_cdf_MLFTS_int(data_F = fdata_F[1:(31+ij),], data_M = fdata_M[1:(31+ij),],
                                           aux_variable = fdata_common[1:(31+ij),], fh = horizon, fmethod = "ets",
                                           method_ncomp = way_ncomp, no_boot = boot_number, alpha_level = level_alpha)
        forecast_val_F[,,ij] = dum$mlfts_fore_F_boot
        forecast_val_M[,,ij] = dum$mlfts_fore_M_boot
        rm(ij); rm(dum)
    }

    holdout_val_F = t(matrix(fdata_F[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))
    holdout_val_M = t(matrix(fdata_M[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata_F)))

    int_err_F = interval_score(holdout = holdout_val_F, lb = forecast_val_F[,1,], ub = forecast_val_F[,2,], alpha = level_alpha)
    int_err_M = interval_score(holdout = holdout_val_M, lb = forecast_val_M[,1,], ub = forecast_val_M[,2,], alpha = level_alpha)

    return(list(int_err_F = int_err_F, int_err_M = int_err_M))
}

#####################################################
### interval forecast error for h = 1,...,16 (K = 6)
#####################################################

## alpha_level = 0.2

# EVR

int_fore_national_err_MLFTS_F_EVR = int_fore_national_err_MLFTS_M_EVR = matrix(NA, 16, 3)
for(iw in 1:16)
{
    dum = int_fore_national_cdf_MLFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop,
                                      fdata_common = NULL, horizon = iw, way_ncomp = "EVR",
                                      boot_number = 1000, level_alpha = 0.2)
    int_fore_national_err_MLFTS_F_EVR[iw,] = dum$int_err_F
    int_fore_national_err_MLFTS_M_EVR[iw,]   = dum$int_err_M
    print(iw); rm(iw); rm(dum)
}
colnames(int_fore_national_err_MLFTS_F_EVR) = colnames(int_fore_national_err_MLFTS_M_EVR) = c("score", "CPD", "ECP")
rownames(int_fore_national_err_MLFTS_F_EVR) = rownames(int_fore_national_err_MLFTS_M_EVR) = 1:16

# K = 6

int_fore_national_err_MLFTS_F = int_fore_national_err_MLFTS_M = matrix(NA, 16, 3)
for(iw in 1:16)
{
    dum = int_fore_national_cdf_MLFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop,
                                      fdata_common = NULL, horizon = iw, way_ncomp = "provide",
                                      boot_number = 1000, level_alpha = 0.2)
    int_fore_national_err_MLFTS_F[iw,] = dum$int_err_F
    int_fore_national_err_MLFTS_M[iw,]   = dum$int_err_M
    print(iw); rm(iw); rm(dum)
}
colnames(int_fore_national_err_MLFTS_F) = colnames(int_fore_national_err_MLFTS_M) = c("score", "CPD", "ECP")
rownames(int_fore_national_err_MLFTS_F) = rownames(int_fore_national_err_MLFTS_M) = 1:16

## alpha_level = 0.05

# EVR

int_fore_national_err_MLFTS_F_PI_95_EVR = int_fore_national_err_MLFTS_M_PI_95_EVR = matrix(NA, 16, 3)
for(iw in 1:16)
{
    dum = int_fore_national_cdf_MLFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop,
                                      fdata_common = NULL, horizon = iw, way_ncomp = "EVR",
                                      boot_number = 1000, level_alpha = 0.05)
    int_fore_national_err_MLFTS_F_PI_95_EVR[iw,] = dum$int_err_F
    int_fore_national_err_MLFTS_M_PI_95_EVR[iw,]   = dum$int_err_M
    print(iw); rm(iw); rm(dum)
}
colnames(int_fore_national_err_MLFTS_F_PI_95_EVR) = colnames(int_fore_national_err_MLFTS_M_PI_95_EVR) = c("score", "CPD", "ECP")
rownames(int_fore_national_err_MLFTS_F_PI_95_EVR) = rownames(int_fore_national_err_MLFTS_M_PI_95_EVR) = 1:16

# K = 6

int_fore_national_err_MLFTS_F_PI_95 = int_fore_national_err_MLFTS_M_PI_95 = matrix(NA, 16, 3)
for(iw in 1:16)
{
    dum = int_fore_national_cdf_MLFTS(fdata_F = female_JPN_pop, fdata_M = male_JPN_pop,
                                      fdata_common = NULL, horizon = iw, way_ncomp = "provide", level_alpha = 0.05)
    int_fore_national_err_MLFTS_F_PI_95[iw,] = dum$int_err_F
    int_fore_national_err_MLFTS_M_PI_95[iw,]   = dum$int_err_M
    print(iw); rm(iw); rm(dum)
}
colnames(int_fore_national_err_MLFTS_F_PI_95) = colnames(int_fore_national_err_MLFTS_M_PI_95) = c("score", "CPD", "ECP")
rownames(int_fore_national_err_MLFTS_F_PI_95) = rownames(int_fore_national_err_MLFTS_M_PI_95) = 1:16
