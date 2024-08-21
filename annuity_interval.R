##########################################
# Multilevel functional time series model
##########################################

# data_set: a list of p by n data matrix
# aux_var: an aggregated p by n data matrix
# ncomp_method: method for selecting the number of components
# fh: forecast horizon
# fore_method: univariate time-series forecasting method, "ARIMA" or "ETS"
# B: number of bootstrap samples

# (variability is too small)

MLFTS_model_int_all <- function(data_input, aux_var, ncomp_method, fh, fore_method, B)
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

    forerr = matrix(NA, (n_year - ncomp_aggregate), ncomp_aggregate)
    for(i in 1:(n_year - ncomp_aggregate))
    {
        k = i + (ncomp_aggregate - 1)
        fore = matrix(NA, 1, ncomp_aggregate)
        if(fore_method == "ets")
        {
            for(j in 1:ncomp_aggregate)
            {
                fore[,j] = forecast(ets(ftsm_aggregate$coeff[1:k,j]), h = 1)$mean[1]
            }
        }
        else if(fore_method == "arima")
        {
            for(j in 1:ncomp_aggregate)
            {
                fore[,j] = forecast(auto.arima(ftsm_aggregate$coeff[1:k,j]), h = 1)$mean[1]
            }
        }
        forerr[i, ] = ftsm_aggregate$coeff[k + 1,] - fore
    }

    coef_fore_boot = array(NA, dim = c(1000, ncomp_aggregate, fh))
    for(ik in 1:1000)
    {
        coef_fore_boot[ik,,] = coef_fore + t(forerr[sample(1:nrow(forerr), 50, replace = TRUE),])
    }

    ftsm_aggregate_boot = array(NA, dim = c(1000, 110, fh))
    for(ik in 1:1000)
    {
        ftsm_aggregate_boot[ik,,] = ftsm_aggregate$basis %*% coef_fore_boot[ik,,]
    }
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
        forerr = matrix(NA, (n_year - ncomp_resi[[iw]]), ncomp_resi[[iw]])
        for(i in 1:(n_year - ncomp_resi[[iw]]))
        {
            k = i + (ncomp_resi[[iw]] - 1)
            fore = matrix(NA, 1, ncomp_resi[[iw]])
            if(fore_method == "ets")
            {
                for(j in 1:ncomp_resi[[iw]])
                {
                    fore[,j] = forecast(ets(ftsm_resi[[iw]]$coeff[1:k,j]), h = 1)$mean[1]
                }
            }
            else if(fore_method == "arima")
            {
                for(j in 1:ncomp_resi[[iw]])
                {
                    fore[,j] = forecast(auto.arima(ftsm_resi[[iw]]$coeff[1:k,j]), h = 1)$mean[1]
                }
            }
            else
            {
                warning("Forecasting method can either be ARIMA or ETS.")
            }
            forerr[i, ] = ftsm_resi[[iw]]$coeff[k + 1,] - fore
        }

        coef_fore_boot = array(NA, dim = c(1000, ncomp_aggregate, fh))
        for(ik in 1:1000)
        {
            coef_fore_boot[ik,,] = coef_fore_resi_list[[iw]] + t(forerr[sample(1:nrow(forerr), fh, replace = TRUE),])
        }
        coef_fore_resi_list_boot[[iw]] = coef_fore_boot
    }

    # residual bootstrap samples

    resi_fore_boot = list()
    for(iw in 1:n_pop)
    {
        dum = array(NA, dim = c(1000, 110, 50))
        for(ij in 1:1000)
        {
            dum[ij,,] = ftsm_resi[[iw]]$basis[,1:ncomp_resi[[iw]]] %*% (coef_fore_resi_list_boot[[iw]])[ij,,]
        }
        resi_fore_boot[[iw]] = dum
        rm(iw); rm(dum)
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
        dum = array(NA, dim = c(1000, 110, fh))
        for(ij in 1:1000)
        {
            dum[ij,,] = mean_function_list[[iw]] + ftsm_aggregate_boot[ij,,] + aggregate_residuals_boot[,ij] +
                                              (resi_fore_boot[[iw]])[ij,,] + (resi_residual_boot[[iw]])[,ij]
        }
        final_fore_boot[[iw]] = dum
        rm(iw); rm(dum)
    }
    return(final_fore_boot)
}

# interval forecasts

fore_national_cdf_MLFTS_int <- function(data_F, data_M, fh, fmethod, no_boot,
                                        method_ncomp, alpha_level)
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

    data_comb = array(NA, dim = c((ncol(data_F) - 1), nrow(data_F), 2))
    data_comb[,,1] = t(data_cumsum_logit_F)
    data_comb[,,2] = t(data_cumsum_logit_M)

    # implementing the multilevel functional data model

    dum = MLFTS_model_int_all(data_input = data_comb, aux_var = data_common, ncomp_method = method_ncomp,
                          fh = fh, fore_method = fmethod, B = no_boot)
    data_cumsum_logit_F_fore_boot = dum[[1]]
    data_cumsum_logit_M_fore_boot = dum[[2]]
    rm(dum)

    data_cumsum_logit_F_fore_add_boot = data_cumsum_logit_M_fore_add_boot = array(NA, dim = c(1000, 111, 50))
    for(ik in 1:50)
    {
        data_cumsum_logit_F_fore_add_boot[,,ik] = cbind(invlogit(data_cumsum_logit_F_fore_boot[,,ik]), rep(1, no_boot))
        data_cumsum_logit_M_fore_add_boot[,,ik] = cbind(invlogit(data_cumsum_logit_M_fore_boot[,,ik]), rep(1, no_boot))
    }


    data_cumsum_logit_F_fore_add_diff_boot = data_cumsum_logit_M_fore_add_diff_boot = array(NA, dim = c(1000, 111, 50))
    for(iw in 1:no_boot)
    {
        for(ij in 1:50)
        {
            data_cumsum_logit_F_fore_add_diff_boot[iw,,ij] = c(data_cumsum_logit_F_fore_add_boot[iw,1,ij], diff(data_cumsum_logit_F_fore_add_boot[iw,,ij]))
            data_cumsum_logit_M_fore_add_diff_boot[iw,,ij] = c(data_cumsum_logit_M_fore_add_boot[iw,1,ij], diff(data_cumsum_logit_M_fore_add_boot[iw,,ij]))
        }
        rm(iw)
    }
    return(list(boot_F = data_cumsum_logit_F_fore_add_diff_boot, boot_M = data_cumsum_logit_M_fore_add_diff_boot))
}

boot_F_mat = aperm(boot_F, c(2, 1, 3)) * 10^5
boot_M_mat = aperm(boot_M, c(2, 1, 3)) * 10^5

######################################################
# Work out lx on the basis of forecast death count dx
######################################################

lx_female_boot = array(NA, dim = c(111, 50, 1000))
for(iw in 1:1000)
{
    for(ij in 1:50)
    {
        for(ik in 1:111)
        {
            lx_female_boot[ik,ij,iw] = 10^5 - sum(boot_F_mat[1:ik,iw,ij])
        }
    }
}

######################################################
# Work out survival probability px based on dx and lx
######################################################

px_female_boot = array(NA, dim = c(111, 50, 1000))
for(iw in 1:1000)
{
    for(ij in 1:50)
    {
        for(ik in 1:111)
        {
            px_female_boot[,ij,iw] = 1 - round(boot_F_mat[,iw,ij]/c(10^5, lx_female_boot[1:110,ij,iw]),4)
        }
    }
    print(iw)
}

annuities_female_boot = array(NA, dim = c(length(maturities), length(ages), 1000))
for(iu in 1:1000)
{
    for(ij in 1:length(maturities))
    {
        for(iw in 1:length(ages))
        {
            annuities_female_boot[ij,iw,iu] = try(AnnuityPrice_point(y.predict = diag(px_female_boot[61:110,,iu]), age = ages[iw],
                                                  maturity = maturities[ij], inRate = 0.25/100), silent = TRUE)
        }
    }
    print(iu)
}

female_boot_annuities = array(NA, dim = c(length(ages), length(maturities), 1000))
for(iw in 1:1000)
{
    female_boot_annuities[,,iw] = t(matrix(as.numeric(annuities_female_boot[,,iw]), nrow = length(maturities)))
}

female_annuities_lb = round(apply(female_boot_annuities, c(1,2), quantile, 0.025, na.rm = TRUE), 3)
female_annuities_ub = round(apply(female_boot_annuities, c(1,2), quantile, 0.975, na.rm = TRUE), 3)
rownames(female_annuities_lb) = rownames(female_annuities_ub) = ages

xtable(cbind(female_annuities_lb[,1], female_annuities_ub[,1],
             female_annuities_lb[,2], female_annuities_ub[,2],
             female_annuities_lb[,3], female_annuities_ub[,3],
             female_annuities_lb[,4], female_annuities_ub[,4],
             female_annuities_lb[,5], female_annuities_ub[,5],
             female_annuities_lb[,6], female_annuities_ub[,6]), digits = 3)

rm(lx_female_boot); rm(px_female_boot); rm(annuities_female_boot)
