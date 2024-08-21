##########################################
# Multilevel functional time series model
##########################################

# data_set: a list of p by n data matrix
# aux_var: an aggregated p by n data matrix
# ncomp_method: method for selecting the number of components
# fh: forecast horizon
# fore_method: univariate time-series forecasting method, "ARIMA" or "ETS"

MLFTS_model_all <- function(data_input, aux_var, ncomp_method, fh, fore_method)
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
        final_fore[[iw]] = mean_function_list[[iw]] + (aggregate_fore + resi_fore[[iw]])
        rm(iw)
    }
    return(final_fore)
}

# cohort-based point forecasts

fore_national_cdf_MLFTS_out <- function(data_F, data_M, fh, fmethod, method_ncomp)
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

    dum = MLFTS_model_all(data_input = data_comb, aux_var = data_common, ncomp_method = method_ncomp,
                      fh = fh, fore_method = fmethod)
    data_cumsum_logit_F_fore = dum[[1]]
    data_cumsum_logit_M_fore = dum[[2]]
    rm(dum)

    data_cumsum_logit_F_fore_add = data_cumsum_logit_M_fore_add = matrix(NA, ncol(data_F), fh)
    for(ik in 1:fh)
    {
        data_cumsum_logit_F_fore_add[,ik] = c(invlogit(data_cumsum_logit_F_fore[,ik]), 1)
        data_cumsum_logit_M_fore_add[,ik] = c(invlogit(data_cumsum_logit_M_fore[,ik]), 1)
    }

    data_cumsum_logit_F_fore_add_diff = data_cumsum_logit_M_fore_add_diff = matrix(NA, ncol(data_F), fh)
    for(ik in 1:fh)
    {
        data_cumsum_logit_F_fore_add_diff[,ik] = c(data_cumsum_logit_F_fore_add[1,ik], diff(data_cumsum_logit_F_fore_add[,ik]))
        data_cumsum_logit_M_fore_add_diff[,ik] = c(data_cumsum_logit_M_fore_add[1,ik], diff(data_cumsum_logit_M_fore_add[,ik]))
    }
    colnames(data_cumsum_logit_F_fore_add_diff) = colnames(data_cumsum_logit_M_fore_add_diff) = 2023:2072
    return(list(dx_F = data_cumsum_logit_F_fore_add_diff, dx_M = data_cumsum_logit_M_fore_add_diff))
}

dx_fore = fore_national_cdf_MLFTS_out(data_F = female_JPN_pop, data_M = male_JPN_pop, fh = 50,
                                      fmethod = "ets", method_ncomp = "provide")
Japan_F_fore = dx_fore$dx_F * 10^5
Japan_M_fore = dx_fore$dx_M * 10^5

# save figures

savefig("Fig_7a", width = 12, height = 10, toplines = 0.5, type = "png")
plot(fts(0:110, Japan_F_fore), xlab = "Age", ylab = "Life-table death counts",
     main = "Japan: female data (2023-2072)")
dev.off()

savefig("Fig_7b", width = 12, height = 10, toplines = 0.5, type = "png")
plot(fts(0:110, Japan_M_fore), xlab = "Age", ylab = "", main = "Japan: male data (2023-2072)")
dev.off()

# Work out lx

lx_female = lx_male = matrix(NA, 111, 50)
for(ij in 1:50)
{
    for(ik in 1:111)
    {
        lx_female[ik,ij] = 10^5 - sum(Japan_F_fore[1:ik, ij])
        lx_male[ik,ij]   = 10^5 - sum(Japan_M_fore[1:ik, ij])
    }
}

# Work out survival probability px

px_female = px_male = matrix(NA, 111, 50)
for(ij in 1:50)
{
    for(ik in 1:111)
    {
        px_female[,ij] = 1 - round(Japan_F_fore[,ij]/c(10^5, lx_female[1:110,ij]), 4)
        px_male[,ij]   = 1 - round(Japan_M_fore[,ij]/c(10^5, lx_male[1:110,ij]), 4)
    }
}

# Annuity price

AnnuityPrice_point <- function(y.predict, age, maturity, inRate)
{
    if(age < 60 || age > 110) return(print("age is outside range"))
    if(age + maturity > 111) return(print("NA"))

    surv_prob = vector("numeric", maturity)
    for(iw in 1:maturity)
    {
        surv_prob[iw] = y.predict[age - 60 + iw]
    }
    surv_curve = cumprod(surv_prob)

    annuity.price = 0
    for(iw in 1:maturity)
    {
        annuity.price = annuity.price + exp(-inRate * iw) * surv_curve[iw]
    }
    return(round(annuity.price, 4))
}

#############################
# define ages and maturities
#############################

ages = seq(60, 105, by = 5)
maturities = seq(5, 30, by = 5)

#############################
# intererst rate eta = 0.25%
#############################

annuities_female = annuities_male = matrix(NA, length(maturities), length(ages))
for(ij in 1:length(maturities))
{
    for(iw in 1:length(ages))
    {
        annuities_female[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_female[61:110,]), age = ages[iw],
                                                     maturity = maturities[ij], inRate = 0.25/100), silent = TRUE)

        annuities_male[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_male[61:110,]), age = ages[iw],
                                                       maturity = maturities[ij], inRate = 0.25/100), silent = TRUE)
    }
}

annuities_female_mat = t(annuities_female)
annuities_male_mat = t(annuities_male)
rownames(annuities_female_mat) = rownames(annuities_male_mat) = ages
colnames(annuities_female_mat) = colnames(annuities_male_mat) = maturities

annuities_female_mat_numeric = matrix(as.numeric(annuities_female_mat), length(ages),)
annuities_male_mat_numeric   = matrix(as.numeric(annuities_male_mat),   length(ages),)
annuities_both_mat_numeric = cbind(annuities_female_mat_numeric, annuities_male_mat_numeric)
colnames(annuities_both_mat_numeric) = c(maturities, maturities)
rownames(annuities_both_mat_numeric) = ages

xtable(annuities_female_mat_numeric, digits = 3)
xtable(annuities_male_mat_numeric, digits = 3)
xtable(annuities_both_mat_numeric, digits = 3)

###########################
# interest rate eta = 0.1%
###########################

annuities_female_eta_0.001 = annuities_male_eta_0.001 = matrix(NA, length(maturities), length(ages))
for(ij in 1:length(maturities))
{
    for(iw in 1:length(ages))
    {
        annuities_female_eta_0.001[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_female[61:110,]), age = ages[iw],
                                                         maturity = maturities[ij], inRate = 0.1/100), silent = TRUE)

        annuities_male_eta_0.001[ij,iw] = try(AnnuityPrice_point(y.predict = diag(px_male[61:110,]), age = ages[iw],
                                                       maturity = maturities[ij], inRate = 0.1/100), silent = TRUE)
    }
}

annuities_female_mat_eta_0.001 = t(annuities_female_eta_0.001)
annuities_male_mat_eta_0.001 = t(annuities_male_eta_0.001)
rownames(annuities_female_mat_eta_0.001) = rownames(annuities_male_mat_eta_0.001) = ages
colnames(annuities_female_mat_eta_0.001) = colnames(annuities_male_mat_eta_0.001) = maturities

annuities_female_mat_numeric_eta_0.001 = matrix(as.numeric(annuities_female_mat_eta_0.001), length(ages),)
annuities_male_mat_numeric_eta_0.001   = matrix(as.numeric(annuities_male_mat_eta_0.001),   length(ages),)
annuities_both_mat_numeric_eta_0.001 = cbind(annuities_female_mat_numeric_eta_0.001, annuities_male_mat_numeric_eta_0.001)
colnames(annuities_both_mat_numeric_eta_0.001) = c(maturities, maturities)
rownames(annuities_both_mat_numeric_eta_0.001) = ages

xtable(annuities_female_mat_numeric_eta_0.001, digits = 3)
xtable(annuities_male_mat_numeric_eta_0.001, digits = 3)
xtable(annuities_both_mat_numeric_eta_0.001, digits = 3)


rm(lx_female); rm(px_female); rm(annuities_female)
