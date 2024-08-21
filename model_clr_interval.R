##########################
# interval forecast error
##########################

# fdata: (n x p) data matrix
# horizon: forecast horizon
# alpha_val: level of significance
# CLR_ncomp_selection: way of selecting number of components
# no_boot: number of bootstrap samples

interval_fore_clr <- function(fdata, horizon, alpha_val, CLR_ncomp_selection, no_boot)
{
    forecast_val_lb = forecast_val_ub = matrix(NA, ncol(fdata), (17 - horizon))
    for(ij in 1:(17 - horizon))
    {
        dum = clr_fun(fdata = fdata[1:(31+ij),], ncomp_selection = CLR_ncomp_selection,
                      fh = horizon, fore_method = "ets",
                      object_interest = "interval", B = no_boot, alpha = alpha_val)$PI
        forecast_val_lb[,ij] = dum[1,,horizon]
        forecast_val_ub[,ij] = dum[2,,horizon]
        rm(ij)
    }
    holdout_val = t(matrix(fdata[(32 + horizon):48,], length((32 + horizon):48), ncol(fdata)))
    return(interval_score(holdout = holdout_val, lb = forecast_val_lb, ub = forecast_val_ub,
                          alpha = alpha_val))
}

############################################
## interval forecast error for h = 1,...,16
############################################

# ncomp = "EVR"

interval_fore_national_EVR_clr_female_PI_80 = interval_fore_national_EVR_clr_male_PI_80 =
interval_fore_national_EVR_clr_female_PI_95 = interval_fore_national_EVR_clr_male_PI_95 = matrix(NA, 16, 3)
for(iw in 1:16)
{
    # alpha = 0.2

    interval_fore_national_EVR_clr_female_PI_80[iw,] = interval_fore_clr(fdata = female_JPN_pop,
                                                                         horizon = iw, alpha_val = 0.2, CLR_ncomp_selection = "EVR", no_boot = 399)

    interval_fore_national_EVR_clr_male_PI_80[iw,] = interval_fore_clr(fdata = male_JPN_pop,
                                                                       horizon = iw, alpha_val = 0.2, CLR_ncomp_selection = "EVR", no_boot = 399)

    # alpha = 0.05

    interval_fore_national_EVR_clr_female_PI_95[iw,] = interval_fore_clr(fdata = female_JPN_pop,
                                                                         horizon = iw, alpha_val = 0.05, CLR_ncomp_selection = "EVR", no_boot = 399)

    interval_fore_national_EVR_clr_male_PI_95[iw,] = interval_fore_clr(fdata = male_JPN_pop,
                                                                       horizon = iw, alpha_val = 0.05, CLR_ncomp_selection = "EVR", no_boot = 399)
    print(iw); rm(iw)
}

colnames(interval_fore_national_EVR_clr_female_PI_80) = colnames(interval_fore_national_EVR_clr_male_PI_80) =
colnames(interval_fore_national_EVR_clr_female_PI_95) = colnames(interval_fore_national_EVR_clr_male_PI_95) = c("score", "CPD", "ECP")

# ncomp = 6

interval_fore_national_fixed_clr_female_PI_80 = interval_fore_national_fixed_clr_male_PI_80 =
interval_fore_national_fixed_clr_female_PI_95 = interval_fore_national_fixed_clr_male_PI_95 = matrix(NA, 16, 3)
for(iw in 1:16)
{
    # alpha = 0.2

    interval_fore_national_fixed_clr_female_PI_80[iw,] = interval_fore_clr(fdata = female_JPN_pop,
                                                                           horizon = iw, alpha_val = 0.2,
                                                                           CLR_ncomp_selection = "fixed", no_boot = 399)

    interval_fore_national_fixed_clr_male_PI_80[iw,] = interval_fore_clr(fdata = male_JPN_pop,
                                                                         horizon = iw, alpha_val = 0.2,
                                                                         CLR_ncomp_selection = "fixed", no_boot = 399)

    # alpha = 0.05

    interval_fore_national_fixed_clr_female_PI_95[iw,] = interval_fore_clr(fdata = female_JPN_pop,
                                                                           horizon = iw, alpha_val = 0.05,
                                                                           CLR_ncomp_selection = "fixed", no_boot = 399)

    interval_fore_national_fixed_clr_male_PI_95[iw,] = interval_fore_clr(fdata = male_JPN_pop,
                                                                         horizon = iw, alpha_val = 0.05,
                                                                         CLR_ncomp_selection = "fixed", no_boot = 399)
    print(iw); rm(iw)
}

colnames(interval_fore_national_fixed_clr_female_PI_80) = colnames(interval_fore_national_fixed_clr_male_PI_80) =
colnames(interval_fore_national_fixed_clr_female_PI_95) = colnames(interval_fore_national_fixed_clr_male_PI_95) = c("score", "CPD", "ECP")

rownames(interval_fore_national_fixed_clr_female_PI_80) = rownames(interval_fore_national_fixed_clr_male_PI_80) =
rownames(interval_fore_national_fixed_clr_female_PI_95) = rownames(interval_fore_national_fixed_clr_male_PI_95) = 1:16
