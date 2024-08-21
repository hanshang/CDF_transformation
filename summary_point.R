###########
## summary
###########

# UFTS

round(colMeans(point_fore_national_err_female), 4) # 0.0118 0.0032
round(colMeans(point_fore_national_err_male), 4)   # 0.0044 0.0011

round(colMeans(point_fore_national_err_female_EVR), 4) # 0.0127 0.0035
round(colMeans(point_fore_national_err_male_EVR), 4)   # 0.0045 0.0011

xtable(rbind(point_fore_national_err_female, colMeans(point_fore_national_err_female)), digits = 4)
xtable(rbind(point_fore_national_err_male,   colMeans(point_fore_national_err_male)),   digits = 4)

# MFTS

round(colMeans(point_fore_national_err_MFTS_F), 4) # 0.0072 0.0019
round(colMeans(point_fore_national_err_MFTS_M), 4) # 0.0041 0.0010

round(colMeans(point_fore_national_err_MFTS_F_EVR), 4) # 0.0097 0.0027
round(colMeans(point_fore_national_err_MFTS_M_EVR), 4) # 0.0041 0.0011

xtable(rbind(point_fore_national_err_MFTS_F, colMeans(point_fore_national_err_MFTS_F)), digits = 4)
xtable(rbind(point_fore_national_err_MFTS_M, colMeans(point_fore_national_err_MFTS_M)), digits = 4)

# MLFTS

round(colMeans(point_fore_national_err_MLFTS_F), 4) # 0.0070 0.0019
round(colMeans(point_fore_national_err_MLFTS_M), 4) # 0.0033 0.0009

round(colMeans(point_fore_national_err_MLFTS_F_EVR), 4) # 0.0045 0.0012
round(colMeans(point_fore_national_err_MLFTS_M_EVR), 4) # 0.0038 0.0010

xtable(rbind(point_fore_national_err_MLFTS_F, colMeans(point_fore_national_err_MLFTS_F)), digits = 4)
xtable(rbind(point_fore_national_err_MLFTS_M, colMeans(point_fore_national_err_MLFTS_M)), digits = 4)

#################
# xtable summary
#################

##########
## errors
##########

# EVR

point_fore_national_err_F_EVR = cbind(point_fore_national_err_female_EVR,
                                      point_fore_national_err_MFTS_F_EVR,
                                      point_fore_national_err_MLFTS_F_EVR,
                                      point_fore_national_clr_EVR_err_female)

point_fore_national_err_M_EVR = cbind(point_fore_national_err_male_EVR,
                                      point_fore_national_err_MFTS_M_EVR,
                                      point_fore_national_err_MLFTS_M_EVR,
                                      point_fore_national_clr_EVR_err_male)

# K = 6

point_fore_national_err_F = cbind(point_fore_national_err_female,
                                  point_fore_national_err_MFTS_F,
                                  point_fore_national_err_MLFTS_F,
                                  point_fore_national_clr_fixed_err_female)

point_fore_national_err_M = cbind(point_fore_national_err_male,
                                  point_fore_national_err_MFTS_M,
                                  point_fore_national_err_MLFTS_M,
                                  point_fore_national_clr_fixed_err_male)

point_fore_national_err_F_mean = rbind(point_fore_national_err_F, colMeans(point_fore_national_err_F))
point_fore_national_err_M_mean = rbind(point_fore_national_err_M, colMeans(point_fore_national_err_M))

colnames(point_fore_national_err_F_mean) = colnames(point_fore_national_err_M_mean) = rep(c("KLD", "JSD"), 5)
rownames(point_fore_national_err_F_mean) = rownames(point_fore_national_err_M_mean) = c(1:16, "Mean")

xtable(point_fore_national_err_F_mean * 100, digits = 4)
xtable(point_fore_national_err_M_mean * 100, digits = 4)

#############
## which min
#############

# EVR

apply(point_fore_national_err_F_EVR[,seq(1, 8, by = 2)], 1, which.min)
apply(point_fore_national_err_F_EVR[,seq(2, 8, by = 2)], 1, which.min)

apply(point_fore_national_err_M_EVR[,seq(1, 8, by = 2)], 1, which.min)
apply(point_fore_national_err_M_EVR[,seq(2, 8, by = 2)], 1, which.min)

# K = 6

apply(point_fore_national_err_F[,seq(1, 10, by = 2)], 1, which.min)
apply(point_fore_national_err_F[,seq(2, 10, by = 2)], 1, which.min)

apply(point_fore_national_err_M[,seq(1, 10, by = 2)], 1, which.min)
apply(point_fore_national_err_M[,seq(2, 10, by = 2)], 1, which.min)

#####
# e0
#####

e0_err_K6_female = cbind(point_fore_national_err_female_lt,
                         point_fore_national_err_MFTS_F_lt,
                         point_fore_national_err_MLFTS_F_lt,
                         point_fore_national_clr_fixed_err_female_lt)

e0_err_K6_male = cbind(point_fore_national_err_male_lt,
                       point_fore_national_err_MFTS_M_lt,
                       point_fore_national_err_MLFTS_M_lt,
                       point_fore_national_clr_fixed_err_male_lt)

e0_err_EVR_female = cbind(point_fore_national_err_female_EVR_lt,
                          point_fore_national_err_MFTS_F_EVR_lt,
                          point_fore_national_err_MLFTS_F_EVR_lt,
                          point_fore_national_clr_EVR_err_female_lt)

e0_err_EVR_male = cbind(point_fore_national_err_male_EVR_lt,
                        point_fore_national_err_MFTS_M_EVR_lt,
                        point_fore_national_err_MLFTS_M_EVR_lt,
                        point_fore_national_clr_EVR_err_male_lt)

require(xtable)
xtable(rbind(e0_err_K6_female, colMeans(e0_err_K6_female)), digits = 4)
xtable(rbind(e0_err_K6_male,   colMeans(e0_err_K6_male)),   digits = 4)

xtable(rbind(e0_err_EVR_female, colMeans(e0_err_EVR_female)), digits = 4)
xtable(rbind(e0_err_EVR_male,   colMeans(e0_err_EVR_male)),   digits = 4)

############
# which.min
############

apply(e0_err_K6_female[,seq(1,8,by=2)], 1, which.min)
apply(e0_err_K6_female[,seq(2,8,by=2)], 1, which.min)

apply(e0_err_K6_male[,seq(1,8,by=2)], 1, which.min)
apply(e0_err_K6_male[,seq(2,8,by=2)], 1, which.min)

apply(e0_err_EVR_female[,seq(1,8,by=2)], 1, which.min)
apply(e0_err_EVR_female[,seq(2,8,by=2)], 1, which.min)

apply(e0_err_EVR_male[,seq(1,8,by=2)], 1, which.min)
apply(e0_err_EVR_male[,seq(2,8,by=2)], 1, which.min)
