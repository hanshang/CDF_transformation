###################
# summary interval
###################

# UFTS (K = 6)

round(colMeans(interval_fore_national_err_female_PI_80), 4) # 0.0046 0.0682 0.8682
round(colMeans(interval_fore_national_err_male_PI_80), 4)   # 0.0036 0.0734 0.8562

round(colMeans(interval_fore_national_err_female_PI_95), 4) # 0.0061 0.0095 0.9546
round(colMeans(interval_fore_national_err_male_PI_95), 4)   # 0.0049 0.0205 0.9679

# UFTS (EVR)

round(colMeans(interval_fore_national_err_female_PI_80_EVR), 4) # 0.0060 0.1819 0.9819
round(colMeans(interval_fore_national_err_male_PI_80_EVR), 4)   # 0.0053 0.1965 0.9965

round(colMeans(interval_fore_national_err_female_PI_95_EVR), 4) # 0.0094 0.0498 0.9998
round(colMeans(interval_fore_national_err_male_PI_95_EVR), 4)   # 0.0082 0.0494 0.9994

# CLR (K = 6)

round(colMeans(interval_fore_national_fixed_clr_female_PI_80), 4) # 0.0045 0.0959 0.8959
round(colMeans(interval_fore_national_fixed_clr_male_PI_80), 4)   # 0.0031 0.0904 0.8904

round(colMeans(interval_fore_national_fixed_clr_female_PI_95), 4) # 0.0061 0.0245 0.9611
round(colMeans(interval_fore_national_fixed_clr_male_PI_95), 4)   # 0.0043 0.0201 0.9701

# CLR (EVR)

round(colMeans(interval_fore_national_EVR_clr_female_PI_80), 4) # 0.0047 0.1112 0.9112
round(colMeans(interval_fore_national_EVR_clr_male_PI_80), 4)   # 0.0035 0.0888 0.8888

round(colMeans(interval_fore_national_EVR_clr_female_PI_95), 4) # 0.0076 0.0301 0.9725
round(colMeans(interval_fore_national_EVR_clr_male_PI_95), 4)   # 0.0054 0.0269 0.9759

# MFTS (K = 6)

round(colMeans(interval_fore_national_err_MFTS_F), 4)       # 0.0051 0.1322 0.9322
round(colMeans(interval_fore_national_err_MFTS_M), 4)       # 0.0035 0.1179 0.9179

round(colMeans(interval_fore_national_err_MFTS_F_PI_95), 4) # 0.0070 0.0301 0.9801
round(colMeans(interval_fore_national_err_MFTS_M_PI_95), 4) # 0.0047 0.0407 0.9907

# MFTS (EVR)

round(colMeans(interval_fore_national_err_MFTS_F_EVR), 4)       # 0.0050 0.1309 0.9309
round(colMeans(interval_fore_national_err_MFTS_M_EVR), 4)       # 0.0035 0.1117 0.9117

round(colMeans(interval_fore_national_err_MFTS_F_PI_95_EVR), 4) # 0.0071 0.0305 0.9805
round(colMeans(interval_fore_national_err_MFTS_M_PI_95_EVR), 4) # 0.0047 0.0414 0.9914

# MLFTS (K = 6)

round(colMeans(int_fore_national_err_MLFTS_F), 4) # 0.0046 0.0310 0.8124
round(colMeans(int_fore_national_err_MLFTS_M), 4) # 0.0032 0.0540 0.8500

round(colMeans(int_fore_national_err_MLFTS_F_PI_95), 4) # 0.0061 0.0376 0.9142
round(colMeans(int_fore_national_err_MLFTS_M_PI_95), 4) # 0.0039 0.0159 0.9522

##########
# summary
##########

### 80% nominal coverage probability

## EVR

int_summary_female_EVR_PI_80 = cbind(interval_fore_national_err_female_PI_80_EVR[,1:2],
                                     interval_fore_national_err_MFTS_F_EVR[,1:2],
                                     int_fore_national_err_MLFTS_F_EVR[,1:2],
                                     interval_fore_national_EVR_clr_female_PI_80[,1:2])

int_summary_male_EVR_PI_80 = cbind(interval_fore_national_err_male_PI_80_EVR[,1:2],
                                   interval_fore_national_err_MFTS_M_EVR[,1:2],
                                   int_fore_national_err_MLFTS_M_EVR[,1:2],
                                   interval_fore_national_EVR_clr_male_PI_80[,1:2])

xtable(rbind(int_summary_female_EVR_PI_80, colMeans(int_summary_female_EVR_PI_80)), digits = 4)
xtable(rbind(int_summary_male_EVR_PI_80,   colMeans(int_summary_male_EVR_PI_80)), digits = 4)

## K = 6

int_summary_female_PI_80 = cbind(interval_fore_national_err_female_PI_80[,1:2],
                                 interval_fore_national_err_MFTS_F[,1:2],
                                 int_fore_national_err_MLFTS_F[,1:2],
                                 interval_fore_national_fixed_clr_female_PI_80[,1:2])

int_summary_male_PI_80 = cbind(interval_fore_national_err_male_PI_80[,1:2],
                               interval_fore_national_err_MFTS_M[,1:2],
                               int_fore_national_err_MLFTS_M[,1:2],
                               interval_fore_national_fixed_clr_male_PI_80[,1:2])

xtable(rbind(int_summary_female_PI_80, colMeans(int_summary_female_PI_80)), digits = 4)
xtable(rbind(int_summary_male_PI_80, colMeans(int_summary_male_PI_80)), digits = 4)

#############
## which min
#############

## EVR

apply(int_summary_female_EVR_PI_80[,seq(1, 8, by = 2)], 1, which.min)
apply(int_summary_female_EVR_PI_80[,seq(2, 8, by = 2)], 1, which.min)

apply(int_summary_male_EVR_PI_80[,seq(1, 8, by = 2)], 1, which.min)
apply(int_summary_male_EVR_PI_80[,seq(2, 8, by = 2)], 1, which.min)

## K = 6

apply(int_summary_female_PI_80[,seq(1, 10, by = 2)], 1, which.min)
apply(int_summary_female_PI_80[,seq(2, 10, by = 2)], 1, which.min)

apply(int_summary_male_PI_80[,seq(1, 10, by = 2)], 1, which.min)
apply(int_summary_male_PI_80[,seq(2, 10, by = 2)], 1, which.min)

#####################################
### 95% nominal coverage probability
#####################################

## EVR

int_summary_female_PI_95_EVR = cbind(interval_fore_national_err_female_PI_95_EVR[,1:2],
                                     interval_fore_national_err_MFTS_F_PI_95_EVR[,1:2],
                                     int_fore_national_err_MLFTS_F_PI_95_EVR[,1:2],
                                     interval_fore_national_EVR_clr_female_PI_95[,1:2])

int_summary_male_PI_95_EVR = cbind(interval_fore_national_err_male_PI_95_EVR[,1:2],
                                   interval_fore_national_err_MFTS_M_PI_95_EVR[,1:2],
                                   int_fore_national_err_MLFTS_M_PI_95_EVR[,1:2],
                                   interval_fore_national_EVR_clr_male_PI_95[,1:2])
colnames(int_summary_female_PI_95_EVR) = colnames(int_summary_male_PI_95_EVR) = rep(c("score", "CPD"), 4)

xtable(rbind(int_summary_female_PI_95_EVR, colMeans(int_summary_female_PI_95_EVR)), digits = 4)
xtable(rbind(int_summary_male_PI_95_EVR, colMeans(int_summary_male_PI_95_EVR)), digits = 4)

## K = 6

int_summary_female_PI_95 = cbind(interval_fore_national_err_female_PI_95[,1:2],
                                 interval_fore_national_err_MFTS_F_PI_95[,1:2],
                                 int_fore_national_err_MLFTS_F_PI_95[,1:2],
                                 interval_fore_national_fixed_clr_female_PI_95[,1:2])

int_summary_male_PI_95 = cbind(interval_fore_national_err_male_PI_95[,1:2],
                                 interval_fore_national_err_MFTS_M_PI_95[,1:2],
                                 int_fore_national_err_MLFTS_M_PI_95[,1:2],
                                 interval_fore_national_fixed_clr_male_PI_95[,1:2])

xtable(rbind(int_summary_female_PI_95, colMeans(int_summary_female_PI_95)), digits = 4)
xtable(rbind(int_summary_male_PI_95, colMeans(int_summary_male_PI_95)), digits = 4)

#############
## which min
#############

# EVR

apply(int_summary_female_PI_95_EVR[,seq(1, 8, by = 2)], 1, which.min)
apply(int_summary_female_PI_95_EVR[,seq(2, 8, by = 2)], 1, which.min)

apply(int_summary_male_PI_95_EVR[,seq(1, 8, by = 2)], 1, which.min)
apply(int_summary_male_PI_95_EVR[,seq(2, 8, by = 2)], 1, which.min)

# K = 6

apply(int_summary_female_PI_95[,seq(1, 8, by = 2)], 1, which.min)
apply(int_summary_female_PI_95[,seq(2, 8, by = 2)], 1, which.min)

apply(int_summary_male_PI_95[,seq(1, 8, by = 2)], 1, which.min)
apply(int_summary_male_PI_95[,seq(2, 8, by = 2)], 1, which.min)
