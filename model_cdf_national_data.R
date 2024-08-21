##############################
# install and load R packages
##############################

install.packages(c("ftsa", "LaplacesDemon", "flexmix", "psych", "easyCODA", "doMC", "MortalityLaws", "DescTools", "xtable"))

require(ftsa)
require(LaplacesDemon)
require(flexmix)
require(psych)
require(easyCODA)
require(doMC)
require(MortalityLaws)
require(DescTools)
require(xtable)

# change the working directory

setwd("/Users/hanlinshang/Dropbox/Todos/HDFTS_CDF/data/national_data")

# read data set

years = 1975:2022
n_year = length(years)
ages = 0:110
n_age = length(ages)

# age distribution of deaths

female_JPN_qx = t(matrix(read.table("JPN_female_lt.txt", header = TRUE)$qx, n_age, n_year))
male_JPN_qx   = t(matrix(read.table("JPN_male_lt.txt",   header = TRUE)$qx, n_age, n_year))
total_JPN_qx  = t(matrix(read.table("JPN_total_lt.txt",  header = TRUE)$qx, n_age, n_year))

# age-specific life expectancy

female_JPN_ex = t(matrix(read.table("JPN_female_lt.txt", header = TRUE)$ex, n_age, n_year))
male_JPN_ex   = t(matrix(read.table("JPN_male_lt.txt",   header = TRUE)$ex, n_age, n_year))
total_JPN_ex  = t(matrix(read.table("JPN_total_lt.txt",  header = TRUE)$ex, n_age, n_year))

############
# rescaling
############

female_JPN_pop = male_JPN_pop = total_JPN_pop = matrix(NA, nrow(female_JPN_qx), ncol(female_JPN_qx))
for(ij in 1:nrow(female_JPN_qx))
{
    # set radix
    start_pop_female = start_pop_male = start_pop_total = 1
    for(ik in 1:ncol(female_JPN_qx))
    {
        female_JPN_pop[ij,ik] = female_JPN_qx[ij,ik] * start_pop_female
        start_pop_female = start_pop_female - female_JPN_pop[ij,ik]

        male_JPN_pop[ij,ik] = male_JPN_qx[ij,ik] * start_pop_male
        start_pop_male = start_pop_male - male_JPN_pop[ij,ik]

        total_JPN_pop[ij,ik] = total_JPN_qx[ij,ik] * start_pop_total
        start_pop_total = start_pop_total - total_JPN_pop[ij,ik]
    }
}

# save figures of the rainbow plots

savefig("Fig_1a", width = 12, height = 10, toplines = 0.5, type = "png")
plot(fts(ages, t(female_JPN_pop)), xlab = "Age", ylab = "dx", main = "Japan: female data (1975-2022)")
dev.off()

savefig("Fig_1b", width = 12, height = 10, toplines = 0.5, type = "png")
plot(fts(ages, t(male_JPN_pop)), xlab = "Age", ylab = "", main = "Japan: male data (1975-2022)")
dev.off()

####################
# Gini coefficients
####################

Gini_female_JPN_pop = Gini_male_JPN_pop = vector("numeric", nrow(female_JPN_pop))
for(ik in 1:nrow(female_JPN_pop))
{
    Gini_female_JPN_pop[ik] = Gini(female_JPN_pop[ik,])
    Gini_male_JPN_pop[ik] = Gini(male_JPN_pop[ik,])
}

savefig("Fig_2", width = 12, height = 10, toplines = 1, type = "png")
plot(years, Gini_female_JPN_pop, xlab = "Year", ylab = "Gini coefficient", type = "l",
     main = "Japanese data (1975-2022)", ylim = c(0.65, 0.74))
lines(years, Gini_male_JPN_pop,   xlab = "Year", ylab = "", type = "l",
      ylim = c(0.65, 0.74), lty = 2, col = 2)
legend("topleft", c("Female", "Male"), lty = c(1, 2), col = c(1, 2), cex = 0.8)
dev.off()

#############################
# From PDF to CDF via cumsum
#############################

female_JPN_pop_cumsum = male_JPN_pop_cumsum = matrix(NA, nrow(female_JPN_qx), ncol(female_JPN_qx))
for(ij in 1:nrow(female_JPN_qx))
{
    female_JPN_pop_cumsum[ij,] = cumsum(female_JPN_pop[ij,])
    male_JPN_pop_cumsum[ij,]   = cumsum(male_JPN_pop[ij,])
    print(ij); rm(ij)
}

savefig("Fig_3a", width = 12, height = 10, toplines = 0.5, type = "png")
plot(fts(ages, t(female_JPN_pop_cumsum)), xlab = "Age", ylab = "CDF", main = "Japan: female data (1975-2022)")
dev.off()

savefig("Fig_3b", width = 12, height = 10, toplines = 0.5, type = "png")
plot(fts(ages, t(male_JPN_pop_cumsum)), xlab = "Age", ylab = "", main = "Japan: male data (1975-2022)")
dev.off()

##########################################################################
# From CDF to PDF via diff (check if the reconstruction == original data)
##########################################################################

female_JPN_pop_recon = male_JPN_pop_recon = matrix(NA, nrow(female_JPN_qx), ncol(female_JPN_qx))
for(ij in 1:nrow(female_JPN_qx))
{
    female_JPN_pop_recon[ij,] = c(female_JPN_pop_cumsum[ij,1], diff(female_JPN_pop_cumsum[ij,]))
    male_JPN_pop_recon[ij,]   = c(male_JPN_pop_cumsum[ij,1],   diff(male_JPN_pop_cumsum[ij,]))
    print(ij); rm(ij)
}

all(round(female_JPN_pop, 7) == round(female_JPN_pop_recon, 7)) # TRUE
all(round(male_JPN_pop, 7)   == round(male_JPN_pop_recon, 7))   # TRUE

# model a time series of CDF via logistic transformation
# (note: we leave out the last column as it is always 1)

female_JPN_pop_cumsum_logit = male_JPN_pop_cumsum_logit = matrix(NA, nrow(female_JPN_qx), (ncol(female_JPN_qx) - 1))
for(ij in 1:nrow(female_JPN_qx))
{
    female_JPN_pop_cumsum_logit[ij,] = logit(female_JPN_pop_cumsum[ij, 1:110])
    male_JPN_pop_cumsum_logit[ij,]   = logit(male_JPN_pop_cumsum[ij, 1:110])
    rm(ij)
}
colnames(female_JPN_pop_cumsum_logit) = colnames(male_JPN_pop_cumsum_logit) = ages[2:111]
rownames(female_JPN_pop_cumsum_logit) = rownames(male_JPN_pop_cumsum_logit) = years

savefig("Fig_4a", width = 12, height = 10, toplines = 0.5, type = "png")
plot(fts(ages[1:110], t(female_JPN_pop_cumsum_logit)), xlab = "Age", ylab = "Logit transformation of CDF",
         main = "Japan: female data (1975-2022)")
dev.off()

savefig("Fig_4b", width = 12, height = 10, toplines = 0.5, type = "png")
plot(fts(ages[1:110], t(male_JPN_pop_cumsum_logit)), xlab = "Age", ylab = "",
     main = "Japan: female data (1975-2022)")
dev.off()

################################################################
# produce 20-steps-ahead forecasts (add the last columns of 1s)
################################################################

female_JPN_pop_cumsum_logit_fore = forecast(ftsm(fts(ages[1:110], t(female_JPN_pop_cumsum_logit))), h = 20, method = "ets")
female_JPN_pop_cumsum_logit_fore_add = rbind(invlogit(female_JPN_pop_cumsum_logit_fore$mean$y), rep(1, 20))

male_JPN_pop_cumsum_logit_fore = forecast(ftsm(fts(ages[1:110], t(male_JPN_pop_cumsum_logit))), h = 20, method = "ets")
male_JPN_pop_cumsum_logit_fore_add = rbind(invlogit(male_JPN_pop_cumsum_logit_fore$mean$y), rep(1, 20))

savefig("Fig_5a", width = 12, height = 10, toplines = 0.5, type = "png")
plot(fts(ages, t(female_JPN_pop_cumsum)), xlab = "Age", ylab = "CDF", main = "Japan: female data",
     colorchoice = "gray")
lines(fts(ages, female_JPN_pop_cumsum_logit_fore_add), xlab = "Age")
dev.off()

savefig("Fig_5b", width = 12, height = 10, toplines = 0.5, type = "png")
plot(fts(ages, t(male_JPN_pop_cumsum)), xlab = "Age", ylab = "", main = "Japan: male data",
     colorchoice = "gray")
lines(fts(ages, male_JPN_pop_cumsum_logit_fore_add), xlab = "Age")
dev.off()

#####################################################
# transform back from CDF forecasts to PDF forecasts
#####################################################

female_JPN_pop_cumsum_logit_fore_add_diff = male_JPN_pop_cumsum_logit_fore_add_diff = matrix(NA, ncol(female_JPN_qx), 20)
for(ij in 1:20)
{
    female_JPN_pop_cumsum_logit_fore_add_diff[,ij] = c(female_JPN_pop_cumsum_logit_fore_add[1,ij], diff(female_JPN_pop_cumsum_logit_fore_add[,ij]))
    male_JPN_pop_cumsum_logit_fore_add_diff[,ij]   = c(male_JPN_pop_cumsum_logit_fore_add[1,ij],   diff(male_JPN_pop_cumsum_logit_fore_add[,ij]))
    print(ij); rm(ij)
}

#####################################
# plot historical data and forecasts
#####################################

savefig("Fig_6a", width = 12, height = 10, toplines = 0.5, type = "png")
plot(fts(ages, t(female_JPN_pop)), main = "Japan: female data", xlab = "Age", ylab = "Life-table death counts", colorchoice = "gray")
lines(fts(ages, female_JPN_pop_cumsum_logit_fore_add_diff), xlab = "Age", ylab = "Life-table death counts")
dev.off()

savefig("Fig_6b", width = 12, height = 10, toplines = 0.5, type = "png")
plot(fts(ages, t(male_JPN_pop)), main = "Japan: male data", xlab = "Age", ylab = "", colorchoice = "gray")
lines(fts(ages, male_JPN_pop_cumsum_logit_fore_add_diff), xlab = "Age", ylab = "")
dev.off()
