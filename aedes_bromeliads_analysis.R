#############################
### DATA ANALYSIS
#############################

# Required package
library(glmmTMB)


#############################
### DATA PREPARATION
#############################

# Import data
dataset = read.table("https://raw.githubusercontent.com/buenoas/bromeliad/main/mosquito_data.txt",
                     header = TRUE)

# Breeding-site type
# Ovitraps are the reference category.
dataset$trap = factor(
  dataset$trap,
  levels = c("ovitrap", "bromeliad"))

# Sampling location and week
dataset$location = factor(dataset$location)
dataset$week = factor(dataset$week)


#################################
### 1. AE. AEGYPTI OCCURRENCE ###
#################################

# Convert abundance to a binary occurrence variable.
# A value of 1 indicates the presence of at least one
# emerged Ae. aegypti individual in the sample.

dataset$presence = as.integer(dataset$Aedes_aegypti_count > 0)

# Binomial GLMM for Ae. aegypti occurrence.
#
# Trap type is included as a fixed effect.
# Location and week are included as random intercepts.

m_occ = glmmTMB(
  presence ~
    trap +
    (1 | location) +
    (1 | week),
  family = binomial,
  data = dataset)

summary(m_occ)


################################
### 2. AE. AEGYPTI ABUNDANCE ###
################################

# Negative binomial GLMM for the number of emerged
# Ae. aegypti individuals.
#
# Trap type is included as a fixed effect.
# Location and week are included as random intercepts.
# Dispersion is allowed to vary between trap types.

m_abund = glmmTMB(
  Aedes_aegypti_count ~
    trap +
    (1 | location) +
    (1 | week),
  dispformula = ~ trap,
  family = nbinom2,
  data = dataset)

summary(m_abund)


################################
### 3. WATER PHYSICOCHEMICAL ###
###    PARAMETERS            ###
################################

# Temperature
# Gaussian GLMM with trap type as a fixed effect
# and location and week as random intercepts.

m_temp = glmmTMB(
  temperature ~
    trap +
    (1 | location) +
    (1 | week),
  family = gaussian(),
  data = dataset)

summary(m_temp)


# Total dissolved solids (TDS)
# Gamma GLMM with a log link.
# Trap type is included as a fixed effect.
# Location and week are included as random intercepts.

m_TDS = glmmTMB(
  TDS ~
    trap +
    (1 | location) +
    (1 | week),
  family = Gamma(link = "log"),
  data = dataset)

summary(m_TDS)


# pH
# Gaussian GLMM with trap type as a fixed effect
# and location and week as random intercepts.

m_pH = glmmTMB(
  pH ~
    trap +
    (1 | location) +
    (1 | week),
  family = gaussian(),
  data = dataset)

summary(m_pH)


###########################
### 4. EFFECT ESTIMATES ###
###########################

# Extract the effect of bromeliads relative to ovitraps.
# Ovitraps are the reference category.


# Ae. aegypti occurrence

coefs_occ = summary(m_occ)$coefficients$cond

beta_occ = coefs_occ["trapbromeliad", "Estimate"]
se_occ = coefs_occ["trapbromeliad", "Std. Error"]
z_occ = coefs_occ["trapbromeliad", "z value"]
p_occ = coefs_occ["trapbromeliad", "Pr(>|z|)"]

odds_ratio_occ = exp(beta_occ)


# Ae. aegypti abundance

coefs_abund = summary(m_abund)$coefficients$cond

beta_abund = coefs_abund["trapbromeliad", "Estimate"]
se_abund = coefs_abund["trapbromeliad", "Std. Error"]
z_abund = coefs_abund["trapbromeliad", "z value"]
p_abund = coefs_abund["trapbromeliad", "Pr(>|z|)"]

abund_ratio = exp(beta_abund)


# Temperature

coefs_temp = summary(m_temp)$coefficients$cond

beta_temp = coefs_temp["trapbromeliad", "Estimate"]
se_temp = coefs_temp["trapbromeliad", "Std. Error"]
z_temp = coefs_temp["trapbromeliad", "z value"]
p_temp = coefs_temp["trapbromeliad", "Pr(>|z|)"]


# Total dissolved solids (TDS)

coefs_TDS = summary(m_TDS)$coefficients$cond

beta_TDS = coefs_TDS["trapbromeliad", "Estimate"]
se_TDS = coefs_TDS["trapbromeliad", "Std. Error"]
z_TDS = coefs_TDS["trapbromeliad", "z value"]
p_TDS = coefs_TDS["trapbromeliad", "Pr(>|z|)"]


# pH

coefs_pH = summary(m_pH)$coefficients$cond

beta_pH = coefs_pH["trapbromeliad", "Estimate"]
se_pH = coefs_pH["trapbromeliad", "Std. Error"]
z_pH = coefs_pH["trapbromeliad", "z value"]
p_pH = coefs_pH["trapbromeliad", "Pr(>|z|)"]


###############################
### 5. CONFIDENCE INTERVALS ###
###############################

# 95% confidence intervals

confint(m_occ)
confint(m_abund)
confint(m_temp)
confint(m_TDS)
confint(m_pH)


#######################
### 6. EFFECT SIZES ###
#######################

# Odds ratio for Ae. aegypti occurrence:
# values > 1 indicate higher odds of occurrence
# in bromeliads relative to ovitraps.

odds_ratio_occ


# Abundance ratio for Ae. aegypti:
# values < 1 indicate lower abundance
# in bromeliads relative to ovitraps.

abund_ratio
