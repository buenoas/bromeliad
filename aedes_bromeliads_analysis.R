#####################
### DATA ANALYSIS ###
#####################

# Required package
library(glmmTMB)
library(ggplot2)
library(gridExtra)


########################
### DATA PREPARATION ###
########################

# Import data
dataset = read.table("https://raw.githubusercontent.com/buenoas/bromeliad/main/mosquito_data.txt", header = TRUE)

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


##########################
### DATA VISUALIZATION ###
##########################

#################################
### 1. AE. AEGYPTI OCCURRENCE ###
#################################

# Format model estimates

beta.glmm.occ = round(beta_occ, 2)
z.glmm.occ = round(z_occ, 2)

p.glmm.occ = ifelse(
  p_occ < 0.001,
  "< 0.001",
  paste0("= ", sprintf("%.3f", p_occ)))


# Organize the data

df = data.frame(
  location = dataset$location,
  week = dataset$week,
  trap = dataset$trap,
  presence = dataset$presence)

df = df[order(df$trap, df$week, df$presence),]
df = df[order(df$location, decreasing = TRUE),]

df$location = factor(df$location, levels = unique(df$location))
df$week = factor(df$week, levels = unique(df$week))

# Change presence values in bromeliads for visualization only

df[which(df$trap == "bromeliad" & df$presence == 1), "presence"] = 2


# Occurrence plot

g_occ =
  ggplot(data = df,
         aes(x = week, y = location,
             fill = as.factor(presence))) +
  
  geom_tile(colour = "white", linewidth = 1) +
  
  facet_wrap(~ trap,
             labeller = labeller(
               trap = c(
                 ovitrap = "Ovitrap",
                 bromeliad = "Bromeliad"))) +
  
  labs(x = "Week", y = "Site",
       subtitle = bquote(
         "Binomial GLMM:" ~
           italic(beta) * " = " * .(beta.glmm.occ) * ", " *
           italic(z) * " = " * .(z.glmm.occ) * ", " *
           italic(p) * " " * .(p.glmm.occ))) +
  
  scale_fill_manual(
    values = c(
      "0" = "grey90",
      "1" = "#414142",
      "2" = "#008E00")) +
  
  scale_y_discrete(
    labels = c("J", "I", "H", "G", "F",
               "E", "D", "C", "B", "A")) +
  
  theme_minimal(base_size = 16) +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 12),
        axis.text = element_text(size = 14, colour = "black"),
        panel.grid = element_blank(),
        axis.title = element_text(colour = "black"))

g_occ


ggsave(g_occ, filename = "graph_occurrence.png", dpi = 600,
       width = 16, height = 16, units = "cm")


################################
### 2. AE. AEGYPTI ABUNDANCE ###
################################

# Format model estimates

beta.glmm.abund = round(beta_abund, 2)
z.glmm.abund = round(z_abund, 2)

p.glmm.abund = ifelse(
  p_abund < 0.001,
  "< 0.001",
  paste0("= ", sprintf("%.3f", p_abund)))


# Abundance plot

g_ab =
  ggplot(data = dataset,
         aes(x = trap, y = Aedes_aegypti_count,
             colour = trap, fill = trap)) +
  
  labs(x = "Breeding site type",
       y = expression("Number of emerged mosquitoes"),
       subtitle = bquote(
         "Negative binomial GLMM:" ~
           italic(beta) * " = " * .(beta.glmm.abund) * ", " *
           italic(z) * " = " * .(z.glmm.abund) * ", " *
           italic(p) * " " * .(p.glmm.abund))) + 
  
  geom_boxplot(outlier.shape = NA, fill = NA,
               width = 0.1, position = position_nudge(x = -0.3)) +
  
  geom_point(position = position_dodge2(0.3),
             shape = 21, size = 2.5, alpha = 0.5) +
  
  stat_summary(
    fun = mean, geom = "point", shape = 23,
    size = 2.5, position = position_nudge(x = -0.3)) +
  
  scale_colour_manual(
    values = c("ovitrap" = "#414142",
               "bromeliad" = "#008E00")) +
  
  scale_fill_manual(
    values = c("ovitrap" = "#414142",
               "bromeliad" = "#008E00")) +
  
  scale_x_discrete(
    labels = c(bromeliad = "Bromeliad",
               ovitrap = "Ovitrap")) +
  
  theme_classic(base_size = 16) +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.line = element_line(linewidth = 1/2),
        axis.ticks = element_line(linewidth = 1/2))

g_ab


ggsave(g_ab, filename = "graph_abundance.png", dpi = 600,
       width = 16, height = 16, units = "cm")


###################################
### 3. ENVIRONMENTAL CONDITIONS ###
###################################

# Temperature

g_temp =
  ggplot(data = dataset,
         aes(x = trap, y = temperature,
             colour = trap, fill = trap)) +
  
  labs(x = "Breeding site type",
       y = "Temperature (°C)",
       tag = "(a)",
       subtitle = bquote(
         "Gaussian GLMM:" ~
           italic(beta) * " = " * .(round(beta_temp, 2)) * ", " *
           italic(z) * " = " * .(round(z_temp, 2)) * ", " *
           italic(p) * " " * .(ifelse(
             p_temp < 0.001,
             "< 0.001",
             paste0("= ", sprintf("%.3f", p_temp)))))) + 
  
  geom_boxplot(outlier.shape = NA, fill = NA,
               width = 0.1, position = position_nudge(x = -0.3)) +
  
  geom_point(position = position_dodge2(0.3),
             shape = 21, size = 2.5, alpha = 0.5) +
  
  stat_summary(
    fun = mean, geom = "point", shape = 23,
    size = 2.5, position = position_nudge(x = -0.3)) +
  
  scale_colour_manual(
    values = c("ovitrap" = "#414142",
               "bromeliad" = "#008E00")) +
  
  scale_fill_manual(
    values = c("ovitrap" = "#414142",
               "bromeliad" = "#008E00")) +
  
  scale_x_discrete(
    labels = c(bromeliad = "Bromeliad",
               ovitrap = "Ovitrap")) +
  
  scale_y_continuous(
    limits = c(20, 35)) +
  
  theme_classic(base_size = 16) +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.line = element_line(linewidth = 1/2),
        axis.ticks = element_line(linewidth = 1/2))

g_temp

ggsave(g_temp, filename = "graph_temperature.png", dpi = 600,
       width = 16, height = 16, units = "cm")


# Total dissolved solids (TDS)

g_TDS =
  ggplot(data = dataset,
         aes(x = trap, y = TDS,
             colour = trap, fill = trap)) +
  
  labs(x = "Breeding site type",
       y = "Total dissolved solids (ppm)",
       tag = "(b)",
       subtitle = bquote(
         "Gamma GLMM:" ~
           italic(beta) * " = " * .(round(beta_TDS, 2)) * ", " *
           italic(z) * " = " * .(round(z_TDS, 2)) * ", " *
           italic(p) * " " * .(ifelse(
             p_TDS < 0.001,
             "< 0.001",
             paste0("= ", sprintf("%.3f", p_TDS)))))) + 
  
  geom_boxplot(outlier.shape = NA, fill = NA,
               width = 0.1, position = position_nudge(x = -0.3)) +
  
  geom_point(position = position_dodge2(0.3),
             shape = 21, size = 2.5, alpha = 0.5) +
  
  stat_summary(
    fun = mean, geom = "point", shape = 23,
    size = 2.5, position = position_nudge(x = -0.3)) +
  
  scale_colour_manual(
    values = c("ovitrap" = "#414142",
               "bromeliad" = "#008E00")) +
  
  scale_fill_manual(
    values = c("ovitrap" = "#414142",
               "bromeliad" = "#008E00")) +
  
  scale_x_discrete(
    labels = c(bromeliad = "Bromeliad",
               ovitrap = "Ovitrap")) +
  
  scale_y_continuous(
    limits = c(0, 400)) +
  
  theme_classic(base_size = 16) +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.line = element_line(linewidth = 1/2),
        axis.ticks = element_line(linewidth = 1/2))

g_TDS

ggsave(g_TDS, filename = "graph_TDS.png", dpi = 600,
       width = 16, height = 16, units = "cm")


# pH

g_pH =
  ggplot(data = dataset,
         aes(x = trap, y = pH,
             colour = trap, fill = trap)) +
  
  labs(x = "Breeding site type",
       y = "pH",
       tag = "(c)",
       subtitle = bquote(
         "Gaussian GLMM:" ~
           italic(beta) * " = " * .(round(beta_pH, 2)) * ", " *
           italic(z) * " = " * .(round(z_pH, 2)) * ", " *
           italic(p) * " " * .(ifelse(
             p_pH < 0.001,
             "< 0.001",
             paste0("= ", sprintf("%.3f", p_pH)))))) + 
  
  geom_boxplot(outlier.shape = NA, fill = NA,
               width = 0.1, position = position_nudge(x = -0.3)) +
  
  geom_point(position = position_dodge2(0.3),
             shape = 21, size = 2.5, alpha = 0.5) +
  
  stat_summary(
    fun = mean, geom = "point", shape = 23,
    size = 2.5, position = position_nudge(x = -0.3)) +
  
  scale_colour_manual(
    values = c("ovitrap" = "#414142",
               "bromeliad" = "#008E00")) +
  
  scale_fill_manual(
    values = c("ovitrap" = "#414142",
               "bromeliad" = "#008E00")) +
  
  scale_x_discrete(
    labels = c(bromeliad = "Bromeliad",
               ovitrap = "Ovitrap")) +
  
  theme_classic(base_size = 16) +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.line = element_line(linewidth = 1/2),
        axis.ticks = element_line(linewidth = 1/2))

g_pH


ggsave(g_pH, filename = "graph_pH.png", dpi = 600,
       width = 16, height = 16, units = "cm")


# Salva os tres graficos
ggsave(grid.arrange(g_temp, g_TDS, g_pH, ncol = 3),
       filename = "environmental_factors.png",
       dpi = 600, width = 16*3, height = 16, units = "cm")
