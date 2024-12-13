#Main analysis

# Main analysis restricted to patients for whom
#this dx is their first ever C-code cancer diagnosis

# 1. Unstaged = EARLY analysis -------------------------------------------------------------------
## Data cleaning 

pacman::p_load(tidyverse, remotes, glmmTMB, ggeffects, merTools, sjPlot, reshape2, 
               ggeffects, performance, here, tidymodels, multilevelmod, Metrics,
               table1, kableExtra, knitr, rmarkdown, NDRSAfunctions)

events <- readRDS(paste0(Sys.getenv("MAIHDA_data"), "/colorectal_av2022.Rds"))

#restrict to each patient's first ever malignant cancer diagnosis
firstcan <- events |> filter(first_cancer == 1) |>
  filter(diagnosisyear < 2020) 


events <- firstcan |>
  mutate(eth = as.factor(case_when(ethnicity ==  "0" ~ "White",
                                   ethnicity ==  "A" ~ "White",
                                   ethnicity ==  "B" ~ "White",
                                   grepl("C", ethnicity) ~ "White", 
                                   ethnicity ==  "D" ~ "Other",
                                   ethnicity ==  "E" ~ "Other",
                                   ethnicity ==  "F" ~ "Other",
                                   ethnicity ==  "G" ~ "Other",
                                   ethnicity ==  "H" ~ "Asian",
                                   ethnicity ==  "J" ~ "Asian",
                                   ethnicity ==  "K" ~ "Asian",
                                   ethnicity ==  "L" ~ "Asian",
                                   ethnicity ==  "M" ~ "Black",
                                   ethnicity ==  "N" ~ "Black",
                                   ethnicity ==  "P" ~ "Black",
                                   ethnicity ==  "R" ~ "Asian",
                                   ethnicity ==  "S" ~ "Other",
                                   ethnicity ==  "8" ~ "Other",
                                   ethnicity ==  "X" ~ NA,
                                   ethnicity ==  "Z" ~ NA,
                                   is.na(ethnicity) ~ NA,
                                   .default = "Other"))) |>
  mutate(gender = as.factor(case_when(gender == "1" ~ "male",
                                      gender == "2" ~ "female",
                                      .default = NA))) |>
  mutate(agegrp = as.factor(ifelse(age < 40, 1, 
                                   ifelse(age>39 & age <80, floor((age/10)-3), 5)))) |>
  mutate(agegrp = factor(agegrp, labels = c("0 to 49","50 to 59",
                                            "60 to 69","70 to 79","80 plus"))) |>
  mutate(stage = case_when(grepl("1|2", stage_best) ~ 0,
                           grepl("3|4", stage_best) ~ 1,
                           grepl("X", stage_best) ~ 1, #X indicates unstaged but stageable
                           .default = NA)) |>
  mutate(reg = as.factor(canalliance_2024_name),
         imd = as.factor(quintile_2019),
         table_stage = factor(stage, levels = c(0,1), labels = c("Early","Advanced"))) |>
  drop_na(stage, gender, agegrp, eth)

## a) Table 1----
label(events$eth) <- "Ethnic group"
label(events$gender) <- "Gender"
label(events$agegrp) <- "5-Year Age Band"
units(events$agegrp) <- "years"
label(events$imd) <- "Deprivation level"
label(events$table_stage) <- "Stage at diagnosis"

length(events$stage_best[events$stage_best == "X"])/nrow(events) * 100
t1 <- table1::table1(~ agegrp + gender + eth + imd | table_stage, data = events)
kable(t1) 
write.csv(t1, paste0(Sys.getenv("MAIHDA_out"), "/Table1_SA2.csv"), row.names = FALSE)

#sensible default levels
events <- within(events, eth <- relevel(eth, ref = "White"))
events <- within(events, agegrp <- relevel(agegrp, ref = 4))

# create the strata
events$stratum <- as.factor(paste(events$gender, events$agegrp, 
                                  events$eth, events$imd, sep = "_"))
#add counts of tumours per stratum
events <- events  |>
  group_by(stratum) |>
  mutate(strataN = n()) |>
  arrange(stratum)

#Fit the two-level logistic regression with no covariates
model1_log_colo <- glmmTMB::glmmTMB(stage ~ (1|stratum), data=events, family=binomial)

#predict the fitted linear predictor (on the probability scale)
events$m2Axbu <- predict(model1_log_colo, type = "response")

#predict the linear predictor for the fixed portion only
events$m2Axb <- predict (model1_log_colo, type = "response", re.form = NA)

#Fit a 2-level regression with covariates ("additive main effects model") 
model2_log_colo <- glmmTMB::glmmTMB(stage ~ eth + agegrp + gender + imd + 
                                      (1|stratum), data=events, family=binomial)
#get the L2 variance with more decimal places than the table gives
##b) L2----
l2_late <- as.data.frame(VarCorr(model2_log_colo)[["cond"]][["stratum"]]) 

#Get the estimates as Odds ratios (and SEs on the odds scale)
tab_model(model2_log_colo, digits.re=8, show.se=T,
          file = paste0(Sys.getenv("MAIHDA_out"), "/Model2_late.doc"),
          show.reflvl = TRUE)

#predict the fitted linear predictor (on the probability scale)
events$m2Bmfit <- predict(model2_log_colo, type = "response")

# Predict the linear predictor for the fixed portion of the model only
events$m2BmF <- predict(model2_log_colo, type="response", re.form=NA)

#predict the SE of the fitted linear predictor
events$m2Bmse <- predict(model2_log_colo, type="response", se.fit=TRUE)$se.fit

#create confidence intervals for the fitted linear predictor
events$m2Bmupr <- events$m2Bmfit + 1.96*events$m2Bmse
events$m2Bmlwr <- events$m2Bmfit - 1.96*events$m2Bmse

#predict the fitted linear predictor, and SE, on the logit 
# scale 
events$m2BmfitL <- predict(model2_log_colo, type="link")
events$m2BmseL <- predict(model2_log_colo, type="link", se.fit=TRUE)$se.fit

# Collapse the data down to a stratum-level dataset 
stratum_level <- 
  events |>
  group_by(gender, eth, imd, agegrp, stratum, strataN, m2Bmfit, m2BmF,m2Bmupr,
           m2Bmlwr, m2BmfitL, m2BmseL) |>
  summarise(stage = mean(as.numeric(stage), na.rm = T))

# convert the outcome from a proportion to a percentage
stratum_level$stage <- stratum_level$stage*100

# Table 2: Calculate stratum-level descriptive statistics

# generate binary indicators for whether each stratum has more than X 
# individuals
stratum_level$n100plus <- ifelse(stratum_level$strataN>=100, 1,0)
stratum_level$n50plus <- ifelse(stratum_level$strataN>=50, 1,0)
stratum_level$n30plus <- ifelse(stratum_level$strataN>=30, 1,0)
stratum_level$n20plus <- ifelse(stratum_level$strataN>=20, 1,0)
stratum_level$n10plus <- ifelse(stratum_level$strataN>=10, 1,0)
stratum_level$nlessthan10 <- ifelse(stratum_level$strataN<10, 1,0)

##c) Table 2 ----

table2 <- data.frame("Sample Size Per Stratum" = c("100 or more","50 or more","30 or more","20 or more","10 or more","Less than 10"),
                     "Number of Strata" = c(sum(stratum_level$n100plus), sum(stratum_level$n50plus), sum(stratum_level$n30plus),
                                            sum(stratum_level$n20plus), sum(stratum_level$n10plus), sum(stratum_level$nlessthan10)),
                     "% of Strata" = c((sum(stratum_level$n100plus)/200)*100, (sum(stratum_level$n50plus)/200)*100, 
                                       (sum(stratum_level$n30plus)/200)*100,(sum(stratum_level$n20plus)/200)*100, (sum(stratum_level$n10plus)/200)*100, (sum(stratum_level$nlessthan10)/200)*100),
                     check.names = FALSE)


kable(table2)
write.csv(table2, file = paste0(Sys.getenv("MAIHDA_out"), "/Table_2_SA2.csv"))

# Get the estimates as Odds ratios (and SEs on the odds scale)
##d) Models----
tab_model(model1_log_colo, model2_log_colo, show.se=T,
          file = paste0(Sys.getenv("MAIHDA_out"), "/Model_1and2_SA2.doc"),
          show.reflvl = TRUE)

# Calculate the area under the receiver operating characteristic (ROC) curve
# for model 2a - based on intercept and stratum random effects
AUC2A <- auc(events$stage, events$m2Axbu)

#for model 2A - based on only the fixed portion of the model
AUC2AF <- auc(events$stage, events$m2Axb)

#for model 2b - based on intercept, main effects, and stratum random effects
AUC2B <- auc(events$stage, events$m2Bmfit)
#for model 2b - based on the fixed portion of the model (main effects)
AUC2BF <- auc(events$stage, events$m2BmF) 

#output the AUC calculations
##e) AUC ----
AUC2A
AUC2AF
AUC2B
AUC2BF

# Panel B - histogram of the observed stratum means
##f) histogram ----
png(paste0(Sys.getenv("MAIHDA_out"), "/Obs_stratummeans_SA2.png"),
    res = 72)
ggplot(stratum_level, aes(x=stage)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=2, 
                 boundary=22) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("Percentage of strata") +
  xlab("Percentage of Advanced Stage Colorectal Cancer Diagnoses") +
  geom_vline(aes(xintercept=mean(stage))) + 
  annotate("text", x=55, y=0.2, size= 4,label=paste0("Grand Mean = ", round(mean(as.numeric(events$stage))*100, 1), "%")) +
  theme_bw() +
  labs(title = "Histogram of Percentage of Advanced Stage Colorectal Cancer Diagnoses") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        title = element_text(size = 12))
dev.off()

# Panel B

# Rank the predicted stratum probabilities
stratum_level <- stratum_level |>
  ungroup() |>
  mutate(rank2=rank(m2Bmfit))

# convert probabilities to percentages
stratum_level$m2Bmfit <- stratum_level$m2Bmfit * 100
stratum_level$m2Bmupr <- stratum_level$m2Bmupr * 100
stratum_level$m2Bmlwr <- stratum_level$m2Bmlwr * 100

# Plot the caterpillar plot of the predicted stratum means
##g) predicted percent----
png(paste0(Sys.getenv("MAIHDA_out"), "/predicted_percentage_SA2.png"),
    res = 72)
stratum_level |>
  filter(nlessthan10 == 0) |>
  ggplot(aes(y=m2Bmfit, x=rank2)) +
  geom_point() +
  geom_pointrange(aes(ymin=m2Bmlwr, ymax=m2Bmupr)) +
  ylab("Predicted Percent Advanced Stage Diagnoses") +
  xlab("Stratum Rank") +
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        title = element_text(size = 12))
dev.off()

##h) both ----

lowest <- stratum_level |> ungroup() |> slice_min(m2Bmfit, n=6) |>
  dplyr::select("Rank" = rank2, "Gender" = gender, "Ethnic Group" = eth, 
                "IMD Quintile" = imd, "Age Group" = agegrp,
                "n" = strataN, "Predicted proportion of late stage diagnoses" = m2Bmfit,
                "Lower 95% CI" = m2Bmlwr, "Upper 95% CI" = m2Bmupr) |>
  mutate("Predicted proportion of late stage diagnoses" = round(`Predicted proportion of late stage diagnoses`, 2))

highest <- top <- stratum_level |> ungroup() |> slice_max(m2Bmfit, n=6) |>
  dplyr::select("Rank" = rank2, "Gender" = gender, "Ethnic Group" = eth, 
                "IMD Quintile" = imd, "Age Group" = agegrp,
                "n" = strataN, "Predicted proportion of late stage diagnoses" = m2Bmfit,
                "Lower 95% CI" = m2Bmlwr, "Upper 95% CI" = m2Bmupr) |>
  mutate("Predicted proportion of late stage diagnoses" = round(`Predicted proportion of late stage diagnoses`, 2))

both <- bind_rows(lowest, highest)
write.csv(both, paste0(Sys.getenv("MAIHDA_out"), "/both_SA2.csv"), row.names = F)
kable(both) 


# from section 7 of https://www.rpubs.com/Hailstone/388481


# caterpillar plot

# get random effects
random_effects <- ranef(model2_log_colo)
random_intercept <- random_effects$cond


# get variances
random_effect_var <- TMB::sdreport(model2_log_colo$obj, getJointPrecision=TRUE)
random_effect_sd <- sqrt(random_effect_var$diag.cov.random)

caterpillar_data <- data.frame(
  "intercepts"=random_intercept$stratum$`(Intercept)`,
  "sd"=random_effect_sd,
  "stratum"=factor(row.names(random_effects$cond$stratum))
)

# calc confidence interval
caterpillar_data$ucl <- caterpillar_data$intercepts + (caterpillar_data$sd * 1.96)
caterpillar_data$lcl <- caterpillar_data$intercepts - (caterpillar_data$sd * 1.96)

# categorise for colour coding in plot

caterpillar_data$category <- ifelse(caterpillar_data$lcl > 0, "High",
                                    ifelse(caterpillar_data$ucl < 0, "Low",
                                           "Average"))

# reorder the ccg names factor
caterpillar_data$stratum <- fct_reorder(caterpillar_data$stratum, caterpillar_data$intercepts)

# add quintiles
caterpillar_data$intercepts_quintile <- fct_rev(factor(ntile(caterpillar_data$intercepts, 5)))

# Difference in Predicted Percent Advanced stage diagnosis due to Interactions
##i) catepillar logit -----
png(paste0(Sys.getenv("MAIHDA_out"), "/Caterpillar_logit_SA2.png"),
    res = 72)
ggplot(caterpillar_data,
       aes(stratum,
           intercepts,
           colour=category)) +
  geom_hline(yintercept=0) +
  geom_point(size=4)  +
  geom_errorbar(aes(ymin=lcl, ymax=ucl)) +
  scale_colour_manual(values=c("grey", "#0571b0","#ca0020")) +
  xlab("Stratum Rank") +
  ylab("Difference in Predicted Percent Advanced Stage Colorectal Cancer Diagnosis due to Interactions") +
  theme(axis.text = element_text(size= 12),
        axis.text.x = element_blank(),
        legend.position = "none",
        axis.title= element_text(size = 12)) 
dev.off()

#significant strata based on logit
sig <- caterpillar_data |> 
  filter((ucl>0 & lcl >0|(ucl<0 & lcl < 0)))
write.csv(sig, paste0(Sys.getenv("MAIHDA_out"), "/significant_strata_logit_SA2.csv"))

stratumsim <- rbind(stratum_level, 
                    stratum_level[rep(1:nrow(stratum_level),999),])

set.seed(354612)

stratumsim$m2Bpsim <- 100*invlogit(stratumsim$m2BmfitL + 
                                     rnorm(nrow(stratumsim), mean=0, sd=stratumsim$m2BmseL))
stratumsim$m2BpAsim <- 100*(stratumsim$m2BmF)
stratumsim$m2BpBsim <- stratumsim$m2Bpsim - stratumsim$m2BpAsim
stratumsim <- stratumsim[order(stratumsim$stratum),]

stratumsim2 <- stratumsim |>
  group_by(stratum) |>
  summarise(mean=mean(m2BpBsim), std=sd(m2BpBsim)) |>
  mutate(rank=rank(mean)) |>
  mutate(hi=(mean + 1.96*std),
         lo=(mean - 1.96*std))
stratumsim2$category <- ifelse(stratumsim2$lo > 0, "High",
                               ifelse(stratumsim2$hi < 0, "Low",
                                      "Average"))

# plot the caterpillar plot of the predicted stratum percentage differences
##j) catepillar probability----
png(paste0(Sys.getenv("MAIHDA_out"), "/Caterpillar_probability_SA2.png"),
    res = 72)
ggplot(stratumsim2, aes(x=rank, y=mean, colour = category)) +
  geom_hline(yintercept=0, color="red", linewidth=1) +
  geom_point(size=4) +
  geom_pointrange(aes(ymin=lo, ymax=hi)) + 
  scale_colour_manual(values=c("grey", "#0571b0","#ca0020")) +
  xlab("Stratum Rank") +
  ylab("Difference in Predicted Percent Advanced Stage Colorectal\nCancer Diagnosis due to Interactions") +
  theme_bw() +
  theme(axis.text =element_text(size = 12),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 10))
dev.off()

sig <- stratumsim2 |> 
  filter((hi>0 & lo >0|(hi<0 & lo < 0)))
write.csv(sig, paste0(Sys.getenv("MAIHDA_out"), "/significant_strata_prob_SA2.csv"))
