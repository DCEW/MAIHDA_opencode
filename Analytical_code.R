#Analytical code#

#Load relevant packages
pacman::p_load(tidyverse, remotes, glmmTMB, ggeffects, merTools, sjPlot, reshape2, 
               ggeffects, performance, here, tidymodels, multilevelmod, Metrics,
               table1, kableExtra, NDRSAfunctions)

## 1. EXTRACT CANCER STAGE DATA -----
con <- createConnection(sid = Sys.getenv("sid"), port = 1525, 
                        username = Sys.getenv("analysis_username"))

#stage at diagnosis at National level
SELECTQUERY <- "with stage2021 as
(
    select  /*+ USE_HASH(avt spu)*/
    avt.tumourid
		, case when to_char(avt.morph_icd10_o2) = 'X' then '8000' else avt.morph_icd10_o2 end as morph_fix

        , CASE
             when spu.stage_pi_detail = 'Y' and avt.stage_best_system = 'Binet' and substr(avt.stage_best,1,1) in ('A','B') then 'Stage 1, 2'
             when spu.stage_pi_detail = 'Y' and avt.stage_best_system = 'Binet' and substr(avt.stage_best,1,1) in ('C') then 'Stage 3, 4'
			 when spu.stage_pi_detail = 'Y' and avt.stage_best_system = 'INRG' and substr(avt.stage_best,1,1) in ('L') then 'Stage 1, 2'
             when spu.stage_pi_detail = 'Y' and avt.stage_best_system = 'INRG' and substr(avt.stage_best,1,1) in ('M') then 'Stage 3, 4'
             when spu.stage_pi_detail = 'Y' and avt.stage_best_system = 'Wilms' and substr(avt.stage_best,1,1) in ('5') then 'Stage 3, 4'
             when spu.stage_pi_detail = 'Y' and avt.stage_best_system = 'Chang' and substr(avt.stage_best,1,2) in ('M0') then 'Stage 1, 2'
             when spu.stage_pi_detail = 'Y' and avt.stage_best_system = 'Chang' and substr(avt.stage_best,1,2) in ('M1', 'M2', 'M3', 'M4') then 'Stage 3, 4'
             when spu.stage_pi_detail = 'Y' and substr(avt.stage_best,1,1) in ('1','2') then 'Stage 1, 2'
             when spu.stage_pi_detail = 'Y' and (substr(avt.stage_best,1,1) in ('3','4'))  then  'Stage 3, 4'
             else 'Error'
             end AS early_advanced_stage

    from av2021.at_tumour_england@casref01 avt
    left join analysisncr.at_tumour_england@cas2404 spu
        on avt.tumourid = spu.tumourid

    where avt.diagnosisyear between 2013 and 2021
        and substr(avt.SITE_ICD10R4_O2_FROM2013,1,1)= 'C' and substr(avt.SITE_ICD10R4_O2_FROM2013,1,3)<> 'C44'
        and (not (substr(nvl(avt.stage_best,-1),1,1)='0'))
		and spu.stage_pi_detail is not null
)

, first_tumour as (
select patientid, tumourid
, rank() over(partition by patientid order by diagnosisdatebest, tumourid) as ranked
 from av2021.at_tumour_england@casref01
 where ctry_code = 'E'
    and statusofregistration = 'F'
    and dedup_flag = 1
    and age between 0 and 200
    and gender in ('1','2')
    order by patientid, diagnosisdatebest, tumourid
   )
    
    
select /*+ USE_HASH(geo nspl) USE_HASH(nspl lookup) USE_HASH(geo imd15) USE_HASH(geo imd19)*/
distinct avt.tumourid
, avt.age
, avt.diagnosisyear
, avt.gender
, avt.SITE_ICD10R4_O2_FROM2013
, morph_fix
, lookup.cal23cd
, lookup.cal23nm
, s21.stage_pi
, avt.stage_best
, s21.early_advanced_stage
, avt.ethnicity
, imd15.imd15_quintile_lsoas  as quintile_2015
, imd19.imd19_quintile_lsoas as quintile_2019

, case when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C34' then 'Lung'
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C15' or avt.SITE_ICD10R4_O2_FROM2013 = 'C160' then 'Oesophagus'
        when avt.SITE_ICD10R4_O2_FROM2013 in ('C161', 'C162', 'C163', 'C164', 'C165', 'C166', 'C167', 'C168', 'C169') then 'Stomach'     --Do not report on sites with <70% completeness for the latest 3 years of data available. For 2019-2021: cervix, stomach, and thyroid
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C18' then 'Colon'
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) in ('C19', 'C20') then 'Rectum'
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C25' then 'Pancreas'
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C43' then 'Melanoma of skin'
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C50' and gender = 2 then 'Breast'
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C53' and gender = 2 then 'Cervix'     --Do not report on sites with <70% completeness for the latest 3 years of data available. For 2019-2021: cervix, stomach, and thyroid
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) in ('C54', 'C55') and gender = 2 then 'Uterus'
        when (substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) in ('C56') or avt.SITE_ICD10R4_O2_FROM2013 in ('C570', 'C571', 'C573', 'C572', 'C574', 'C575', 'C576')
            or (substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) in ('C48') and (morph_fix not in (8693, 8800, 8801, 8802, 8803, 8804, 8805, 8806, 8810, 8963,
                    8990, 8991, 9040, 9041, 9042, 9043, 9044, 9490, 9500)
                    and (morph_fix not between 8811 and 8921)
                    and (morph_fix not between 9120 and 9373)
                    and (morph_fix not between 9530 and 9582)
                    )))
            and gender = 2 then 'Ovary'
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C61' and gender = 1 then 'Prostate'
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C62' and gender = 1 then 'Testis'
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C64' then 'Kidney'
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C67' then 'Bladder'
        when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C81' then 'Hodgkin lymphoma'
       when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C73' then 'Thyroid'     --Do not report on sites with <70% completeness for the latest 3 years of data available. For 2019-2021: cervix, stomach, and thyroid
		when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) = 'C32' or avt.SITE_ICD10R4_O2_FROM2013 in ('C101') then 'Larynx'

	when diagnosisyear <= 2018 and (substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) in ('C01', 'C09') or avt.SITE_ICD10R4_O2_FROM2013 in ('C051', 'C052', 'C100', 'C102', 'C103', 'C104', 'C105', 'C106', 'C107', 'C108', 'C109')) then 'Oropharynx'
    when diagnosisyear > 2018 and (substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) in ('C01', 'C09') or avt.SITE_ICD10R4_O2_FROM2013 in ('C051', 'C052', 'C100', 'C102', 'C103', 'C104', 'C105', 'C106', 'C107', 'C108', 'C109', 'C024')) then 'Oropharynx'

	when diagnosisyear <= 2018 and (substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) in ('C02', 'C03', 'C04', 'C06') or avt.SITE_ICD10R4_O2_FROM2013 in ('C050', 'C003', 'C004', 'C005')) then 'Oral cavity'
    when diagnosisyear > 2018 and (substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) in ('C03', 'C04', 'C06') or avt.SITE_ICD10R4_O2_FROM2013 in ('C050', 'C003', 'C004', 'C005', 'C020', 'C021', 'C022', 'C023', 'C028', 'C029')) then 'Oral cavity'

	when substr(avt.SITE_ICD10R4_O2_FROM2013,1,3) in ('C83', 'C84', 'C85', 'C82','C86') then 'NHL'
		else 'Other'
		end as stage_cancer

, case when diagnosisyear >= 2019 and avt.SITE_ICD10R4_O2_3char_FROM2013 in ('C00', 'C01', 'C02', 'C03', 'C04', 'C05', 'C06', 'C07', 'C08', 'C09', 'C10', 'C11', 'C12', 'C13', 'C14', 'C32') then 1
		when diagnosisyear >= 2019 and avt.SITE_ICD10R4_O2_FROM2013 in ('C300', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316', 'C317', 'C318', 'C319') then 1
		when diagnosisyear >= 2018 and (avt.SITE_ICD10R4_O2_3char_FROM2013 not in ('C00', 'C01', 'C02', 'C03', 'C04', 'C05', 'C06', 'C07', 'C08', 'C09', 'C10', 'C11', 'C12', 'C13', 'C14', 'C32')
			or avt.SITE_ICD10R4_O2_FROM2013 not in ('C300', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316', 'C317', 'C318', 'C319')) then 1
		else 0
		end as tnm_flag

from analysisncr.AT_TUMOUR_ENGLAND@cas2309.encore.nhs.uk avt
inner join first_tumour ft on avt.tumourid = ft.tumourid
inner join STAGE2021 S21
    ON AVT.TUMOURID = S21.TUMOURID
inner join analysisncr.at_geography_england@cas2309.encore.nhs.uk geo
    on avt.tumourid=geo.tumourid
inner join alt_nspl@casref01 nspl
    on nspl.pcd = geo.postcode
inner join lsoa21_sicbl_icb_ca_la@casref01 lookup
    on lookup.lsoa21cd = nspl.lsoa21
inner join imd.imd2015_equal_lsoas@casref01  imd15
    on geo.lsoa11_code=imd15.lsoa11_code
inner join imd.imd2019_equal_lsoas@casref01  imd19
    on geo.lsoa11_code=imd19.lsoa11_code

where
    avt.ctry_code = 'E'
    and statusofregistration = 'F'
    and dedup_flag = 1
    and age between 0 and 200
    and gender in ('1','2')
    and early_advanced_stage<> 'Error'
"

events <- dbGetQueryOracle(con, SELECTQUERY,  rowlimit = NA)

events |>
  group_by(EARLY_ADVANCED_STAGE, DIAGNOSISYEAR) |>
  summarise(freq = n())

saveRDS(events, file = "Data/Stage_data_MAIHDA.Rds")

## 2. DATA CLEANING ----

events <- events |>
  mutate(eth = as.factor(case_when(ETHNICITY ==  "0" ~ "White",
                                   ETHNICITY ==  "A" ~ "White",
                                   ETHNICITY ==  "B" ~ "White",
                                   grepl("C", ETHNICITY) ~ "White", 
                                   ETHNICITY ==  "D" ~ "Other",
                                   ETHNICITY ==  "E" ~ "Other",
                                   ETHNICITY ==  "F" ~ "Other",
                                   ETHNICITY ==  "G" ~ "Other",
                                   ETHNICITY ==  "H" ~ "Asian",
                                   ETHNICITY ==  "J" ~ "Asian",
                                   ETHNICITY ==  "K" ~ "Asian",
                                   ETHNICITY ==  "L" ~ "Asian",
                                   ETHNICITY ==  "M" ~ "Black",
                                   ETHNICITY ==  "N" ~ "Black",
                                   ETHNICITY ==  "P" ~ "Black",
                                   ETHNICITY ==  "R" ~ "Asian",
                                   ETHNICITY ==  "S" ~ "Other",
                                   ETHNICITY ==  "8" ~ "Other",
                                   ETHNICITY ==  "X" ~ "Other",
                                   ETHNICITY ==  "Z" ~ "Other",
                                   .default = "Other"))) |>
  mutate(gender = as.factor(case_when(GENDER == "1" ~ "male",
                                      GENDER == "2" ~ "female",
                                      .default = NA))) |>
  mutate(agegrp = as.factor(ifelse(AGE < 40, 1, 
                                   ifelse(AGE>39 & AGE <80, floor((AGE/10)-3), 
                                          5)))) |>
  mutate(agegrp = factor(agegrp, labels = c("0 to 49","50 to 59",
                                            "60 to 69","70 to 79","80 plus"))) |>
  mutate(stage_ch = EARLY_ADVANCED_STAGE,
         stage = case_when(stage_ch == "Stage 1, 2" ~ 0,
                           stage_ch == "Stage 3, 4" ~ 1,
                           .default = NA)) |>
  mutate(reg = as.factor(CAL23NM),
         imd = as.factor(QUINTILE_2019),
         imd2 = case_when(imd == "1 - most deprived"| imd == "2" ~ 1, 
                          .default = 0)) |> 
  drop_na(stage, gender, agegrp, eth) |>
  filter(DIAGNOSISYEAR < 2020) 


#restrict to colorectal cancer only
colo <- events |> filter(SITE_ICD10R4_O2_FROM2013 %in% c("C18","C19","C20"))

label(colo$eth) <- "Ethnic group"
label(colo$gender) <- "Gender"
label(colo$agegrp) <- "5-Year Age Band"
units(colo$agegrp) <- "years"
label(colo$imd) <- "Deprivation level"
label(colo$stage_ch) <- "Stage at diagnosis"

t1 <- table1::table1(~ agegrp + gender + eth + imd | stage_ch, data = colo )
kable(t1)

#make sure the referent categories are the most populous 
colo <- within(colo, eth <- relevel(eth, ref = "White"))
colo <- within(colo, agegrp <- relevel(agegrp, ref = 4))

# create the strata
colo$stratum <- as.factor(paste(colo$gender, colo$agegrp, 
                                colo$eth, colo$imd, sep = "_"))

colo <- colo  |>
  group_by(stratum) |>
  mutate(strataN = n()) |>
  arrange(stratum)

## 3. MODELLING ----
model1_log_colo <- glmmTMB::glmmTMB(stage ~ (1|stratum), data=colo, 
                                    family=binomial)

#predict the fitted linear predictor (on the probability scale)
colo$m2Axbu <- predict(model1_log_colo, type = "response")

#predict the linear predictor for the fixed portion only
colo$m2Axb <- predict (model1_log_colo, type = "response", re.form = NA)

#Fit a 2-level regression with covariates ("additive main effects model") 
model2_log_colo <- glmmTMB::glmmTMB(stage ~ eth + agegrp + gender + imd + 
                                      (1|stratum), data=colo, family=binomial)

#Get the estimates as Odds ratios (and SEs on the odds scale)
tab_model(model2_log_colo, digits.re=8, show.se=T,
          file = paste0(Sys.getenv("MAIHDA_out"), "/Model2.doc"),
          show.reflvl = TRUE)

#predict the fitted linear predictor (on the probability scale)
colo$m2Bmfit <- predict(model2_log_colo, type = "response")

# Predict the linear predictor for the fixed portion of the model only
colo$m2BmF <- predict(model2_log_colo, type="response", re.form=NA)

#predict the SE of the fitted linear predictor
colo$m2Bmse <- predict(model2_log_colo, type="response", se.fit=TRUE)$se.fit

#create confidence intervals for the fitted linear predictor
colo$m2Bmupr <- colo$m2Bmfit + 1.96*colo$m2Bmse
colo$m2Bmlwr <- colo$m2Bmfit - 1.96*colo$m2Bmse

#predict the fitted linear predictor, and SE, on the logit 
# scale 
colo$m2BmfitL <- predict(model2_log_colo, type="link")
colo$m2BmseL <- predict(model2_log_colo, type="link", se.fit=TRUE)$se.fit

# Collapse the data down to a stratum-level dataset 
stratum_level <- 
  colo |>
  group_by(gender, eth, imd, agegrp, stratum, strataN, m2Bmfit, m2BmF,m2Bmupr,
           m2Bmlwr, m2BmfitL, m2BmseL) |>
  summarise(stage = mean(stage, na.rm = T))

# convert the outcome from a proportion to a percentage
stratum_level$stage <- stratum_level$stage*100

# generate binary indicators for whether each stratum has more than X individuals
stratum_level$n100plus <- ifelse(stratum_level$strataN>=100, 1,0)
stratum_level$n50plus <- ifelse(stratum_level$strataN>=50, 1,0)
stratum_level$n30plus <- ifelse(stratum_level$strataN>=30, 1,0)
stratum_level$n20plus <- ifelse(stratum_level$strataN>=20, 1,0)
stratum_level$n10plus <- ifelse(stratum_level$strataN>=10, 1,0)
stratum_level$nlessthan10 <- ifelse(stratum_level$strataN<10, 1,0)

table2 <- 
  data.frame("Sample Size Per Stratum" = 
               c("100 or more","50 or more","30 or more","20 or more",
                 "10 or more","Less than 10"),
                     "Number of Strata" = 
               c(sum(stratum_level$n100plus), sum(stratum_level$n50plus), 
                 sum(stratum_level$n30plus), sum(stratum_level$n20plus), 
                 sum(stratum_level$n10plus), sum(stratum_level$nlessthan10)),
                     "% of Strata" = c((sum(stratum_level$n100plus)/200)*100, 
                                       (sum(stratum_level$n50plus)/200)*100, 
                                       (sum(stratum_level$n30plus)/200)*100,
                                       (sum(stratum_level$n20plus)/200)*100, 
                                       (sum(stratum_level$n10plus)/200)*100, 
                                       (sum(stratum_level$nlessthan10)/200)*100),
                     check.names = FALSE)

kable(table2)

# Get the estimates as Odds ratios (and SEs on the odds scale)
tab_model(model1_log_colo, model2_log_colo, show.se=T,
          file = paste0(Sys.getenv("MAIHDA_out"), "/Model1_2.doc"),
          show.reflvl = TRUE)

# Partitioning Coefficients (VPC)
tab_model(model1_log_colo, model2_log_colo, p.style="stars",
          file = paste0(Sys.getenv("MAIHDA_out"), "/Models_VPC.doc"),
          show.reflvl = TRUE)

# Calculate the area under the receiver operating characteristic (ROC) curve
# for model 2a - based on intercept and stratum random effects
AUC2A <- auc(colo$stage, colo$m2Axbu)

#for model 2A - based on only the fixed portion of the model
AUC2AF <- auc(colo$stage, colo$m2Axb)

#for model 2b - based on intercept, main effects, and stratum random effects
AUC2B <- auc(colo$stage, colo$m2Bmfit)
#for model 2b - based on the fixed portion of the model (main effects)
AUC2BF <- auc(colo$stage, colo$m2BmF) 

figure1 <- ggplot(stratum_level, aes(x=stage)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth=2, 
                 boundary=22) + 
  scale_y_continuous(labels=scales::percent) +
  ylab("Percentage of strata") +
  xlab("Percentage of late stage diagnoses") +
  geom_vline(aes(xintercept=mean(stratum_level$stage))) + #this is the grand mean
  annotate("text", x=45, y=0.2, label=paste0("Grand Mean = ", 
                                             round(mean(colo$stage)*100, 1), "%")) +
  theme_bw() +
  labs(title = "Histogram of percentage of late stage diagnoses")

# Rank the predicted stratum probabilities
stratum_level <- stratum_level |>
  ungroup() |>
  mutate(rank2=rank(m2Bmfit))

# convert probabilities to percentages
stratum_level$m2Bmfit <- stratum_level$m2Bmfit * 100
stratum_level$m2Bmupr <- stratum_level$m2Bmupr * 100
stratum_level$m2Bmlwr <- stratum_level$m2Bmlwr * 100

figure2 <- ggplot(stratum_level, aes(y=m2Bmfit, x=rank2)) +
  geom_point() +
  geom_pointrange(aes(ymin=m2Bmlwr, ymax=m2Bmupr)) +
  ylab("Predicted Percent Late Stage Diagnoses") +
  xlab("Stratum Rank") +
  theme_bw()

#obtain the highest and lowest predicted late stage
lowest <- stratum_level |> ungroup() |> slice_min(m2Bmfit, n=6) |>
  dplyr::select("Rank" = rank2, "Gender" = gender, "Ethnic Group" = eth, 
                "IMD Quintile" = imd, "Age Group" = agegrp,
                "n" = strataN, "Predicted proportion of late stage diagnoses" = 
                  m2Bmfit,
                "Lower 95% CI" = m2Bmlwr, "Upper 95% CI" = m2Bmupr) |>
  mutate("Predicted proportion of late stage diagnoses" = 
           round(`Predicted proportion of late stage diagnoses`, 2))

highest <- top <- stratum_level |> ungroup() |> slice_max(m2Bmfit, n=6) |>
  dplyr::select("Rank" = rank2, "Gender" = gender, "Ethnic Group" = eth, 
                "IMD Quintile" = imd, "Age Group" = agegrp,
                "n" = strataN, "Predicted proportion of late stage diagnoses" = 
                  m2Bmfit,
                "Lower 95% CI" = m2Bmlwr, "Upper 95% CI" = m2Bmupr) |>
  mutate("Predicted proportion of late stage diagnoses" 
         = round(`Predicted proportion of late stage diagnoses`, 2))

both <- bind_rows(lowest, highest)
kable(both) 