#data extract


pacman::p_load(tidyverse, remotes, glmmTMB, ggeffects, merTools, sjPlot, reshape2, 
               ggeffects, performance, here, tidymodels, multilevelmod, Metrics,
               table1, kableExtra, knitr, rmarkdown, NDRSAfunctions)

# Extract -----------------------------------------

source("ns_stage.R")

#switch here to recreate the dataset from CAS if needed but load from local if not.
if (!file.exists(paste0(Sys.getenv("MAIHDA_data"), "/colorectal_av2022.Rds"))){
  
  con <- createConnection(username = Sys.getenv("analysis_username"))

  SELECTQUERY <- "with t1 as (
select
	t.tumourid
	, t.age
	, t.diagnosisyear
	, t.gender
	, s.ndrs_main
	, s.ndrs_detailed
	, t.site_icd10
	, g.canalliance_2024_code
	, g.canalliance_2024_name
	, t.stage_best
	, t.stage_best_system
	, t.ethnicity
	, imd.imd2015_quintiles as quintile_2015
	, imd.imd2019_quintiles as quintile_2019
    
from av2022.at_tumour_england@casref01 t
    left outer join av2022.at_site_england@casref01 s on t.tumourid = s.tumourid
    left outer join av2022.at_transformation_haem_england@casref01 th on t.tumourid = th.transformed_id
    left outer join av2022.at_geography_england@casref01 g on t.tumourid = g.tumourid
    left outer join LSOA21_IMD_15_19_OHID@casref01 imd on g.lsoa21_code = imd.lsoa21cd
where 
	t.ctry_code = 'E' -- England residents using country code 
	and t.statusofregistration = 'F' -- Finalised cases 
	and t.age between 0 and 200 -- 
	and t.gender in ('1','2') -- Known gender
	-- Years of interest
	and t.diagnosisyear between 2013 and 2022
	-- Site restrictions

	and s.ndrs_main in ('Colorectal')
	    --- Exclude records in at_transformation_haem_england, as these are (not-primary) haematological neoplasms but they are not excluded above as ndrs_main is null 
	and th.transformed_id is null 

)

, first_tumour as (
select patientid, tumourid
, rank() over(partition by patientid order by diagnosisdatebest, tumourid) as ranked_all_cancer
from av2022.at_tumour_england@casref01
where ctry_code = 'E'
    and statusofregistration = 'F'
    and dedup_flag = 1
    and age between 0 and 200
    and gender in ('1','2')
    and site_icd10  like 'C%'
    order by patientid, diagnosisdatebest, tumourid
   )
   
, first_colorectal as (
select t.patientid, t.tumourid
, rank() over(partition by t.patientid order by t.diagnosisdatebest, t.tumourid) as ranked_colorectal
from av2022.at_tumour_england@casref01 t
left outer join av2022.at_site_england@casref01 s on t.tumourid = s.tumourid
where ctry_code = 'E'
    and statusofregistration = 'F'
    and dedup_flag = 1
    and age between 0 and 200
    and gender in ('1','2')
    and (site_icd10_3char in ('C18', 'C19', 'C20') 
      or ndrs_main = 'Colorectal')
    order by patientid, diagnosisdatebest, tumourid
   )   

select t1.*
, case when ft.ranked_all_cancer = 1 then 1
      else 0 
      end as first_cancer
, case when ct.ranked_colorectal = 1 then 1 
      else 0
      end as first_colorectal
from t1
left outer join first_tumour ft on t1.tumourid = ft.tumourid
left outer join first_colorectal ct on t1.tumourid = ct.tumourid

  "
events <- dbGetQueryOracle(con, SELECTQUERY,  rowlimit = NA)
names(events) <- stringr::str_to_lower(names(events))
events <- events |>
  ns_stage(
    input_stage_best_field = "stage_best",
    input_stage_best_system_field = "stage_best_system",
    ndrs_detailed_field = "ndrs_detailed"
  ) 

saveRDS(events, 
        file = paste0(Sys.getenv("MAIHDA_data"), "/colorectal_av2022.Rds"))

} else{
  
  events <- readRDS(paste0(Sys.getenv("MAIHDA_data"), "/colorectal_av2022.Rds"))}

