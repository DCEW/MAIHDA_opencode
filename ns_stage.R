# Libraries and dependencies ----
library(dplyr)
library(stringr)

#' # Function to assign early and advanced stage
#' @param data A dataframe with stage_best and stage_best_system variables
#' @param input_stage_best_field A named field with the existing stage
#' @param input_stage_best_system_field A named field with the staging system
#' @param main_field_name A desired name for the main stage field
#' @param detailed_field_name A desired name for the detailed stage field
#' @param ndrs_detailed_field A named field with the NDRS detailed field
#' @return A data frame with the additional columns
#' @examples

ns_stage <- function(data, input_stage_best_field, input_stage_best_system_field, ndrs_detailed_field, main_field_name = "stage_summary", detailed_field_name = "stage_detailed"){
  data <- data |> dplyr::mutate(
    !!detailed_field_name := dplyr::case_when(
      (!!sym(ndrs_detailed_field) %in% c("BCC", "cSCC")) ~ "Do not include in stage analysis - NMSC",
      (!!sym(input_stage_best_system_field) == "Binet" & stringr::str_sub(!!sym(input_stage_best_field), 1, 1) %in% c("A", "B")) ~ "Staged - other early",
      (!!sym(input_stage_best_system_field) == "Binet" & stringr::str_sub(!!sym(input_stage_best_field), 1, 1) %in% c("C")) ~ "Staged - other advanced",
      (!!sym(input_stage_best_system_field) == "Chang" & stringr::str_sub(!!sym(input_stage_best_field), 1, 2) %in% c("M0")) ~ "Staged - other early",
      (!!sym(input_stage_best_system_field) == "Chang" & stringr::str_sub(!!sym(input_stage_best_field), 1, 2) %in% c("M1", "M2", "M3", "M4")) ~ "Staged - other advanced",
      (!!sym(input_stage_best_system_field) == "INRG" & stringr::str_sub(!!sym(input_stage_best_field), 1, 1) %in% c("L")) ~ "Staged - other early",
      (!!sym(input_stage_best_system_field) == "INRG" & stringr::str_sub(!!sym(input_stage_best_field), 1, 1) %in% c("M")) ~ "Staged - other advanced",
      (!!sym(input_stage_best_system_field) == "ISS" & stringr::str_sub(!!sym(input_stage_best_field), 1, 1) %in% c("1", "2")) ~ "Staged - other early",
      (!!sym(input_stage_best_system_field) == "ISS" & stringr::str_sub(!!sym(input_stage_best_field), 1, 1) %in% c("3")) ~ "Staged - other advanced",
      (!!sym(input_stage_best_system_field) == "Wilms" & stringr::str_sub(!!sym(input_stage_best_field), 1, 1) %in% c("5")) ~ "4",
      (stringr::str_sub(!!sym(input_stage_best_field), 1, 1) %in% c("1", "2", "3", "4")) ~ stringr::str_sub(!!sym(input_stage_best_field), 1, 1),
      (!!sym(input_stage_best_field) %in% c("X")) ~ "Missing",
      (!!sym(input_stage_best_field) %in% c("U") | is.na(!!sym(input_stage_best_field))) ~ "Unstageable",
      .default = "Error"
    ),
    !!main_field_name := dplyr::case_when(
      (!!sym(detailed_field_name) %in% c("1", "2", "Staged - other early")) ~ "Stage 1 & 2",
      (!!sym(detailed_field_name) %in% c("3", "4", "Staged - other advanced")) ~ "Stage 3 & 4",
      .default = !!sym(detailed_field_name)
    )
  )
  return(data)
}


