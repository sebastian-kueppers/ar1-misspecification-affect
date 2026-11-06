
# required packages
library(tidyverse)
library(dplyr)

# create folder structure
# dir.create("data/raw/")
# dir.create("data/processed/")
# dir.create("inputs/coded/", recursive = T)
# dir.create("results/")

# Preparation of input files ===================================================
# Based on the request summary file, create two separate files for relevant
# ESM variables and person-level variables.

###### @Seb: You don not need to run this code section again.

# prepare the input sheet
inputs <- read.csv("data/raw/summary/data_reports_VRXLEPA17Y_2026-06-11_request_summary_S6LEUis.csv")

# remove certain variables
inputs <- inputs[inputs$assessment_type != "ESM2",]
inputs <- inputs[inputs$assessment_type != "Follow_up",]
inputs <- inputs[inputs$assessment_type != "Follow_up2",]
inputs <- inputs[!grepl(inputs$constructs, pattern = "REGULATION"),]
inputs <- inputs[!grepl(inputs$constructs, pattern = "SITUATION_CONTEXT"),]
inputs <- inputs[!grepl(inputs$constructs, pattern = "COPING"),]
inputs <- inputs[!grepl(inputs$constructs, pattern = "APPRAISAL"),]
inputs <- inputs[!grepl(inputs$constructs, pattern = "OPENNESS"),]
inputs <- inputs[!grepl(inputs$constructs, pattern = "EXTRAVERSION"),]
inputs <- inputs[!grepl(inputs$constructs, pattern = "AGREEABLENESS"),]
inputs <- inputs[!grepl(inputs$constructs, pattern = "CONSCIENTIOUSNESS"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "_RT_"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "change"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "GOAL"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "AFF"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "STRATEGY"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "EVENT"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "EFFORT"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "MOTIVE"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "_TS_"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "GRID"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "survey"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "Survey"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "SURVEY"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "Day"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "day"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "DAY"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "number"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "Date"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "STUDY"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "TMZ"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "Time"),]
inputs <- inputs[!grepl(inputs$variable, pattern = "NO_"),]

inputs <- inputs[!(inputs$variable %in% c(
  "APP_ES", "ABAND_ES", "ANHE_ES", "CONC_ES", "ENER_ES", "IMPUL_ES",
  "PAIN_ES","REL_QUAL_ES","SELFH_ES","SOCEX_ES", "COUPLE_ID_BL","BEEP_ES", "ESM_BREAK_ES",
  "MONTH_ES", "YEAR_ES", "CONNECT_ES", "EXP_COND_ES", "INDEBT_ES",
  "STUDY_PHASE_ES","TRUST_ES","GENDER_BL", "INCPT_ES", "SOCINT_ES", "TSST_TIME_ES", "AMUSE_ES",
  "AWE_ES", "CURIOUS_ES","COMP_ES","CLARITY_ES", "EFFICACY_ES", "CHANGE_BHVR_ES"

)),]


# prepare columns for coding
inputs$scale_min = NA
inputs$scale_max = NA
scales = strsplit(inputs$scale_levels, split = ",")
for(i in 1:nrow(inputs)){
  inputs$scale_min[i] = as.numeric(scales[[i]][1])
  inputs$scale_max[i] = ifelse(length(scales[[i]])!=0, as.numeric(scales[[i]][length(scales[[i]])]), NA)
}
inputs$recode_item <- 0


# auto-choose construct
inputs$construct <- NA
## well-being
inputs$construct <- ifelse(startsWith(inputs$variable, prefix = "BFI"), "NEUR", inputs$construct)
inputs$construct <- ifelse(startsWith(inputs$variable, prefix = "bfi"), "NEUR", inputs$construct)
inputs$construct <- ifelse(startsWith(inputs$variable, prefix = "NEO"), "NEUR", inputs$construct)
inputs$construct <- ifelse(startsWith(inputs$variable, prefix = "TIPI"), "NEUR", inputs$construct)
inputs$construct <- ifelse(startsWith(inputs$variable, prefix = "DASS"), "DEP", inputs$construct)
inputs$construct <- ifelse(startsWith(inputs$variable, prefix = "CESD"), "DEP", inputs$construct)
inputs$construct <- ifelse(startsWith(inputs$variable, prefix = "SWL"), "SWL", inputs$construct)



## momentary affect
inputs$construct <- ifelse(grepl(inputs$constructs, pattern = "NEGATIVE_AFFECT"), "NEG", inputs$construct)
inputs$construct <- ifelse(grepl(inputs$constructs, pattern = "POSITIVE_AFFECT"), "POS", inputs$construct)

## momentary affect - assign specific construct per emotion variable
# POS/NEG assignment from constructs column (kept for btw-level variables)
inputs$construct <- ifelse(grepl(inputs$constructs, pattern = "NEGATIVE_AFFECT"), "NEG", inputs$construct)
inputs$construct <- ifelse(grepl(inputs$constructs, pattern = "POSITIVE_AFFECT"), "POS", inputs$construct)

# assign specific emotion constructs by variable name
# anger
inputs$construct <- ifelse(inputs$variable == "ANG_ES", "ANG", inputs$construct)
# anxiety
inputs$construct <- ifelse(inputs$variable == "ANX_ES", "ANX", inputs$construct)
# calm
inputs$construct <- ifelse(inputs$variable == "CALM_ES", "CALM", inputs$construct)
# cheerful
inputs$construct <- ifelse(inputs$variable == "CHEER_ES", "CHEER", inputs$construct)
# self-esteem / self-confidence (CONF, SELFE are the same construct)
inputs$construct <- ifelse(inputs$variable %in% c("CONF_ES", "SELFE_ES"), "SELFE", inputs$construct)
# contempt
inputs$construct <- ifelse(inputs$variable == "CONTEMPT_ES", "CONTEMPT", inputs$construct)
# content
inputs$construct <- ifelse(inputs$variable == "CONT_ES", "CONT", inputs$construct)
# depressed
inputs$construct <- ifelse(inputs$variable == "DEP_ES", "DEP", inputs$construct)
# disappointed
inputs$construct <- ifelse(inputs$variable == "DISAP_ES", "DISAP", inputs$construct)
# dull
inputs$construct <- ifelse(inputs$variable == "DULL_ES", "DULL", inputs$construct)
# elated
inputs$construct <- ifelse(inputs$variable == "ELATED_ES", "ELATED", inputs$construct)
# embarrassed
inputs$construct <- ifelse(inputs$variable == "EMB_ES", "EMB", inputs$construct)
# energetic
inputs$construct <- ifelse(inputs$variable == "ENERG_ES", "ENERG", inputs$construct)
# enthusiastic
inputs$construct <- ifelse(inputs$variable == "ENTHU_ES", "ENTHU", inputs$construct)
# euphoric
inputs$construct <- ifelse(inputs$variable == "EUPH_ES", "EUPH", inputs$construct)
# excited
inputs$construct <- ifelse(inputs$variable == "EXCIT_ES", "EXCIT", inputs$construct)
# fear
inputs$construct <- ifelse(inputs$variable == "FEAR_ES", "FEAR", inputs$construct)
# grateful
inputs$construct <- ifelse(inputs$variable == "GRATE_ES", "GRATE", inputs$construct)
# guilty (GUILT and GUIL are the same construct)
inputs$construct <- ifelse(inputs$variable %in% c("GUILT_ES", "GUIL_ES"), "GUILT", inputs$construct)
# happy
inputs$construct <- ifelse(inputs$variable == "HAP_ES", "HAP", inputs$construct)
# hopeful
inputs$construct <- ifelse(inputs$variable == "HOPE_ES", "HOPE", inputs$construct)
# irritated
inputs$construct <- ifelse(inputs$variable == "IRR_ES", "IRR", inputs$construct)
# lonely
inputs$construct <- ifelse(inputs$variable == "LONE_ES", "LONE", inputs$construct)
# loving
inputs$construct <- ifelse(inputs$variable == "LOV_ES", "LOV", inputs$construct)
# peaceful
inputs$construct <- ifelse(inputs$variable == "PEACE_ES", "PEACE", inputs$construct)
# proud
inputs$construct <- ifelse(inputs$variable == "PROUD_ES", "PROUD", inputs$construct)
# relieved
inputs$construct <- ifelse(inputs$variable == "RELIEF_ES", "RELIEF", inputs$construct)
# restless
inputs$construct <- ifelse(inputs$variable == "RESTL_ES", "RESTL", inputs$construct)
# relaxed
inputs$construct <- ifelse(inputs$variable == "RLX_ES", "RLX", inputs$construct)
# sad
inputs$construct <- ifelse(inputs$variable == "SAD_ES", "SAD", inputs$construct)
# satisfied
inputs$construct <- ifelse(inputs$variable == "SAT_ES", "SAT", inputs$construct)
# self-conscious
inputs$construct <- ifelse(inputs$variable == "SELFC_ES", "SELFC", inputs$construct)
# ashamed
inputs$construct <- ifelse(inputs$variable == "SHAM_ES", "SHAM", inputs$construct)
# stressed
inputs$construct <- ifelse(inputs$variable == "STR_ES", "STR", inputs$construct)
# tired
inputs$construct <- ifelse(inputs$variable == "TIR_ES", NA, inputs$construct)
# success - not in preregistered item pool
inputs$construct <- ifelse(inputs$variable == "SUCCESS_ES", NA, inputs$construct)
# event valence items - not affect items per se
inputs$construct <- ifelse(inputs$variable %in% c("NEG_EXTENT_ES", "POS_EXTENT_ES"), NA, inputs$construct)

# assign affect_type (POS/NEG) based on construct
inputs$affect_type <- NA
inputs$affect_type <- ifelse(inputs$construct %in% c(
  "ANG", "ANX", "CONTEMPT", "DEP", "DISAP", "DULL", "EMB",
  "FEAR", "GUILT", "IRR", "LONE", "RESTL", "SAD", "SELFC",
  "SHAM", "STR"
), "NEG", inputs$affect_type)
inputs$affect_type <- ifelse(inputs$construct %in% c(
  "CALM", "CHEER", "CONT", "ELATED", "ENERG", "ENTHU", "EUPH",
  "EXCIT", "GRATE", "HAP", "HOPE", "LOV", "PEACE", "PROUD",
  "RELIEF", "RLX", "SAT", "SELFE"
), "POS", inputs$affect_type)


# show remaining variables
inputs$variable[is.na(inputs$construct)]

#
inputs <- inputs %>% arrange(dataset, assessment_type, variable)

# store file to add inputs
writexl::write_xlsx(inputs, path = "inputs/input_file_all.xlsx")

# ESM only
inputsESM <- inputs %>%
  filter(assessment_type == "ESM") %>%
  arrange(dataset, construct, variable)

# MULTICHOICE items are primarily emotion regulation strategy variables, not affect items - remove
inputsESM %>% filter(response_type == "MULTICHOICE") %>% select(dataset_presentation, variable, var_info_eng) %>% as.data.frame() %>% print()
inputsESM <- inputsESM %>% filter(response_type != "MULTICHOICE")

# RADIO item is a participant type variable, not an affect item - remove
inputsESM %>% filter(response_type == "RADIO") %>% select(dataset_presentation, variable, var_info_eng) %>% as.data.frame() %>% print()
inputsESM <- inputsESM %>% filter(response_type != "RADIO")

# SYSTEM items: ANG_ES in Leuven couples 2016 is a valid affect item (anger) miscoded as SYSTEM - keep it
# The others (TRIGGER_ES, WAVE_ES, toofast_count_ES) are system/admin variables - remove
inputsESM %>% filter(response_type == "SYSTEM") %>% select(dataset_presentation, variable, var_info_eng, scale_min, scale_max) %>% as.data.frame() %>% print()
inputsESM <- inputsESM %>%
  filter(!(response_type == "SYSTEM" & variable != "ANG_ES"))

# TIME items are timestamp variables, not affect items - remove
inputsESM %>% filter(response_type == "TIME") %>% select(dataset_presentation, variable, var_info_eng) %>% as.data.frame() %>% print()
inputsESM <- inputsESM %>% filter(response_type != "TIME")

# store file
writexl::write_xlsx(inputsESM, path = "inputs/input_file_ESM.xlsx")

# between-level variables
inputsbtw <- inputs %>%
  filter(assessment_type != "ESM") %>%
  arrange(construct, dataset, variable)

inputsbtw$var_info_eng[inputs$variable == "CESD_1_BL"]


# store file to add inputs
writexl::write_xlsx(inputsbtw, path = "inputs/input_file_BTW.xlsx")




# ==============================================================================

# Processing of timestamps

# @Seb: You can start here, once you have stored the data under 'data/raw/'

# ==============================================================================

# load datasets
files <- Sys.glob("data/raw/*.csv")

# process time stamps =====

df.list <- lapply(1:length(files), function(x){

  # load file
  df <- read.csv(files[x])

  # print number for debugging
  print(df$dataset[1])

  # use only wave 1 data
  if(df$dataset[1] == "Leuven 3-wave longitudinal study"){
    df <- df[df$WAVE_ES==1,]
  }

  # running integer person ID
  df$ID = as.numeric(factor(df$UUID))

  # ensure UUID is stored as character for consistency across datasets
  df$UUID <- as.character(df$UUID)

  # process time stamps =====

  # Some datasets need to be treated differently because the time format is not
  # consistent ...


  ## first, bring date into right format
  if(df$dataset[1] %in% c("Emotion Regulation Effort (2020)",
                          "Social Processes Under Daily Situations (SPUDS)")){


    if(df$dataset[1] == "Emotion Regulation Effort (2020)"){
      df$date <- parse_date(df$Date_Local, format = "%Y/%m/%d")
      df <- df[!is.na(df$date),]
      df <- df[df$Time_Local != "",]
      datetime <- paste0(df$Date_Local, " ", df$Time_Local)
      df$datetime <- parse_date_time(datetime, "%Y/%m/%d %H:%M:%S", exact = TRUE)
    }

    if(df$dataset[1] == "Social Processes Under Daily Situations (SPUDS)"){
      # replace missing entries with scheduled times
      scheduled <- str_split(df$SCHEDULED_TS_ES, pattern = " ")
      scheduled.date <- unlist(lapply(scheduled, function(x){
        date_raw = str_split(x[[1]], pattern = "/")
        date_form = paste0("20",date_raw[[1]][3],"/",
                           str_pad(date_raw[[1]][1], width = 2, pad = 0),"/",
                           str_pad(date_raw[[1]][2], width = 2, pad = 0))
      }))
      scheduled.time <- unlist(lapply(scheduled, function(x){time_raw = paste0(x[[2]],":00")}))
      df$Date_Local <- ifelse(df$Date_Local == "", scheduled.date, df$Date_Local)
      df$Time_Local <- ifelse(df$Time_Local == "", scheduled.time, df$Time_Local)
      df$date = parse_date(df$Date_Local, format = "%Y/%m/%d")
      datetime = paste0(df$Date_Local, " ", df$Time_Local)
      df$datetime = parse_date_time(c(datetime), c("%Y/%m/%d %H:%M:%S"), exact = TRUE)
    }


  } else if(df$dataset[1] == "Goals, Emotions, Moods, and Strategies (GEMS) 2024"){
    # Date_Local is already YYYY-MM-DD, just normalize separators to YYYY/MM/DD
    df$Date_Local <- gsub("-", "/", df$Date_Local)
    df$date <- parse_date(df$Date_Local, format = "%Y/%m/%d")
    df <- df[!is.na(df$date),]
    df <- df[df$Time_Local != "",]
    datetime <- paste0(df$Date_Local, " ", df$Time_Local)
    df$datetime <- parse_date_time(datetime, "%Y/%m/%d %H:%M:%S", exact = TRUE)

  } else {

    if(df$dataset[1] == "Assessing, Understanding, and Regulating Affect (AURA)"){
      # replace missing entries with scheduled times
      scheduled <- str_split(df$SCHEDULED_TS_ES, pattern = " ")
      scheduled.date <- unlist(lapply(scheduled, function(x){
        if(x[[1]] == "") return("")
        date_raw <- str_split(x[[1]], pattern = "/")[[1]]
        paste0("20", str_pad(date_raw[3], 2, pad = "0"), "/",
               str_pad(date_raw[1], 2, pad = "0"), "/",
               str_pad(date_raw[2], 2, pad = "0"))
      }))
      scheduled.time <- unlist(lapply(scheduled, function(x){
        ifelse(length(x) > 1, paste0(x[[2]], ":00"), "")
      }))
      df$Date_Local <- ifelse(df$Date_Local == "", scheduled.date, df$Date_Local)
      df$Time_Local <- ifelse(df$Time_Local == "", scheduled.time, df$Time_Local)
    }

    if(df$dataset[1] == "Assessing, Understanding, and Regulating Affect (AURA)"){
      df$date <- parse_date(df$Date_Local, format = "%Y/%m/%d")
      df <- df[!is.na(df$date),]
      df <- df[df$Time_Local != "",]
      datetime <- paste0(df$Date_Local, " ", df$Time_Local)
      df$datetime <- parse_date_time(datetime, "%Y/%m/%d %H:%M:%S", exact = TRUE)
    } else {
      df$date <- parse_date(df$Date_Local, format = "%d/%m/%Y")
      df <- df[!is.na(df$date),]
      df <- df[df$Time_Local != "",]
      datetime <- paste0(df$Date_Local, " ", df$Time_Local)
      df$datetime <- parse_date_time(datetime, "%d/%m/%Y %H:%M:%S", exact = TRUE)
    }


  }

  ## check for strictly increasing continuous time
  ### First, mark time points with same timestamp
  df <- df %>%
    mutate(row_no = 1:n()) %>%
    group_by(ID, datetime) %>%
    mutate(n_same_time = n())

  # count missings on ESM variables and select rows with less missing entries
  df_exl <- df %>% select(row_no, n_same_time, ID, datetime, ends_with("_ES")) %>%
    filter(n_same_time > 1)

  if(nrow(df_exl)>0){
    df_exl$n_miss <- rowSums(is.na(df_exl))
    df_exl <- df_exl %>%
      arrange(ID, datetime, n_miss) %>%
      group_by(ID, datetime) %>%
      mutate(dupl_row = 1:n()) %>%
      filter(dupl_row != 1)

    # exclude the duplicate rows
    df <- df[!(df$row_no %in% df_exl$row_no),]
  }

  ## now create running beep counter variables
  df = df %>%
    group_by(ID) %>%
    arrange(datetime) %>%
    mutate(
      dayno = as.numeric(date - (min(date)-1)),
      cont_time = as.numeric(difftime(datetime, min(datetime), units = "hours")),
      beep = 1:n()) %>%
    group_by(ID, date) %>%
    arrange(beep) %>%
    mutate(beep_day = 1:n(),
           beep.interval = ifelse(beep_day > 1,
                                  as.numeric(difftime(datetime, lag(datetime), units = "hours")),
                                  NA)) %>%
    ungroup() %>%
    arrange(ID, beep)

  ## for one dataset there are still bogus duplicate scores with 1sec time difference
  if(df$dataset[1] == "Leuven couples study 2014"){
    df <- df %>% filter(beep.interval > 0.001 | is.na(beep.interval))

    df = df %>%
      group_by(ID) %>%
      arrange(datetime) %>%
      mutate(
        dayno = as.numeric(date - (min(date)-1)),
        cont_time = as.numeric(difftime(datetime, min(datetime), units = "hours")),
        beep = 1:n()) %>%
      group_by(ID, date) %>%
      arrange(beep) %>%
      mutate(beep_day = 1:n(),
             beep.interval = ifelse(beep_day > 1,
                                    as.numeric(difftime(datetime, lag(datetime), units = "hours")),
                                    NA)) %>%
      ungroup() %>%
      arrange(ID, beep)

  }

  df

})


# check that continuous time is strictly increasing
for(k in 1:length(df.list)){
print(df.list[[k]] %>%
    filter(!duplicated(.)) %>%
    group_by(ID, cont_time) %>%
    mutate(n_time = n()) %>%
    filter(n_time > 1) %>%
    dim(.))

}

### looks fine now









# ==============================================================================
# Composite Scores of Well-Being
# ==============================================================================

##### ATTENTION: ONLY RUN THIS SECTION ONLY ONCE!

# load coding inputs
# inputsbtw <- readxl::read_xlsx("inputs/coded/input_file_BTW.xlsx")
inputsbtw <- readxl::read_xlsx("inputs/input_file_BTW.xlsx")

# Some datasets have wrong entries for scales, e.g., min = max = 123:
inputsbtw %>% filter(dataset_presentation == "Assessing, Understanding, and Regulating Affect (AURA)", construct == "DEP") %>% slice(1) %>% select(variable, scale_min, scale_max) %>% as.data.frame()

# we extract the min and max vals from the scale_levels columns
parse_scale_levels <- function(x) {
  # if there are letters, it is unparseable
  if (grepl("[A-Za-z]", x)) return(c(NA_real_, NA_real_))
  # if no separator present, treat each digit as a scale point
  if(!grepl("[,. ]", x)){
    nums <- as.numeric(strsplit(x, "")[[1]])
  } else {
    nums <- suppressWarnings(as.numeric(unlist(regmatches(x, gregexpr("-?[0-9]+", x)))))
  }
  nums <- nums[!is.na(nums)]
  if(length(nums) == 0) return(c(NA, NA))
  c(min(nums), max(nums))
}

parsed <- t(sapply(inputsbtw$scale_levels, parse_scale_levels))
inputsbtw$scale_min <- parsed[, 1]
inputsbtw$scale_max <- parsed[, 2]

# ERE: constructs column is empty so all 60 BFI2 items were assigned NEUR.
# Per BFI-2 scoring key, NEUR items are 4,9,14,19,24,29,34,39,44,49,54,59.
# However, items 4,9,24,29,44,49 are absent from the CSV (data export issue).
# Only 6 NEUR items are available: 14,19,34,39,54,59 — none need reverse-scoring.
cat("NEUR items missing from CSV:", paste(paste0("BFI2_", c(4,9,24,29,44,49), "_BL"), collapse = ", "), "\n")

# which BFI items are for NEUR?
neuro_items_all <- paste0("BFI2_", c(4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59), "_BL")


inputsbtw <- inputsbtw %>%
  filter(!(dataset_presentation == "Emotion Regulation Effort (2020)" &
             grepl("BFI2", variable) &
             !(variable %in% neuro_items_all))) %>%
  mutate(
    construct = ifelse(
      dataset_presentation == "Emotion Regulation Effort (2020)" &
        variable %in% neuro_items_all,
      "NEUR", construct
    )
  )

# Leuven couples study has the actual item response text in scale
inputsbtw %>%
  filter(dataset_presentation == "Leuven couples study 2016", construct == "DEP") %>%
  select(variable, scale_levels) %>% head(3) %>% as.data.frame()

# apply min/max scale manually
inputsbtw <- inputsbtw %>%
  mutate(
    scale_min = ifelse(dataset_presentation == "Leuven couples study 2016" & construct == "DEP", 0, scale_min),
    scale_max = ifelse(dataset_presentation == "Leuven couples study 2016" & construct == "DEP", 3, scale_max)
  )

# SWLS_8_BL in Leuven 2012 uses a 1-4 scale ("All in all, would you say you:")
# which is a different instrument than the standard SWLS (1-7). Dropping it.
inputsbtw %>% filter(dataset_presentation == "Leuven emotions in daily life 2012", variable == "SWLS_8_BL") %>% select(variable, var_info_eng, scale_levels) %>% as.data.frame()
inputsbtw <- inputsbtw %>%
  filter(!(dataset_presentation == "Leuven emotions in daily life 2012" & variable == "SWLS_8_BL"))

# BFI_31_BL in GEMS 2024 is an attention check item ("Please select 'Disagree strongly'")
# with zero variance - it must be excluded from the NEUR composite
inputsbtw %>% filter(dataset_presentation == "Goals, Emotions, Moods, and Strategies (GEMS) 2024", variable == "BFI_31_BL") %>% select(variable, var_info_eng) %>% as.data.frame()
inputsbtw <- inputsbtw %>%
  filter(!(dataset_presentation == "Goals, Emotions, Moods, and Strategies (GEMS) 2024" & variable == "BFI_31_BL"))


# store omega results
omegas <- data.frame()

# reverse-keyed items
reverse_items <- c(
  paste0("BFI2_", c(4, 9, 24, 29, 44, 49), "_BL"),
  "TIPI_9_BL",
  paste0("NEO_",  c(1, 16, 31, 46), "_BL"),
  paste0("CESD_", c(4, 8, 12, 16), "_BL")
)

for(k in 1:length(df.list)){
  # get study
  study = unique(df.list[[k]]$dataset)
  print(paste0(k, ": ", study))
  inpstudy = inputsbtw %>% filter(dataset_presentation == study) %>% filter(!is.na(construct))
  # get outcomes
  scales = unique(inpstudy$construct)
  for(i in 1:length(scales)){
    print(scales[i])
    # get items
    items = inpstudy$variable[inpstudy$construct == scales[i]]
    # drop items with no column in this dataset (e.g. 6 BFI2 NEUR items absent from ERE export)
    absent = setdiff(items, names(df.list[[k]]))
    if(length(absent) > 0)
      warning(sprintf("%s / %s: dropped items absent from data: %s",
                      study, scales[i], paste(absent, collapse = ", ")))
    items = items[items %in% names(df.list[[k]])]

    # reverse-keyed items present in THIS study & construct
    items_rec = items[items %in% reverse_items]
    # for some studies responses are given as character string
    if((study == "Leuven couples study 2014" | study == "Leuven couples study 2016") & scales[i] == "DEP"){
      for(j in items){
        df.list[[k]][,j][df.list[[k]][,j] == "zelden of nooit (minder dan 1 dag)"] <- "0"
        df.list[[k]][,j][df.list[[k]][,j] == "soms of weinig (1 tot 2 dagen)"] <- "1"
        df.list[[k]][,j][df.list[[k]][,j] == "regelmatig (3 tot 4 dagen)"] <- "2"
        df.list[[k]][,j][df.list[[k]][,j] == "meestal of altijd (5 tot 7 dagen)"] <- "3"
        df.list[[k]][,j] <- as.numeric(unlist(df.list[[k]][,j]))
      }
    }
    # reverse-code flagged items, each on its OWN response scale
    if(length(items_rec) > 0){
      for(v in items_rec){
        s_min <- inpstudy$scale_min[inpstudy$variable == v]
        s_max <- inpstudy$scale_max[inpstudy$variable == v]
        if(length(s_min) != 1 || is.na(s_min) || is.na(s_max)){
          warning(sprintf("Skipping reverse of %s in %s: missing scale endpoints", v, study))
          next
        }
        df.list[[k]][[v]] <- s_min + s_max - df.list[[k]][[v]]
      }
    }
    # calculate omega
    sub = df.list[[k]] %>% select_at(c("UUID",items)) %>% distinct() %>% select_at(c(items))
    omeg <- psych::omega(sub, nfactors = 1)
    omegas <- rbind(omegas,
                    data.frame(
                      "study" = study,
                      "outcome" = scales[i],
                      "N items" = length(items),
                      "alpha" = omeg$alpha,
                      "omega" = omeg$omega.tot
                    ))
    # calculate composite scores
    df.list[[k]][,scales[i]] <- rowMeans(df.list[[k]][,items], na.rm = T)
  }
}


# save results
save(omegas, file = "results/WB_reliaties.Rdata")



# ==============================================================================

# Rescaling of ESM items

# ==============================================================================

# load coding inputs
# inputsESM <- readxl::read_xlsx("inputs/coded/input_file_ESM.xlsx")
inputsESM <- readxl::read_xlsx("inputs/input_file_ESM.xlsx")

for(k in 1:length(df.list)){

  # get study
  study = unique(df.list[[k]]$dataset)
  print(paste0(k, ": ", study))
  inpstudy = inputsESM %>% filter(dataset_presentation == study) %>% filter(!is.na(construct))

  # get ESM items by study
  items = inpstudy$variable
  s_min = unique(inpstudy$scale_min)
  s_max = unique(inpstudy$scale_max)

  # range of -2.5 to 2.5
  for(j in items){
    df.list[[k]][,j] <- (df.list[[k]][,j] - s_min) / (s_max-s_min) * 5 - 2.5
  }
}



# ==============================================================================

### Apply Exclusion Criteria

# ==============================================================================


excl <- data.frame()F
excl.df <- list()

for(k in 1:length(df.list)){

  # get study
  study = unique(df.list[[k]]$dataset)
  inpstudy = inputsESM %>% filter(dataset_presentation == study) %>% filter(!is.na(construct))

  # get ESM items by study
  items = inpstudy$variable

  # exclusion criteria:
  backup <- df.list[[k]]

  # 1. less than 25 surveys completed =========================================
  df.list[[k]][,"complete"] = ifelse(!is.na(rowMeans(df.list[[k]][,items])),1,0)

  # get cases
  excl_1 <- df.list[[k]] %>% group_by(ID) %>%
    summarise(n_completed = sum(complete)) %>% filter(n_completed <= 25)

  # apply criterion
  df.list[[k]] <- df.list[[k]] %>% filter(!(ID %in% excl_1$ID))

  # 2. IVAR = 0 on any ESM indicator(s) =======================================
  excl_2 <- df.list[[k]] %>%
    select_at(c("ID", items)) %>%
    pivot_longer(cols = 2:ncol(.)) %>%
    group_by(ID, name) %>%
    summarise(IVAR = var(value, na.rm = T)) %>%
    filter(IVAR == 0) %>%
    select(ID) %>% distinct()

  # apply criterion
  df.list[[k]] <- df.list[[k]] %>% filter(!(ID %in% excl_2$ID))

  # 3. less than three response categories used ================================
  excl_3 <- df.list[[k]] %>%
    select_at(c("ID", items)) %>%
    pivot_longer(cols = 2:ncol(.)) %>%
    group_by(ID, name) %>%
    summarise(unique_cat = length(unique(value[!is.na(value)]))) %>%
    filter(unique_cat < 3) %>%
    select(ID) %>% distinct()

  # apply criterion
  df.list[[k]] <- df.list[[k]] %>% filter(!(ID %in% excl_3$ID))

  # store excluded cases
  excl.df[[k]] <- backup %>% filter(ID %in% c(excl_1$ID, excl_2$ID, excl_3$ID))

  # get summary
  excl = rbind(
    excl,
    data.frame(
      "study"      = study,
      "excl_TP25"  = ifelse(nrow(excl_1) == 0, 0, nrow(excl_1)),
      "excl_IVAR0" = ifelse(nrow(excl_2) == 0, 0, nrow(excl_2)),
      "excl_2cat"  = ifelse(nrow(excl_3) == 0, 0, nrow(excl_3))
    )
  )
}

save(excl, file = "results/exclusion_summary.Rdata")
save(excl.df, file = "results/excluded_subjects_data.Rdata")

# get sample sizes
n_final <- sapply(df.list, function(d) n_distinct(d$ID))
n_excl  <- sapply(excl.df, function(d) n_distinct(d$ID))
n_start <- n_final + n_excl

sum(n_start) # 2390
sum(n_final) # 2108


# get n due exlcusions
colSums(excl[, c("excl_TP25", "excl_IVAR0", "excl_2cat")])
# excl_TP25 excl_IVAR0  excl_2cat
# 96         92         94

# ==============================================================================
#
# ESM Preprocessing
#
# ==============================================================================

# within-person center each ESM item
for(k in 1:length(df.list)){
  study <- unique(df.list[[k]]$dataset)
  inpstudy <- inputsESM %>% filter(dataset_presentation == study) %>% filter(!is.na(construct))
  items <- inpstudy$variable

  df.list[[k]] <- df.list[[k]] %>%
    group_by(ID) %>%
    arrange(beep) %>%
    mutate(across(
      all_of(items),
      ~ . - mean(., na.rm = TRUE),
      .names = "{.col}_c"
    )) %>%
    ungroup() %>%
    arrange(ID, beep)
}

df.list <- lapply(df.list, function(df){

  df %>%
    group_by(ID) %>%
    arrange(beep) %>%
    group_modify(function(person_df, key){

      result <- list()

      # insert overnight lags
      for(i in 1:nrow(person_df)){
        if(i > 1 && person_df$dayno[i] != person_df$dayno[i - 1]){
          na_row <- person_df[i, ]
          # set everything to NA except ID, UUID, dataset
          na_row[, setdiff(names(na_row), c("ID", "UUID", "dataset"))] <- NA
          result[[length(result) + 1]] <- na_row
        }
        result[[length(result) + 1]] <- person_df[i, ]
      }

      bind_rows(result)
    }) %>%
    ungroup()
})


# save preprocessed data
save(df.list, file = "data/processed/data_processed.Rdata")







