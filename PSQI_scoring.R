# Pittsburgh Sleep Quality Index (PSQI) - Long Form Scoring
# It has 21 items ranging from sleep quality, experience, and amount.  
# The PSQI can be classified as continuous and can also be cutoff 0-4 is Good sleep quality; 5-21 is Poor sleep quality. 

# -----------------------------------------------------------------------------------------
# COMPONENT AND GLOBAL SCORE CALCULATION
# -----------------------------------------------------------------------------------------

# Import data
secondary = read.csv("Z:/Analysis/Kye/PSQI Analysis/PSQI_raw.csv")

library(hms)
library(dplyr)
library(lubridate)

# Handle missing values
secondary[secondary == ""] <- NA
secondary[secondary == "na"] <- NA

# Set factor variables
secondary$subject_number <- as.factor(secondary$subject_number)
secondary$session_number <- as.factor(secondary$session_number)
secondary$ip <- as.factor(secondary$ip)
secondary$arm <- as.factor(secondary$arm)
secondary$preorpost <- as.factor(secondary$preorpost)

# Recode times
secondary$component4_2_2 <- parse_hm(secondary$component4_2_2)
secondary$component4_2_1 <- parse_hm(secondary$component4_2_1) 

# Continuous numeric columns
secondary$component2_3 <- as.numeric(secondary$component2_3)
secondary$component4_1 <- as.numeric(secondary$component4_1)

cols_to_adjust <- c(
  "component2_2and3and4", "component5", "component5_v2", "component5_v3",
  "component5_v4", "component5_v5", "component5_v6", "component5_v7",
  "component5_v8", "component5_v9", "component1", "component6", 
  "component7", "component7_9"
)

# Subtract 1 from all selected columns (survey uses 1-4 Likert scale but responses are scored 0-3)
secondary[cols_to_adjust] <- secondary[cols_to_adjust] - 1

# Component 1: SUBJECTIVE SLEEP QUALITY
  # "component1" column already holds the final component 1 score

# Component 2: SLEEP LATENCY
  # Sleep latency in minutes (component2_3) allocated 0-3 score (≤15=0, 16-20=1, 31-60=2, >60=3)
  # Add this score to Q5a (component2_2and3and4), giving a score from 0-6
  # Collapse into 0-3 score (0=0, 1-2=1, 3-4=2, 5-6=3)
secondary <- secondary %>%
  mutate(component2_3_grp = case_when(
    component2_3 <= 15 ~ 0,
    component2_3 > 15 & component2_3 <= 30 ~ 1,
    component2_3 > 30 & component2_3 <= 60 ~ 2,
    component2_3 > 60 ~ 3,
    FALSE ~ NA_real_  # In case the sum is outside these ranges
  )
  )

secondary <- secondary %>%
  mutate(component2_sum = component2_2and3and4 + component2_3_grp)


secondary <- secondary %>%
  mutate(component2_total = case_when(
    component2_sum <= 0 ~ 0,
    component2_sum >= 1 & component2_sum <= 2 ~ 1,
    component2_sum >= 3 & component2_sum <= 4 ~ 2,
    component2_sum >= 5 ~ 3,
    FALSE ~ NA_real_  # In case the sum is outside these ranges
  )
  )

#Component 3: DURATION OF SLEEP
  # Duration of sleep (component4_1) assigned score from 0-3 (>7=0, 6-7=1, 5-6=2, <5=3)
secondary <- secondary %>%
  mutate(component4_1_grouped = case_when(
    component4_1 >= 7 ~ 0,
    component4_1 >= 6 & component4_1 < 7 ~ 1,
    component4_1 >= 5 & component4_1 < 6 ~ 2,
    component4_1 < 5 ~ 3,
    FALSE ~ NA_real_  # In case the sum is outside these ranges
  )
  )

#Component 4: SLEEP EFFICIENCY
secondary <- secondary %>%
  mutate(
    # Convert bedtime and wake time hh:mm format to numeric hours
    component4_2_2_hours = hour(component4_2_2) + minute(component4_2_2) / 60,
    component4_2_1_hours = hour(component4_2_1) + minute(component4_2_1) / 60,
    
    # Calculate total hours in bed:
    # If bed time is later than wake time (e.g. 23:00 vs 07:00), then sleep period = (24 - bed time) + wake time;
    # Otherwise (bed time after midnight), sleep period = wake time - bed time.
    total_hours_in_bed = ifelse(
      component4_2_2_hours > component4_2_1_hours,
      (24 - component4_2_2_hours) + component4_2_1_hours,
      component4_2_1_hours - component4_2_2_hours
    ),
    
    # If total_hours_in_bed is NA or 0, set it to NA to avoid division errors
    total_hours_in_bed = ifelse(is.na(total_hours_in_bed) | total_hours_in_bed == 0, NA, total_hours_in_bed),
    
    # Calculate sleep efficiency: (# hours slept / # hours in bed) x 100%
    sleep_efficiency = ifelse(!is.na(total_hours_in_bed) & total_hours_in_bed != 0, 
                              (component4_1 / total_hours_in_bed) * 100, 
                              NA),
    
    # Assign Component 4 score based on sleep efficiency
    component4_score = case_when(
      sleep_efficiency > 85  ~ 0,
      sleep_efficiency >= 75 ~ 1,
      sleep_efficiency >= 65 ~ 2,
      sleep_efficiency < 65  ~ 3,
      FALSE ~ NA_real_  # Assign NA if sleep efficiency is NA
    )
  )

#Component 5: SLEEP DISTURBANCE

secondary <- secondary %>%
  mutate(
    # Calculate component5_sum using rowSums with na.rm = FALSE so that sessions with incomplete responses are excluded from analysis
    component5_sum = rowSums(select(., component5, component5_v2, component5_v3, component5_v4, 
                                    component5_v5, component5_v6, component5_v7, component5_v8, component5_v9), 
                             na.rm = FALSE),
    # Create the grouping variable based on the calculated component5_sum
    component5_total = case_when(
      component5_sum == 0 ~ 0,
      component5_sum >= 1 & component5_sum <= 9 ~ 1,
      component5_sum > 9 & component5_sum <= 18 ~ 2,
      component5_sum > 18 & component5_sum <= 27 ~ 3,
      FALSE ~ NA_real_  # In case the sum is outside these ranges
    )
  )

# Component 6: USE OF SLEEP MEDICATION
  # "component6" column already holds the final component 6 score

# Component 7: DAYTIME DYSFUNCTION
secondary <- secondary %>%
  mutate(component7_sum = component7 + component7_9)

secondary <- secondary %>%
  mutate(component7_total = case_when(
    component7_sum == 0 ~ 0,
    component7_sum >= 1 & component7_sum <= 2 ~ 1,
    component7_sum >= 3 & component7_sum <= 4 ~ 2,
    component7_sum >= 5 & component7_sum <= 6 ~ 3,
    FALSE ~ NA_real_  # Default case for unexpected values
  ))

# PSQI Global SCORE
secondary <- secondary %>% 
  rowwise() %>% 
  mutate(PSQI_Component1 = sum(c_across(c("component1")), na.rm = FALSE), #Subjective Sleep Quality
         PSQI_Component2 = sum(c_across(c("component2_total")), na.rm = FALSE), #Sleep Latency 
         PSQI_Component3 = sum(c_across(c("component4_1_grouped")), na.rm = FALSE), #Sleep Duration
         PSQI_Component4 = sum(c_across(c("component4_score")), na.rm = FALSE), #Habitual Sleep Efficiency
         PSQI_Component5 = sum(c_across(c("component5_total")), na.rm = FALSE), #Sleep Disturbances
         PSQI_Component6 = sum(c_across(c("component6")), na.rm = FALSE), #Use of Sleep Medication
         PSQI_Component7 = sum(c_across(c("component7_total")), na.rm = FALSE), # Daytime Dysfunction
         PSQI_Global_Score = sum(c_across(c("PSQI_Component1", "PSQI_Component2", "PSQI_Component3", 
                                            "PSQI_Component4", "PSQI_Component5", "PSQI_Component6", 
                                            "PSQI_Component7")), na.rm = FALSE)) %>%  #Sum of all scores (0 to 21), PSQI >5 indicates poor sleep quality
  ungroup()

# Export PSQI component and global scores
secondary %>%
  select(subject_number, session_number, preorpost, ip, arm,
         PSQI_Component1, PSQI_Component2, PSQI_Component3,
         PSQI_Component4, PSQI_Component5, PSQI_Component6,
         PSQI_Component7, PSQI_Global_Score) %>%
  write.csv("Z:/Analysis/Kye/PSQI Analysis/PSQI_scored.csv", row.names = FALSE)
