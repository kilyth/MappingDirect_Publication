#######################################################
### Deep Brain Stimulation: When to go directional.
### Data Cleaning
### 2022 Katrin Petermann
#######################################################

cleandata <- function(dd){
  # cleaning up the original input file
  # subsetting the data to raw
  
  dd <- subset(dd, select = c("Patient.Name",
                              "Date.of.birth",
                              "Gender",
                              "Disease_onset",
                              "rater",
                              "date.of.Mapping",
                              "date.of.surgery",
                              "ledd_pre",
                              "ledd_6m",
                              "ledd_1y",
                              "MU3_pre_off",
                              "MU3_pre_on",
                              "MU3_1y_off",
                              "MU3_1y_on",
                              "Target..STN.1..GPi.2.",
                              "Side_C1..right.1..left.2.",
                              "Rig_base_C1.8",
                              "Rig_base_C9.16",
                              "Rig_best_C1.8",
                              "Rig_best_C9.16",
                              "Rig_best_C1",
                              "Rig_best_C2",
                              "Rig_best_C3",
                              "Rig_best_C4",
                              "Rig_best_C5",
                              "Rig_best_C6",
                              "Rig_best_C7",
                              "Rig_best_C8",
                              "Rig_best_C234",
                              "Rig_best_C567",
                              "Rig_best_C9",
                              "Rig_best_C10",
                              "Rig_best_C11",
                              "Rig_best_C12",
                              "Rig_best_C13",
                              "Rig_best_C14",
                              "Rig_best_C15",
                              "Rig_best_C16",
                              "Rig_best_C101112",
                              "Rig_best_C131415",
                              "Rig_mA_C1",
                              "Rig_mA_C2",
                              "Rig_mA_C3",
                              "Rig_mA_C4",
                              "Rig_mA_C5",
                              "Rig_mA_C6",
                              "Rig_mA_C7",
                              "Rig_mA_C8",
                              "Rig_mA_C234",
                              "Rig_mA_C567",
                              "Rig_mA_C9",
                              "Rig_mA_C10",
                              "Rig_mA_C11",
                              "Rig_mA_C12",
                              "Rig_mA_C13",
                              "Rig_mA_C14",
                              "Rig_mA_C15",
                              "Rig_mA_C16",
                              "Rig_mA_C101112",
                              "Rig_mA_C131415",
                              "Rig_Fr",
                              "Rig_pw",
                              "Side_eff_mA_C1",
                              "Side_eff_mA_C2",
                              "Side_eff_mA_C3",
                              "Side_eff_mA_C4",
                              "Side_eff_mA_C5",
                              "Side_eff_mA_C6",
                              "Side_eff_mA_C7",
                              "Side_eff_mA_C8",
                              "Side_eff_mA_C234",
                              "Side_eff_mA_C567",
                              "Side_eff_mA_C9",
                              "Side_eff_mA_C10",
                              "Side_eff_mA_C11",
                              "Side_eff_mA_C12",
                              "Side_eff_mA_C13",
                              "Side_eff_mA_C14",
                              "Side_eff_mA_C15",
                              "Side_eff_mA_C16",
                              "Side_eff_mA_C101112",
                              "Side_eff_mA_C131415",
                              "date_last_visit",
                              "steering_last_visit",
                              "steering_past",
                              "stim_mA_l",
                              "stim_pw_l",
                              "stim_mA_r",
                              "stim_pw_r",
                              "stim_fq"))
  
  # adjusting column names
  colnames(dd) <- c("patient",
                    "DOB",
                    "gender",
                    "onset",
                    "rater",
                    "mapping.date",
                    "surgery.date",
                    "ledd.pre",
                    "ledd.6m",
                    "ledd.1y",
                    "MU3.pre.off",
                    "MU3.pre.on",
                    "MU3.1y.off",
                    "MU3.1y.on",
                    "target",
                    "Side.C1",
                    "Rig.base.C1-8",
                    "Rig.base.C9-16",
                    "Rig.best.C1-8",
                    "Rig.best.C9-16",
                    "Rig.best.C1",
                    "Rig.best.C2",
                    "Rig.best.C3",
                    "Rig.best.C4",
                    "Rig.best.C5",
                    "Rig.best.C6",
                    "Rig.best.C7",
                    "Rig.best.C8",
                    "Rig.best.C234",
                    "Rig.best.C567",
                    "Rig.best.C9",
                    "Rig.best.C10",
                    "Rig.best.C11",
                    "Rig.best.C12",
                    "Rig.best.C13",
                    "Rig.best.C14",
                    "Rig.best.C15",
                    "Rig.best.C16",
                    "Rig.best.C101112",
                    "Rig.best.C131415",
                    "Rig.mA.C1",
                    "Rig.mA.C2",
                    "Rig.mA.C3",
                    "Rig.mA.C4",
                    "Rig.mA.C5",
                    "Rig.mA.C6",
                    "Rig.mA.C7",
                    "Rig.mA.C8",
                    "Rig.mA.C234",
                    "Rig.mA.C567",
                    "Rig.mA.C9",
                    "Rig.mA.C10",
                    "Rig.mA.C11",
                    "Rig.mA.C12",
                    "Rig.mA.C13",
                    "Rig.mA.C14",
                    "Rig.mA.C15",
                    "Rig.mA.C16",
                    "Rig.mA.C101112",
                    "Rig.mA.C131415",
                    "Rig.Fr",
                    "Rig.pw",
                    "Side.eff.mA.C1",
                    "Side.eff.mA.C2",
                    "Side.eff.mA.C3",
                    "Side.eff.mA.C4",
                    "Side.eff.mA.C5",
                    "Side.eff.mA.C6",
                    "Side.eff.mA.C7",
                    "Side.eff.mA.C8",
                    "Side.eff.mA.C234",
                    "Side.eff.mA.C567",
                    "Side.eff.mA.C9",
                    "Side.eff.mA.C10",
                    "Side.eff.mA.C11",
                    "Side.eff.mA.C12",
                    "Side.eff.mA.C13",
                    "Side.eff.mA.C14",
                    "Side.eff.mA.C15",
                    "Side.eff.mA.C16",
                    "Side.eff.mA.C101112",
                    "Side.eff.mA.C131415",
                    "date.last.visit",
                    "steering.last.visit",
                    "steering.past",
                    "stim.mA.l",
                    "stim.pw.l",
                    "stim.mA.r",
                    "stim.pw.r",
                    "stim.fq")
  
  ## Adjust date formats
  dd$DOB <- as.POSIXct(dd$DOB, format = "%d.%m.%Y")
  dd$mapping.date <- as.POSIXct(dd$mapping.date, format = "%d.%m.%Y")
  dd$date.last.visit <- as.POSIXct(dd$date.last.visit, format = "%d.%m.%Y")
  dd$surgery.date <- as.POSIXct(dd$surgery.date, format = "%d.%m.%Y")
  dd$onset <- as.POSIXct(as.character(dd$onset), format = "%Y")
  
  ## Adjust data formats
  dd$target <- factor(dd$target, levels = c(1, 2), labels = c("STN", "GPi"))
  dd$Side.C1 <- factor(dd$Side.C1, levels = c(1, 2), labels = c("right", "left"))
  dd$steering.last.visit <- as.logical(dd$steering.last.visit)
  dd$steering.past <- as.logical(dd$steering.past)
  
  ## round effect threshold to 0.5mA (very few exceptions)
  cols <- which(colnames(dd) == "Rig.mA.C1") : which(colnames(dd) == "Rig.mA.C131415")
  dd[, cols] <- round(dd[, cols]/0.5)*0.5
  
  ## set max effect threshold to 8mA (very few exceptions)
  dd[, cols] <- sapply(dd[, cols], function(x){pmin(8, x, na.rm = FALSE)})
  
  ## round side effect thresholds to 0.5mA (very few exceptions)
  cols <- which(colnames(dd) == "Side.eff.mA.C1") : which(colnames(dd) == "Side.eff.mA.C131415")
  dd[, cols] <- round(dd[, cols]/0.5)*0.5
  
  ## set max side effect threshold to 8mA (very few exceptions)
  dd[, cols] <- sapply(dd[, cols], function(x){pmin(8, x, na.rm = FALSE)})
  
  return(dd)
  rm(t)
}