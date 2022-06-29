# MappingDirect
# Data derivation
# 2022 Katrin Petermann

derivedata <- function(dd){
  # calculating derivations for best contact, therapeutic window and efficacy
  
  ## calculate missing time values
  dd$time.after.surgery <- difftime(dd$mapping.date, dd$surgery.date, units = "weeks")
  dd$time.after.surgery <- as.numeric(dd$time.after.surgery)
  dd$time.last.visit <- difftime(dd$date.last.visit, dd$surgery.date, units = "weeks")
  dd$age.at.surgery <- as.numeric(format(as.Date(dd$surgery.date, format="%d/%m/%Y"),"%Y")) - as.numeric(format(as.Date(dd$DOB, format="%d/%m/%Y"),"%Y"))
  dd$duration <- difftime(dd$surgery.date, dd$onset, units = "days")/365.25
  dd$duration <- as.numeric(dd$duration)
  dd$rater <- factor(dd$rater, levels = c("KP", "FK", "ID", "CL", "SF"), ordered = TRUE)
  dd$mu3_improvement_pre <- 100 * (dd$MU3.pre.off - dd$MU3.pre.on)/dd$MU3.pre.off
  dd$mu3_improvement_1y <- 100 * (dd$MU3.1y.off - dd$MU3.1y.on)/dd$MU3.1y.off
  
  #*****************************************************************************
  #*****************************************************************************
  # calculating therapeutic window for each contact
  #*****************************************************************************
  #*****************************************************************************
  
  ## TW is set to zero if negative
  
  # left
  dd$TW.C1 <- ifelse(dd$Rig.mA.C1 == 0 & !is.na(dd$Rig.mA.C1), NA, 
                     pmax(0, dd$Side.eff.mA.C1-dd$Rig.mA.C1, na.rm = TRUE))
  dd$TW.C2 <- ifelse(dd$Rig.mA.C2 == 0 & !is.na(dd$Rig.mA.C2), NA, 
                     pmax(0, dd$Side.eff.mA.C2-dd$Rig.mA.C2, na.rm = TRUE))
  dd$TW.C3 <- ifelse(dd$Rig.mA.C3 == 0 & !is.na(dd$Rig.mA.C3), NA, 
                     pmax(0, dd$Side.eff.mA.C3-dd$Rig.mA.C3, na.rm = TRUE))
  dd$TW.C4 <- ifelse(dd$Rig.mA.C4 == 0 & !is.na(dd$Rig.mA.C4), NA, 
                     pmax(0, dd$Side.eff.mA.C4-dd$Rig.mA.C4, na.rm = TRUE))
  dd$TW.C5 <- ifelse(dd$Rig.mA.C5 == 0 & !is.na(dd$Rig.mA.C5), NA, 
                     pmax(0, dd$Side.eff.mA.C5-dd$Rig.mA.C5, na.rm = TRUE))
  dd$TW.C6 <- ifelse(dd$Rig.mA.C6 == 0 & !is.na(dd$Rig.mA.C6), NA, 
                     pmax(0, dd$Side.eff.mA.C6-dd$Rig.mA.C6, na.rm = TRUE))
  dd$TW.C7 <- ifelse(dd$Rig.mA.C7 == 0 & !is.na(dd$Rig.mA.C7), NA, 
                     pmax(0, dd$Side.eff.mA.C7-dd$Rig.mA.C7, na.rm 
                          = TRUE))
  dd$TW.C8 <- ifelse(dd$Rig.mA.C8 == 0 & !is.na(dd$Rig.mA.C8), NA, pmax(0, dd$Side.eff.mA.C8-dd$Rig.mA.C8, na.rm = TRUE))
  dd$TW.C234 <- ifelse(dd$Rig.mA.C234 == 0 & !is.na(dd$Rig.mA.C234), NA, 
                       pmax(0, dd$Side.eff.mA.C234-dd$Rig.mA.C234, na.rm = TRUE))
  dd$TW.C567 <- ifelse(dd$Rig.mA.C567 == 0 & !is.na(dd$Rig.mA.C567), NA, 
                       pmax(0, dd$Side.eff.mA.C567-dd$Rig.mA.C567, na.rm = TRUE))
  #right
  dd$TW.C9 <- ifelse(dd$Rig.mA.C9 == 0 & !is.na(dd$Rig.mA.C9), NA, 
                     pmax(0, dd$Side.eff.mA.C9-dd$Rig.mA.C9, na.rm = TRUE))
  dd$TW.C10 <- ifelse(dd$Rig.mA.C10 == 0 & !is.na(dd$Rig.mA.C10), NA, 
                      pmax(0, dd$Side.eff.mA.C10-dd$Rig.mA.C10, na.rm = TRUE))
  dd$TW.C11 <- ifelse(dd$Rig.mA.C11 == 0 & !is.na(dd$Rig.mA.C11), NA, 
                      pmax(0, dd$Side.eff.mA.C11-dd$Rig.mA.C11, na.rm = TRUE))
  dd$TW.C12 <- ifelse(dd$Rig.mA.C12 == 0 & !is.na(dd$Rig.mA.C12), NA, 
                      pmax(0, dd$Side.eff.mA.C12-dd$Rig.mA.C12, na.rm = TRUE))
  dd$TW.C13 <- ifelse(dd$Rig.mA.C13 == 0 & !is.na(dd$Rig.mA.C13), NA, 
                      pmax(0, dd$Side.eff.mA.C13-dd$Rig.mA.C13, na.rm = TRUE))
  dd$TW.C14 <- ifelse(dd$Rig.mA.C14 == 0 & !is.na(dd$Rig.mA.C14), NA, 
                      pmax(0, dd$Side.eff.mA.C14-dd$Rig.mA.C14, na.rm = TRUE))
  dd$TW.C15 <- ifelse(dd$Rig.mA.C15 == 0 & !is.na(dd$Rig.mA.C15), NA, 
                      pmax(0, dd$Side.eff.mA.C15-dd$Rig.mA.C15, na.rm = TRUE))
  dd$TW.C16 <- ifelse(dd$Rig.mA.C16 == 0 & !is.na(dd$Rig.mA.C16), NA, 
                      pmax(0, dd$Side.eff.mA.C16-dd$Rig.mA.C16, na.rm = TRUE))
  dd$TW.C101112 <- ifelse(dd$Rig.mA.C101112 == 0 & !is.na(dd$Rig.mA.C101112), NA, 
                          pmax(0, dd$Side.eff.mA.C101112-dd$Rig.mA.C101112, na.rm = TRUE))
  dd$TW.C131415 <- ifelse(dd$Rig.mA.C131415 == 0 & !is.na(dd$Rig.mA.C131415), NA, 
                          pmax(0, dd$Side.eff.mA.C131415-dd$Rig.mA.C131415, na.rm = TRUE))
  
  
  #*****************************************************************************
  #*****************************************************************************
  # defining best level for therapeutic window
  #*****************************************************************************
  #*****************************************************************************
  
  #left
  
  bl.tw.l <- function(dd){
    if(is.na(dd$TW.C234) & !is.na(dd$TW.C567)){ return(c("C567", "b"))} # if C234 is missing
    if(!is.na(dd$TW.C234) & is.na(dd$TW.C567)){ return(c("C234", "a"))} # if C567 is missing
    if(is.na(dd$TW.C234) & is.na(dd$TW.C567)){ return(c(NA, NA))} # if both levels are missing
    if(dd$TW.C234 > dd$TW.C567){
        return(c("C234", "a"))
    } else {
      if(dd$TW.C234 < dd$TW.C567){
          return(c("C567", "b"))
      } else { ## C234 == C567
        if(sum(c(dd$TW.C2, dd$TW.C3, dd$TW.C4), na.rm = TRUE) > sum(c(dd$TW.C5, dd$TW.C6, dd$TW.C7), na.rm = TRUE)){
            return(c("C234", "c"))
        } else {
          if(sum(c(dd$TW.C2, dd$TW.C3, dd$TW.C4), na.rm = TRUE) < sum(c(dd$TW.C5, dd$TW.C6, dd$TW.C7), na.rm = TRUE)){
              return(c("C567", "d"))
          } else { return(c(NA, NA))}
        }
      }
    }
  }
  
  for(i in 1:nrow(dd)){
      tmp <- bl.tw.l(dd[i, ])
      dd$`BLTW_C1-8`[i] <- tmp[1]
      dd$BLTW1_dt[i] <- tmp[2]
  }
  
  
  #right
  bl.tw.r <- function(dd){
    if(is.na(dd$TW.C101112) & !is.na(dd$TW.C131415)){ return(c("C131415", "b"))}
    if(!is.na(dd$TW.C101112) & is.na(dd$TW.C131415)){ return(c("C101112", "a"))}
    if(is.na(dd$TW.C101112) & is.na(dd$TW.C131415)){ return(c(NA, NA))}
    if(dd$TW.C101112 > dd$TW.C131415){
        return(c("C101112", "a"))
    } else {
      if(dd$TW.C101112 < dd$TW.C131415){
          return(c("C131415", "b"))
      } else {
        if(sum(c(dd$TW.C10, dd$TW.C11, dd$TW.C12), na.rm = TRUE) > sum(c(dd$TW.C13, dd$TW.C14, dd$TW.C15), na.rm = TRUE)){
            return(c("C101112", "c"))
        } else {
          if(sum(c(dd$TW.C10, dd$TW.C11, dd$TW.C12), na.rm = TRUE) < sum(c(dd$TW.C13, dd$TW.C14, dd$TW.C15), na.rm = TRUE)){
              return(c("C131415", "d"))
          } else { return(c(NA, NA))}
        }
      }
    }
  }
  
  for(i in 1:nrow(dd)){
      tmp <- bl.tw.r(dd[i, ])
    dd$`BLTW_C9-16`[i] <- tmp[1]
    dd$BLTW9_dt[i] <- tmp[2]
  }
  
  #*****************************************************************************
  #*****************************************************************************
  # defining best level for effect threshold
  #*****************************************************************************
  #*****************************************************************************
  
  #left
  bl.et.l <- function(dd){
    if(is.na(dd$Rig.mA.C234) & !is.na(dd$Rig.mA.C567)){ return("C567")}
    if(!is.na(dd$Rig.mA.C234) & is.na(dd$Rig.mA.C567)){ return("C234")}
    if(is.na(dd$Rig.mA.C234) & is.na(dd$Rig.mA.C567)){ return(NA)}
    if(dd$Rig.mA.C234 < dd$Rig.mA.C567){
      return("C234")
    } else {
      if(dd$Rig.mA.C234 > dd$Rig.mA.C567){
        return("C567")
      } else {
        if(sum(c(dd$Rig.mA.C2, dd$Rig.mA.C3, dd$Rig.mA.C4), na.rm = TRUE) < sum(c(dd$Rig.mA.C5, dd$Rig.mA.C6, dd$Rig.mA.C7), na.rm = TRUE)){
          return("C234")
        } else {
          if(sum(c(dd$Rig.mA.C2, dd$Rig.mA.C3, dd$Rig.mA.C4), na.rm = TRUE) > sum(c(dd$Rig.mA.C5, dd$Rig.mA.C6, dd$Rig.mA.C7), na.rm = TRUE)){
            return("C567")
          } else { return(NA)}
        }
      }
    }
  }
  
  for(i in 1:nrow(dd)){
    dd$`BLET_C1-8`[i] <- bl.et.l(dd[i, ])
  }
  
  
  #right
  bl.et.r <- function(dd){
    if(is.na(dd$Rig.mA.C101112) & !is.na(dd$Rig.mA.C131415)){ return("C131415")}
    if(!is.na(dd$Rig.mA.C101112) & is.na(dd$Rig.mA.C131415)){ return("C101112")}
    if(is.na(dd$Rig.mA.C101112) & is.na(dd$Rig.mA.C131415)){ return(NA)}
    if(dd$Rig.mA.C101112 < dd$Rig.mA.C131415){
      return("C101112")
    } else {
      if(dd$Rig.mA.C101112 > dd$Rig.mA.C131415){
        return("C131415")
      } else {
        if(sum(c(dd$Rig.mA.C10, dd$Rig.mA.C11, dd$Rig.mA.C12), na.rm = TRUE) < sum(c(dd$Rig.mA.C13, dd$Rig.mA.C14, dd$Rig.mA.C15), na.rm = TRUE)){
          return("C101112")
        } else {
          if(sum(c(dd$Rig.mA.C10, dd$Rig.mA.C11, dd$Rig.mA.C12), na.rm = TRUE) > sum(c(dd$Rig.mA.C13, dd$Rig.mA.C14, dd$Rig.mA.C15), na.rm = TRUE)){
            return("C131415")
          } else { return(NA)}
        }
      }
    }
  }
  
  for(i in 1:nrow(dd)){
    dd$`BLET_C9-16`[i] <- bl.et.r(dd[i, ])
  }
  
  
  #*****************************************************************************
  #*****************************************************************************
  # cross checking BL when NA between therapeutic window and effect threshold
  #*****************************************************************************
  #*****************************************************************************
  
  for(i in 1:nrow(dd)){
      if(is.na(dd$`BLTW_C1-8`[i])){
          dd$`BLTW_C1-8`[i] <- dd$`BLET_C1-8`[i]
          if(is.na(dd$`BLTW_C1-8`[i])) { 
              dd$BLTW1_dt[i] <- "g" 
              next
          }
          if(dd$`BLTW_C1-8`[i] == "C234") { 
              dd$BLTW1_dt[i] <- "e" 
              next
          }
          if(dd$`BLTW_C1-8`[i] == "C567") { 
              dd$BLTW1_dt[i] <- "f" 
          }
      }
  }
  for(i in 1:nrow(dd)){
      if(is.na(dd$`BLTW_C9-16`[i])){
          dd$`BLTW_C9-16`[i] <- dd$`BLET_C9-16`[i]
          if(is.na(dd$`BLTW_C9-16`[i])) { 
              dd$BLTW9_dt[i] <- "g"
              next
          }
          if(dd$`BLTW_C9-16`[i] == "C101112") { 
              dd$BLTW9_dt[i] <- "e" 
              next
          }
          if(dd$`BLTW_C9-16`[i] == "C131415") { 
              dd$BLTW9_dt[i] <- "f"
          }
      }
  }
  
  #*****************************************************************************
  #*****************************************************************************
  # defining contact ranking for therapeutic window
  #*****************************************************************************
  #*****************************************************************************
  
  #left
  for(i in 1:nrow(dd)){
    if(is.na(dd$`BLTW_C1-8`[i])){
      dd$`TW_BL_C1-8`[i] <- NA
      dd$`TW_BL-I_C1-8`[i] <- NA
      dd$`TW_BL-II_C1-8`[i] <- NA
      dd$`TW_BL-III_C1-8`[i] <- NA
      dd$`TW_WL_C1-8`[i] <- NA
      dd$`TW_WL-I_C1-8`[i] <- NA
      dd$`TW_WL-II_C1-8`[i] <- NA
      dd$`TW_WL-III_C1-8`[i] <- NA
      dd$`TW_BL_ET_C1-8`[i] <- NA
      dd$`TW_BL-I_ET_C1-8`[i] <- NA
      dd$`TW_BL-II_ET_C1-8`[i] <- NA
      dd$`TW_BL-III_ET_C1-8`[i] <- NA
      dd$`TW_WL_ET_C1-8`[i] <- NA
      dd$`TW_WL-I_ET_C1-8`[i] <- NA
      dd$`TW_WL-II_ET_C1-8`[i] <- NA
      dd$`TW_WL-III_ET_C1-8`[i] <- NA
      dd$`TW_BL_SET_C1-8`[i] <- NA
      dd$`TW_BL-I_SET_C1-8`[i] <- NA
      dd$`TW_BL-II_SET_C1-8`[i] <- NA
      dd$`TW_BL-III_SET_C1-8`[i] <- NA
      dd$`TW_WL_SET_C1-8`[i] <- NA
      dd$`TW_WL-I_SET_C1-8`[i] <- NA
      dd$`TW_WL-II_SET_C1-8`[i] <- NA
      dd$`TW_WL-III_SET_C1-8`[i] <- NA
      
    } else {
      if(dd$`BLTW_C1-8`[i] == "C234"){
        dd$`TW_BL_C1-8`[i] <- dd$TW.C234[i]
        dd$`TW_BL_ET_C1-8`[i] <- dd$Rig.mA.C234[i]
        dd$`TW_BL_SET_C1-8`[i] <- dd$Side.eff.mA.C234[i]
        
        csort <- sort(c(dd$TW.C2[i], dd$TW.C3[i], dd$TW.C4[i]), decreasing = TRUE, na.last = TRUE, index.return = TRUE)
        
        # define TW ET and SET for best level
        
        if(csort$ix[1] == 1){
          dd$`TW_BL-I_C1-8`[i] <- dd$TW.C2[i]
          dd$`TW_BL-I_ET_C1-8`[i] <- dd$Rig.mA.C2[i]
          dd$`TW_BL-I_SET_C1-8`[i] <- dd$Side.eff.mA.C2[i]
          if(csort$ix[2] == 2){
            dd$`TW_BL-II_C1-8`[i] <- dd$TW.C3[i]
            dd$`TW_BL-II_ET_C1-8`[i] <- dd$Rig.mA.C3[i]
            dd$`TW_BL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C3[i]
            dd$`TW_BL-III_C1-8`[i] <- dd$TW.C4[i]
            dd$`TW_BL-III_ET_C1-8`[i] <- dd$Rig.mA.C4[i]
            dd$`TW_BL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C4[i]
          }
          if(csort$ix[2] == 3){
            dd$`TW_BL-II_C1-8`[i] <- dd$TW.C4[i]
            dd$`TW_BL-II_ET_C1-8`[i] <- dd$Rig.mA.C4[i]
            dd$`TW_BL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C4[i]
            dd$`TW_BL-III_C1-8`[i] <- dd$TW.C3[i]
            dd$`TW_BL-III_ET_C1-8`[i] <- dd$Rig.mA.C3[i]
            dd$`TW_BL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C3[i]
          }
        }
        
        if(csort$ix[1] == 2){
          dd$`TW_BL-I_C1-8`[i] <- dd$TW.C3[i]
          dd$`TW_BL-I_ET_C1-8`[i] <- dd$Rig.mA.C3[i]
          dd$`TW_BL-I_SET_C1-8`[i] <- dd$Side.eff.mA.C3[i]
          if(csort$ix[2] == 1){
            dd$`TW_BL-II_C1-8`[i] <- dd$TW.C2[i]
            dd$`TW_BL-II_ET_C1-8`[i] <- dd$Rig.mA.C2[i]
            dd$`TW_BL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C2[i]
            dd$`TW_BL-III_C1-8`[i] <- dd$TW.C4[i]
            dd$`TW_BL-III_ET_C1-8`[i] <- dd$Rig.mA.C4[i]
            dd$`TW_BL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C4[i]
          }
          if(csort$ix[2] == 3){
            dd$`TW_BL-II_C1-8`[i] <- dd$TW.C4[i]
            dd$`TW_BL-II_ET_C1-8`[i] <- dd$Rig.mA.C4[i]
            dd$`TW_BL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C4[i]
            dd$`TW_BL-III_C1-8`[i] <- dd$TW.C2[i]
            dd$`TW_BL-III_ET_C1-8`[i] <- dd$Rig.mA.C2[i]
            dd$`TW_BL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C2[i]
          }
        }
        
        if(csort$ix[1] == 3){
          dd$`TW_BL-I_C1-8`[i] <- dd$TW.C4[i]
          dd$`TW_BL-I_ET_C1-8`[i] <- dd$Rig.mA.C4[i]
          dd$`TW_BL-I_SET_C1-8`[i] <- dd$Side.eff.mA.C4[i]
          if(csort$ix[2] == 1){
            dd$`TW_BL-II_C1-8`[i] <- dd$TW.C2[i]
            dd$`TW_BL-II_ET_C1-8`[i] <- dd$Rig.mA.C2[i]
            dd$`TW_BL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C2[i]
            dd$`TW_BL-III_C1-8`[i] <- dd$TW.C3[i]
            dd$`TW_BL-III_ET_C1-8`[i] <- dd$Rig.mA.C3[i]
            dd$`TW_BL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C3[i]
          }
          if(csort$ix[2] == 2){
            dd$`TW_BL-II_C1-8`[i] <- dd$TW.C3[i]
            dd$`TW_BL-II_ET_C1-8`[i] <- dd$Rig.mA.C3[i]
            dd$`TW_BL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C3[i]
            dd$`TW_BL-III_C1-8`[i] <- dd$TW.C2[i]
            dd$`TW_BL-III_ET_C1-8`[i] <- dd$Rig.mA.C2[i]
            dd$`TW_BL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C2[i]
          }
        }
             
        dd$`TW_WL_C1-8`[i] <- dd$TW.C567[i]
        dd$`TW_WL_ET_C1-8`[i] <- dd$Rig.mA.C567[i]
        dd$`TW_WL_SET_C1-8`[i] <- dd$Side.eff.mA.C567[i]
        
        
        csort <- sort(c(dd$TW.C5[i], dd$TW.C6[i], dd$TW.C7[i]), decreasing = TRUE, na.last = TRUE, index.return = TRUE)
        
        # define TW ET and SET for worst level
        
        if(csort$ix[1] == 1){
          dd$`TW_WL-I_C1-8`[i] <- dd$TW.C5[i]
          dd$`TW_WL-I_ET_C1-8`[i] <- dd$Rig.mA.C5[i]
          dd$`TW_WL-I_SET_C1-8`[i] <- dd$Side.eff.mA.C5[i]
          if(csort$ix[2] == 2){
            dd$`TW_WL-II_C1-8`[i] <- dd$TW.C6[i]
            dd$`TW_WL-II_ET_C1-8`[i] <- dd$Rig.mA.C6[i]
            dd$`TW_WL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C6[i]
            dd$`TW_WL-III_C1-8`[i] <- dd$TW.C7[i]
            dd$`TW_WL-III_ET_C1-8`[i] <- dd$Rig.mA.C7[i]
            dd$`TW_WL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C7[i]
          }
          if(csort$ix[2] == 3){
            dd$`TW_WL-II_C1-8`[i] <- dd$TW.C7[i]
            dd$`TW_WL-II_ET_C1-8`[i] <- dd$Rig.mA.C7[i]
            dd$`TW_WL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C7[i]
            dd$`TW_WL-III_C1-8`[i] <- dd$TW.C6[i]
            dd$`TW_WL-III_ET_C1-8`[i] <- dd$Rig.mA.C6[i]
            dd$`TW_WL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C6[i]
          }
        }
        
        if(csort$ix[1] == 2){
          dd$`TW_WL-I_C1-8`[i] <- dd$TW.C6[i]
          dd$`TW_WL-I_ET_C1-8`[i] <- dd$Rig.mA.C6[i]
          dd$`TW_WL-I_SET_C1-8`[i] <- dd$Side.eff.mA.C6[i]
          if(csort$ix[2] == 1){
            dd$`TW_WL-II_C1-8`[i] <- dd$TW.C5[i]
            dd$`TW_WL-II_ET_C1-8`[i] <- dd$Rig.mA.C5[i]
            dd$`TW_WL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C5[i]
            dd$`TW_WL-III_C1-8`[i] <- dd$TW.C7[i]
            dd$`TW_WL-III_ET_C1-8`[i] <- dd$Rig.mA.C7[i]
            dd$`TW_WL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C7[i]
          }
          if(csort$ix[2] == 3){
            dd$`TW_WL-II_C1-8`[i] <- dd$TW.C7[i]
            dd$`TW_WL-II_ET_C1-8`[i] <- dd$Rig.mA.C7[i]
            dd$`TW_WL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C7[i]
            dd$`TW_WL-III_C1-8`[i] <- dd$TW.C5[i]
            dd$`TW_WL-III_ET_C1-8`[i] <- dd$Rig.mA.C5[i]
            dd$`TW_WL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C5[i]
          }
        }
        
        if(csort$ix[1] == 3){
          dd$`TW_WL-I_C1-8`[i] <- dd$TW.C7[i]
          dd$`TW_WL-I_ET_C1-8`[i] <- dd$Rig.mA.C7[i]
          dd$`TW_WL-I_SET_C1-8`[i] <- dd$Side.eff.mA.C7[i]
          if(csort$ix[2] == 1){
            dd$`TW_WL-II_C1-8`[i] <- dd$TW.C5[i]
            dd$`TW_WL-II_ET_C1-8`[i] <- dd$Rig.mA.C5[i]
            dd$`TW_WL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C5[i]
            dd$`TW_WL-III_C1-8`[i] <- dd$TW.C6[i]
            dd$`TW_WL-III_ET_C1-8`[i] <- dd$Rig.mA.C6[i]
            dd$`TW_WL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C6[i]
          }
          if(csort$ix[2] == 2){
            dd$`TW_WL-II_C1-8`[i] <- dd$TW.C6[i]
            dd$`TW_WL-II_ET_C1-8`[i] <- dd$Rig.mA.C6[i]
            dd$`TW_WL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C6[i]
            dd$`TW_WL-III_C1-8`[i] <- dd$TW.C5[i]
            dd$`TW_WL-III_ET_C1-8`[i] <- dd$Rig.mA.C5[i]
            dd$`TW_WL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C5[i]
          }
        }
        
      } else {
        if(dd$`BLTW_C1-8`[i] == "C567"){
          
          dd$`TW_BL_C1-8`[i] <- dd$TW.C567[i]
          dd$`TW_BL_ET_C1-8`[i] <- dd$Rig.mA.C567[i]
          dd$`TW_BL_SET_C1-8`[i] <- dd$Side.eff.mA.C567[i]
          
          csort <- sort(c(dd$TW.C5[i], dd$TW.C6[i], dd$TW.C7[i]), decreasing = TRUE, na.last = TRUE, index.return = TRUE)
          
          # define TW ET and SET for worst level
          
          if(csort$ix[1] == 1){
            dd$`TW_BL-I_C1-8`[i] <- dd$TW.C5[i]
            dd$`TW_BL-I_ET_C1-8`[i] <- dd$Rig.mA.C5[i]
            dd$`TW_BL-I_SET_C1-8`[i] <- dd$Side.eff.mA.C5[i]
            if(csort$ix[2] == 2){
              dd$`TW_BL-II_C1-8`[i] <- dd$TW.C6[i]
              dd$`TW_BL-II_ET_C1-8`[i] <- dd$Rig.mA.C6[i]
              dd$`TW_BL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C6[i]
              dd$`TW_BL-III_C1-8`[i] <- dd$TW.C7[i]
              dd$`TW_BL-III_ET_C1-8`[i] <- dd$Rig.mA.C7[i]
              dd$`TW_BL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C7[i]
            }
            if(csort$ix[2] == 3){
              dd$`TW_BL-II_C1-8`[i] <- dd$TW.C7[i]
              dd$`TW_BL-II_ET_C1-8`[i] <- dd$Rig.mA.C7[i]
              dd$`TW_BL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C7[i]
              dd$`TW_BL-III_C1-8`[i] <- dd$TW.C6[i]
              dd$`TW_BL-III_ET_C1-8`[i] <- dd$Rig.mA.C6[i]
              dd$`TW_BL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C6[i]
            }
          }
          
          if(csort$ix[1] == 2){
            dd$`TW_BL-I_C1-8`[i] <- dd$TW.C6[i]
            dd$`TW_BL-I_ET_C1-8`[i] <- dd$Rig.mA.C6[i]
            dd$`TW_BL-I_SET_C1-8`[i] <- dd$Side.eff.mA.C6[i]
            if(csort$ix[2] == 1){
              dd$`TW_BL-II_C1-8`[i] <- dd$TW.C5[i]
              dd$`TW_BL-II_ET_C1-8`[i] <- dd$Rig.mA.C5[i]
              dd$`TW_BL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C5[i]
              dd$`TW_BL-III_C1-8`[i] <- dd$TW.C7[i]
              dd$`TW_BL-III_ET_C1-8`[i] <- dd$Rig.mA.C7[i]
              dd$`TW_BL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C7[i]
            }
            if(csort$ix[2] == 3){
              dd$`TW_BL-II_C1-8`[i] <- dd$TW.C7[i]
              dd$`TW_BL-II_ET_C1-8`[i] <- dd$Rig.mA.C7[i]
              dd$`TW_BL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C7[i]
              dd$`TW_BL-III_C1-8`[i] <- dd$TW.C5[i]
              dd$`TW_BL-III_ET_C1-8`[i] <- dd$Rig.mA.C5[i]
              dd$`TW_BL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C5[i]
            }
          }
          
          if(csort$ix[1] == 3){
            dd$`TW_BL-I_C1-8`[i] <- dd$TW.C7[i]
            dd$`TW_BL-I_ET_C1-8`[i] <- dd$Rig.mA.C7[i]
            dd$`TW_BL-I_SET_C1-8`[i] <- dd$Side.eff.mA.C7[i]
            if(csort$ix[2] == 1){
              dd$`TW_BL-II_C1-8`[i] <- dd$TW.C5[i]
              dd$`TW_BL-II_ET_C1-8`[i] <- dd$Rig.mA.C5[i]
              dd$`TW_BL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C5[i]
              dd$`TW_BL-III_C1-8`[i] <- dd$TW.C6[i]
              dd$`TW_BL-III_ET_C1-8`[i] <- dd$Rig.mA.C6[i]
              dd$`TW_BL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C6[i]
            }
            if(csort$ix[2] == 2){
              dd$`TW_BL-II_C1-8`[i] <- dd$TW.C6[i]
              dd$`TW_BL-II_ET_C1-8`[i] <- dd$Rig.mA.C6[i]
              dd$`TW_BL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C6[i]
              dd$`TW_BL-III_C1-8`[i] <- dd$TW.C5[i]
              dd$`TW_BL-III_ET_C1-8`[i] <- dd$Rig.mA.C5[i]
              dd$`TW_BL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C5[i]
            }
          }
          
            
          dd$`TW_WL_C1-8`[i] <- dd$TW.C234[i]
          dd$`TW_WL_ET_C1-8`[i] <- dd$Rig.mA.C234[i]
          dd$`TW_WL_SET_C1-8`[i] <- dd$Side.eff.mA.C234[i]
          
          if(csort$ix[1] == 1){
            dd$`TW_WL-I_C1-8`[i] <- dd$TW.C2[i]
            dd$`TW_WL-I_ET_C1-8`[i] <- dd$Rig.mA.C2[i]
            dd$`TW_WL-I_SET_C1-8`[i] <- dd$Side.eff.mA.C2[i]
            if(csort$ix[2] == 2){
              dd$`TW_WL-II_C1-8`[i] <- dd$TW.C3[i]
              dd$`TW_WL-II_ET_C1-8`[i] <- dd$Rig.mA.C3[i]
              dd$`TW_WL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C3[i]
              dd$`TW_WL-III_C1-8`[i] <- dd$TW.C4[i]
              dd$`TW_WL-III_ET_C1-8`[i] <- dd$Rig.mA.C4[i]
              dd$`TW_WL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C4[i]
            }
            if(csort$ix[2] == 3){
              dd$`TW_WL-II_C1-8`[i] <- dd$TW.C4[i]
              dd$`TW_WL-II_ET_C1-8`[i] <- dd$Rig.mA.C4[i]
              dd$`TW_WL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C4[i]
              dd$`TW_WL-III_C1-8`[i] <- dd$TW.C3[i]
              dd$`TW_WL-III_ET_C1-8`[i] <- dd$Rig.mA.C3[i]
              dd$`TW_WL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C3[i]
            }
          }
          
          if(csort$ix[1] == 2){
            dd$`TW_WL-I_C1-8`[i] <- dd$TW.C3[i]
            dd$`TW_WL-I_ET_C1-8`[i] <- dd$Rig.mA.C3[i]
            dd$`TW_WL-I_SET_C1-8`[i] <- dd$Side.eff.mA.C3[i]
            if(csort$ix[2] == 1){
              dd$`TW_WL-II_C1-8`[i] <- dd$TW.C2[i]
              dd$`TW_WL-II_ET_C1-8`[i] <- dd$Rig.mA.C2[i]
              dd$`TW_WL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C2[i]
              dd$`TW_WL-III_C1-8`[i] <- dd$TW.C4[i]
              dd$`TW_WL-III_ET_C1-8`[i] <- dd$Rig.mA.C4[i]
              dd$`TW_WL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C4[i]
            }
            if(csort$ix[2] == 3){
              dd$`TW_WL-II_C1-8`[i] <- dd$TW.C4[i]
              dd$`TW_WL-II_ET_C1-8`[i] <- dd$Rig.mA.C4[i]
              dd$`TW_WL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C4[i]
              dd$`TW_WL-III_C1-8`[i] <- dd$TW.C2[i]
              dd$`TW_WL-III_ET_C1-8`[i] <- dd$Rig.mA.C2[i]
              dd$`TW_WL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C2[i]
            }
          }
          
          if(csort$ix[1] == 3){
            dd$`TW_WL-I_C1-8`[i] <- dd$TW.C4[i]
            dd$`TW_WL-I_ET_C1-8`[i] <- dd$Rig.mA.C4[i]
            dd$`TW_WL-I_SET_C1-8`[i] <- dd$Side.eff.mA.C4[i]
            if(csort$ix[2] == 1){
              dd$`TW_WL-II_C1-8`[i] <- dd$TW.C2[i]
              dd$`TW_WL-II_ET_C1-8`[i] <- dd$Rig.mA.C2[i]
              dd$`TW_WL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C2[i]
              dd$`TW_WL-III_C1-8`[i] <- dd$TW.C3[i]
              dd$`TW_WL-III_ET_C1-8`[i] <- dd$Rig.mA.C3[i]
              dd$`TW_WL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C3[i]
            }
            if(csort$ix[2] == 2){
              dd$`TW_WL-II_C1-8`[i] <- dd$TW.C3[i]
              dd$`TW_WL-II_ET_C1-8`[i] <- dd$Rig.mA.C3[i]
              dd$`TW_WL-II_SET_C1-8`[i] <- dd$Side.eff.mA.C3[i]
              dd$`TW_WL-III_C1-8`[i] <- dd$TW.C2[i]
              dd$`TW_WL-III_ET_C1-8`[i] <- dd$Rig.mA.C2[i]
              dd$`TW_WL-III_SET_C1-8`[i] <- dd$Side.eff.mA.C2[i]
            }
          }
          

        } else { cat("ERROR in derivations contact ranking TW left")}
      }
    }
  }
  
  #right
  for(i in 1:nrow(dd)){
    if(is.na(dd$`BLTW_C9-16`[i])){
      dd$`TW_BL_C9-16`[i] <- NA
      dd$`TW_BL-I_C9-16`[i] <- NA
      dd$`TW_BL-II_C9-16`[i] <- NA
      dd$`TW_BL-III_C9-16`[i] <- NA
      dd$`TW_WL_C9-16`[i] <- NA
      dd$`TW_WL-I_C9-16`[i] <- NA
      dd$`TW_WL-II_C9-16`[i] <- NA
      dd$`TW_WL-III_C9-16`[i] <- NA
      dd$`TW_BL_ET_C9-16`[i] <- NA
      dd$`TW_BL-I_ET_C9-16`[i] <- NA
      dd$`TW_BL-II_ET_C9-16`[i] <- NA
      dd$`TW_BL-III_ET_C9-16`[i] <- NA
      dd$`TW_WL_ET_C9-16`[i] <- NA
      dd$`TW_WL-I_ET_C9-16`[i] <- NA
      dd$`TW_WL-II_ET_C9-16`[i] <- NA
      dd$`TW_WL-III_ET_C9-16`[i] <- NA
      dd$`TW_BL_SET_C9-16`[i] <- NA
      dd$`TW_BL-I_SET_C9-16`[i] <- NA
      dd$`TW_BL-II_SET_C9-16`[i] <- NA
      dd$`TW_BL-III_SET_C9-16`[i] <- NA
      dd$`TW_WL_SET_C9-16`[i] <- NA
      dd$`TW_WL-I_SET_C9-16`[i] <- NA
      dd$`TW_WL-II_SET_C9-16`[i] <- NA
      dd$`TW_WL-III_SET_C9-16`[i] <- NA
      
    } else {
      if(dd$`BLTW_C9-16`[i] == "C101112"){
        dd$`TW_BL_C9-16`[i] <- dd$TW.C101112[i]
        dd$`TW_BL_ET_C9-16`[i] <- dd$Rig.mA.C101112[i]
        dd$`TW_BL_SET_C9-16`[i] <- dd$Side.eff.mA.C101112[i]
        
        csort <- sort(c(dd$TW.C10[i], dd$TW.C11[i], dd$TW.C12[i]), decreasing = TRUE, na.last = TRUE, index.return = TRUE)
        
        # define TW ET and SET for best level
        
        if(csort$ix[1] == 1){
          dd$`TW_BL-I_C9-16`[i] <- dd$TW.C10[i]
          dd$`TW_BL-I_ET_C9-16`[i] <- dd$Rig.mA.C10[i]
          dd$`TW_BL-I_SET_C9-16`[i] <- dd$Side.eff.mA.C10[i]
          if(csort$ix[2] == 2){
            dd$`TW_BL-II_C9-16`[i] <- dd$TW.C11[i]
            dd$`TW_BL-II_ET_C9-16`[i] <- dd$Rig.mA.C11[i]
            dd$`TW_BL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C11[i]
            dd$`TW_BL-III_C9-16`[i] <- dd$TW.C12[i]
            dd$`TW_BL-III_ET_C9-16`[i] <- dd$Rig.mA.C12[i]
            dd$`TW_BL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C12[i]
          }
          if(csort$ix[2] == 3){
            dd$`TW_BL-II_C9-16`[i] <- dd$TW.C12[i]
            dd$`TW_BL-II_ET_C9-16`[i] <- dd$Rig.mA.C12[i]
            dd$`TW_BL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C12[i]
            dd$`TW_BL-III_C9-16`[i] <- dd$TW.C11[i]
            dd$`TW_BL-III_ET_C9-16`[i] <- dd$Rig.mA.C11[i]
            dd$`TW_BL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C11[i]
          }
        }
        
        if(csort$ix[1] == 2){
          dd$`TW_BL-I_C9-16`[i] <- dd$TW.C11[i]
          dd$`TW_BL-I_ET_C9-16`[i] <- dd$Rig.mA.C11[i]
          dd$`TW_BL-I_SET_C9-16`[i] <- dd$Side.eff.mA.C11[i]
          if(csort$ix[2] == 1){
            dd$`TW_BL-II_C9-16`[i] <- dd$TW.C10[i]
            dd$`TW_BL-II_ET_C9-16`[i] <- dd$Rig.mA.C10[i]
            dd$`TW_BL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C10[i]
            dd$`TW_BL-III_C9-16`[i] <- dd$TW.C12[i]
            dd$`TW_BL-III_ET_C9-16`[i] <- dd$Rig.mA.C12[i]
            dd$`TW_BL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C12[i]
          }
          if(csort$ix[2] == 3){
            dd$`TW_BL-II_C9-16`[i] <- dd$TW.C12[i]
            dd$`TW_BL-II_ET_C9-16`[i] <- dd$Rig.mA.C12[i]
            dd$`TW_BL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C12[i]
            dd$`TW_BL-III_C9-16`[i] <- dd$TW.C10[i]
            dd$`TW_BL-III_ET_C9-16`[i] <- dd$Rig.mA.C10[i]
            dd$`TW_BL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C10[i]
          }
        }
        
        if(csort$ix[1] == 3){
          dd$`TW_BL-I_C9-16`[i] <- dd$TW.C12[i]
          dd$`TW_BL-I_ET_C9-16`[i] <- dd$Rig.mA.C12[i]
          dd$`TW_BL-I_SET_C9-16`[i] <- dd$Side.eff.mA.C12[i]
          if(csort$ix[2] == 1){
            dd$`TW_BL-II_C9-16`[i] <- dd$TW.C10[i]
            dd$`TW_BL-II_ET_C9-16`[i] <- dd$Rig.mA.C10[i]
            dd$`TW_BL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C10[i]
            dd$`TW_BL-III_C9-16`[i] <- dd$TW.C11[i]
            dd$`TW_BL-III_ET_C9-16`[i] <- dd$Rig.mA.C11[i]
            dd$`TW_BL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C11[i]
          }
          if(csort$ix[2] == 2){
            dd$`TW_BL-II_C9-16`[i] <- dd$TW.C11[i]
            dd$`TW_BL-II_ET_C9-16`[i] <- dd$Rig.mA.C11[i]
            dd$`TW_BL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C11[i]
            dd$`TW_BL-III_C9-16`[i] <- dd$TW.C10[i]
            dd$`TW_BL-III_ET_C9-16`[i] <- dd$Rig.mA.C10[i]
            dd$`TW_BL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C10[i]
          }
        }
        
        dd$`TW_WL_C9-16`[i] <- dd$TW.C131415[i]
        dd$`TW_WL_ET_C9-16`[i] <- dd$Rig.mA.C131415[i]
        dd$`TW_WL_SET_C9-16`[i] <- dd$Side.eff.mA.C131415[i]
        
        
        csort <- sort(c(dd$TW.C13[i], dd$TW.C14[i], dd$TW.C15[i]), decreasing = TRUE, na.last = TRUE, index.return = TRUE)
        
        # define TW ET and SET for worst level
        
        if(csort$ix[1] == 1){
          dd$`TW_WL-I_C9-16`[i] <- dd$TW.C13[i]
          dd$`TW_WL-I_ET_C9-16`[i] <- dd$Rig.mA.C13[i]
          dd$`TW_WL-I_SET_C9-16`[i] <- dd$Side.eff.mA.C13[i]
          if(csort$ix[2] == 2){
            dd$`TW_WL-II_C9-16`[i] <- dd$TW.C14[i]
            dd$`TW_WL-II_ET_C9-16`[i] <- dd$Rig.mA.C14[i]
            dd$`TW_WL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C14[i]
            dd$`TW_WL-III_C9-16`[i] <- dd$TW.C15[i]
            dd$`TW_WL-III_ET_C9-16`[i] <- dd$Rig.mA.C15[i]
            dd$`TW_WL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C15[i]
          }
          if(csort$ix[2] == 3){
            dd$`TW_WL-II_C9-16`[i] <- dd$TW.C15[i]
            dd$`TW_WL-II_ET_C9-16`[i] <- dd$Rig.mA.C15[i]
            dd$`TW_WL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C15[i]
            dd$`TW_WL-III_C9-16`[i] <- dd$TW.C14[i]
            dd$`TW_WL-III_ET_C9-16`[i] <- dd$Rig.mA.C14[i]
            dd$`TW_WL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C14[i]
          }
        }
        
        if(csort$ix[1] == 2){
          dd$`TW_WL-I_C9-16`[i] <- dd$TW.C14[i]
          dd$`TW_WL-I_ET_C9-16`[i] <- dd$Rig.mA.C14[i]
          dd$`TW_WL-I_SET_C9-16`[i] <- dd$Side.eff.mA.C14[i]
          if(csort$ix[2] == 1){
            dd$`TW_WL-II_C9-16`[i] <- dd$TW.C13[i]
            dd$`TW_WL-II_ET_C9-16`[i] <- dd$Rig.mA.C13[i]
            dd$`TW_WL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C13[i]
            dd$`TW_WL-III_C9-16`[i] <- dd$TW.C15[i]
            dd$`TW_WL-III_ET_C9-16`[i] <- dd$Rig.mA.C15[i]
            dd$`TW_WL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C15[i]
          }
          if(csort$ix[2] == 3){
            dd$`TW_WL-II_C9-16`[i] <- dd$TW.C15[i]
            dd$`TW_WL-II_ET_C9-16`[i] <- dd$Rig.mA.C15[i]
            dd$`TW_WL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C15[i]
            dd$`TW_WL-III_C9-16`[i] <- dd$TW.C13[i]
            dd$`TW_WL-III_ET_C9-16`[i] <- dd$Rig.mA.C13[i]
            dd$`TW_WL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C13[i]
          }
        }
        
        if(csort$ix[1] == 3){
          dd$`TW_WL-I_C9-16`[i] <- dd$TW.C15[i]
          dd$`TW_WL-I_ET_C9-16`[i] <- dd$Rig.mA.C15[i]
          dd$`TW_WL-I_SET_C9-16`[i] <- dd$Side.eff.mA.C15[i]
          if(csort$ix[2] == 1){
            dd$`TW_WL-II_C9-16`[i] <- dd$TW.C13[i]
            dd$`TW_WL-II_ET_C9-16`[i] <- dd$Rig.mA.C13[i]
            dd$`TW_WL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C13[i]
            dd$`TW_WL-III_C9-16`[i] <- dd$TW.C14[i]
            dd$`TW_WL-III_ET_C9-16`[i] <- dd$Rig.mA.C14[i]
            dd$`TW_WL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C14[i]
          }
          if(csort$ix[2] == 2){
            dd$`TW_WL-II_C9-16`[i] <- dd$TW.C14[i]
            dd$`TW_WL-II_ET_C9-16`[i] <- dd$Rig.mA.C14[i]
            dd$`TW_WL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C14[i]
            dd$`TW_WL-III_C9-16`[i] <- dd$TW.C13[i]
            dd$`TW_WL-III_ET_C9-16`[i] <- dd$Rig.mA.C13[i]
            dd$`TW_WL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C13[i]
          }
        }
        
      } else {
        if(dd$`BLTW_C9-16`[i] == "C131415"){
          
          dd$`TW_BL_C9-16`[i] <- dd$TW.C131415[i]
          dd$`TW_BL_ET_C9-16`[i] <- dd$Rig.mA.C131415[i]
          dd$`TW_BL_SET_C9-16`[i] <- dd$Side.eff.mA.C131415[i]
          
          csort <- sort(c(dd$TW.C13[i], dd$TW.C14[i], dd$TW.C15[i]), decreasing = TRUE, na.last = TRUE, index.return = TRUE)
          
          # define TW ET and SET for worst level
          
          if(csort$ix[1] == 1){
            dd$`TW_BL-I_C9-16`[i] <- dd$TW.C13[i]
            dd$`TW_BL-I_ET_C9-16`[i] <- dd$Rig.mA.C13[i]
            dd$`TW_BL-I_SET_C9-16`[i] <- dd$Side.eff.mA.C13[i]
            if(csort$ix[2] == 2){
              dd$`TW_BL-II_C9-16`[i] <- dd$TW.C14[i]
              dd$`TW_BL-II_ET_C9-16`[i] <- dd$Rig.mA.C14[i]
              dd$`TW_BL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C14[i]
              dd$`TW_BL-III_C9-16`[i] <- dd$TW.C15[i]
              dd$`TW_BL-III_ET_C9-16`[i] <- dd$Rig.mA.C15[i]
              dd$`TW_BL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C15[i]
            }
            if(csort$ix[2] == 3){
              dd$`TW_BL-II_C9-16`[i] <- dd$TW.C15[i]
              dd$`TW_BL-II_ET_C9-16`[i] <- dd$Rig.mA.C15[i]
              dd$`TW_BL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C15[i]
              dd$`TW_BL-III_C9-16`[i] <- dd$TW.C14[i]
              dd$`TW_BL-III_ET_C9-16`[i] <- dd$Rig.mA.C14[i]
              dd$`TW_BL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C14[i]
            }
          }
          
          if(csort$ix[1] == 2){
            dd$`TW_BL-I_C9-16`[i] <- dd$TW.C14[i]
            dd$`TW_BL-I_ET_C9-16`[i] <- dd$Rig.mA.C14[i]
            dd$`TW_BL-I_SET_C9-16`[i] <- dd$Side.eff.mA.C14[i]
            if(csort$ix[2] == 1){
              dd$`TW_BL-II_C9-16`[i] <- dd$TW.C13[i]
              dd$`TW_BL-II_ET_C9-16`[i] <- dd$Rig.mA.C13[i]
              dd$`TW_BL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C13[i]
              dd$`TW_BL-III_C9-16`[i] <- dd$TW.C15[i]
              dd$`TW_BL-III_ET_C9-16`[i] <- dd$Rig.mA.C15[i]
              dd$`TW_BL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C15[i]
            }
            if(csort$ix[2] == 3){
              dd$`TW_BL-II_C9-16`[i] <- dd$TW.C15[i]
              dd$`TW_BL-II_ET_C9-16`[i] <- dd$Rig.mA.C15[i]
              dd$`TW_BL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C15[i]
              dd$`TW_BL-III_C9-16`[i] <- dd$TW.C13[i]
              dd$`TW_BL-III_ET_C9-16`[i] <- dd$Rig.mA.C13[i]
              dd$`TW_BL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C13[i]
            }
          }
          
          if(csort$ix[1] == 3){
            dd$`TW_BL-I_C9-16`[i] <- dd$TW.C15[i]
            dd$`TW_BL-I_ET_C9-16`[i] <- dd$Rig.mA.C15[i]
            dd$`TW_BL-I_SET_C9-16`[i] <- dd$Side.eff.mA.C15[i]
            if(csort$ix[2] == 1){
              dd$`TW_BL-II_C9-16`[i] <- dd$TW.C13[i]
              dd$`TW_BL-II_ET_C9-16`[i] <- dd$Rig.mA.C13[i]
              dd$`TW_BL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C13[i]
              dd$`TW_BL-III_C9-16`[i] <- dd$TW.C14[i]
              dd$`TW_BL-III_ET_C9-16`[i] <- dd$Rig.mA.C14[i]
              dd$`TW_BL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C14[i]
            }
            if(csort$ix[2] == 2){
              dd$`TW_BL-II_C9-16`[i] <- dd$TW.C14[i]
              dd$`TW_BL-II_ET_C9-16`[i] <- dd$Rig.mA.C14[i]
              dd$`TW_BL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C14[i]
              dd$`TW_BL-III_C9-16`[i] <- dd$TW.C13[i]
              dd$`TW_BL-III_ET_C9-16`[i] <- dd$Rig.mA.C13[i]
              dd$`TW_BL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C13[i]
            }
          }
          
          
          dd$`TW_WL_C9-16`[i] <- dd$TW.C101112[i]
          dd$`TW_WL_ET_C9-16`[i] <- dd$Rig.mA.C101112[i]
          dd$`TW_WL_SET_C9-16`[i] <- dd$Side.eff.mA.C101112[i]
          
          if(csort$ix[1] == 1){
            dd$`TW_WL-I_C9-16`[i] <- dd$TW.C10[i]
            dd$`TW_WL-I_ET_C9-16`[i] <- dd$Rig.mA.C10[i]
            dd$`TW_WL-I_SET_C9-16`[i] <- dd$Side.eff.mA.C10[i]
            if(csort$ix[2] == 2){
              dd$`TW_WL-II_C9-16`[i] <- dd$TW.C11[i]
              dd$`TW_WL-II_ET_C9-16`[i] <- dd$Rig.mA.C11[i]
              dd$`TW_WL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C11[i]
              dd$`TW_WL-III_C9-16`[i] <- dd$TW.C12[i]
              dd$`TW_WL-III_ET_C9-16`[i] <- dd$Rig.mA.C12[i]
              dd$`TW_WL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C12[i]
            }
            if(csort$ix[2] == 3){
              dd$`TW_WL-II_C9-16`[i] <- dd$TW.C12[i]
              dd$`TW_WL-II_ET_C9-16`[i] <- dd$Rig.mA.C12[i]
              dd$`TW_WL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C12[i]
              dd$`TW_WL-III_C9-16`[i] <- dd$TW.C11[i]
              dd$`TW_WL-III_ET_C9-16`[i] <- dd$Rig.mA.C11[i]
              dd$`TW_WL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C11[i]
            }
          }
          
          if(csort$ix[1] == 2){
            dd$`TW_WL-I_C9-16`[i] <- dd$TW.C11[i]
            dd$`TW_WL-I_ET_C9-16`[i] <- dd$Rig.mA.C11[i]
            dd$`TW_WL-I_SET_C9-16`[i] <- dd$Side.eff.mA.C11[i]
            if(csort$ix[2] == 1){
              dd$`TW_WL-II_C9-16`[i] <- dd$TW.C10[i]
              dd$`TW_WL-II_ET_C9-16`[i] <- dd$Rig.mA.C10[i]
              dd$`TW_WL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C10[i]
              dd$`TW_WL-III_C9-16`[i] <- dd$TW.C12[i]
              dd$`TW_WL-III_ET_C9-16`[i] <- dd$Rig.mA.C12[i]
              dd$`TW_WL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C12[i]
            }
            if(csort$ix[2] == 3){
              dd$`TW_WL-II_C9-16`[i] <- dd$TW.C12[i]
              dd$`TW_WL-II_ET_C9-16`[i] <- dd$Rig.mA.C12[i]
              dd$`TW_WL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C12[i]
              dd$`TW_WL-III_C9-16`[i] <- dd$TW.C10[i]
              dd$`TW_WL-III_ET_C9-16`[i] <- dd$Rig.mA.C10[i]
              dd$`TW_WL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C10[i]
            }
          }
          
          if(csort$ix[1] == 3){
            dd$`TW_WL-I_C9-16`[i] <- dd$TW.C12[i]
            dd$`TW_WL-I_ET_C9-16`[i] <- dd$Rig.mA.C12[i]
            dd$`TW_WL-I_SET_C9-16`[i] <- dd$Side.eff.mA.C12[i]
            if(csort$ix[2] == 1){
              dd$`TW_WL-II_C9-16`[i] <- dd$TW.C10[i]
              dd$`TW_WL-II_ET_C9-16`[i] <- dd$Rig.mA.C10[i]
              dd$`TW_WL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C10[i]
              dd$`TW_WL-III_C9-16`[i] <- dd$TW.C11[i]
              dd$`TW_WL-III_ET_C9-16`[i] <- dd$Rig.mA.C11[i]
              dd$`TW_WL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C11[i]
            }
            if(csort$ix[2] == 2){
              dd$`TW_WL-II_C9-16`[i] <- dd$TW.C11[i]
              dd$`TW_WL-II_ET_C9-16`[i] <- dd$Rig.mA.C11[i]
              dd$`TW_WL-II_SET_C9-16`[i] <- dd$Side.eff.mA.C11[i]
              dd$`TW_WL-III_C9-16`[i] <- dd$TW.C10[i]
              dd$`TW_WL-III_ET_C9-16`[i] <- dd$Rig.mA.C10[i]
              dd$`TW_WL-III_SET_C9-16`[i] <- dd$Side.eff.mA.C10[i]
            }
          }
          
        } else { cat("ERROR in derivations contact ranking TW left")}
      }
    }
  }
  
  rm(csort, i)

#------------------------------------------------

  dd$directional <- 0

  for(i in 1:nrow(dd)){
    if(!is.na(dd$`TW_BL_C1-8`[i]) & !is.na(dd$`TW_BL-I_C1-8`[i])){
      if(dd$`TW_BL_C1-8`[i] < dd$`TW_BL-I_C1-8`[i]){
        dd$directional[i] <- dd$directional[i] + 1
      }
    }
    if(!is.na(dd$`TW_BL_C9-16`[i]) & !is.na(dd$`TW_BL-I_C9-16`[i])){
      if(dd$`TW_BL_C9-16`[i] < dd$`TW_BL-I_C9-16`[i]){
        dd$directional[i] <- dd$directional[i] + 1
      }
    }
  }
  
  dd$directional <- as.factor(dd$directional)
  dd$dirlog <- as.logical(dd$directional)
  
  return(dd)
}

makeLongData <- function(dd){
  # creating long format table
  SETcols <- c("TW_BL_SET_C1-8", "TW_BL-I_SET_C1-8", "TW_BL-II_SET_C1-8", "TW_BL-III_SET_C1-8",
               "TW_WL_SET_C1-8", "TW_WL-I_SET_C1-8", "TW_WL-II_SET_C1-8", "TW_WL-III_SET_C1-8",
               "TW_BL_SET_C9-16", "TW_BL-I_SET_C9-16", "TW_BL-II_SET_C9-16", "TW_BL-III_SET_C9-16",
               "TW_WL_SET_C9-16", "TW_WL-I_SET_C9-16", "TW_WL-II_SET_C9-16", "TW_WL-III_SET_C9-16")
  
  ETcols <- c("TW_BL_ET_C1-8", "TW_BL-I_ET_C1-8", "TW_BL-II_ET_C1-8", "TW_BL-III_ET_C1-8",
              "TW_WL_ET_C1-8", "TW_WL-I_ET_C1-8", "TW_WL-II_ET_C1-8", "TW_WL-III_ET_C1-8",
              "TW_BL_ET_C9-16", "TW_BL-I_ET_C9-16", "TW_BL-II_ET_C9-16", "TW_BL-III_ET_C9-16",
              "TW_WL_ET_C9-16", "TW_WL-I_ET_C9-16", "TW_WL-II_ET_C9-16", "TW_WL-III_ET_C9-16")
  
  TWcols <- c("TW_BL_C1-8", "TW_BL-I_C1-8", "TW_BL-II_C1-8", "TW_BL-III_C1-8",
              "TW_WL_C1-8", "TW_WL-I_C1-8", "TW_WL-II_C1-8", "TW_WL-III_C1-8",
              "TW_BL_C9-16", "TW_BL-I_C9-16", "TW_BL-II_C9-16", "TW_BL-III_C9-16",
              "TW_WL_C9-16", "TW_WL-I_C9-16", "TW_WL-II_C9-16", "TW_WL-III_C9-16")
  
  
  onePatientLong <- expand.grid(direction = c("ring level", "best \n direction", "second best \n direction", "worst \n direction"), 
                                level = c("best Level", "worst Level"),
                                side = c("1-8", "9-16"))
  longData <- data.frame(patient = factor(),
                         side = factor(),
                         level = factor(),
                         direction = factor(),
                         SET = double(),
                         ET = double(),
                         TW = double())
  
  vec <- rep(c(1, 5, 9, 13), each = 4)
  patients <- dd$patient
  
  testinf <- function(x){
    if(is.infinite(x)){
      return(NA)
    } else {
      return(x)
    }
  }
  
  for(i in patients){
    
    onePatientLong$patient <-  i
    whichRow <- which(dd$patient == i)
    
    onePatientLong$SET <- as.double(dd[whichRow, SETcols])
    onePatientLong$ET <- as.double(dd[whichRow, ETcols])
    onePatientLong$TW <- as.double(dd[whichRow, TWcols])
    onePatientLong$diffSET <- as.double(dd[whichRow, SETcols] - dd[whichRow, SETcols[vec]])
    onePatientLong$diffET <- as.double(dd[whichRow, ETcols] - dd[whichRow, ETcols[vec]])
    onePatientLong$diffTW <- as.double(dd[whichRow, TWcols] - dd[whichRow, TWcols[vec]])
    
    longData <- rbind(longData, onePatientLong)
  }
  longData$hemisphere <- interaction(longData$patient, longData$side)
  return(longData)
}

