setwd("~/Documents/GitHub.nosync/NationalTrialEnv/")
library(ggbiplot)
library(tidyverse)
library(GWASpoly)
source(file = "newReadGWASpoly.fn.R")
library(fuzzyjoin)
library(ggcorrplot)
only_letters <- function(x) { gsub("^([[:alpha:]]*).*$","\\1",x) }
library(cowplot)
conflicted::conflict_prefer_all("dplyr", quiet = T)

geno <- read_csv("ncpt_geno_2023-08-19.csv", show_col_types = F) %>% 
  mutate(chrom = factor(gsub("chr0?", "", chrom), levels = as.character(1:12)))

names(geno) <- str_replace(names(geno), "-", ".")
names(geno) <- str_replace(names(geno), " ", ".")
names(geno) <- str_replace(names(geno), "/", ".")

field.files <- list.files("./ncpt_field/", full.names = TRUE)
ncptList <- lapply(field.files, function(i){read_csv(i, show_col_types = F)})
chipFieldData <- bind_rows(ncptList)
rm(list = c("ncptList", "field.files"))

field_data_all <- chipFieldData %>% 
  filter(!is.na(variety), variety != "Trial Information", variety != "Variety") %>%
  select(variety, trial, yield_total) %>% 
  mutate(variety = str_replace(variety, "-", "."),
         variety = str_replace(variety, " ", "."),
         variety = str_replace(variety, "/", "."),
         yield_total = as.numeric(yield_total)) %>%
  separate(trial, c("year", NA, "location"), " ") %>%
  filter(yield_total > 0)

full.model <- lme4::lmer(yield_total ~ (year) + (location) + (year:location) +
                           (1|variety) + (1|variety:year) + (1|variety:location) +
                           (1|variety:year:location), 
                         data = field_data_all)

sum.full.model <- summary(full.model)

varComps <- tibble(term = c(names(unlist(sum.full.model$varcor)), "residuals"), 
                   component = c(unlist(sum.full.model$varcor), sqrt(sum.full.model$sigma)), 
                   pve = component/sum(component)*100)

field_data_all <- filter(field_data_all, variety %in% names(geno))

field_data_sel <- field_data_all %>%
  group_by(variety) %>%
  mutate(year.min = min(year)) %>%
  filter(year == year.min) %>%
  ungroup() %>%
  filter(!duplicated(variety)) %>%
  mutate(code = only_letters(variety))

programCodes <- read_csv("programCodes.csv")

selSites <- field_data_sel %>% 
  regex_left_join(programCodes, by = "code") %>% 
  filter(!is.na(sel.site)) %>%
  select(variety, year, sel.site) %>% ungroup()

precip <- read_csv("precip_final.csv", show_col_types = F) %>% 
  select(Feb:Sep, state, year) %>% 
  filter(year %in% 2007:2021) %>%
  pivot_longer(cols = Feb:Sep, names_to = "month", values_to = "precip") %>%
  filter(ifelse(state == "NorthCarolina", month %in% c("Mar", "Apr", "May", "Jun"), 
                ifelse(state %in% c("Wisconsin", "Michigan", "NewYork", "Colorado", "Texas", "Maine"), 
                       month %in% c("May", "Jun", "Jul", "Aug"),
                       ifelse(state %in% c("Idaho", "Oregon"),
                              month %in% c("Apr", "May", "Jun", "Jul"),
                              ifelse(state == "NorthDakota", month %in% c("Jun", "Jul", "Aug", "Sep"), NA))))) %>%
  mutate(precip = precip*2.54, month = factor(month, levels = month.abb))

maxTemp <- read_csv("max_temp_final.csv", show_col_types = F) %>% 
  select(Feb:Sep, state, year) %>%
  filter(year %in% 2007:2021) %>%
  pivot_longer(cols = Feb:Sep, names_to = "month", values_to = "maxTemp") %>%
  filter(ifelse(state == "NorthCarolina", month %in% c("Mar", "Apr", "May", "Jun"), 
                ifelse(state %in% c("Wisconsin", "Michigan", "NewYork", "Colorado", "Texas", "Maine"), 
                       month %in% c("May", "Jun", "Jul", "Aug"),
                       ifelse(state %in% c("Idaho", "Oregon"),
                              month %in% c("Apr", "May", "Jun", "Jul"),
                              ifelse(state == "NorthDakota", month %in% c("Jun", "Jul", "Aug", "Sep"), NA))))) %>%
  dplyr::mutate(maxTemp = (maxTemp - 32) * (5/9),
                month = factor(month, levels = month.abb))

minTemp <- read_csv("min_temp_final.csv", show_col_types = F) %>%
  select(Feb:Sep, state, year) %>%
  filter(year %in% 2007:2021) %>%
  pivot_longer(cols = Feb:Sep, names_to = "month", values_to = "minTemp") %>%
  filter(ifelse(state == "NorthCarolina", month %in% c("Mar", "Apr", "May", "Jun"), 
                ifelse(state %in% c("Wisconsin", "Michigan", "NewYork", "Colorado", "Texas", "Maine"), 
                       month %in% c("May", "Jun", "Jul", "Aug"),
                       ifelse(state %in% c("Idaho", "Oregon"),
                              month %in% c("Apr", "May", "Jun", "Jul"),
                              ifelse(state == "NorthDakota", month %in% c("Jun", "Jul", "Aug", "Sep"), NA))))) %>%
  dplyr::mutate(minTemp = (minTemp - 32) * (5/9),
                month = factor(month, levels = month.abb))

selSiteEnvs <- left_join(left_join(maxTemp, minTemp), 
                                  precip) %>% 
  group_by(state, year) %>% 
  dplyr::summarise(precip = sum(precip), maxTemp = mean(maxTemp), 
                   minTemp = mean(minTemp)) %>%  
  ungroup() %>%
  rename(sel.site = state)


envs.3.years <- tibble(sel.site = "State", year = 0, precip = 0, maxTemp = 0, minTemp = 0)
for(i in 2010:2022){
  start = i-3
  end = i-1
  envs.3.years <- add_row(envs.3.years, selSiteEnvs %>% 
                            group_by(sel.site) %>% 
                            filter(year %in% start:end) %>% 
                            summarise(precip = mean(precip),
                                             maxTemp = mean(maxTemp),
                                             minTemp = mean(minTemp), year = i))
}

envs.3.years <- filter(envs.3.years, sel.site != "State") %>% ungroup()


envs.selSite <- selSites %>% ungroup() %>% 
  mutate(year = as.numeric(year)) %>%
  left_join(envs.3.years, by = c("sel.site", "year")) %>%
  select(variety, year, precip, minTemp, maxTemp)

geno.sel <- select(geno, marker, chrom, bp, envs.selSite$variety)

params <- set.params(MAF = 0.01, n.PC = 5)
selData <- read.GWASpoly.object(ploidy = 4, phenos = envs.selSite, genos = geno.sel, format = "numeric", n.traits = 4)
selData <- set.K(data = selData, n.core = 4, LOCO = T)

models = c("additive")
selGWAS <- GWASpoly(selData, models = models, params = params, n.core = 4)
selThresh <- set.threshold(selGWAS, method = "M.eff", level = 0.05)

selectTrialMonths <- function(x) {
  trait = substitute(x)
  x = x %>%
  select(Feb:Sep, Location, Year) %>%
  pivot_longer(cols = Feb:Sep, names_to = "month", values_to = "trait") %>%
  filter(ifelse(Location %in% c("NC", "MO"), 
                month %in% c("Mar", "Apr", "May", "Jun"), 
                ifelse(Location %in% c("WI", "MI", "NY","TX"),
                       month %in% c("May", "Jun", "Jul", "Aug"),
                       ifelse(Location %in% c("OR"),
                              month %in% c("Apr", "May", "Jun", "Jul"),
                              ifelse(Location %in% c("CA", "FL"),
                                     month %in% c("Feb", "Mar", "Apr", "May"),
                                     ifelse(Location == "ND",
                                            month %in% c("Jun", "Jul", "Aug", "Sep"),
                                            NA)))))) %>%
  group_by(Location, Year) %>%
  dplyr::summarise(trait= mean(trait)) %>%
  group_by(Year) %>%
  unite(Year, Location, remove = T, col = "trial") %>% ungroup()
  
  if(gsub("precip", "", trait) == ""){
    x <- mutate(x, trait = 4 * trait * 2.54) # times 4 to give total rainfall not average over 4 months
  }else{
    x <- mutate(x, trait = (trait - 32)*5/9)
  }
  
  colnames(x) <- c("trial",trait)
  return(x)
}

precip <- read_csv("nationalTrialPrecipitation.csv", show_col_types = F)
maxTemp <- read_csv("newEnvs/maxTempTrialUpdated.csv", show_col_types = F)
minTemp <- read_csv("newEnvs/minTempTrialUpdated.csv", show_col_types = F)

envTrials <- selectTrialMonths(precip) %>% 
  left_join(selectTrialMonths(minTemp)) %>% 
  left_join(selectTrialMonths(maxTemp))

envTrialsObvs <- field_data_all %>%
  mutate(yield_total = as.numeric(yield_total)) %>%
  unite(year, location, col = "trial") %>%
  left_join(envTrials, by = "trial") %>% 
  group_by(trial) %>%
  filter(yield_total > 0, !is.na(yield_total)) %>%
  mutate(normalized_yield = log(yield_total/mean(yield_total, na.rm = T))) %>%
  ungroup()

envYield.lm <- envTrialsObvs %>% 
  group_by(variety) %>% 
  do({
    mod.precip = lm(normalized_yield ~ precip, data = .)
    mod.minTemp = lm(normalized_yield ~ minTemp, data = .)
    mod.maxTemp = lm(normalized_yield ~ maxTemp, data = .)
    data.frame(precip = coef(mod.precip)[2],
               minTemp = coef(mod.minTemp)[2],
               maxTemp = coef(mod.maxTemp)[2])
  }) %>% filter(!is.na(precip)) %>% ungroup()

geno.slope <- select(geno, marker, chrom, bp, envYield.lm$variety)

slopeData <- read.GWASpoly.object(ploidy = 4, phenos = envYield.lm, genos = geno.slope, format = "numeric", n.traits = 3)
slopeData <- set.K(data = slopeData, n.core = 4, LOCO = T)

models = c("additive")

slopeGWAS <- GWASpoly(slopeData, models = models, params = params, n.core = 4)
slopeThresh <- set.threshold(slopeGWAS, method = "M.eff", level = 0.05)

gps <- tibble(sel.site = c("NorthCarolina", "Texas", "Michigan", "NorthDakota", "NewYork",
                           "Wisconsin", "Colorado", "Oregon", "Maine", "Idaho"),
              long = c(-76.65000, -102.74410,  -85.17680,  -97.62810,  -76.34481,  
                       -89.52318, -106.14472, -119.28399,  -68.00620, -111.27726),
              lat = c(35.83330, 35.97130, 43.35050, 48.54690,42.51215, 
                      44.13358, 37.70517, 45.81683, 46.65390, 43.85560)) %>%
  left_join(selSites, by = "sel.site") %>% select(variety, lat, long)

gpsData <- read.GWASpoly.object(4, gps, genos = geno.sel, format = "numeric", n.traits = 2)
gpsData.k <- set.K(gpsData, LOCO = T)

gpsDataGWAS <- GWASpoly(gpsData.k,
                        models=c("additive"),
                        params=params,
                        n.core = 4)

gpsThresh <- set.threshold(gpsDataGWAS, method = "M.eff")

table(selSites$sel.site)/nrow(selSites) # select only states that make up >=10% of the total genotypes

state_binary <- selSites %>% mutate(NY = ifelse(sel.site == "NewYork", 1, 0),
                                    ME = ifelse(sel.site == "Maine", 1, 0),
                                    MI = ifelse(sel.site == "Michigan", 1, 0),
                                    WI = ifelse(sel.site == "Wisconsin", 1, 0)) %>%
                             select(variety, NY:WI)

stateData <- read.GWASpoly.object(4, state_binary, geno.sel, format = "numeric", n.traits = 4)
stateData <- set.K(stateData, LOCO = T)

stateScan <- GWASpoly(stateData,
                      models = "additive",
                      params = params,
                      n.core = 4)

stateThresh <- set.threshold(stateScan, method = "M.eff", n.core = 4)

yearQTL <- get.QTL(data=selThresh,trait = "year", bp.window = 5e6)
minQTL <- get.QTL(data=selThresh,trait = "minTemp", bp.window = 5e6)
precipQTL <- get.QTL(data=selThresh,trait = "precip", bp.window = 5e6)
minRegQTL <- get.QTL(data=slopeThresh,trait = "minTemp", bp.window = 5e6)
maxRegQTL <- get.QTL(data=slopeThresh,trait = "maxTemp", bp.window = 5e6)
latQTL <- get.QTL(data=gpsThresh,trait = "lat", bp.window = 5e6)
longQTL <- get.QTL(data=gpsThresh,trait = "long", bp.window = 5e6)
meQTL <- get.QTL(data=stateThresh,trait = "ME", bp.window = 5e6)
miQTL <- get.QTL(data=stateThresh,trait = "MI", bp.window = 5e6)
nyQTL <- get.QTL(data=stateThresh,trait = "NY", bp.window = 5e6)
wiQTL <- get.QTL(data=stateThresh,trait = "WI", bp.window = 5e6)

yearFit <- fit.QTL(data=selThresh,trait = "year",
        qtl=yearQTL[,c("Marker","Model")])

minFit <- fit.QTL(data=selThresh,trait = "minTemp",
        qtl=minQTL[,c("Marker","Model")])

precipFit <- fit.QTL(data=selThresh,trait = "precip",
        qtl=precipQTL[,c("Marker","Model")])

minRegFit <- fit.QTL(data=slopeThresh,trait = "minTemp",
        qtl=minRegQTL[,c("Marker","Model")])

maxRegFit <- fit.QTL(data=slopeThresh,trait = "maxTemp",
        qtl=maxRegQTL[,c("Marker","Model")])

latFit <- fit.QTL(data=gpsThresh,trait = "lat",
        qtl=latQTL[,c("Marker","Model")])

longFit <- fit.QTL(data=gpsThresh,trait = "long",
        qtl=longQTL[,c("Marker","Model")])

meFit <- fit.QTL(data=stateThresh,trait = "ME",
        qtl=meQTL[,c("Marker","Model")])

miFit <- fit.QTL(data=stateThresh,trait = "MI",
        qtl=miQTL[,c("Marker","Model")])

nyFit <- fit.QTL(data=stateThresh,trait = "NY",
        qtl=nyQTL[,c("Marker","Model")])

wiFit <- fit.QTL(data=stateThresh,trait = "WI",
        qtl=wiQTL[,c("Marker","Model")])

yearFit <- yearFit %>% 
  select(-pval, -Model) %>% 
  mutate(R2 = round(R2, digits = 4))

minFit <- minFit %>% 
  select(-pval, -Model) %>% 
  mutate(R2 = round(R2, digits = 4))

precipFit <- precipFit %>% 
  select(-pval, -Model) %>% 
  mutate(R2 = round(R2, digits = 4))

minRegFit <- minRegFit %>% 
  select(-pval, -Model) %>% 
  mutate(R2 = round(R2, digits = 4))

maxRegFit <- maxRegFit %>% 
  select(-pval, -Model) %>% 
  mutate(R2 = round(R2, digits = 4))

latFit <- latFit %>% 
  select(-pval, -Model) %>% 
  mutate(R2 = round(R2, digits = 4))

longFit <- longFit %>% 
  select(-pval, -Model) %>% 
  mutate(R2 = round(R2, digits = 4))

meFit <- meFit %>% 
  select(-pval, -Model) %>% 
  mutate(R2 = round(R2, digits = 4))

miFit <- miFit %>% 
  select(-pval, -Model) %>% 
  mutate(R2 = round(R2, digits = 4))

nyFit <- nyFit %>% 
  select(-pval, -Model) %>% 
  mutate(R2 = round(R2, digits = 4))

wiFit <- wiFit %>% 
  select(-pval, -Model) %>% 
  mutate(R2 = round(R2, digits = 4))

sum(minFit$R2)
sum(precipFit$R2)
sum(minRegFit$R2)
sum(maxRegFit$R2)
sum(yearFit$R2)
sum(latFit$R2)
sum(longFit$R2)
sum(meFit$R2)
sum(miFit$R2)
sum(nyFit$R2)
sum(wiFit$R2)

minFit1 <- minQTL %>% mutate(R2 = minFit$R2) %>% select(Trait, Marker, Chrom, Position, Score, Effect, R2)
precipFit1 <- precipQTL %>% mutate(R2 = precipFit$R2) %>% select(Trait, Marker, Chrom, Position, Score, Effect, R2)
minRegFit1 <- minRegQTL %>% mutate(R2 = minRegFit$R2) %>% select(Trait, Marker, Chrom, Position, Score, Effect, R2)
maxRegFit1 <- maxRegQTL %>% mutate(R2 = maxRegFit$R2) %>% select(Trait, Marker, Chrom, Position, Score, Effect, R2)
yearFit1 <- yearQTL %>% mutate(R2 = yearFit$R2) %>% select(Trait, Marker, Chrom, Position, Score, Effect, R2)
latFit1 <- latQTL %>% mutate(R2 = latFit$R2) %>% select(Trait, Marker, Chrom, Position, Score, Effect, R2)
longFit1 <- longQTL %>% mutate(R2 = longFit$R2) %>% select(Trait, Marker, Chrom, Position, Score, Effect, R2)
meFit1 <- meQTL %>% mutate(R2 = meFit$R2) %>% select(Trait, Marker, Chrom, Position, Score, Effect, R2)
miFit1 <- miQTL %>% mutate(R2 = miFit$R2) %>% select(Trait, Marker, Chrom, Position, Score, Effect, R2)
nyFit1 <- nyQTL %>% mutate(R2 = nyFit$R2) %>% select(Trait, Marker, Chrom, Position, Score, Effect, R2)
wiFit1 <- wiQTL %>% mutate(R2 = wiFit$R2) %>% select(Trait, Marker, Chrom, Position, Score, Effect, R2)

rbind(meFit1, miFit1, nyFit1, wiFit1, minFit1, precipFit1, minRegFit1, maxRegFit1, latFit1, longFit1, yearFit1) %>% arrange(Chrom, Position)
