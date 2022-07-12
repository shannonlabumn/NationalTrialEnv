setwd("~/NationalTrialEnv")
library(tidyverse); library(GWASpoly)
source(file = "ascer.fn.R")
source(file = "newReadGWASpoly.fn.R")
chipPhenoSelected <- read_csv("chipPhenoSelected.csv")
chipGenoSelected <- read_csv("chipGenoSelected.csv")
wisconsinChips <- chipPhenoSelected %>% 
  select(variety, sel.site, trial) %>%
  mutate(sel.site = ifelse(sel.site == "Wisconsin", 1, 0))
params <- set.params(geno.freq  = 1 - (10/605))
data <- read.GWASpoly.object(ploidy = 4, phenos = wisconsinChips, genos = chipGenoSelected, n.traits = 1, format = "numeric")
data.original <- set.K(data, LOCO=FALSE)
data.original.WI <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
data.LOCO <- set.K(data, LOCO=TRUE)
data.LOCO.scan.WI <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
idahoChips <- chipPhenoSelected %>% 
  select(variety, sel.site, trial) %>%
  mutate(sel.site = ifelse(sel.site == "Idaho", 1, 0))
data <- read.GWASpoly.object(ploidy = 4, phenos = idahoChips, genos = chipGenoSelected, n.traits = 1, format = "numeric")
data.original <- set.K(data, LOCO=FALSE)
data.original.ID <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
data.LOCO <- set.K(data, LOCO=TRUE)
data.LOCO.scan.ID <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
michiganChips <- chipPhenoSelected %>% 
  select(variety, sel.site, trial) %>%
  mutate(sel.site = ifelse(sel.site == "Michigan", 1, 0))
data <- read.GWASpoly.object(ploidy = 4, phenos = michiganChips, genos = chipGenoSelected, n.traits = 1, format = "numeric")
data.original <- set.K(data, LOCO=FALSE)
data.original.MI <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
data.LOCO <- set.K(data, LOCO=TRUE)
data.LOCO.scan.MI <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
northDakotaChips <- chipPhenoSelected %>% 
  select(variety, sel.site, trial) %>%
  mutate(sel.site = ifelse(sel.site == "NorthDakota", 1, 0))
data <- read.GWASpoly.object(ploidy = 4, phenos = northDakotaChips, genos = chipGenoSelected, n.traits = 1, format = "numeric")
data.original <- set.K(data, LOCO=FALSE)
data.original.ND <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
data.LOCO <- set.K(data, LOCO=TRUE)
data.LOCO.scan.ND <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
newYorkChips <- chipPhenoSelected %>% 
  select(variety, sel.site, trial) %>%
  mutate(sel.site = ifelse(sel.site == "NewYork", 1, 0))
data <- read.GWASpoly.object(ploidy = 4, phenos = newYorkChips, genos = chipGenoSelected, n.traits = 1, format = "numeric")
data.original <- set.K(data, LOCO=FALSE)
data.original.NY <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
data.LOCO <- set.K(data, LOCO=TRUE)
data.LOCO.scan.NY <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
coloradoChips <- chipPhenoSelected %>% 
  select(variety, sel.site, trial) %>%
  mutate(sel.site = ifelse(sel.site == "Colorado", 1, 0))
data <- read.GWASpoly.object(ploidy = 4, phenos = coloradoChips, genos = chipGenoSelected, n.traits = 1, format = "numeric")
data.original <- set.K(data, LOCO=FALSE)
data.original.CO <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
data.LOCO <- set.K(data, LOCO=TRUE)
data.LOCO.scan.CO <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
maineChips <- chipPhenoSelected %>% 
  select(variety, sel.site, trial) %>%
  mutate(sel.site = ifelse(sel.site == "Maine", 1, 0))
data <- read.GWASpoly.object(ploidy = 4, phenos = maineChips, genos = chipGenoSelected, n.traits = 1, format = "numeric")
data.original <- set.K(data, LOCO=FALSE)
data.original.MA <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
data.LOCO <- set.K(data, LOCO=TRUE)
data.LOCO.scan.MA <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
texasChips <- chipPhenoSelected %>% 
  select(variety, sel.site, trial) %>%
  mutate(sel.site = ifelse(sel.site == "Texas", 1, 0))
data <- read.GWASpoly.object(ploidy = 4, phenos = texasChips, genos = chipGenoSelected, n.traits = 1, format = "numeric")
data.original <- set.K(data, LOCO=FALSE)
data.original.TX <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
data.LOCO <- set.K(data, LOCO=TRUE)
data.LOCO.scan.TX <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
northCarolinaChips <- chipPhenoSelected %>% 
  select(variety, sel.site, trial) %>%
  mutate(sel.site = ifelse(sel.site == "North Carolina", 1, 0))
data <- read.GWASpoly.object(ploidy = 4, phenos = northCarolinaChips, genos = chipGenoSelected, n.traits = 1, format = "numeric")
data.original <- set.K(data, LOCO=FALSE)
params <- set.params(geno.freq  = 1 - (10/605))
data.original.NC <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
data.LOCO <- set.K(data, LOCO=TRUE)
data.LOCO.scan.NC <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
oregonChips <- chipPhenoSelected %>% 
  select(variety, sel.site, trial) %>%
  mutate(sel.site = ifelse(sel.site == "Oregon", 1, 0))
data <- read.GWASpoly.object(ploidy = 4, phenos = oregonChips, genos = chipGenoSelected, n.traits = 1, format = "numeric")
data.original <- set.K(data, LOCO=FALSE)
data.original.OR <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
data.LOCO <- set.K(data, LOCO=TRUE)
data.LOCO.scan.OR <- GWASpoly(data.original,models=c("additive", "1-dom", "2-dom"),traits=c("sel.site"),params=params)
rm(list = c("data", "data.LOCO", "chipPhenoSelected", "chipGenoSelected", "wisconsinChips", "idahoChips", "michiganChips", "northDakotaChips", "newYorkChips", "coloradoChips", "maineChips", "texasChips", "northCarolinaChips", "oregonChips"))
save(list = ls(all.names = TRUE), file = "selsiteGWAS.RData")