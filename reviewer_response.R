sigMarkers <- c("ST4.03ch06_57692614", "ST4.03ch09_411487", "solcap_snp_c1_11105",
                "ST4.03ch02_14148961", "ST4.03ch06_57244824", "ST4.03ch08_1934839")


geno.t <- filter(geno.sel, marker %in% sigMarkers) %>% t()

geno.df <- as.data.frame(geno.t)
names(geno.df) <- sigMarkers
geno.df$variety <- rownames(geno.df)
geno.df <- geno.df[-c(1:3),]
markerYear <- left_join(selSites, geno.df)
selSiteMarkerYear <- markerYear %>%
  mutate(across(ST4.03ch06_57692614:ST4.03ch08_1934839, as.numeric)) %>%
  unite("selYear", sel.site:year, remove = F) %>%
  group_by(sel.site, year, selYear) %>%
  summarise(across(ST4.03ch06_57692614:ST4.03ch08_1934839, ~mean(.x))) %>%
  arrange(sel.site, year)

envs.sel.avg <- envs.3.years %>%
  arrange(sel.site, year) %>%
  unite("selYear", sel.site:year, remove = F) %>%
  filter(selYear %in% selSiteMarkerYear$selYear)

cor(selSiteMarkerYear$ST4.03ch06_57692614, envs.sel.avg$minTemp)
cor(selSiteMarkerYear$ST4.03ch09_411487, envs.sel.avg$minTemp)
cor(selSiteMarkerYear$solcap_snp_c1_11105, envs.sel.avg$minTemp)
cor(selSiteMarkerYear$ST4.03ch02_14148961, envs.sel.avg$precip)
cor(selSiteMarkerYear$ST4.03ch06_57244824, envs.sel.avg$precip)
cor(selSiteMarkerYear$ST4.03ch08_1934839, envs.sel.avg$precip)

states <- c("CO", "ID", "ME", "MI", "NC",
            "ND", "NY", "OR", "TX", "WI")

lat <- c(37.71, 42.95, 46.65, 43.35, 35.87,
          48.53, 42.43, 45.81, 36.08, 44.12)

long <- c(-106.14, -112.83, -68.01, -85.18, -76.65,
          -97.62, -76.39, -119.28, -102.60, -89.54)

elevation <- c(2339, 1341.0, 138.6, 281.9, 5.1,
               272, 309.4, 140, 1213, 331)

cor(elevation, lat)
cor(elevation, long)

envs.3.years %>%
  group_by(sel.site) %>%
  summarise(precip = mean(precip),
            minTemp = mean(minTemp),
            maxTemp = mean(maxTemp)) %>%
  arrange(sel.site) %>%
  mutate(elevation = elevation) %>%
  summarise(precip.long = cor(precip, long),
            minTemp.long = cor(minTemp, long),
            maxTemp.long = cor(maxTemp, long),
            elevation.long = cor(elevation, long),
            precip.lat = cor(precip, lat),
            minTemp.lat = cor(minTemp, lat),
            maxTemp.lat = cor(maxTemp, lat),
            elevation.lat = cor(elevation, lat))


soil.ph <- c(8, 7.7, 5.1, 5.8, 5.5, 7.9, 5.3, 7.2, 7.4, 5.2)
cor(soil.ph, long)

field_data1 <- mutate(field_data_all, year = as.numeric(year)) %>%
  arrange(yield_total) %>%
  group_by(variety, location, year) %>%
  summarise(yield_total = mean(yield_total, na.rm = T))

mgPerHa <- read_csv(file = "~/Documents/GitHub.nosync/NationalTrialEnv/mgPerHa.csv") %>%
  mutate(variety = str_replace(variety, "-", "."),
         variety = str_replace(variety, " ", "."),
         variety = str_replace(variety, "/", ".")) %>%
  group_by(variety, location, year) %>%
  summarise(yield.standard = mean(yield.standard, na.rm = T))

field_data1 %>%
  right_join(mgPerHa, by = c("variety", "year", "location")) %>%
  mutate(sqMeters = round((yield_total/yield.standard*10), digits = 1)) %>%
  view()

field_data1 %>%
  left_join(mgPerHa, by = c("variety", "year", "location")) %>%
  mutate(sqMeters = round((yield_total/yield.standard*10), digits = 1)) %>%
  filter(!is.na(sqMeters)) %>%
  group_by(location, year) %>%
  summarise(sqMeters = round(mean(sqMeters, na.rm = T), digits = 1),
            n = n()) %>%
  mutate(sqMeters = ifelse(location == "ND", 3.9, sqMeters)) %>% view()

field_data1 %>%
  left_join(mgPerHa, by = c("variety", "year", "location")) %>%
  mutate(sqMeters = round((yield_total/yield.standard*10), digits = 1)) %>%
  filter(!is.na(sqMeters),
         location == "CA") %>%
  view()

field_data_all %>%
  group_by(year, location) %>%
  filter(yield_total > 0 , location == "FL") %>%
  ggplot(aes(x = yield_total, group = year, fill = factor(year), color= factor(year))) +
  ggridges::geom_density_ridges(aes(y = factor(year)), show.legend = F)

field_data_all %>%
  group_by(year, location) %>%
  filter(yield_total > 0 , location == "OR") %>%
  ggplot(aes(x = yield_total, group = year, fill = factor(year), color= factor(year))) +
  ggridges::geom_density_ridges(aes(y = factor(year)), show.legend = F)

field_data_all %>%
  group_by(year, location) %>%
  filter(yield_total > 0 , location == "NC") %>%
  ggplot(aes(x = yield_total, group = year, fill = factor(year), color= factor(year))) +
  ggridges::geom_density_ridges(aes(y = factor(year)), show.legend = F)

field_data_all %>%
  group_by(year, location) %>%
  filter(yield_total > 0 , location == "TX") %>%
  ggplot(aes(x = yield_total, group = year, fill = factor(year), color= factor(year))) +
  ggridges::geom_density_ridges(aes(y = factor(year)), show.legend = F)

chipFieldData %>%
  filter(!is.na(variety), variety != "Trial Information", variety != "Variety") %>%
  select(variety, trial, yield_total, yield_total_cwt) %>%
  mutate(variety = str_replace(variety, "-", "."),
         variety = str_replace(variety, " ", "."),
         variety = str_replace(variety, "/", "."),
         yield_total = as.numeric(yield_total),
         yield_total_cwt = as.numeric(yield_total_cwt)) %>%
  separate(trial, c("year", NA, "location"), " ") %>%
  filter(yield_total > 0,
         !is.na(yield_total_cwt),
         yield_total != 76.2) %>%
  mutate(m2 = round((yield_total/(yield_total_cwt/7.966))*10, digits = 2)) %>%
  filter(!is.na(m2)) %>%
  group_by(year, location) %>%
  summarise(yield_max = max(yield_total),
            yield_max_cwt = max(yield_total_cwt),
            mean = median(m2),
            sd = sd(m2),
            min = min(m2),
            max = max(m2)) %>%
  select(location, year, mean, sd, min, max) %>%
  arrange(location, year) %>% view()

chipFieldData %>%
  filter(!is.na(variety), variety != "Trial Information", variety != "Variety") %>%
  select(variety, trial, yield_total, yield_total_cwt) %>%
  mutate(variety = str_replace(variety, "-", "."),
         variety = str_replace(variety, " ", "."),
         variety = str_replace(variety, "/", "."),
         yield_total = as.numeric(yield_total),
         yield_total_cwt = as.numeric(yield_total_cwt)) %>%
  separate(trial, c("year", NA, "location"), " ") %>%
  filter(yield_total > 0, !is.na(yield_total_cwt)) %>%
  mutate(m2 = round((yield_total/(yield_total_cwt/7.966))*10, digits = 2)) %>%
  filter(location == "MI") %>% arrange(year, m2) %>% view()

field_data1 <- field_data_all %>%
  mutate(year = as.numeric(year),
         yield_total = ifelse(location == "CA" & year %in% c(2014, 2015, 2017), yield_total/2.2, yield_total),
         yield_total = ifelse(location == "FL" & year == 2010, yield_total/2.2, yield_total),
         ) %>%
  filter(location != "MD")

OR_corrected <- read_csv("phenotypes/OR_2020_corrected_phenotypes.csv", show_col_types = F) %>%
  mutate(year = 2020) %>%
  select(variety, location, year, yield_total) %>%
  group_by(variety) %>% arrange(yield_total) %>%
  mutate(rep = 1:n())

field_data2 <- field_data1 %>% group_by(variety, location, year) %>%
  arrange(yield_total) %>% mutate(rep = 1:n()) %>%
  ungroup() %>%
  rows_update(OR_corrected, by = c("location", "year", "variety", "rep"),
                            unmatched = "ignore")

sqMeters <- tibble(location = c("CA", "FL", "MI", "NC", "TX", "ND", "OR", "WI", "NY", "MO"),
                   area =      c(2.5,  2.8,  3.3, 3.7,  2.4,  4.4,  2.8,  4.2,  2.8,  3.3))

field_data3 <- field_data2 %>%
  left_join(sqMeters, by = "location") %>%
  mutate(yield_per_ha = (10*yield_total)/area)

yield_per_ha <- field_data3 %>%
  filter(!is.na(yield_total), yield_total>0,
         year %in% 2014:2022,
         location != "MD") %>%
  ggplot(aes(x = fct_reorder(location, yield_per_ha, median), y = yield_per_ha)) +
  geom_boxplot() +
  theme_cowplot(12) +
  xlab("Trial Location") +
  ylab(expression(Yield~(tonnes%*%ha^-1)))

yield_per_plot <- field_data3 %>%
  filter(!is.na(yield_total), yield_total>0,
         year %in% 2014:2020) %>%
  ggplot(aes(x = fct_reorder(location, yield_per_ha, median), y = yield_total)) +
  geom_boxplot() +
  theme_cowplot(12) +
  xlab("Trial Location") +
  ylab(expression(Yield~(kg %*% plot^-1)))

both_yield_plots <- ggpubr::ggarrange(yield_per_ha, yield_per_plot, labels = c("A", "B"))

ggsave("both_yield_plots.png", plot = both_yield_plots, units = "in", width = 6.5, height = 4)
