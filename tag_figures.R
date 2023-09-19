## IMPORT DATA FILES FROM NCPT_envGWAS.R

gxe_plot <- field_data_all %>% 
  filter(!is.na(yield_total), yield_total>0) %>%
  ggplot(aes(x = fct_reorder(location, yield_total, median), y = yield_total)) + 
  geom_boxplot() + 
  theme_cowplot(12) +
  xlab("Trial Location") +
  ylab(expression(Yield~(kg)))

ggsave("gxe_plot.png", plot = gxe_plot, units = "in", width = 6.5, height = 4)

precip.named <- "Precipitation"
names(precip.named) <- "precip"

precipMan <- manhattan.plot(selThresh, trait = "precip", models = "additive") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black")) + 
  facet_grid(~trait, labeller = labeller(trait = precip.named))

minMan <- manhattan.plot(selThresh, trait = "minTemp", models = "additive") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black"))

maxMan <- manhattan.plot(selThresh, trait = "maxTemp", models = "additive") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black"))

precip.facet <- "Precipitation Regression"
names(precip.facet) <- "precip"

precipYldMan <- manhattan.plot(slopeThresh, traits = "precip", models = "additive") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black")) + 
  facet_grid(~trait, labeller = labeller(trait = precip.facet))

minTemp.facet <- "minTemp Regression"
names(minTemp.facet) <- "minTemp"

minYldMan <- manhattan.plot(slopeThresh, traits = "minTemp", models = "additive") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black")) + facet_grid(~trait, labeller = labeller(trait = minTemp.facet))

maxTemp.facet <- "maxTemp Regression"
names(maxTemp.facet) <- "maxTemp"

maxYldMan <- manhattan.plot(slopeThresh, traits = "maxTemp", models = "additive") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black")) + facet_grid(~trait, labeller = labeller(trait = maxTemp.facet))

yearMan <- manhattan.plot(selThresh, models = "additive", traits = "year") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black"))

trialEnvs <- envTrials %>% rename(Precipitation = "precip")

trial.cor <- ggcorrplot(cor(select(trialEnvs, -trial)), hc.order =TRUE,type ="lower",
                        outline.color ="white", lab =TRUE) + theme_cowplot(12) + xlab("") + ylab("") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

selSiteRawEnvs <- selSiteEnvs %>% rename(Precipitation = "precip") %>% select(Precipitation, minTemp, maxTemp)

sel.site.cor <- ggcorrplot(cor(selSiteRawEnvs), hc.order =TRUE,type ="lower", 
                           outline.color ="white", lab =TRUE) + theme_cowplot(12) + xlab("") + ylab("") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

envGWASplots <- plot_grid(minMan, maxMan, precipMan, sel.site.cor, labels = c("A", "B", "C", "D"), ncol = 2)
slopeGWASplots <- plot_grid(minYldMan, maxYldMan, precipYldMan, trial.cor, labels = c("A", "B", "C", "D"), ncol = 2)

ggsave("envGWASplots.png", envGWASplots, width = 7.5, height = 4, units = "in")
ggsave("slopeGWASplots.png", slopeGWASplots, width = 7.5, height = 4, units = "in")
ggsave("yearMan.png", yearMan, width = 3.25, height = 3, units = "in")

## Supplementary Figures

me.facet <- "Maine"
names(me.facet) <- "ME"

meMan<- manhattan.plot(stateThresh, models = "additive", traits = "ME") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black")) + facet_grid(~trait, labeller = labeller(trait = me.facet))

mi.facet <- "Michigan"
names(mi.facet) <- "MI"

miMan<- manhattan.plot(stateThresh, models = "additive", traits = "MI") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black")) + facet_grid(~trait, labeller = labeller(trait = mi.facet))

ny.facet <- "New York"
names(ny.facet) <- "NY"

nyMan<- manhattan.plot(stateThresh, models = "additive", traits = "NY") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black")) + facet_grid(~trait, labeller = labeller(trait = ny.facet))

wi.facet <- "Wisconsin"
names(wi.facet) <- "WI"

wiMan <- manhattan.plot(stateThresh, models = "additive", traits = "WI") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black")) + 
  facet_grid(~trait, labeller = labeller(trait = wi.facet))


stateGWASplots <- plot_grid(meMan, miMan, nyMan, wiMan, 
                            labels = c("A", "B", "C", "D"), ncol = 2)

lat.facet <- "Latitude"
names(lat.facet) <- "lat"

long.facet <- "Longitude"
names(long.facet) <- "long"

latMan <- manhattan.plot(gpsThresh, models = "additive", traits = "lat") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black")) + 
  facet_grid(~trait, labeller = labeller(trait = lat.facet))

longMan <- manhattan.plot(gpsThresh, models = "additive", traits = "long") + 
  theme_cowplot(12) + 
  theme(legend.position = "none", strip.background = element_rect(linewidth = 1, color = "black")) + 
  facet_grid(~trait, labeller = labeller(trait = long.facet))

gpsGWASplots <- plot_grid(latMan, longMan, labels = c("A", "B"), ncol = 2)

ggsave("stateGWASplots.png", plot = stateGWASplots, units = "in", width = 6.5, height = 9)

ggsave("gpsGWASplots.png", plot = gpsGWASplots, units = "in", width = 6.5, height = 4)

t.geno <- t(select(geno.sel, -marker, -chrom, -bp))
colnames(t.geno) <- geno.sel$marker
pca <- prcomp(t.geno)

pca.plot <- ggbiplot::ggbiplot(pca, var.axes = F, labels.size = 0, groups = selSites$sel.site, ellipse = T, size = 5) + 
  labs(colour = "Breeding program") + theme_cowplot(12)


ggsave("/Users/husainagha/Documents/Experimental/NationalTrialEnv/pca_plot.png", pca.plot, device = "png", width = 7.5, height = 8, units = "in")


