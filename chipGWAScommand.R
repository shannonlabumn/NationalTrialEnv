setwd("~/NationalEnv")
library(tidyverse); library(GWASpoly)

data <- read.GWASpoly(ploidy = 4, 
                      pheno.file = "chipPheno12k.csv", 
                      geno.file = "chipGeno12k.csv", 
                      format = "numeric", 
                      n.traits = 1)

data.original <- set.K(data, 
                       LOCO=FALSE,
                       n.core = 16)

params <- set.params(geno.freq  = 1 - 10/334,
                     fixed      = c("program", "trial"),
                     fixed.type = c("factor", "factor"))
Sys.time()

data.original.scan <- GWASpoly(data.original,
                               models=c("additive"),
                               traits=c("yield_total"),
                               params=params,
                               n.core = 16)

Sys.time()

qq.plot(data.original.scan,trait="yield_total") + ggtitle(label="Original")
data2 <- set.threshold(data.original.scan,method="Bonferroni",level=0.05)

p <- manhattan.plot(data2,traits="yield_total") +
      theme(axis.text.x = element_text(angle=90,vjust=0.5))

ggsave("chipManhattan.png", 
       p, 
       width = 583,
       height = 483,
       units = "px")
