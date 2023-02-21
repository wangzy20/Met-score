library(dplyr)
library(IOBR)

#read data
load("Rdata/raw.Rdata")
sur <- cbind(phe[11:18], expr)

#cox
res <- batch_surv(
  pdata    = sur,
  variable = colnames(sur)[9:ncol(sur)],
  time     = "PFS",
  status   = "PFSstatus"
)

res$FDR <- p.adjust(res$P, method = "fdr")
res$adj.P <-  p.adjust(res$P, method = "BH")
#write.csv(res,"Rdata/raw_res_all.csv")
res <- filter(res, FDR < 0.05)
res_f <- filter(res, FDR < 0.05)
ID <- filter(res, FDR < 0.05) %>% select("ID") %>% unlist()


expr <- expr[, ID]
#multi-cox
dat_cox <- cbind(expr, phe)
head(dat_cox)
## prepare variale for Surv(), paste gene names
multivariate <- paste(sort(colnames(expr)), collapse = '+')

## 1.2 multi-cox----
dat_cox$time <- dat_cox$PFS
dat_cox$event <- dat_cox$PFSstatus
s <-  paste0(' Surv(time, event) ~  ', multivariate)
model <- coxph(as.formula(s), data = dat_cox)

library(survminer)
library(survival)
ggforest(
  model,
  data = dat_cox,
  main = "Hazard ratio",
  cpositions = c(0.06, 0.22, 0.4),
  fontsize = 1.0,
  refLabel = "1",
  noDigits = 4
)

set.seed(1234)
library(My.stepwise)
My.stepwise.coxph(
  Time = 'time',
  Status = 'event',
  variable.list = colnames(expr),
  data = dat_cox
)

#final variance:MEDN0554 MEDP1294 MADP0391 MEDP1710 MEDP1234 MADP0036 MADP0090 MEDP1916 MEDN1406 
#VIF:
#MEDN0554 MEDP1294 MADP0391 MEDP1710 MEDP1234 MADP0036 MADP0090 
#2.276538 1.256070 1.049747 1.054619 1.024008 1.250590 1.623014 
#MEDP1916 MEDN1406 
#3.931309 2.375514 
