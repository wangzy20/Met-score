rm(list=ls())
gc()

library(data.table)
library(dplyr)
library(readxl)
library(IOBR)

expr <- readxl::read_excel("data/expr data.xlsx") %>%
  column_to_rownames("Index") %>%
  select(-c(1,2,3)) %>%
  t() %>%
  data.frame()
  
phe <- read_excel("data/metadata.xlsx")

#corrdinate the data
phe <- phe %>% 
  mutate(sample_id = paste0("X",phe$标本号))

phe <- phe[phe$sample_id %in% colnames(expr),]
phe <- data.frame(phe) %>%
  column_to_rownames('sample_id')
table(rownames(phe) %in% colnames(expr)) == 320

#transform
expr <- t(expr) %>% as.data.frame()
expr <- expr[rownames(phe),]
identical(rownames(expr), rownames(phe))
expr=as.data.frame(lapply(expr,as.numeric))

#PCA
library(tidydr)
x <- dr(data = expr, fun = prcomp)
autoplot(x)

#scale
expr=apply(expr, 2,scale) %>% data.frame
rownames(expr) <- rownames(phe)

save(expr, phe , file  = "Rdata/raw.Rdata")
