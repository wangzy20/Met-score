library(dplyr)
library(survival)
library(survminer)
library(timeROC)
library(RColorBrewer)
library(ggstatsplot)
#读取风险分数数据
trainscore <- read.csv("data/train_riskscore.csv")
testscore <- read.csv("data/test_riskscore .csv.csv")
phe_train <- readxl::read_excel("data/train dataset.xlsx") %>% data.frame()
phe_train$sampleid <- paste0("X",phe_train$标本号)
rownames(phe_train) <- phe_train$sampleid

phe_train <- phe_train[trainscore$X,]
identical(rownames(phe_train),trainscore$X)
phe_train$riskscore <- trainscore$Risk_score

#1 Met-score and different variance 
dat_cox <- phe_train[,-c(1:3,9:14,27,28)]
colnames(dat_cox)[c(6,7)] <- c("event","time")
multivariate <- paste(colnames(dat_cox)[-c(6,7,1,2)], collapse = '+')
s <-  paste0(' Surv(time, event) ~  ', multivariate )
model <- coxph(as.formula(s), data = dat_cox )
ggforest(model, data = dat_cox, 
         main = "Hazard ratio", 
         cpositions = c(0.06, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)


#2.tradition model
model_clinic <- coxph(Surv(time, event) ~ STAGE + sex + EBVDNA4000为cutoff值, data = dat_cox)
riskscore_train <- predict(model_clinic, type = 'risk', data = dat_cox)
dat_cox$riskscore_train <- riskscore_train



ggforest(model_clinic, data = dat_cox, 
         main = "Hazard ratio", 
         cpositions = c(0.06, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)


dat_cox$time <- dat_cox$time/12
## need 3 cols，time、event and risk scores
result1 <- with(dat_cox, timeROC(T=time,
                                delta=event,
                                marker=riskscore_train,
                                cause=1,
                                times = c(1,3,5),
                                iid = TRUE))
#identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
data1 = data.frame(fpr = as.numeric(result1$FP),
                  tpr = as.numeric(result1$TP),
                  time = rep(as.factor(c(1,3,5)),each = nrow(result1$TP)))
#3.mixed model
model_conbine <- coxph(Surv(time, event) ~ STAGE+ sex+ EBVDNA4000为cutoff值 + riskscore, data = dat_cox)

riskscore2_train <- predict(model_conbine, type = 'risk', data = dat_cox)
dat_cox$riskscore2_train <- riskscore2_train

result2 <- with(dat_cox, timeROC(T=time,
                                delta=event,
                                marker=riskscore2_train,
                                cause=1,
                                times = c(1,3,5),
                                iid = TRUE))
#identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
data2 = data.frame(fpr = as.numeric(result2$FP),
                  tpr = as.numeric(result2$TP),
                  time = rep(as.factor(c(1,3,5)),each = nrow(result2$TP)))


#4.test 
dat_cox2 <- readxl::read_excel("data/test dataset.xlsx") %>% data.frame()
rownames(dat_cox2) <- paste0("X",dat_cox2$标本号)
dat_cox2 <- dat_cox2[testscore$X,]
dat_cox2 <- dat_cox2[,-c(1:3,9:14,27)]
dat_cox2$riskscore <- testscore$Risk_score
colnames(dat_cox2)[c(6,7)] <- c("event","time")

riskscore_test <- predict(model_clinic, type = "risk", newdata = dat_cox2)
dat_cox2$riskscore_test <- riskscore_test

ggforest(model_clinic, data = dat_cox2, 
         main = "Hazard ratio", 
         cpositions = c(0.06, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)



dat_cox2$time <- dat_cox2$time/12
result3 <- with(dat_cox2, timeROC(T=time,
                                delta=event,
                                marker=riskscore_test,
                                cause=1,
                                times = c(1,3,5),
                                iid = TRUE))
data3 = data.frame(fpr = as.numeric(result3$FP),
                   tpr = as.numeric(result3$TP),
                   time = rep(as.factor(c(1,3,5)),each = nrow(result3$TP)))

riskscore2_test <- predict(model_conbine, type = "risk", newdata = dat_cox2)
dat_cox2$riskscore2_test <- riskscore2_test
result4 <- with(dat_cox2, timeROC(T=time,
                                 delta=event,
                                 marker=riskscore2_test,
                                 cause=1,
                                 times = c(1,3,5),
                                 iid = TRUE))

data4 = data.frame(fpr = as.numeric(result4$FP),
                   tpr = as.numeric(result4$TP),
                   time = rep(as.factor(c(1,3,5)),each = nrow(result4$TP)))

data1$model <- "train_clinic"
data2$model <- "train_combine"
data3$model <- "test_clinic"
data4$model <- "test_combine"

data <- rbind(data1,data2)
data <- rbind(data3,data4)
data_5 <- dplyr::filter(data,time == 5)


model_names = c("train_clinic","train_combine")
model_names = c('test_clinic',"test_combine")
auc <- c(result1$AUC[3],result2$AUC[3])
auc<- c(result3$AUC[3],result4$AUC[3])

ggplot() + 
  geom_line(data = data_5, 
              aes(x = fpr, y = tpr, color = model), 
              size = 1,
              se = FALSE) + 
  scale_color_manual(name = NULL,
                     values = brewer.pal(4,"Set1"),
                     labels = paste0("AUC of ",model_names ," model : ",
                                     format(round(auc,2),nsmall = 2)))+
  geom_line(aes(x = c(0, 1), y = c(0,1)), 
            color = "grey",
            linetype = 'dotdash')+
  theme_ggstatsplot()+
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.line = element_line(linetype = 1),
        panel.grid = element_blank(),
        legend.background = element_rect(linetype = 2, size = 0.2, colour = "black"),
        legend.position = c(0.665,0.135))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()

x <- compare(result3,result4,adjusted = T)
round(x$p_values_AUC,5)

