rm(list = ls())
load("Rdata/raw.Rdata")

library(data.table)
library(dplyr)
library(survminer)
library(survival)
library(loose.rock)
library(futile.logger)
library(glmSparseNet)
library(ggrisk)
library(pheatmap)
library(ggplot2)
library(ggsci)
library(ggstatsplot)
library(rms)
library(IOBR)

phe_new <- select(phe, sample_id,PFS, PFSstatus)
colnames(phe_new) <- c("sample_id","time", 'event')
#seed = 123
train_test <- SplitTrainTest(x = expr, y = phe_new, train_ratio = .7, 
                             type = "survival", seed = 123)
train.x = train_test$train.x
train.y <- train_test$train.y
test.x = train_test$test.x
test.y <- train_test$test.y


build_cox_model <- function(xdata,ydata,pro){
  n <- apply(xdata,2,scale) %>%
    as.data.frame()
  rownames(n) <- rownames(xdata)
  boxplot(n)
  dat_cox <- cbind(ydata,n)
  head(dat_cox)
  ## prepare variale for Surv(), paste gene names
  
  formula = Surv(time, event) ~ MEDN0554 + MEDP1710 + MEDP1294 + 
    MEDN1224 + MEDN1475 + MEDP1234 + MEDP1962 + MADP0090 + MEDN1406
  attach(dat_cox)
  model <- coxph(formula = formula, data = dat_cox )
  summary(model, data = dat_cox)
  risk_score <- predict(model, type = 'risk', data = dat_cox)
  dat_cox$Risk_score <- risk_score 
  
  # 2.  visualization----
  ## 2.1 forest----
  library(survminer)
  options(scipen=1)
  ggforest(model, data = dat_cox, 
           main = "Hazard ratio", 
           cpositions = c(0.06, 0.22, 0.4), 
           fontsize = 1.0, 
           refLabel = "1", noDigits = 4)
  ggsave(paste0(pro,'_multicox_forest.pdf'), height = 5, width = 13)
  
  library(timeROC)
  attach(dat_cox)
  head(dat_cox)
  new_dat <- dat_cox[, c('event', 'time','Risk_score')]
  new_dat$time = new_dat$time/12
  ## need 3 cols，time、event and risk scores
  result <- with(new_dat, timeROC(T=time,
                                  delta=event,
                                  marker=Risk_score,
                                  cause=1,
                                  times = c(1, 3, 5),
                                  iid = TRUE))
  #identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
  data = data.frame(fpr = as.numeric(result$FP),
                    tpr = as.numeric(result$TP),
                    time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))
  
  ggplot() + 
    geom_smooth(data = data, 
                aes(x = fpr, y = tpr, color = time), 
                size = 1,
                method = "loess",
                se = FALSE) + 
    scale_color_manual(name = NULL,
                       values = c("#92C5DE", "#F4A582", "#66C2A5"),
                       labels = paste0("AUC of ",c(1,3,5),"-year survival: ",
                                       format(round(result$AUC,2),nsmall = 2)))+
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
  ggsave(paste0(pro,'_multicox_smooth_ROC.pdf'),width = 8,height = 8)
  
  ## 3.4 km----
  attach(new_dat)
  head(new_dat)
  risk_level <- as.factor(ifelse(new_dat$Risk_score > median(new_dat$Risk_score),"High","Low"))
  new_dat$Risk_level <- risk_level

  sfit <- survfit(Surv(time, event)~Risk_level, data=new_dat)
  sfit
  summary(sfit)
  #注意，这里event要用数值，而且TRUE不能写成T
  ## more complicate figures.
  survp=ggsurvplot(
    sfit,                     # survfit object with calculated statistics.
    legend.title = 'Risk level', 
    legend = "top",#图例位置
    legend.labs = c('High', 'Low'),
    pval = TRUE, #在图上添加log rank检验的p值
    risk.table = TRUE, 
    risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
    xlab = "Time in years", #x轴标题
    # xlim = c(0, 10), #展示x轴的范围
    size = 1.5, #线条大小
    ggtheme = theme_ggstatsplot(),
    palette="nejm", #配色
    data = new_dat
  )
  pdf(paste0(pro,'_multicox_KM.pdf'), onefile = F)
  print(survp)
  dev.off()
  
  save(new_dat,dat_cox,file = paste0(pro, "model_summary.Rdata"))
  return(model)
  
  fit <-  cph(formula, dat_cox[,-1])
  library(ggrisk)
  # save in pdf----
  pdf(paste0(pro,'_riskscore.pdf'), onefile = F)
  ggrisk(fit,
         color.A = c(low = "#0B5E9D", high = "#EA5205"),
         color.B = c(code.0 = "#0B5E9D", code.1 = "#EA5205"),
         color.C = c(low = "#0B5E9D", median = "white", high = "#EA5205"))
  dev.off() 
  save(model, dat_cox, new_dat, file = paste0(pro,'_multicox_model.Rdata'))
}

model <- build_cox_model(xdata = train.x, ydata = train.y,pro = "train_123")

x <- read.csv("data/test_riskscore.csv")
x <-x[match(rownames(xdata),x$X),]

check_model<- function(xdata, ydata, model, pro){
  n2 <- apply(xdata,2,scale) %>%
    as.data.frame()
  rownames(n2) <- rownames(xdata)
  boxplot(n2)
  dat_cox2 <- cbind(x,n2)
  head(dat_cox2)
  risk_score2 <- predict(model, type = 'risk', newdata = dat_cox2)
  dat_cox2$Risk_score <- x$Risk_score
  
  library(timeROC)
  new_dat2 <- dat_cox2[, c('event', 'time','Risk_score')]
  new_dat2$time <- new_dat2$time/12
  result <- with(new_dat2, timeROC(T=time,
                                   delta=event,
                                   marker=Risk_score,
                                   cause=1,
                                   times = c(1, 3, 5),
                                   iid = TRUE))
  #identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
  data = data.frame(fpr = as.numeric(result$FP),
                    tpr = as.numeric(result$TP),
                    time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))
  
  ## 将geom_line()改为geom_smooth(method = "loess")
  ## 在数学中有种算法叫“样条插补法”，这种方法可以获得过点的平滑曲线
  ggplot() + 
    geom_smooth(data = data, 
                aes(x = fpr, y = tpr, color = time), 
                size = 1,
                method = "loess",
                se = FALSE) + 
    scale_color_manual(name = NULL,
                       values = c("#92C5DE", "#F4A582", "#66C2A5"),
                       labels = paste0("AUC of ",c(1,3,5),"-year survival: ",
                                       format(round(result$AUC,2),nsmall = 2)))+
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
  ggsave(paste0(pro,'_multicox_smooth_ROC.pdf'),width = 8,height = 8)
  
  risk_level <- as.factor(ifelse(new_dat2$Risk_score > median(new_dat2$Risk_score),"High","Low"))
  new_dat2$Risk_level <- risk_level
  
  
  sfit <- survfit(Surv(time, event)~Risk_level, data=new_dat2)
  sfit
  summary(sfit)
  #注意，这里event要用数值，而且TRUE不能写成T
  ## more complicate figures.
  survp=ggsurvplot(
    sfit,                     # survfit object with calculated statistics.
    legend.title = 'Risk level', 
    legend = "top",#图例位置
    legend.labs = c('High', 'Low'),
    pval = TRUE, #在图上添加log rank检验的p值
    risk.table = TRUE, 
    risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
    xlab = "Time in years", #x轴标题
    # xlim = c(0, 10), #展示x轴的范围
    size = 1.5, #线条大小
    ggtheme = theme_ggstatsplot(),
    palette="nejm", #配色
    data = new_dat2
  )
  pdf(paste0(pro,'_multicox_KM.pdf'), onefile = F)
  print(survp)
  dev.off()
  
  save(new_dat2,dat_cox2,model,file = paste0(pro, "model_summary.Rdata"))
}
check_model(test.x,test.y,model, pro = "test_123")
