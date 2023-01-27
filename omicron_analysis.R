setwd("C:/marta - covid/")

## color order - red;yellow;lightblue;purple;aqua;orange;darkblue;magenta
coul=c("#E8421D","#FFD423","#11B9E8","#981FFF","#4AF072","#FC6B13","#2B7BD6","#EC2DF2")

library(adegenet)
library(survival)
library(MASS)
library(VGAM)
library(rms)
library(coxme)
library(mice)
library(VIM)
library(lattice)
library(psych)
library(readxl)
library(dplyr)
library(survminer)
library(rpart)
library(party)
library(rpart.plot)
library(caret)
library(ggpubr)
library(rstatix)

data.sero <- read_excel("HCW_cases_omi_ricardo_17jan2022_2.xlsx", sheet = "test")

data <- data.sero %>%
  mutate(recent_vacc = as.factor(recent_vacc), prior_exposure = as.factor(prior_exposure), symp_visit=as.factor(symp_visit), booster = as.factor(booster)) %>%
  mutate(ever_positive = as.factor(ever_positive), serology_result = as.factor(serology_result), vaccine_new = as.factor(vaccine_new)) %>%
  as.data.frame()

####################################    ALL SAMPLES    ##################################################################################

data.f<-data[-which(is.na(data$serology_result)),]
data.f<-data.f[,-10]
data.f$time_since_inf[data.f$time_since_inf>4000]<-5000

training.samples <- data.f$ever_positive %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- data.f[training.samples, ]
test.data <- data.f[-training.samples, ]
model <- train(
  ever_positive ~., data = train.data, method = "ctree2",
  trControl = trainControl("cv", number = 1000),
  tuneGrid = expand.grid(maxdepth = 5, mincriterion = 0.9)
)
plot(model$finalModel)

# Make predictions on the test data
predicted.classes <- model %>% predict(test.data)
# Compute model accuracy rate on test data
mean(predicted.classes == test.data$ever_positive)


####################################     SYMPTOMATIC ONLY     ##################################################################################
# data.f2<-data[-which(is.na(data$serology_result) | data$symp_visit!=1 | data$vaccine_new==2),]
data.f2<-data[-which(is.na(data$serology_result) | data$symp_visit!=1 ) ,]
# data.f2<-data.f2[-which(data.f2$ever_positive==1 & data.f2$symptoms_less_2==0),]
data.f2<-data.f2[,c(2:10,15)]
# data.f2$vaccine_new[data.f2$vaccine_new!=1]<-0
# data.f2$time_since_inf[data.f2$time_since_inf>4000]<-5000
colnames(data.f2)[2]<-"Previous_SARS_CoV2_NAAT_confirmed_infection"
colnames(data.f2)[7]<-"Anti_spike_IgG_binding_antibody_units"

training.samples2 <- data.f2$ever_positive %>%
  createDataPartition(p = 0.9, list = FALSE)
train.data2  <- data.f2[training.samples2, ]
test.data2 <- data.f2[-training.samples2, ]
control = ctree_control(mincriterion = 1-(.05*ncol(train.data2)-1))
model2 <- train(
  ever_positive ~ Previous_SARS_CoV2_NAAT_confirmed_infection+recent_vacc+booster+Anti_spike_IgG_binding_antibody_units, data = train.data2, method = "ctree2",
  trControl = trainControl(method="repeatedcv", repeats=25),
  tuneGrid = expand.grid(maxdepth = 5, mincriterion = 0.75)
)
plot(model2$finalModel,inner_panel=innerWeights,
          terminal_panel=node_barplot2,
          tp_args = list(ylines = c(2, 4))) # this arg. modifies the spacing between barplots
p1<-grid.grab()
p1


# Make predictions on the test data
predicted.classes2 <- model2 %>% predict(test.data2)
# Compute model accuracy rate on test data
mean(predicted.classes2 == test.data2$ever_positive)
table(predicted.classes2,test.data2$ever_positive)

library(scales)
options(scipen=5)
p2<-ggplot(data.f2, aes(x=Anti_spike_IgG_binding_antibody_units, color=ever_positive)) +
  geom_density(size=2, show.legend = F)+
  stat_density(aes(x=Anti_spike_IgG_binding_antibody_units, color=ever_positive),geom="line",position="identity",size=2)+
  xlab("Anti-spike IgG binding antibody units")+ylab("Density")+
  scale_x_continuous(label=comma)
p2<-p2+theme_classic()+geom_vline(xintercept = 1550)+theme(legend.position = "top", text = element_text(size=16))+
  scale_color_manual(values = c("grey75","red3"),name="NAAT reactive during the study period")
p2


data.f2<-data[-which(is.na(data$serology_result) | data$symp_visit!=1 | data$vaccine_new==3),]
data.f2<-data.f2[,c(2:10,15)]
colnames(data.f2)[2]<-"Previous_SARS_CoV2_NAAT_confirmed_infection"
colnames(data.f2)[7]<-"Anti_spike_IgG_binding_antibody_units"

data.f2$time_since_inf2<-data.f2$time_since_inf
w1<-which(data.f2$time_since_inf<180 & data.f2$time_since_inf>=14)
w2<-which(data.f2$time_since_inf>=180 & data.f2$time_since_inf<360)
w3<-which(data.f2$time_since_inf>=360 | data.f2$time_since_inf<14)
data.f2$time_since_inf2[w1]<-1
data.f2$time_since_inf2[w2]<-2
data.f2$time_since_inf2[w3]<-3
data.f2$time_since_inf2<-as.factor(data.f2$time_since_inf2)
levels(data.f2$recent_vacc) <- c("< 6 months ago", "between 6 and 12 months ago", "over 1 year ago")
levels(data.f2$vaccine_new) <- c("No Vaccine", "JJ", "Pf", "JJ+Pf")
levels(data.f2$time_since_inf2) <- c("< 6 months ago", "between 6 and 12 months ago", "over 1 year ago or none")
# levels(data.f2$booster) <- c("single dose", "booster dose")
ww<-which(data.f2$vaccine_new=="Pf" & data.f2$booster==0)
data.f2<-data.f2[-ww,]
data.f2$Anti_spike_IgG_binding_antibody_units[data.f2$Anti_spike_IgG_binding_antibody_units<2]<-1

comps<- list( c("0", "1") )
library(viridis)
p<-ggplot(data.f2, aes(x=Previous_SARS_CoV2_NAAT_confirmed_infection, y=Anti_spike_IgG_binding_antibody_units, fill=vaccine_new:booster)) +
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(dotsize = .35, binaxis='y', stackdir='center', stackratio = .7, position=position_dodge(0.8))
p<-p+theme_classic()+theme(legend.position = "top",text = element_text(size=16))+ylab("Anti-spike IgG binding antibody units")+xlab("Previous SARS-CoV-2 NAAT confirmed infection")+
  scale_fill_viridis_d(labels = c("No vaccine", "Single dose Ad26.CoV.2", "Double dose Ad26.CoV.2", "Double dose BNT162b2"),
                    # labels=c("< 6 months ago", "between 6 and 12 months ago", "over 1 year ago or none"),
                    name = "Vaccination status")+
  scale_y_continuous(trans='log10',label=comma_format(accuracy=1))+
stat_compare_means(comparisons = comps)+stat_compare_means(label.y = 4.75)
p+guides( fill=guide_legend(nrow = 2))
  

tiff("Abs.tiff",width = 50, height = 35, units = "cm",res = 300)
ggarrange(p1,
ggarrange(p2+geom_vline(xintercept = 1550) , p+guides(fill = guide_legend(nrow = 2)) , 
          labels = c("B", "C"),
          ncol = 2, nrow = 1),nrow=2, labels="A")
dev.off()


stat.test <- data.f2 %>%
  group_by(prior_exposure) %>%
  t_test(spike_concentration ~ vaccine_new) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test


bxp <- ggboxplot(
  data.f2, x = "prior_exposure", y = "spike_concentration", 
  color = "vaccine_new", palette = c("#00AFBB", "#E7B800","#00AFB1")
)

stat.test <- stat.test %>%
  add_xy_position(x = "prior_exposure", dodge = 0.8)
bxp + stat_pvalue_manual(step.group.by = "booster",
  stat.test,  label = "p", tip.length = 0
)
