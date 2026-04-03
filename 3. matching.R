rm(list=setdiff(ls(), lsf.str()))   #removing everything but functions

foo <- function(x)
{for( i in x ){
  #  require returns TRUE invisibly if it was able to load package
  if( ! require( i , character.only = TRUE ) ){
    #  If package was not able to be loaded then re-install
    install.packages( i , dependencies = TRUE )
    #  Load package after installing
    require( i , character.only = TRUE ) }}}

foo(c("ggplot2" , "reshape2", "zoo", "hrbrthemes", "GGally", "ggpubr", "esquisse", "ggthemes",  
      "viridis",  "data.table" , "logr", "janitor", "readr", "lubridate",  "ggalluvial", 
      "naniar", "stringdist",  "purrr", "stringr", "mondate", "xlsx", "lemon", 
      "readr", "readxl", "dplyr", "tidyr", "kableExtra", "directlabels", "tableone",
      "ggrepel", "gridExtra", "lattice", "zoo", "DT", "grid", "plotly", "MatchIt", "broom", "gridExtra", "optmatch"))


setwd("C:/Users/rsodhi.ic/Box/Ridh_workspace/CFL/Connect for Life - Manuscript Submission/data")

# loading the main data
d1<-read_csv("C:/Users/rsodhi.ic/Box/Ridh_workspace/CFL/Connect for Life - Manuscript Submission/data/output/analytical.data.csv", guess_max=20000)[-1]


table(d1$test.control, useNA="ifany")
table(d1$male, useNA="ifany")
table(d1$outcome, useNA="ifany")

## 1.2 Difference-in-means: pre-treatment covariates

d1.cov <- c('age', 'male', 'cbnaat.num', 'conn.follow.ups', 'ep')
d1%>%
  group_by(cfl) %>%
  select(one_of(d1.cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))

d1.cov <- c('age', 'male', 'cbnaat.num', 'conn.follow.ups', 'ep', 'clinic.name', 'fdc.num')

# score estimation----
d1.ps <- glm(cfl ~   male + age.class+ cbnaat.num + fdc.num+   ep + factor(clinic.name) + diag.qtr, 
             family = binomial(), data = d1)
summary(d1.ps)

# saving the results

d1.ps.exp = exp(coef(d1.ps)) 

library(texreg)
library(stargazer)
stargazer(d1.ps, 
          title = "Impact of CfL: Average Treatment Effect on Matched dataset",
          coef=list(d1.ps.exp), 
          type = "text",
          model.numbers = F,
          omit.table.layout = "n",
          digits = 3, 
          column.labels = c("Model: Propensity Score using Logistic Regression"),
          dep.var.caption  = "Dependent Variable: CfL",
          single.row=TRUE,
          dep.var.labels.include = FALSE,
          ci = TRUE,
          out="C:/Users/rsodhi.ic/Box/Ridh_workspace/CFL/Connect for Life - Manuscript Submission/data/results/results.propensity.html") 




#estimation
tidy(d1.ps, exponentiate=TRUE, conf.int = TRUE)


# predicting----

pred.d1 <- data.frame(pr.score = predict(d1.ps, type = "response"),
                      cfl = d1.ps$model$cfl)
head(pred.d1)


predicted_prob <- predict(d1.ps, d1, type = "response")

prediction <- as.integer(predicted_prob > 0.4)
confusion_mat <- addmargins(table(d1$cfl, prediction))
confusion_mat

names(dimnames(confusion_mat)) <- c("True status", "Prediction")
colnames(confusion_mat) <- c("Fail", "Success", "Total")
rownames(confusion_mat) <- c("Fail", "Success", "Total")
confusion_mat

# Accuracy

accuracy=(confusion_mat[1,1]+confusion_mat[2,2])/count(d1)
accuracy

TYPE.1.Error.Specificity=(confusion_mat[1,2])/(confusion_mat[1,1]+confusion_mat[1,2])
TYPE.1.Error.Specificity
# Of the total observations predicted to not have the treatment, 26% are incorrect. 

true.p.rate=(confusion_mat[2,2])/confusion_mat[2,3]
true.p.rate #(Sensitivity; tp/(fn+tp))



# After estimating the propensity score, it is useful to plot histograms of 
#the estimated propensity scores by treatment status:
  
labs <- paste("Actual status: ", c("CFL", "No CFL"))
pred.d1%>%
  mutate(cfl = ifelse(cfl == 1, labs[1], labs[2])) %>%
  ggplot(aes(x = pr.score)) +
  geom_histogram(color = "white", bins = 30) +
  facet_wrap(~cfl) +
  xlab("Probability of getting CFL") +
  theme_bw()


# prob of score 

t1<-pred.d1%>%
  group_by(cfl)%>%
  dplyr::summarise(n=n(), 
                   mean.score=mean(pr.score), 
                   med.score=median(pr.score))
clipr::write_clip(t1)

t1



### 3.1 Nearest Neighbor Matching

d1.cov <- c('age', 'male', 'cbnaat.num', 'conn.follow.ups', 'ep', 'clinic.name', 'fdc.num')

d1.nomiss <- d1%>%  # MatchIt does not allow missing values
  select(patient.id, outcome, treatment.outcome, age, cfl, ep, clinic.name, one_of(d1.cov), fdc.num, cfl, age.class, diag.qtr)%>%
  mutate(clinic.name=as.factor(clinic.name))


# Match1. (Full) Matching dataset without weights and withought replacement---------

mod.match<- matchit(cfl ~ age.class + male + cbnaat.num + fdc.num + clinic.name + diag.qtr + ep, 
                    method = "full", data = d1.nomiss,  link="logit", 
                    distance = "glm", caliper = .2,
                    weights = "weights", 
                    #ratio=2, replace=F,
                    exact = ~male+ ep+ cbnaat.num + fdc.num)

table(d1$treatment.outcome)
d1.m<- match.data(mod.match) 
dim(d1.m)
summary(mod.match, standardize = T)

table(d1.m$cfl)
table( d1.m$cfl, d1.m$outcome)

table(d1.m$treatment.outcome)

prop.table(table( d1.m$cfl, d1.m$outcome))

write.csv(d1.m, "output/cfl.matched.full.csv")




# Match 2 (nearest neighbor)-----------
mod.match<- matchit(cfl ~ age.class + male + cbnaat.num + fdc.num + clinic.name + diag.qtr + ep, 
                    method = "nearest", data = d1.nomiss,  link="logit", 
                    distance = "glm", caliper = .2,
                    weights = "weights", 
                    #ratio=2, replace=F,
                    exact = ~male+ ep+ cbnaat.num)

table(d1$treatment.outcome)
d1.m.1<- match.data(mod.match) 
dim(d1.m.1)
summary(mod.match, standardize = T)



write.csv(d1.m.1, 
          "output/cfl.matched.nearest.csv")



# Useful plots----------------

# Plot 1-----------

library(cobalt)
#install.packages("cobalt")
k1<-d1.m%>%
  group_by(cfl)%>%
  dplyr::summarise(n=n(), 
                   mean.male=mean(male), 
                   mean.age=mean(age))
k1
par(mar=c(2,2,2,2))

plot(mod.match, type="hist")

# plot2--------
b3<-bal.plot(mod.match, var.name='male', which='both', grid=TRUE, type="ecdf")
b4<-bal.plot(mod.match, var.name='cbnaat.num', which='both', grid=TRUE, type="ecdf")

grid.arrange(b3, b4)



# loveplot-----------

love.plot(bal.tab(mod.match, m.threshold=0.1), 
          stat="mean.diffs", 
          # grid=F,
          stars="raw", 
          # line = T,
          position = "right",
          sample.names = c("Observational", "Matched"))

# checking diff-----------
table_match2 <- CreateTableOne(vars = d1.cov,strata = "cfl",data = d1.m,test = FALSE)
print(table_match2, smd = TRUE)
k2<-d1.m%>%
  group_by(cfl)%>%
  dplyr::summarise(n=n(), 
                   mean.follow=mean(conn.follow.ups), 
                   median.follow=median(conn.follow.ups),
                   mean.male=mean(male))
k2




# We can get some information about how successful the matching was using summary(mod_match) and plot(mod_match). 

t0<-summary(mod.match)
t1<-summary(mod.match)

p1<-as.data.frame(t1$sum.matched, digits = 4, 
                  align = 'c', 
                  caption = 'Table 3: Summary of balance for matched data')
colnames(p1)<-paste0(colnames(p1), "99")

p2<-as.data.frame(t0$sum.all, digits = 4, 
                  align = 'c', 
                  caption = 'Table 3: Summary of balance for matched data')
#p2.1<-p2[-1, ]
p3<-bind_cols(p2[1:4], p1[1:4])%>%
  mutate(diff=`Means Treated`-`Means Control`, 
         diff.1=`Means Treated99`-`Means Control99`, )%>%
  select(`Means Treated`, `Means Control`, diff, `Std. Mean Diff.`, `Var. Ratio`, 
         `Means Treated99`, `Means Control99`, diff.1, `Std. Mean Diff.99`, `Var. Ratio99`, )

colnames(p3)
clipr::write_clip(p3)  
p3



#Matched data summary----------
match.summary <- d1.m%>%select(`Males` = male, 
                               `Age Category` = age.class, 
                               `Age` = age, 
                               `Access to FDCs`=fdc.num, 
                               `Xpert Testing`=cbnaat.num, 
                               `Extra Pulmonary` = ep,
                               `Follow Ups` = conn.follow.ups,
                               `Clinic`=clinic.name, 
                               `CFL` = cfl)
(match.summary%>% tbl_summary(by = `CFL`)%>%add_p())


d1.sel.1<-d1.m%>%
  select(clinic.name, male, age, 
         CFL=cfl, `Successful follow Ups`=conn.follow.ups, 
         `Xpert Testing for Diagnosis` = cbnaat.num, 
         `Average Adherence`=avg.adh.calc, 
         `Total Doses` = tot.doses, 
         `Total Doses taken`= tot.taken, 
         `Extra Pulmonary`= ep)%>%
  ungroup()
