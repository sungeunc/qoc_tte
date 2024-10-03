require(tidyr)
library(naniar)
require(stringr)
require(lubridate)
library(AER)
library(MASS)
library(optmatch)
library(dplyr)
library(tidyverse)
library(dynamichazard)
library(MatchIt)
library(survival)
library(survminer)
library(dynamichazard)
library(cobalt)

data.mod # data that tracks individuals based on the treatment status (adhered vs non-adhered) for tfl 

#### for topical fluoride application

match_data <- matchit(formula = tfl2 ~  
                        age_visit +
                        sex + race + insure_base +
                        proc_eval+proc_seal+proc_cra+ annual_visit+
                        smoke,
                      data = data.mod[which(data.mod$tstart==0),],
                      method = "cem", 
                      estimand = "ATE")

dn.tfl_matched <- match.data(match_data)
ids <- dn.tfl_matched$studyID

match.col<-dplyr::select(dn.tfl_matched, c("studyID","weights","subclass"))
match.col$matched<-"Matched"
# 
d_match<- merge(data.mod,match.col,by="studyID", all.x=TRUE)
d_match$matched<-ifelse(is.na(d_match$matched),"Unmatched",d_match$matched)
d_match$weights<-ifelse(is.na(d_match$weights),1,d_match$weights)

# Save IDs that were matched
ids <- dn.tfl_matched$studyID

## CEM data extraction
match.col<-dplyr::select(dn.tfl_matched, c("studyID","weights","subclass"))
match.col$matched<-"Matched"

d_match<- merge(data.mod,match.col,by="studyID", all.x=TRUE)

# survival curve
fit.km <- survfit(Surv(tstart, tstop, status) ~ tfl2+matched , id=studyID,
                           conf.type="log-log",#weights=weights,
                           data = d_match)#d_match[d_match$matched=="Matched",])


plot <- ggsurvplot(
  fit.km, 
  fun = "event", # plot cumulative incidence
  conf.int = TRUE, # include confidence intervals
  censor = FALSE, # don't include tick marks for events/censorings
  xlim = c(0,24),
  xlab = "Months", # label x-axis
  break.x.by = 2, # display x-axis in 2-week bins
  surv.scale = "percent", # show y-axis in %
  ylab = "Cumulative Incidence (%)", # label y-axis
  ylim = c(0,0.7), # set limits of y-axis
  legend = c(0.15, 0.8), # set legend position
  legend.labs = c("Matched:Untreated", "Unmatched:Untreated","Matched:Treated", "Unmatched:Treated"),
  legend.title = "", 
  title = "Topical Fluoride (TFL-CH)",
  palette = c("#D55E00","#E69F00", "#0072B2", "#56B4E9"),
  linetype = c("solid","dashed", "solid","dashed"),
  ggtheme=custom_theme()) # set colors


d_match_tte<-d_match[d_match$matched=="Matched",]

#unadjusted survival model by treatment status
unadj.hr  <- coxph(formula = Surv(tstart, tstop, status) ~ tfl2,method="breslow",robust=TRUE,
                   cluster = studyID, data=data.mod) #d_match[d_match$matched=="Matched",])

# adjusted survival model 
cox_fit_A <- coxph(formula = Surv(tstart, tstop, status) ~ tfl2+age_visit+
                     sex+race+insure_base+smoke+baseline_cra+annual_visit+
                     proc_eval+proc_seal,method="breslow",robust=TRUE, cluster=subclass,weights=weights,
                   data=d_match_tte) 



v <- data.frame(old = c("age_visit", "sex", "race", "insure_base","proc_cra", 
                        "proc_eval", "proc_seal","smoke","annual_visit"),
                new = c("Age", "Sex", "Race", 
                        "Insurance","Caries risk", "Oral evaluation", "Sealant","Smoking","Annual visit"))


# Generate covariate balance plot  using cobalt package
love.plot(match_data, 
          binary = "std", 
          title = "Covariate balance plot - Topical fluoride",
          sample.names = c("Pre-match", "Post-match"),
          stats = c("mean.diffs", "ks.statistics"),
          #addl = ~ sex+race+insure_base+proc_eval+proc_seal+annual_visit, 
          var.names = v,
          colors = c("#E69F00","#0072B2"),
          shapes = c("triangle", "circle"),
          line = F, 
          abs = F,
          drop.distance = TRUE,
          grid = TRUE)


#### for sealant 

match_data <- matchit(formula = sealant ~ 
                        #age_visit +
                        sex + race + insure_base +
                        proc_eval+proc_fluor+proc_cra+annual_visit,
                      data = data.mod[which(data.mod$tstart==0),],
                      method = "cem",
                      estimand="ATE")

dn.seal_matched <- match.data(match_data)
ids <- dn.seal_matched$studyID

match.col<-dplyr::select(dn.seal_matched, c("studyID","weights","subclass"))
match.col$matched<-"Matched"
# 
d_match<- merge(data.mod,match.col,by="studyID", all.x=TRUE)
d_match$matched<-ifelse(is.na(d_match$matched),"Unmatched",d_match$matched)
d_match$weights<-ifelse(is.na(d_match$weights),1,d_match$weights)

fit.km  <- survfit(Surv(tstart, tstop, status) ~ sealant+matched, id=studyID, 
                   conf.type="log-log", 
                   data = data.mod, weights=weights)

plot <- ggsurvplot(
  fit.km, 
  fun = "event", # plot cumulative incidence
  conf.int = TRUE, # include confidence intervals
  censor = FALSE, # don't include tick marks for events/censorings
  xlim = c(0,24),
  xlab = "Months", # label x-axis
  break.x.by = 2, # display x-axis in 2-week bins
  surv.scale = "percent", # show y-axis in %
  ylab = "Cumulative Incidence (%)", # label y-axis
  ylim = c(0,0.6), # set limits of y-axis
  legend = c(0.15, 0.9), # set legend position
  legend.labs = c("Matched:Untreated", "Unmatched:Untreated","Matched:Treated", "Unmatched:Treated"),
  legend.title = "", # set legend title
  title = "Sealant (SFM-CH)",
  palette = c("#D55E00","#E69F00", "#0072B2", "#56B4E9"),
  linetype = c("solid","dashed", "solid","dashed"),
  ggtheme=custom_theme()) # set colors

d_match_tte<-d_match[d_match$matched=="Matched",]

# survival model

cox_fit_un <- coxph(formula = Surv(tstart, tstop, status) ~ sealant,method="breslow",robust=TRUE,weights=weights,
                    cluster = studyID, data=d_match[d_match$matched=="Matched",])

cox_fit_A <- coxph(formula = Surv(tstart, tstop, status) ~ age_visit+sealant+ 
                     sex+race+insure_base+baseline_cra+annual_visit+
                     proc_eval+proc_fluor,
                   method="breslow",robust=TRUE, cluster = studyID,weights=weights,
                   data =d_match_tte)


