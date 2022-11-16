library(readr)
library(tidyverse)
library(jsonlite)
library(lubridate)
library(dplyr)
library(digest)
library(skimr)
library(timetk)
library(lme4)
library(performance)
library(patchwork)
library(modelr)
library(lmerTest)
library(rstatix)
library(effectsize)
library(visdat)
library(Hmisc)
library(table1)
library(ggdist)
library(gghalves)
library(gmodels)
library(plotrix)



# 1. load data frame ----

sced1 <- readRDS(file = "sced1.RDS") 
sced <- sced1



# 2. subjects: demographics, tinnitus characteristics & app use (Table 1) ----

#create subset subjects
subjects <- sced %>%
  group_by(case_id) %>%
  arrange(con_days) %>%
  filter(row_number()==1) %>%
  select(1:3,6:8,17:43) %>% 
  arrange(case_id) %>% 
  ungroup()

#demographics
table1(~ age + sex + thi_pre + duration_tin_mon + phq9_score + guef_score, data = subjects)

#tinnitus characteristics
subjects %>% count(education)
subjects %>% count(vertigo)
subjects %>% count(tin_num) 
subjects %>% count(tin_on) 
subjects %>% count(tin_change) 
subjects %>% count(tin_qual) 
subjects %>% count(tin_pitch) 
subjects %>% count(tin_loc) 
subjects %>% count(tin_care) 

#app use
table1(~ n_tinedu + n_son + n_diary, data = subjects)

# 3. Change in THI ----

# create new variable thi_diff, create subset thi, switch to long format for boxplot
subjects %>% mutate(thi_diff = thi_pre - thi_post) -> subjects
thi <- select(subjects, case_id, thi_pre, thi_post)
thi <- gather(thi, thi, value, thi_pre:thi_post)
thi$thi <- recode(thi$thi, 'thi_pre' = "Baseline", 'thi_post' = "Final visit")

#Table 2
table1(~ thi_pre + thi_post + thi_diff, data = subjects)
t.test(subjects$thi_pre, subjects$thi_post, paired = TRUE)
cohens_d(subjects$thi_pre, subjects$thi_post, paired = TRUE)

subjects %>% count(thi_diff >= 7)

#Figure 2
thi %>% group_by(thi) %>% 
  summarise(mean = mean(value), se = std.error(value)) -> se

#a)
ggplot(se, aes(thi, mean, colour = thi)) +
  geom_point(size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.1, size = 1, show.legend = FALSE) +
  ylab("THI") + xlab("Intervention moment") + ylim(0,80) + theme_classic(base_size = 17) +
  scale_colour_manual(values=c("#12737b", "#f2953d"))

thi$case_id <- factor(thi$case_id)

#b)
ggplot(thi, aes(thi, value, group=case_id,shape=case_id, colour = thi)) +
  geom_jitter(size = 2.5,show.legend = FALSE) +
  scale_shape_manual(values=1:nlevels(thi$case_id)) +
  ylab("THI") + xlab("Intervention moment") + theme_classic(base_size = 17) + ylim(0,80) +
  scale_colour_manual(values=c("#12737b", "#f2953d"))

# 4. ema data: visualisation & change in tinnitus symptoms  ----

# first prepare data set: remove duplicates, cut dataset after 12 week intervention, remove cases not missing at random & impute missing values

#create subset app
app <- select(sced, case_id, group_id, questionnaire, age, sex, duration_tin_mon, thi_pre, phq9_score, bfi2_extraversion, bfi2_neg_emotion, guef_score, con_days, base_days, int_days, magic_days, loudness_cur:emotion_tod)
app <- filter(app, questionnaire == "diary")
sced %>% count(questionnaire)

#remove duplicates by days and cases -> 79 duplicates removed
app <- distinct(app, case_id, con_days, .keep_all= TRUE) 

#complete data set with missing values and refill demographics
complete(app, case_id, con_days) -> df_imp
df_imp %>% 
  arrange(case_id, con_days) %>% 
  group_by(case_id) %>% 
  fill(group_id, questionnaire, age, sex, duration_tin_mon, thi_pre, phq9_score, bfi2_extraversion, bfi2_neg_emotion, guef_score, .direction = c("downup")) -> df_imp

#filter dataset by time: 12 weeks of intervention for each group
test1 <- filter(df_imp[which(df_imp$group_id==1),], con_days <= 95)
test2 <- filter(df_imp[which(df_imp$group_id==2),], con_days <= 103)
test3 <- filter(df_imp[which(df_imp$group_id==3),], con_days <= 100)
test4 <- filter(df_imp[which(df_imp$group_id==4),], con_days <= 92)

df_imp <- rbind(test1, test2, test3, test4)

#add intervention variable to dataset by group and baseline length
df_imp$intervention <- 99
df_imp <- df_imp %>% relocate(intervention, .after = group_id)

df_imp$intervention[df_imp$group_id == "1" & df_imp$con_days < 11] <- "0"
df_imp$intervention[df_imp$group_id == "1" & df_imp$con_days >= 11] <- "1"
df_imp$intervention[df_imp$group_id == "2" & df_imp$con_days < 19] <- "0"
df_imp$intervention[df_imp$group_id == "2" & df_imp$con_days >= 19] <- "1"
df_imp$intervention[df_imp$group_id == "3" & df_imp$con_days < 16] <- "0"
df_imp$intervention[df_imp$group_id == "3" & df_imp$con_days >= 16] <- "1"
df_imp$intervention[df_imp$group_id == "4" & df_imp$con_days < 8] <- "0"
df_imp$intervention[df_imp$group_id == "4" & df_imp$con_days >= 8] <- "1"

#vis_dat: estimate missingness at random per case
subset(df_imp, case_id == "03") %>% vis_miss() #not
subset(df_imp, case_id == "04") %>% vis_miss() #missing at random (end day 92)
subset(df_imp, case_id == "05") %>% vis_miss() #missing at random
subset(df_imp, case_id == "06") %>% vis_miss() #not
subset(df_imp, case_id == "07") %>% vis_miss() #missing at random
subset(df_imp, case_id == "08") %>% vis_miss() #not
subset(df_imp, case_id == "10") %>% vis_miss() #missing at random
subset(df_imp, case_id == "11") %>% vis_miss() #missing at random
subset(df_imp, case_id == "12") %>% vis_miss() #missing at random (end day 88)
subset(df_imp, case_id == "13") %>% vis_miss() #missing at random (end day 84)
subset(df_imp, case_id == "14") %>% vis_miss() #missing at random 
subset(df_imp, case_id == "15") %>% vis_miss() #missing at random (end day 67)
subset(df_imp, case_id == "16") %>% vis_miss() #missing at random
subset(df_imp, case_id == "17") %>% vis_miss() #missing at random
subset(df_imp, case_id == "18") %>% vis_miss() #missing at random
subset(df_imp, case_id == "19") %>% vis_miss() #missing at random
subset(df_imp, case_id == "20") %>% vis_miss() #missing at random
subset(df_imp, case_id == "21") %>% vis_miss() #not

#select cases where NAs not missing at random & delete NAs at the end
df <- filter(df_imp, !(case_id %in% c("03","06","08","21")))

df %>% arrange(case_id, con_days) -> df
df[-c(93:100, 570:576, 661:671, 841:868),] -> df

sum(is.na(df$loudness_cur)) #151*10 = 1510 data points to be imputed 

#impute: aregImpute
set.seed(1)
imp_arg <- aregImpute(~ loudness_cur + distress_cur + jaw_tension + neck_tension + 
                        tin_thoughts + distress_tod + loudness_tod + movement_tod + stress_tod + emotion_tod,
                      x = TRUE, data = df, n.impute = 5)
set.seed(1)
df_imp <- impute.transcan(imp_arg, imputation=1, data=df, list.out=TRUE,pr=FALSE, check=FALSE) 

as.data.frame(do.call(cbind,df_imp)) -> imp

glimpse(df)
df[,17:26] <- imp

#S4: EMA subset demographics

ema_subjects <- distinct(df, case_id, .keep_all = TRUE)
table1(~ sex + age + phq9_score + thi_pre + duration_tin_mon + guef_score, data = ema_subjects)

#Figure 3: visualisation of EMA distress & loudness
colours <- c("distress_cur" = "#12737b", "loudness_cur" = "#f2953d")
g1 = subset(df, group_id == 1)
g2 = subset(df, group_id == 2)
g3 = subset(df, group_id == 3)
g4 = subset(df, group_id == 4)

ggplot() +
  geom_rect(data = g1,aes(xmin=-Inf, xmax=11, ymin=-Inf, ymax=Inf), fill = "gray85", alpha = 0.05) +
  geom_rect(data = g2,aes(xmin=-Inf, xmax=19, ymin=-Inf, ymax=Inf), fill = "gray85", alpha = 0.05) +
  geom_rect(data = g3,aes(xmin=-Inf, xmax=16, ymin=-Inf, ymax=Inf), fill = "gray85", alpha = 0.05) +
  geom_rect(data = g4,aes(xmin=-Inf, xmax=18, ymin=-Inf, ymax=Inf), fill = "gray85", alpha = 0.05) +
  geom_line(data = df,aes(x = con_days, y = distress_cur, colour = "distress_cur")) +
  geom_line(data = df,aes(x = con_days, y = loudness_cur, colour = "loudness_cur")) +
  facet_wrap(~ case_id, nrow = 6) +
  theme_classic(base_size = 15) +
  labs(x = "Consecutive days",
       y = "Tinnitus symptom strength",
       color = "Legend") +
  scale_colour_manual(values = colours, labels = c("Distress", "Loudness")) + 
  theme(legend.position = c(.9, .05))


# Table 3 & 4: Change of tinnitus symptoms from baseline to end of intervention

#create subset for baseline phase
base <- df %>% group_by(case_id, group_id) %>% filter(intervention == 0)

#create subset for each group and select corresponding time window of intervention phase (last days defined by baseline length per group)
group1 <- filter(df, group_id == 1)
group1 %>% group_by(case_id) %>% arrange(desc(con_days)) %>% filter(row_number() <= 11) -> group1

group2 <- filter(df, group_id == 2)
group2 %>% group_by(case_id) %>% arrange(desc(con_days)) %>% filter(row_number() <= 19) -> group2

group3 <- filter(df, group_id == 3)
group3 %>% group_by(case_id) %>% arrange(desc(con_days)) %>% filter(row_number() <= 16) -> group3

group4 <- filter(df, group_id == 4)
group4 %>% group_by(case_id) %>% arrange(desc(con_days)) %>% filter(row_number() <= 8) -> group4

#df2: combine datasets
rbind(base, group1, group2, group3, group4) -> df2

#summarise to mean and sd
df2 %>% 
  group_by(case_id, intervention) %>%
  summarise(mdist = mean(distress_cur), mloud = mean(loudness_cur), sddist = sd(distress_cur), sdloud = sd(loudness_cur)) -> df2

#convert to wide format, calculate delta (Table 3 & 4)
options(digits=3)

df2 %>% 
  pivot_wider(names_from = "intervention", values_from = c(mdist, sddist, mloud, sdloud)) %>% 
  mutate(ddist = mdist_0 - mdist_1, dloud = mloud_0 - mloud_1) %>% 
  select(case_id, mdist_0, sddist_0, mdist_1, sddist_1, ddist, mloud_0, sdloud_0, mloud_1, sdloud_1, dloud) -> df2

#t-test on a group level for distress and loudness
options(digits=5)
t.test(df2$mdist_0, df2$mdist_1, paired = TRUE)
t.test(df2$mloud_0, df2$mloud_1, paired = TRUE)


# Figure 4: change of correlation for distress and loudness

#correlation at baseline
baseline <- subset(df, intervention == 0)
cor.test(baseline$distress_cur, base$loudness_cur)

#correlation at intervention
tx <- subset(df, intervention == 1)
cor.test(tx$distress_cur, tx$loudness_cur)

#create subsets: intervention phase in weeks 
week1 <- subset(tx, con_days>=(1*7+1) & con_days<=(2*7))
week2 <- subset(tx, con_days>=(2*7+1) & con_days<=(3*7))
week3 <- subset(tx, con_days>=(3*7+1) & con_days<=(4*7))
week4 <- subset(tx, con_days>=(4*7+1) & con_days<=(5*7))
week5 <- subset(tx, con_days>=(5*7+1) & con_days<=(6*7))
week6 <- subset(tx, con_days>=(6*7+1) & con_days<=(7*7))
week7 <- subset(tx, con_days>=(7*7+1) & con_days<=(8*7))
week8 <- subset(tx, con_days>=(8*7+1) & con_days<=(9*7))
week9 <- subset(tx, con_days>=(9*7+1) & con_days<=(10*7))
week10 <- subset(tx, con_days>=(10*7+1) & con_days<=(11*7))
week11 <- subset(tx, con_days>=(11*7+1) & con_days<=(12*7))
week12 <- subset(tx, con_days>=(12*7+1) & con_days<=(13*7))

#correlation at different intervention weeks
cor.test(week1$distress_cur, week1$loudness_cur)
cor.test(week2$distress_cur, week2$loudness_cur)
cor.test(week3$distress_cur, week3$loudness_cur)
cor.test(week4$distress_cur, week4$loudness_cur)
cor.test(week5$distress_cur, week5$loudness_cur)
cor.test(week6$distress_cur, week6$loudness_cur)
cor.test(week7$distress_cur, week7$loudness_cur)
cor.test(week8$distress_cur, week8$loudness_cur)
cor.test(week9$distress_cur, week9$loudness_cur)
cor.test(week10$distress_cur, week10$loudness_cur)
cor.test(week11$distress_cur, week11$loudness_cur)
cor.test(week12$distress_cur, week12$loudness_cur)

#create correlation subset
x <- c(0,1,2,3,4,5,6,7,8,9,10,11,12)
y <- c(0.676, 0.702, 0.715, 0.671, 0.667, 0.646, 0.62, 0.526, 0.509, 0.518, 0.471, 0.349, 0.427)

cor <- data.frame(x,y)

#Figure 5
ggplot(cor) + geom_point(aes(x,y), colour = "#12737b", size = 3) + 
  theme_classic(base_size =16) + labs(x="weeks", y="cor") + 
  ylim(0, 0.75) +
  scale_x_continuous(breaks = seq(0, 12, by = 1)) 
  


# 5. ema data: Mixed effect model (Table 5) ----
set.seed(1)
mod2 <- lmer(distress_cur ~ con_days + 
               intervention +
               age + 
               sex +
               thi_pre +
               phq9_score +
               (1|case_id),
             data = df)
summary(mod2)
performance(mod2)
confint(mod2)
check_model(mod2)

#Train test split method
library(tidymodels)
set.seed(1)
innitial_split <- initial_split(df, prop = 4/5)
training <- training(innitial_split)
testing <- testing(innitial_split)
set.seed(1)
mod2 <- lmer(distress_cur ~ con_days + 
               intervention +
               age + 
               sex +
               thi_pre +
               phq9_score +
               (1|case_id),
             data = training)

test_results <- as.data.frame(predict(mod2, testing) %>% cbind(testing$distress_cur))
rsq_vec(estimate = test_results$`.`, truth = test_results$V2)

ggplot(test_results, aes(x = `.`, y = V2)) +
         geom_point()

# 6. Compare THI with EMA data (Figure 5) ----

#calculate Distress trend for each case (slope from a regression line)
i04 <- subset(df, case_id == "04" & intervention == "1")  
lm(distress_cur ~ con_days, data = i04)

i05 <- subset(df, case_id == "05" & intervention == "1")  
lm(distress_cur ~ con_days, data = i05)

i07 <- subset(df, case_id == "07" & intervention == "1")  
lm(distress_cur ~ con_days, data = i07)

i10 <- subset(df, case_id == "10" & intervention == "1")  
lm(distress_cur ~ con_days, data = i10)

i11 <- subset(df, case_id == "11" & intervention == "1")  
lm(distress_cur ~ con_days, data = i11)

i12 <- subset(df, case_id == "12" & intervention == "1")
lm(distress_cur ~ con_days, data = i12)

i13 <- subset(df, case_id == "13" & intervention == "1")   
lm(distress_cur ~ con_days, data = i13)

i14 <- subset(df, case_id == "14" & intervention == "1")   
lm(distress_cur ~ con_days, data = i14)

i15 <- subset(df, case_id == "15" & intervention == "1")   
lm(distress_cur ~ con_days, data = i15)

i16 <- subset(df, case_id == "16" & intervention == "1")   
lm(distress_cur ~ con_days, data = i16)

i17 <- subset(df, case_id == "17" & intervention == "1")   
lm(distress_cur ~ con_days, data = i17)

i18 <- subset(df, case_id == "18" & intervention == "1")   
lm(distress_cur ~ con_days, data = i18)

i19 <- subset(df, case_id == "19" & intervention == "1")   
lm(distress_cur ~ con_days, data = i19)

i20 <- subset(df, case_id == "20" & intervention == "1")   
lm(distress_cur ~ con_days, data = i20)

#subset df3
df3 <- select(subjects, case_id, thi_diff) %>% 
  filter(case_id %in% c("04","05","07",10,11,12,13,14,15,16,17,18,19,20))

#add Distress trend & Delta Distress (ddist from df2)
slope <- c(0.2,0.4,-0.2,0.1,0.05,-0.2,-0.2,-0.1,0.1,-0.4,-0.1,-0.1,-0.5,0.2)
df3$slope <- slope
df3$ddist <- df2$ddist 

#correlation between THI and EMA change scores
cor.test(df3$slope, df3$thi_diff)
cor.test(df3$ddist, df3$thi_diff)

#Figure 5
lm(thi_diff ~ slope, data = df3)
ggplot(df3, (aes(slope,thi_diff))) + geom_point(color = "#12737b", size = 3) + theme_classic(base_size = 22) + geom_abline(intercept = 8.72, slope = -34.60, color = "#12737b", size = 1) + labs(x = "Distress trend (EMA)", y = "∆THI")

lm(thi_diff ~ ddist, data = df3)
ggplot(df3, (aes(ddist,thi_diff))) + geom_point(color = "#f2953d", size = 3) + theme_classic(base_size = 22) + geom_abline(intercept = 7.78, slope = 0.614, color = "#f2953d", size = 1) + labs(x = "∆Distress (EMA)", y = "∆THI")


