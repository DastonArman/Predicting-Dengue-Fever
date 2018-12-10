#library and source various packages
library(mosaic)
library(sjPlot)
source('http://jgscott.github.io/teaching/r/utils/class_utils.R')

#filter data to only complete cases and with positive values for total case
dengue1 = dengue[complete.cases(dengue) & dengue$total_cases>0,]

### Problem 1
### Can we build a predictive model?

#Figure 1
hist((dengue1$total_cases[dengue1$total_cases < 200]), breaks=40, main='Histogram of Total Cases', xlab='Total Cases', col='light blue')
max(dengue1$total_cases)

#Stepwise Selection
large_model = lm(log(total_cases) ~ ndvi_ne + ndvi_nw + ndvi_se + ndvi_sw + precip_amt_kg_per_m2 + 
                                    specific_humidity_g_per_kg + air_temp_k + min_air_temp_k + 
                                    air_temp_k + precipitation_amt_mm + (city + season + 
                                    avg_temp_k + dew_point_temp_k + 
                                    relative_humidity_percent + 
                                    tdtr_k)^2 + ndvi_ne:city + ndvi_nw:city + ndvi_sw:city + ndvi_se:city, data=dengue1)

step_model = step(large_model)

small_model = lm(log(total_cases) ~ city, data=dengue1)
step_model2 = step(small_model, scope = (log(total_cases) ~ ndvi_ne + ndvi_nw + ndvi_se + ndvi_sw + precip_amt_kg_per_m2 + 
                     specific_humidity_g_per_kg + air_temp_k + min_air_temp_k + 
                     air_temp_k + precipitation_amt_mm + (city + season + 
                     avg_temp_k + dew_point_temp_k + 
                     relative_humidity_percent + 
                     tdtr_k)^2 + ndvi_ne:city + ndvi_nw:city + ndvi_sw:city + ndvi_se:city), data=dengue1)

medium_model = lm(log(total_cases) ~ city + season + ndvi_se + ndvi_sw + dew_point_temp_k + 
                  relative_humidity_percent + tdtr_k, data=dengue1)
step_model3 = step(small_model, scope = (log(total_cases) ~ ndvi_ne + ndvi_nw + ndvi_se + ndvi_sw + precip_amt_kg_per_m2 + 
                                         specific_humidity_g_per_kg + air_temp_k + min_air_temp_k + 
                                         air_temp_k + precipitation_amt_mm + (city + season + 
                                         avg_temp_k + dew_point_temp_k + 
                                         relative_humidity_percent + 
                                         tdtr_k)^2 + ndvi_ne:city + ndvi_nw:city + ndvi_sw:city + ndvi_se:city), data=dengue1)

splits = do(10000)*{
  train_cases = sample(1:nrow(dengue1), size=800)
  dengue_train = dengue1[train_cases,]
  dengue_test = dengue1[-train_cases,]

  step_model_train = update(step_model2, data=dengue_train)
  
  yhat_test_step = predict(step_model2, newdata=dengue_test)
  
  RMSPE_step = sqrt(mean((yhat_test_step - log(dengue_test$total_cases))^2))
  
  c(RMSPE_step)
}

#Appendix B
sjt.lm(step_model2)

#Figure 2
hist(exp(splits$result), breaks=40, main='Histogram of Estimated RMSPEout', col='light blue', xlab='Estimated RMSPEout (cases)')
conf = exp(confint(splits$result, level=.8))
conf
abline(v=conf, col='red', lwd=5)

#Figure 3
par(mfrow=c(1,2))
plot((resid(step_model2))~log(total_cases), data=dengue1, col='blue', main='Residuals vs. log(total_cases)',  xlab='Total Cases (log(cases))', ylab='Residual (log(cases))')
plot(exp(resid(step_model2))~(total_cases), data=dengue1, main='Residuals vs. Total Cases', xlab='Total Cases (cases)', ylab='Residual (cases)', pch=1, col='blue')
par(mfrow=c(1,1))


### Problem 2
### What is the partial slope of "insert precipitation measure" on total cases, holding other factors constant?

#Figure 4: Histogram of Specific Humidity
hist(dengue1$specific_humidity_g_per_kg, breaks=30, col='light blue', xlab='Specific Humidity (g/kg)', main='Histogram of Specifc Humidity')


#Initial model without confounders
slope_model = lm(log(total_cases) ~ specific_humidity_g_per_kg, data=dengue1)
summary(slope_model)
plot(resid(slope_model) ~ log(total_cases), data=dengue1)

#Figure 5: City as a confounder
par(mfrow=c(1,2))
boxplot(log(total_cases) ~ city, data=dengue1, main='Boxplot of log(total_cases) vs. City', col='light blue', xlab='City', ylab='log(total_cases)')
boxplot(specific_humidity_g_per_kg ~ city, data=dengue1, main='Boxplot of Specific Humidity vs. City', col='light blue', xlab='City', ylab='Specific Humidity (g/kg)')
par(mfrow=c(1,1))

#Figure 6: Season as a confounder
par(mfrow=c(1,2))
boxplot(log(total_cases) ~ season, data=dengue1, main='Boxplot of log(total_cases) vs. Season', col='light blue', xlab='Season', ylab='log(total_cases)')
boxplot(specific_humidity_g_per_kg ~ season, data=dengue1, main='Boxplot of Specific Humidity vs. Season', col='light blue', xlab='Season', ylab='Specific Humidity (g/kg)')
par(mfrow=c(1,1))

#Figure 7: City:season as modulator
bwplot(log(total_cases) ~ season | city, data=dengue1, main='Lattice Boxplot of log(total_cases) vs. Season | City', xlab='Season', ylab='log(total_cases)')

#Figure 8: Hypothesis test for city:season as modulator
n=10000
shuffles = do(n)*(
  shuffle_model = lm(log(total_cases) ~ specific_humidity_g_per_kg + city + season + shuffle(city):shuffle(season), data=dengue1)
)

slope_model = lm(log(total_cases) ~ specific_humidity_g_per_kg + city + season + city:season, data=dengue1)
hist((shuffles$r.squared), breaks=30, col='light blue', main='Histogram of R2 from Shuffled Models', xlab='R2')
abline(v=rsquared(slope_model), col='red', lwd=5)
pval = count(shuffles$r.squared > rsquared(slope_model)) / n
pval

#Figure 9: Dew_point_temp_k as confounder
par(mfrow=c(1,2))
plot(log(total_cases) ~ dew_point_temp_k, data=dengue1, col='blue', main='Plot of log(total_cases) vs. Dew Point', xlab='Dew Point (K)')
plot(specific_humidity_g_per_kg ~ dew_point_temp_k, data=dengue1, col='blue', main='Plot of Specific Humidity vs. Dew Point', xlab='Dew Point (K)', ylab='Specific Humidity (g/kg)')
par(mfrow=c(1,1))

#Final slope model and residual plot
slope_model = lm(log(total_cases) ~ specific_humidity_g_per_kg + city + season + city:season + dew_point_temp_k, data=dengue1)
summary(slope_model)
exp(slope_model$coefficients) #Since this is log(total_sales), we need to interpret as percentage change
simple_anova(slope_model)
plot(resid(slope_model) ~ log(total_cases), data=dengue1, main='Plot of Residuals vs. log(total_cases)', ylab='Residual', col='blue')

#Appendix C
sjt.lm(slope_model)

#Figure 10: Bootstrapped distribution for the partial slope
boots = do(10000)*(
  boot_model = lm(log(total_cases) ~ specific_humidity_g_per_kg + city + season + city:season + dew_point_temp_k, data=resample(dengue1))
)

hist(exp(boots$specific_humidity_g_per_kg), main='Histogram of Bootstrapped Coefficients for Specific Humidity', xlab='Specific Humidity (g/kg)', breaks=45, col='light blue')
exp(confint(boots$specific_humidity_g_per_kg))
abline(v=exp(confint(boots$specific_humidity_g_per_kg)), col='red', lwd=5)
