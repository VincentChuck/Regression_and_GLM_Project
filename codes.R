###import the data
setwd("~/Documents/ST300")
body = read.table('body.dat.txt',header = T)
attach(body)

###carry out variable selection using BIC since it suggests a smaller model
#fit all possible covariates
m1 = lm(height~., data=body)
summary(m1)
BIC((m1))
step(m1,direction='both', k=log(507))

###optimal model chosen by stepwise BIC selection
model = lm(height ~ biacromial + pelvic.breadth + knee.diam + 
             ankle.diam + chest.girth + waist.girth + thigh.girth + forearm.girth + 
             calf.girth + weight + gender)

cbind(BIC(m1),BIC(model))
summary(model) #all the fitted covariates have significant t-values at 5% 
               #significant level




###choose variables with significant t-values
#all of the are significant


#check assumptions (NICE-L)
qqnorm(model$resid) ; qqline(model$resid)


#is there a need to do transformation?
#NO !!!!!!!!


###deal with outliers, high leverage points, Cook's D
#outliers
plot(model$fitted,rstandard(model)) ; abline(3,0,lty=3) ; abline(-3,0,lty=3) 
sum(abs(rstandard(model))>2.83) #there are 7 outliers
index = 1:507
o = index[abs(rstandard(model))>2.83]
#o = index[abs(rstandard(model))>3]
no.outliers = lm(height[-o] ~ biacromial[-o] + pelvic.breadth[-o] + knee.diam[-o] + 
                   ankle.diam[-o] + chest.girth[-o] + waist.girth[-o] + 
                   thigh.girth[-o] + forearm.girth[-o] + calf.girth[-o] + 
                   weight[-o] + gender[-o])

plot(no.outliers$fitted,rstandard(no.outliers)) ; abline(3,0,lty=3) ; abline(-3,0,lty=3) 
sum(abs(rstandard(no.outliers))>3)  #no more outliers!!

summary(no.outliers)
cbind(BIC(model),BIC(no.outliers)) #improved BIC and adj.R^2 !!
cbind(summary(model)$adj.r.squared, summary(no.outliers)$adj.r.squared)

qqnorm(no.outliers$resid) ; qqline(no.outliers$resid)

#leverage
sum(hatvalues(model)>2*12/507) #28 high leverage points
index = 1:507
L = index[hatvalues(model)>2*12/507]
no.leverage = lm(height[-L] ~ biacromial[-L] + pelvic.breadth[-L] + knee.diam[-L] + 
                   ankle.diam[-L] + chest.girth[-L] + waist.girth[-L] + 
                   thigh.girth[-L] + forearm.girth[-L] + calf.girth[-L] + 
                   weight[-L] + gender[-L])
sum(hatvalues(no.leverage)>2*12/507)

cbind(BIC(model),BIC(no.leverage)) 
cbind(summary(model)$adj.r.squared, summary(no.leverage)$adj.r.squared)
qqnorm(no.leverage$resid) ; qqline(no.leverage$resid)

#Cook's Distance
sum(cooks.distance(model)>4/507) #29 points with high Cook's Distance
J = index[cooks.distance(model)>4/507]
no.cook = lm(height[-J] ~ biacromial[-J] + pelvic.breadth[-J] + knee.diam[-J] + 
               ankle.diam[-J] + chest.girth[-J] + waist.girth[-J] + 
               thigh.girth[-J] + forearm.girth[-J] + calf.girth[-J] + 
               weight[-J] + gender[-J])
sum(cooks.distance(no.cook)>4/507)

cbind(BIC(model),BIC(no.cook)) 
cbind(summary(model)$adj.r.squared, summary(no.cook)$adj.r.squared)
qqnorm(no.cook$resid) ; qqline(no.cook$resid)



#Leverage and Cook's Distance
#45 data points with high leverage or Cook's Distance
sum(hatvalues(model)>2*12/507 | cooks.distance(model)>4/507) 

#removing points with high leverage and Cook's Disance
index = 1:507
I = index[hatvalues(model)>2*12/507]
J = index[cooks.distance(model)>4/507]
K = sort(union(I,J))
length(K)
model_new = lm(height[-K] ~ biacromial[-K] + pelvic.breadth[-K] + knee.diam[-K] + 
                 ankle.diam[-K] + chest.girth[-K] + waist.girth[-K] + 
                 thigh.girth[-K] + forearm.girth[-K] + calf.girth[-K] + 
                 weight[-K] + gender[-K])

sum(hatvalues(model_new)>2*12/507 & cooks.distance(model_new)>4/507) 

#compare the results
cbind(BIC(model),BIC(model_new)) 
cbind(summary(model)$adj.r.squared, summary(model_new)$adj.r.squared)

#compare the qqplots
par(mfrow = c(1,2))
qqnorm(model$resid, main = "Original") ; qqline(model$resid)
qqnorm(no.leverage$resid, main = "No High Leverage Points") ; qqline(no.leverage$resid)
qqnorm(no.cook$resid, main = "No Points with High Cook's D") ; qqline(no.cook$resid)
qqnorm(model_new$resid, main = "No Influential Points") ; qqline(model_new$resid)

#compare the standardised residuals plot
par(mfrow = c(1,2))
plot(model$fitted,rstandard(model), main = "Original") ; abline(3,0,lty=3) ; abline(-3,0,lty=3) 
plot(no.leverage$fitted,rstandard(no.leverage), main = "No High Leverage Points") ; abline(3,0,lty=3) ; abline(-3,0,lty=3) 
plot(no.cook$fitted,rstandard(no.cook), main = "No Points with High Cook's D") ; abline(3,0,lty=3) ; abline(-3,0,lty=3) 
plot(model_new$fitted,rstandard(model_new), main = "No Points with High Leverage and Cook's D") ; abline(3,0,lty=3) ; abline(-3,0,lty=3) 


###############
#Final Model: #
###############
model_new = lm(height[-K] ~ biacromial[-K] + pelvic.breadth[-K] + knee.diam[-K] + 
                 ankle.diam[-K] + chest.girth[-K] + waist.girth[-K] + 
                 thigh.girth[-K] + forearm.girth[-K] + calf.girth[-K] + 
                 weight[-K] + gender[-K])
summary(model_new)
anova(model_new)



#try removing points with high leverage, cook's d and outliers
L = sort(union(I,J,o))
sum(hatvalues(model)>2*12/507 | cooks.distance(model)>4/507 | 
      abs(rstandard(model))>2.83) 
model_no3 = lm(height[-L] ~ biacromial[-L] + pelvic.breadth[-L] + knee.diam[-L] + 
                 ankle.diam[-L] + chest.girth[-L] + waist.girth[-L] + 
                 thigh.girth[-L] + forearm.girth[-L] + calf.girth[-L] + 
                 weight[-L] + gender[-L])

sum(hatvalues(model_no3)>2*12/507 | cooks.distance(model_no3)>4/507 | 
      abs(rstandard(model_no3))>2.83) 
cbind(BIC(model),BIC(model_no3)) 
cbind(summary(model)$adj.r.squared, summary(model_no3)$adj.r.squared)
qqnorm(model_no3$resid) ; qqline(model_new$resid)

par(mfrow = c(3,2))
plot(model$fitted,rstandard(model), main = "Original") ; abline(3,0,lty=3) ; abline(-3,0,lty=3) 
plot(no.leverage$fitted,rstandard(no.leverage), main = "No High Leverage Points") ; abline(3,0,lty=3) ; abline(-3,0,lty=3) 
plot(no.cook$fitted,rstandard(no.cook), main = "No Points with High Cook's D") ; abline(3,0,lty=3) ; abline(-3,0,lty=3) 
plot(model_new$fitted,rstandard(model_new), main = "No Points with High Leverage and Cook's D") ; abline(3,0,lty=3) ; abline(-3,0,lty=3) 
plot(no.outliers$fitted,rstandard(no.outliers), main = "No Outliers") ; abline(3,0,lty=3) ; abline(-3,0,lty=3) 
#removing points with high leverage and cook's d is enough

#deal with multicollinearity
library(car)
vif(model) #see that weight has a really high vif value
#try dropping weight
noweight = lm(height ~ biacromial + pelvic.breadth + knee.diam + 
              ankle.diam + chest.girth + waist.girth + thigh.girth + forearm.girth + 
              calf.girth + gender)
vif(noweight)
cbind(BIC(m1),BIC(model),BIC(noweight))
cbind(summary(m1)$adj.r.squared, summary(model)$adj.r.squared, 
      summary(noweight)$adj.r.squared)  #dropping weight reduces the adj r^2 too much

nochest.girth = lm(height ~ biacromial + pelvic.breadth + knee.diam + 
            ankle.diam + waist.girth + thigh.girth + forearm.girth + 
            calf.girth + weight + gender) #try removing chest.girth
cbind(summary(m1)$adj.r.squared, summary(model)$adj.r.squared, 
      summary(noweight)$adj.r.squared,summary(nochest.girth)$adj.r.squared) #adj R^2 does not decrease much
vif(nochest.girth)  #but vif of weight is still 19

nochest.girth.and.nowaist.girth = lm(height ~ biacromial + pelvic.breadth + knee.diam + 
            ankle.diam + thigh.girth + forearm.girth + 
            calf.girth + weight + gender) #try removing chest.girth and waist.girth
vif(nochest.girth.and.nowaist.girth)
summary(mod3)
cbind(BIC(m1),BIC(model),BIC(noweight),BIC(nochest.girth),BIC(nochest.girth.and.nowaist.girth))
cbind(summary(m1)$adj.r.squared, summary(model)$adj.r.squared, 
      summary(noweight)$adj.r.squared,summary(nochest.girth)$adj.r.squared,
      summary(nochest.girth.and.nowaist.girth)$adj.r.squared)

#try removing chest.girth and forearm.girth

Height = (205.9880889 + 0.4814908*biacromial + 0.3256312*pelvic.breadth 
        - 0.7208400*knee.diam + 0.7922470*ankle.diam - 0.2627171*chest.girth 
        - 0.6808899*waist.girth - 0.6741001*thigh.girth - 0.9639849*forearm.girth 
        - 0.4847556*calf.girth + 1.3828323*weight + 4.6176846*gender)
summary(Height)