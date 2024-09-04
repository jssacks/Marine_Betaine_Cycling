install.packages("investr")


##Confidence intervals for fits using "investr"

library(tidyverse)
library(investr)
data <- tibble(date = 1:7,
               cases = c(0, 0, 1, 4, 7, 8.5, 8.5))

model <- nls(cases ~ SSlogis(log(date), Asym, xmid, scal), data= data )
new.data <- data.frame(date=seq(1, 10, by = 0.1))
interval <- as_tibble(predFit(model, newdata = new.data, interval = "confidence", level= 0.9)) %>% 
  mutate(date = new.data$date)

p1 <- ggplot(data) +  geom_point(aes(x=date, y=cases),size=2, colour="black") + xlab("Date") + ylab("Cases")  

p1+
  geom_line(data=interval, aes(x = date, y = fit ))+
  geom_ribbon(data=interval, aes(x=date, ymin=lwr, ymax=upr), alpha=0.5, inherit.aes=F, fill="blue")+
  theme_classic()


#Predict 95% confidence intervals using "propogate"
install.packages("propagate")
library("propagate")

x  <- c(25, 25, 10, 10, 5, 5, 2.5, 2.5, 1.25, 1.25)
y <- c(0.0998, 0.0948, 0.076, 0.0724, 0.0557,
       0.0575, 0.0399, 0.0381, 0.017, 0.0253)

m <- nls(y ~ SSmicmen(x, Vm, K), trace = TRUE)

x1 <- seq(0, 25, length = 100)
plot(x, y, xlim = c(0, 25), ylim = c(0, 0.1))
#lines(x, predict(m, data.frame(S = x1)), col = "red")

y.conf <- predictNLS(m, newdata=data.frame(x=x1), interval="confidence", alpha=0.05, nsim=10000)$summary
y.pred <- predictNLS(m, newdata=data.frame(x=x1), interval="prediction", alpha=0.05, nsim=10000)$summary

matlines(x1, y.conf[,c("Sim.2.5%", "Sim.97.5%")], col="red", lty="dashed")
matlines(x1, y.pred[,c("Sim.2.5%", "Sim.97.5%")], col="blue", lty="solid")















