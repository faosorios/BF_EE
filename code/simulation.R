## These simulations were performed on an HP Proliant DL360 Server equipped 
## with Intel Xeon E5-2630 processor at 2.4 GHz and 16 GB of RAM DDR4.
## The total simulation time was 23 hours, 40 minutes, and 3 seconds.

## reading R sources
source("BF.simul.R")

## loading 'alabama' and 'gmm' packages (required to constrained/unconstrained optimization)
library(alabama)
library(gmm)

## creating Table 1... (it will take long time)
out1 <- summary.EE(Nsize = 5000, coef = c(1,10,0.1), sigma = 0.4, alpha = 0.05, trace = FALSE)
out2 <- summary.EE(Nsize = 5000, coef = c(1,5.,0.2), sigma = 0.4, alpha = 0.05, trace = FALSE)
out3 <- summary.EE(Nsize = 5000, coef = c(1,2.,0.5), sigma = 0.4, alpha = 0.05, trace = FALSE)
out4 <- summary.EE(Nsize = 5000, coef = c(1.,1.,1.), sigma = 0.4, alpha = 0.05, trace = FALSE)

table1 <- rbind(out1$output, out2$output, out3$output, out4$output)

print(table1 / 100, digits = 2)
# Output:
#    Wald.A Wald.B  BF.A  BF.B    LM     D
#20   0.420  0.176 0.067 0.064 0.084 0.087
#50   0.282  0.106 0.065 0.059 0.068 0.074
#100  0.197  0.077 0.059 0.059 0.061 0.061
#500  0.104  0.052 0.048 0.048 0.049 0.049
#20   0.277  0.178 0.068 0.065 0.083 0.086
#50   0.171  0.108 0.058 0.057 0.067 0.068
#100  0.127  0.078 0.058 0.058 0.061 0.062
#500  0.070  0.052 0.047 0.047 0.048 0.048
#20   0.145  0.175 0.066 0.067 0.082 0.082
#50   0.096  0.113 0.062 0.057 0.070 0.075
#100  0.078  0.082 0.056 0.056 0.062 0.062
#500  0.049  0.055 0.045 0.045 0.050 0.050
#20   0.140  0.170 0.084 0.070 0.086 0.101
#50   0.095  0.108 0.062 0.062 0.070 0.070
#100  0.074  0.080 0.066 0.066 0.067 0.067
#500  0.055  0.055 0.056 0.056 0.060 0.061

## First row of Figure 1
delta <- seq(0, 2, length = 101)[-1]
out1 <- summary.power(Nsize = 1000, nobs =  20, coef = c(1,10,0.1), delta = delta, sigma = 0.4, alpha = 0.05)
out2 <- summary.power(Nsize = 1000, nobs =  50, coef = c(1,10,0.1), delta = delta, sigma = 0.4, alpha = 0.05)
out3 <- summary.power(Nsize = 1000, nobs = 100, coef = c(1,10,0.1), delta = delta, sigma = 0.4, alpha = 0.05)

## Figure 1.a (1st panel)
plot(delta, out1$power[,5], type = "l", lwd = 3, xlab = "delta", ylab = "Empirical power", ylim = c(0,100))
lines(delta, out1$power[,3], lwd = 3, col = "red", lty = 2)
lines(delta, out1$power[,4], lwd = 3, col = "blue", lty = 3)

## Figure 1.b (2nd panel)
plot(delta, out2$power[,5], type = "l", lwd = 3, xlab = "delta", ylab = "Empirical power", ylim = c(0,100))
lines(delta, out2$power[,3], lwd = 3, col = "red", lty = 2)
lines(delta, out2$power[,4], lwd = 3, col = "blue", lty = 3)

## Figure 1.c (3rd panel)
plot(delta, out3$power[,5], type = "l", lwd = 3, xlab = "delta", ylab = "Empirical power", ylim = c(0,100))
lines(delta, out3$power[,3], lwd = 3, col = "red", lty = 2)
lines(delta, out3$power[,4], lwd = 3, col = "blue", lty = 3)

## and so on..
