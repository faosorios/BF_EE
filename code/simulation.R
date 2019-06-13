## These simulations were performed on an HP Proliant DL360 Server equipped with Intel Xeon E5-2630 processor at 2.4 GHz and 16 GB of RAM DDR4.
## The total simulation time was 4 hours, 18 minutes, and 17 seconds.

## reading R sources
source("BF.simul.R")

## loading 'alabama' and 'gmm' packages (required to constrained/unconstrained optimization)
library(alabama)
library(gmm)

## creating Table 1... (it will take long time)

out1 <- summary.EE(Nsize = 5000, coef = c(1,10,0.1), sigma = 0.4, alpha = 0.05)
out2 <- summary.EE(Nsize = 5000, coef = c(1,5.,0.2), sigma = 0.4, alpha = 0.05)
out3 <- summary.EE(Nsize = 5000, coef = c(1,2.,0.5), sigma = 0.4, alpha = 0.05)
out4 <- summary.EE(Nsize = 5000, coef = c(1.,1.,1.), sigma = 0.4, alpha = 0.05)

table1 <- rbind(out1$output, out2$output, out3$output, out4$output)

print(table1 / 100, digits = 2)
# Output:
#    Wald.A Wald.B  BF.A  BF.B  LM.A  LM.B     D
#20   0.420  0.176 0.067 0.064 0.087 0.084 0.087
#50   0.282  0.106 0.065 0.059 0.074 0.068 0.074
#100  0.197  0.077 0.059 0.059 0.061 0.061 0.061
#500  0.104  0.052 0.048 0.048 0.049 0.049 0.049
#20   0.277  0.178 0.068 0.065 0.086 0.083 0.086
#50   0.171  0.108 0.058 0.057 0.068 0.067 0.068
#100  0.127  0.078 0.058 0.058 0.062 0.061 0.062
#500  0.070  0.052 0.047 0.047 0.048 0.048 0.048
#20   0.145  0.175 0.066 0.067 0.082 0.082 0.082
#50   0.096  0.113 0.062 0.057 0.075 0.070 0.075
#100  0.078  0.082 0.056 0.056 0.062 0.062 0.062
#500  0.049  0.055 0.045 0.045 0.050 0.050 0.050
#20   0.140  0.170 0.084 0.070 0.101 0.086 0.101
#50   0.095  0.108 0.062 0.062 0.070 0.070 0.070
#100  0.074  0.080 0.066 0.066 0.067 0.067 0.067
#500  0.055  0.055 0.056 0.056 0.061 0.060 0.061
