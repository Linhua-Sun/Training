# start R GUI (on Unix, type "R")
# change working directory (File/Change dir)
# load this file (File/Open script)

# '#' starts a commnet line

# R basics

### 1.1.1 as calculator

exp(-2) # press F5 to run this line; press CMD+Enter in Mac

# you can also select a few lines and press F5

rnorm(15) # [1] in the output is the index of the 1st number

?rnorm  # get help for a function

### 1.1.2 assignment

x <- 2 # symbolic variables to represent values

x

# "<-" is a single symbol for "assignment operator"
# spacing around <- is disregarded, but improves readability

# "x < - 2" has different meaning

x + x

# variable names can contain letters, digits, underscore (_), and period (.)
x.1 <- 3

1.x <- 4 # names must not start from digits (or dot followed by digit)

height.1yr <- 3 # use informative names

Height.1yr <- 4 # names are case-sensitive


# avoid the names that are already used by the R system
# c, q, t, C, D, F, I, T, diff, df, pt ...

T # T and F are standard abbreviation for TRUE and FALSE


### 1.1.3 vectorized arithmetic

# R handles entire data vectors as single objects

weight <- c(60, 72, 57, 90, 95, 72) # construct c(...) to define vectors
weight

# there are other ways to read in data

# arithmetic operations of the vectors of the same legnth: element-wise 

height <- c(1.75, 1.80, 1.65, 1.90, 1.74, 1.91)
bmi <- weight / height^2 
# weight in kg, height in m; if using lb and ft, multiply by 4.88

bmi

# arithmetic operations of the vectors of different legnth: the shorter
# vector is recycled

height^2
height^c(2,0,0)
height^c(2,0,0,0) # warning if the two lengths are not multiple

# statistical functions & plotting

t.test(bmi, mu=22.5)

plot(height, weight)

plot(height, weight, pch=2) # change 'plotting character'

# add a line of "normal" weight = 22.5 * height^2

height.new <- seq(from=1.65, to=1.9, by=0.05)
lines(height.new, 22.5 * height.new^2)

# other built-in functions

ls()   # list the objects in the current working space

var(height)   # variance

summary(height)    # summary statistics

curve(x^2 - 10 * x, from=1, to=10)   # plotting a univariate function

curve(expr = sin, 0, 6 * pi)

help.search("optimization")   # search related functions on a topic

RSiteSearch("missing")     # search R-help mailing list and other webpages

rm(height)    # remove an object from working space

# Use the File menu to save and load working space


### 1.2 R language essentials

### 1.2.3 vectors

c("Huey", "Dewey", "Louie")

c('Huey', 'Dewey', 'Louie')  # single quotes are fine, if they are matched


c(T, T, F, T)

bmi > 25


### 1.2.6 functions that create vectors

x <- c(1, 2, 3)   # concatenate
y <- c(10, 20)

c(x, y, 5)

# assign names to the vector elements, useful for dislay

x <- c(red="Huey", blue="Dewey", green="Louie")
x

names(x)  # extracting names

# all elements of a vector have the same type

# concatenating different types will convert to the least 
# "restrictive" types
c(FALSE, 3)  

c(pi, "abc")

c(FALSE, "abc")


# creating a sequence of equidistant numbers, useful in plots

seq(4, 9)

seq(4, 10, 2)

4:9  # for step size 1


# replicating repeat values

oops <- c(7, 9, 13)

rep(oops, 3)  # 2nd argument is a single number

rep(oops, 1:3) # 2nd argument is a vector

rep(oops, 1:2)  # invalid

# useful for group codes

rep(1:2, c(10, 15))  # 10 men first and 15 women next

rep(1:2, each=10)


### 1.2.7 matrices and arrays

x <- 1:12

dim(x) <- c(3, 4)

x    # treats the vector as a 3 by 4 matrix, column-wise


x <- matrix(1:12, nrow=3, byrow=T)  # another way
x

# functions that operate on matrices

rownames(x) <- LETTERS[1:3]  
x

colnames(x) <- 1:4
x

colnames(x)

t(x)  # transposition

# more built-in variables
letters

month.name

month.abb


# glue vectors together

cbind(A=1:4, B=5:8, C=9:12)   # column-wise

rbind(A=1:4, B=5:8, C=9:12)   # row-wise



### 1.2.9 lists

# pre- and postmenstrual energy intake in a group of women

intake.pre <- c(5260,5470,5640,6180,6390,6515,6805,7515,7515,8230,8770)

intake.post <- c(3910,4220,3885,5160,
	5645,4680,5265,5975,6790,6900,7335)
# can use ESC to to stop the waiting sign +

mylist <- list(before=intake.pre, after=intake.post)
mylist

mylist$before   # extracting named components


### 1.2.10 data frames

# create from pre-existing variables

my.frame <- data.frame(intake.pre, intake.post)

my.frame   # data in rows are paired

my.frame$intake.pre   # extract components


### 1.2.11 indexing

intake.pre[5]   # get single element

# can also be used for modifying elements of a vector
intake.pre[5] <- 6390

intake.pre[c(3,5,7)]   # index with a vector

# works if the index vector is a variable
idx <- c(3,5,7)
intake.pre[idx]


intake.pre[1:5]   # get a sequence of elements

intake.pre[-c(3,5,7)]   # negative indexing: all except the 3 elements


### 1.2.12 conditional selection

intake.post[intake.pre > 7000]

intake.pre > 7000   # logical vector of the same length as intake.post

# combine expression

intake.post[intake.pre > 7000 & intake.pre <= 8000]

intake.pre > 7000 & intake.pre <= 8000



### 1.2.13 indexing of data frames

# extract using matrix-like structure

my.frame <- data.frame(intake.pre, intake.post)

my.frame[5,1]  # 5th row, 1st column

my.frame[5,]  # all values for woman no. 5 (all columns for row 5) 

my.frame[2]   # comma above is required


# extract all data for cases that satisfy some criterion
my.frame[my.frame$intake.pre>7000,]   # row names are from the original data frame

my.frame$intake.pre>7000


# check the first few cases
my.frame[1:2,]

head(my.frame)

tail(my.frame)



### flow control



### slide: Fibonacci sequence

Fibonacci <- numeric(12) # allocate vector variable

## Student question 2014.4: why not numeric(0)?

Fibonacci[1] <- Fibonacci[2] <- 1

for (i in 3:12) {
	Fibonacci[i] <- Fibonacci[i - 2] + Fibonacci[i - 1]
}

Fibonacci


### another for() example

x <- seq(0, 1, .05)
plot(x, x, ylab="y", type="l")
for (i in 2:8)
	lines(x, x^i)


### another for() example
# http://www.statmethods.net/graphs/scatterplot.html

install.packages("scatterplot3d") # use a new package from online
# or use menu "Packages/Install package from local zip files".
# Install only once.

library(scatterplot3d)

attach(mtcars)
scatterplot3d(wt,disp,mpg, main="3D Scatterplot")

for (angle in 1:180)
	scatterplot3d(wt,disp,mpg, angle=angle, main="3D Scatterplot")



### slide: if() statement: code example

# Example: a function to add a scattorplot of the data and compute 
# correlations

cor.plot <- function(x, y, plot_it) {
  if (plot_it == TRUE) {
    plot(x, y)
  }
  cor(x, y)
}

cor.plot  # confirm the function is defined

cor.plot(c(2, 5, 7), c(5, 6, 8), FALSE)
cor.plot(c(2, 5, 7), c(5, 6, 8), TRUE)



# any() example
any(1:10 > 5)
any(1:10 > 50)


### one comment about if() format

x <- 2
if (x > 2) { print("x > 2") }
else {print("x <= 2") }   # causing error

# correct form: R finds an incomplete line and will collect 
# all lines before running

if (x > 2) { print("x > 2") 
} else {print("x <= 2") }   # move the brace here

# better form: no confusion for R or you

if (x > 2) { 
	print("x > 2") 
} else {
	print("x <= 2") 
}


### slide: while loop

x.total <- 0
x.count <- 0
while (x.total < 100) {
	x.total <- x.total + runif(1)
	x.count <- x.count + 1
}

list(total=x.total, count=x.count)


### slide: Example: Newton's method for root finding

x0 <- 1

x <- x0
f <- x^3 + 2 * x^2 - 7
tolerance <- 0.000001

while (abs(f) > tolerance) {
  f.prime <- 3 * x^2 + 4 * x
  x <- x - f / f.prime
  f <- x^3 + 2 * x^2 - 7
  print(x)
}

x

# confirm with plot

curve(x^3 + 2 * x^2 - 7, -5, 5)
abline(h = 0, col = "blue")
abline(v = x, col = "red")

# advanced topic: start values as a vector
x0 <- seq(-3, 3, length = 30) 

x <- x0
f <- x^3 + 2 * x^2 - 7
tolerance <- 0.000001

while (max(abs(f)) > tolerance) {   # note the change
  f.prime <- 3 * x^2 + 4 * x
  x <- x - f / f.prime
  f <- x^3 + 2 * x^2 - 7
  #print(x)
}

x




### 2.2.1 plot layout

# override default axis labels, add titles

x <- runif(50, 0, 2)
y <- runif(50, 0, 2)

plot(x, y, main="Main title", sub="subtitle", 
	xlab="x-label", ylab="y-label")

# add points and lines inside the plotting region

text(0.6, 0.6, "text at (0.6, 0.6)")

abline(h=.6, v=.6)   # show the text is centered

# margin coordinates

paste("side", 1:4)
mtext(paste("side", 1:4), side=1:4, line=2, font=2)

paste("line", -1:2)
for (side in 1:4) 
	mtext(paste("line", -1:2), side=side, at=.7, line=-1:2)


### 2.2.2 build a plot from pieces

plot(x, y, type="n", xlab="", ylab="", axes=F)  
  # nothing is drawn, but sets up the plotting region and coordinate systems

points(x, y)
axis(1)
axis(2, at=seq(0.2, 1.8, 0.2))

box()
title(main="Main title", sub="subtitle", xlab="x-label", ylab="y-label")


### 2.2.4 combining plots

x <- rnorm(50)
hist(x, freq=F)
curve(dnorm(x), add=T)  # sometimes the top of the curve may be chopped off

# solution: determine the range of both curve and histgoram

h <- hist(x, plot=F)
h
ylim <- range(0, h$density, dnorm(0))
hist(x, freq=F, ylim=ylim)
curve(dnorm(x), add=T)



### load the ISwR package

install.packages("ISwR")  # only run once

library(ISwR)

help(package=ISwR)

search()   # currently loaded packages



### 4.1 Summary statistics

x <- rnorm(50)

median(x)

quantile(x)

# get other percentiles

percent.vec <- seq(0, 1, 0.1)

quantile(x, percent.vec)


# missing value handling

x2 <- c(rnorm(50), NA, rnorm(50))

mean(x2)  # R will not skip missing values by default

mean(x2, na.rm = T)    # NA remove

length(x2)   # counts NA

sum(!is.na(x2))   # TRUE/FALSE are converted to 1/0 in arithmetic

# summary display

summary(x2)    


# can be applied to data.frame too

juul    # a data set in ISwR, concerning serum IGF-I

?juul

summary(juul)    # don't write as juu1

str(juul)      # display the internal structure of an R object

# change a few variables to factors

juul$sex
juul$sex <- factor(juul$sex, labels=c("M", "F"))
juul$sex

juul$menarche <- factor(juul$menarche, labels=c("No", "Yes"))

juul$tanner <- factor(juul$tanner, labels=c("I", "II", "III", "IV", "V"))

summary(juul)


### display of distributions

### 4.2.1 histograms

hist(x)

hist(x, breaks=5)

# a data set on accident rates by age group
# counts in age groups 0-4, 5-9, 10-15, 16, 17, 18-19, 20-24, 25-59, 60-79

mid.age <- c(2.5, 7.5, 13, 16.5, 17.5, 19, 22.5, 44.5, 70.5)
acc.count <- c(28, 46, 58, 20, 31, 64, 149, 316, 103)

age.acc <- rep(mid.age, acc.count)  # generate pseudo-data

brk <- c(0, 5, 10, 16, 17, 18, 20, 25, 60, 80)
hist(age.acc, breaks=brk)

hist(age.acc)   # for comparison


### 4.2.2. Empirical cumulative distribution

n <- length(x)

plot(sort(x), (1:n)/n, type="s", ylim=c(0,1))

points(sort(x), (1:n)/n)   # for comparison


### 4.2.3 Q-Q plots

qqnorm(x)

qqline(x)   # passes through the first and third quartiles. 

y <- rt(200, df = 5)
qqnorm(y)
qqline(y, col = "red")

qqplot(y, rt(300, df = 5))
qqline(y, col = "blue")

# use density plot to compare
plot(density(x))
points(density(y), type="l", col="blue")

qqplot(x, y)

### Student question: why is a straight line from qqnorm plot indicates
# a normal distribution?

par(mfrow=c(1,1))
x <- rnorm(1000, 10, 20)
qqnorm(x)
qqline(x)

n <- length(x)
normal_sample <- rnorm(n)  # generate a vector of the same length

# compare density
plot(density(x))
points(density(normal_sample), type="l", col="blue")

# standardize x and compare: x - mean / std
x_std <- (x - mean(x) ) / sd(x)
plot(density(x_std))
points(density(normal_sample), type="l", col="blue")

# compare empirical culmulative distribution
plot(sort(x_std), (1:n)/n, type="l", ylim=c(0,1))
points(sort(normal_sample), (1:n)/n, type="l", col="blue")

# compare the two sorted vectors to each other
plot(sort(normal_sample), sort(x_std))

plot(sort(normal_sample), sort(x))  # the same plot except y-axis value

# qqplot() to do the similar
qqplot(normal_sample, x) # two lengths can be different

# qqnorm() with the theoretical quantiles
qqnorm(x)

# now try an non-normal distribution
x <- rchisq(1000, 1)
hist(x)


### 4.2.4 boxplots

par(mfrow=c(1,2))

boxplot(IgM)   # data from ISwR

boxplot(log(IgM))

hist(IgM)      # for comparison
hist(log(IgM))

par(mfrow=c(1,1))   # restore the single plot layout





### save images to file

? pdf    # jpeg(), png()

pdf("pie.pdf")
pie(caff.marital["Single",], main="Single", col=slices.col)
dev.off()   # close file

### slide: 3.1 random sampling

sample(1:40, 5)   # pick 5 numbers at random from 1:40

sample(c("H", "T"), 10, replace=T)  # sampling with replacement

# non-equal probabilities for the outcomes

sample(c("H", "T"), 10, replace=T, prob=c(0.9, 0.1))



### slide: 3.5.1 densities

x <- seq(-4, 4, 0.1)
plot(x, dnorm(x), type="l")

curve(dnorm(x), from=-4, to=4, add=T, col="blue")   # alternative way

# discrete: bionomial distribution, as histogram-like plot

x <- 0:50
plot(x, dbinom(x, size=50, prob=.33), type="h")  


### slide: 3.5.2 cumulative distribution function

1 - pnorm(160, mean=132, sd=13)


### slide: CDF: another example

1 - pbinom(16, size=20, prob=0.5)  # the upper tail probability 

# is this correct?



# the above is prob(17 or more), we need prob(16 or more)

1 - pbinom(15, size=20, prob=0.5)


# two tailed test

(1 - pbinom(15, 20, 0.5)) + pbinom(4, 20, 0.5)


# easier way
?binom.test

binom.test(16, 20, 0.5, alternative="greater")


### slide: 3.5.3 quantiles

# inverse relationship between quantiles and cdf
x <- seq(-4, 4, 0.1)
plot(x, pnorm(x), type="l")

abline(v=1.5, h=pnorm(1.5))

pnorm(1.5)

qnorm(pnorm(1.5))


### slide: qunatiles example

xbar <- 83
sigma <-12
n <- 5
se.mean <- sigma/sqrt(n)  # standard error of the mean
se.mean

xbar + se.mean * qnorm(0.025)
xbar + se.mean * qnorm(0.975)


### slide: 3.5.4 random numbers

rnorm(10)

rnorm(10)

rnorm(10, mean=7, sd=5)

rbinom(10, size=20, prob=.5)

set.seed(565228369)
rnorm(10)

save.image()   # save workspace to ".RData"



### Buffon's needle, Cheng Li's code

# assume lines are horizontal at integers
par(mar=c(2,2,2,2))
plot(0:10, 0:10, type="n", xlab="", ylab="", asp=1)
abline(h=0:10)

# drop a random line with length 1

set.seed(10000) # random seed

one_trial <- function() {

  rand.point.x <- runif(1) * 8 + 1 # random value in [1,9]
  rand.point.y <- runif(1) * 8 + 1 # random value in [1,9]

  rand.angle <- runif(1) * 2 * pi   # random angle in [0,2*Pi]
  end.point.x <- rand.point.x + cos(rand.angle) * 1
  end.point.y <- rand.point.y + sin(rand.angle) * 1

  lines(c(rand.point.x, end.point.x), c(rand.point.y, end.point.y))
  intersected <- floor(rand.point.y) != floor(end.point.y)

  intersected
}

num.trial = 100

#num.trial = 10000

buffon.trial <- replicate(num.trial, one_trial())

prob.intersect <- sum(buffon.trial) / num.trial
pi.estimate = 2 / prob.intersect 

print(paste("my estimate of pi is", pi.estimate, "after", num.trial, "trials"))


## Student question 2014.4: if we use rnorm() instead of runif(), 
## can we still estimate pi?


### 2.1.1 save and load workspace

x <- runif(50, 0, 2)
y <- runif(50, 0, 2)

ls()
save.image()   # save all objects to ".RData" in the current working dir
  # same as menu "File/Save Workspace"

rm(list=ls())  # remove the objects
ls()

load(".RData")  # load workspace, same as menu "File/Load Workspace"
ls()




### Exercise "Three example plots"

par(mfrow=c(2, 2))

# plot 1

?volcano
dim(volcano)
volcano[1:10, 1:5]

z <- 2 * volcano        # Exaggerate the relief
x <- 10 * (1:nrow(z))   # 10 meter spacing (S to N)
y <- 10 * (1:ncol(z))   # 10 meter spacing (E to W)

# Don't draw the grid lines :  border = NA
par(mar=rep(0, 4))

persp(x, y, z, theta = 135, phi = 30, col = "red", scale = FALSE,
      ltheta = -120, shade = 0.75, border = NA, box = FALSE)

mtext("persp()", side=3, line=-2)

# plot 2

par(mar=rep(0.5, 4))
contour(x, y, z, asp=1, labcex=0.35, axes=FALSE)
rect(0, 0, 870, 620)
mtext("contour()", side=3, line=-1.5)

# plot 3

par(mar=c(3, 3, 2, 0.5))

? trees
# Note that  example(trees)  shows more sensible plots!
N <- nrow(trees)
attach(trees)

# Girth is diameter in inches
?symbols
symbols(Height, Volume, circles=Girth/24, inches=FALSE,
        main="", xlab="", ylab="", bg=grey(Girth/max(Girth)))
mtext("symbols()", side=3, line=0.5)


### Interactive spinning 3D Scatterplot

install.packages("rgl") # install only once
library(rgl)

attach(mtcars)
plot3d(wt, disp, mpg, col="red", size=3) 

example(plot3d)  # color as additional dimension
?plot3d  # example() runs the example code in this help page

example(surface3d)

data(package="ISwR")



### Exercise: Graphs

### Task 1, example code of dividing plot regions

N <- 2000
x <- rnorm(N)
op <- par(mar=c(0,0,0,0), oma=c(0,0,0,0)+.1)
layout(matrix(c(1,1,1,2), ncol=1))
  # divide the plot area to 4 vertical layers, figure 1 occupies 
  # the top 3 layers, and figure 2 occupies the bottom 1 layer


## divide the device into two rows and two columns
## allocate figure 1 all of row 1
## allocate figure 2 the intersection of column 2 and row 2
layout(matrix(c(1,1,0,2), 2, 2, byrow = TRUE))
## show the regions that have been allocated to each plot
layout.show(2)

plot(1:10)
plot(1:10)


### Task 2, explain what each line of code does to generate the plot.

?polygon
?rev

op <- par(mfrow=c(1,2), mar=c(3,2,4,2)+.1)
do.it <- function (x, xlab="", ylab="", main="") {
  d <- density(x)
  plot(d, type='l', xlab=xlab, ylab=ylab, main=main)
  q <- quantile(x)
  draw.area <- function (i, col) {
    x <- d$x[i]
    y <- d$y[i]
    polygon( c(x,rev(x)), c(rep(0,length(x)),rev(y)), border=NA, col=col )
  }
  draw.area(d$x <= q[2], 'red')
  draw.area(q[2] <= d$x & d$x <= q[3], 'green')
  draw.area(q[3] <= d$x & d$x <= q[4], 'blue')
  draw.area(d$x >= q[4], 'yellow')
  lines(d, lwd=3)
}
do.it( rnorm(2000), main="Gaussian" )
do.it( rexp(200), main="Exponential" )
par(op)

