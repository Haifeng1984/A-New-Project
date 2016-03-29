#################################################
#Chapter 2 Getting Start (2015/03/12) Haifeng Xu
#################################################
#################################################
#1.1 Workspace
########################
#find out the current working directory
getwd()

#save the long default path before you change it
mine<-getwd()
setwd("c:\\temp")
...
setwd(mine)

#view file names from R
dir("c:\\temp")

#Check for memory assigned to R
memory.size(NA)
memory.limit()

#Set memory
memory.limit(4000)

#List the objects in the current workspace
data(mtcars)
ls()
objects()
object.size(mtcars)
objects("package:base")   #show the names of objects defined in base

#Save the commands history to myfile
savehistory("history_0312")

#Save the workspace to myfile
save.image("workspace_0312")

########################
#1.2 Packages
########################
#The following function shows you where your library is located.
.libPaths()

#The packages installed in R
library()

#The following command tells you which packages are loaded and ready to use
search()

#install command
install.packages("vcd")
update.packages("vcd")

#use packages
library("vcd")
help(package='vcd')
help.start()
install.packages(c("Hmisc","MASS"))
library("Hmisc")
library("MASS")
help(package="Hmisc")
help(package="MASS")

######################################################
#2. R as a Calculator
########################
#2.1 Arithmetic
########################
7+3-5*2
3^2/2
sin(pi/2)
cos(pi/2)
119%/%13
119%%13
15421%%7==0
#Exercise
log(exp(sin(pi/4)^2)*exp(cos(pi/4)^2))

#length
x<-c(1.8,3.14,4,88.169,13)
length(x)
2*x+3
5:1*x+1:5
log(x)

#subsetting
xx<-c("a","b","c","c","d","a")
xx[1]
xx[1:4]
xx[xx>"a"]
u<-xx>"a"
u
xx[u]
x[c(1,4)]
x[-c(2,3,5)]
ones<-rep(1,10)
ones
even<-seq(from=2,to=20,by=2)
even
trend<-1981:2005
trend
c(ones, even)

#Cumulative Sums and Products
x <- c(12,5,13)
cumsum(x)
cumprod(x)
z<-matrix(c(1,5,6, 2,3,13), nrow = 3, ncol = 2, byrow = F)
min(z[,1],z[,2])
pmin(z[,1],z[,2])
nlm(function(x) return(x^2-sin(x)),8)

#Calculus
D(expression(exp(x^2)),"x") # derivative
integrate(function(x) x^2,0,1)

#Matrix Crossproduct
x<-1:5
y<-2*1:5
x*y
t(x)%*%y
crossprod(x,y)
x%*%t(y)
tcrossprod(x,y)
outer(x,y)
A<-matrix(1:6,nrow=2)
A
t(A)
dim(A)
nrow(A)
ncol(A)

#Subsetting a Matrix
B<-A[1,]
C<-A[,1]
D<-A[2,2]
D
class(D)
E<-A[2,2,drop=F]
E
class(E)
F<-A[1, ,drop=F];F
A1<-A[,-2]
A1
det(A1)
eigen(A1)
solve(A1)
A1%*%solve(A1)
diag(4)
diag(2,4,4)
diag(rep(c(1,2),c(3,3)))

B=svd(A1)
B$u%*%diag(B$d)%*%t(B$v)

qr(A1)
cbind(1, A1)
rbind(A1,diag(4,2))
x<-runif(20)
summary(x)
hist(x)

#Exercise
a<-matrix(c(1,1,-1,1),2,2,byrow=T)
b <- c(2,4)
solve(a,b)

#Set Operations
x <- c(1,2,5)
y <- c(5,1,8,9)
union(x,y)
intersect(x,y)
setdiff(x,y)
setdiff(y,x)
setequal(x,y)
setequal(x,c(1,2,5))
2 %in% x
2 %in% y
choose(5,2)

x <- c(13,5,12)
sort(x)
x

x <- c(13,5,12)
x<-sort(x)
########################
#2.2 Complex numbers in R
########################
z<-3.5-8i
Re(z)
Im(z)
Mod(z)
Conj(z)
is.complex(z)
########################
#2.3 Testing for equality with real numbers
########################
x<-sqrt(2)
x*x==2
x*x-2
x<-0.3-0.2
y<-0.1
x==y
identical(x,y)
all.equal(x,y)
########################
#2.4 Three important things
########################
#######################################################
#3. Data structure
########################
attributes(mtcars)
?mtcars
########################
#3.1 Vectors
########################
a<-c(1,2,5,3,6,-2,4)
b<-c("one","two","three")
c<-c(TRUE,TRUE,TRUE,FALSE,TRUE,FALSE)
class(a);class(b);class(c)
a<-c(1,2,5,3,6,-2,4)
a[c(1,3,5)]
a[1:3]
########################
#3.2 Matrices
########################
y<-matrix(1:20,nrow=5,ncol=4)         
y
class(y)
cells<-c(1,26,24,68)
cells
rnames<-c("R1","R2")
cnames<-c("C1","C2")
mymatrix<-matrix(cells,nrow=2,ncol=2,byrow=TRUE,dimnames=list(rnames,cnames))
dim(mymatrix)
mymatrix
y[2, ]
y[,2]
y[1,4]
y[1:2,c(3,4)]
z<-as.vector(y)
z
dim(z)<-c(4,5)
z
########################
#Avoiding Unintended Dimension Reduction
z
r<-z[2,]
r
attributes(z)
attributes(r)
str(z)
str(r)
r <- z[2,, drop=FALSE]
r
dim(r)
########################
#Extended Example: Image Manipulation
library(pixmap)
mtrush1 <- read.pnm("mtrush1.pgm")
mtrush1
plot(mtrush1)
mtrush1@grey[28,88]
locator()
mtrush2 <- mtrush1
mtrush2@grey[84:163,135:177] <- 1
plot(mtrush2)
# adds random noise to img, at the range rows,cols of img; img and the
# return value are both objects of class pixmap; the parameter q
# controls the weight of the noise, with the result being 1-q times the
# original image plus q times the random noise
blurpart <- function(img,rows,cols,q) {
        lrows <- length(rows)
        lcols <- length(cols)
        newimg <- img
        randomnoise <- matrix(nrow=lrows, ncol=lcols,runif(lrows*lcols))
        newimg@grey[rows, cols] <- (1-q) * img@grey[rows, cols] + q * randomnoise
        return(newimg)
}
newimg@grey <- (1-q) * img@grey + q * randomnoise
mtrush3 <- blurpart(mtrush1,84:163,135:177,0.65)
plot(mtrush3)




########################
#3.3 Arrays
########################
z<-(1:24)
dim1<-c("A1","A2")
dim2<-c("B1","B2","B3")
dim3<-c("C1","C2","C3","C4")
z<-array(1:24,c(2,3,4),dimnames=list(dim1, dim2, dim3))
z

#Array Transposition
A<-array(1:24,dim=c(2,3,4));A
B<-aperm(A,c(2,3,1));B
C<-array(1:24,dim=c(3,4,2));C

########################
#3.4 Data frames
########################
patientID<-c(1,2,3,4)
age<-c(25,34,28,52)
diabetes<-c("Type1","Type2","Type1","Type1")
status<-c("Poor","Improved","Excellent","Poor")
patientdata1<-cbind(patientID,age,diabetes,status)
class(patientdata1)
patientdata2<-data.frame(patientID, age, diabetes, status)
class(patientdata2)
patientdata1
patientdata2
class(patientdata1$age)     #error
class(patientdata1$status)  #error
class(patientdata2$age)
class(patientdata2$status)
names(patientdata1)
names(patientdata2)
is.matrix(patientdata1)
is.matrix(patientdata2)
patientdata2[1:2]
patientdata2[c("diabetes","status")]
patientdata2$age

########################
#3.5 Factor
########################
status<-factor(status,ordered=TRUE)
status
str(patientdata2)
summary(patientdata2)
status<-factor(status,order=TRUE,levels=c("Poor","Improved","Excellent"))
status
status2<-c("Poor","Improved","Excellent")
typeof(status2)
status3<-factor(status2,levels=c("Poor","Improved","Excellent"),labels=c("C","B","A"))
status3

#another example
a<-c("A","C","B","C")
b<-as.factor(a)
b[5]<-"D"   #error
c<-as.vector(b)
typeof(c)
c[5]<-"D"
b<-as.factor(c)
b

########################
#3.6 Logical operations
########################
#Exercise:
(TRUE!=FALSE) == !(6==7)
TRUE & c(TRUE,FALSE,FALSE)   
TRUE && c(TRUE,FALSE,FALSE)
TRUE | c(TRUE,FALSE,FALSE)
TRUE || c(TRUE,FALSE,FALSE)
5 > 8 || 6 != 8 && 4 > 3.9
FALSE || TRUE && FALSE
TRUE && 62 < 62 && 44 >= 44
99.99 > 100 || 45 < 7.3 || 4 != 4.0
TRUE && FALSE || 9 >= 4 && 3 < 6
###
isTRUE(6 > 4)

!isTRUE(8 != 5)
!isTRUE(4 < 3)
isTRUE(!TRUE)
isTRUE(3)
isTRUE(NA)
###
identical('Twins', 'twins')

identical(4, 3.1)
!identical(7, 7)
identical('hello', 'Hello')
identical(5 > 4, 3 < 3.1)
###
xor(5 == 6, !FALSE)

xor(!isTRUE(TRUE), 6 > -1)
xor(4 >= 9, 8 != 8.0)
xor(identical(xor, 'xor'), 7 == 7.0)
xor(!!TRUE, !!FALSE)

#Use the which() function to find the indicies of ints that are greater than 7.
ints <- sample(10)
ints
ints>5
which(ints > 7)

########################
#3.7 Type conversions
########################
ab<-c(1,2,3)
is.numeric(ab)
is.vector(ab)
ab<-as.character(ab)
is.numeric(ab)
is.vector(ab)
is.character(ab)
#aaa<-gl(2,5)
#class(aaa)
#mode(aaa)
#typeof(aaa)

########################
#3.8 List
########################
Lst<-list(name="Fred",wife="Mary",child.ages=c(4,7,9))
Lst
Lst[[2]]
Lst[[3]][2]
Lst[["name"]]
Lst$name
#add component
z <- list(a="abc",b=12)
z
z$c <- "sailing" # add a c component
# did c really get added?
z
#You can delete a list component by setting it to NULL.
z$b <- NULL
z
########################
#Partial Matching
x<-list(aardvark=1:5)
x$a
x[["a"]]
x[["a",exact=FALSE]]
########################
findwords <- function(tf) {
        # read in the words from the file, into a vector of mode character    
        txt <- scan(tf,"")   
        wl <- list()    
        for (i in 1:length(txt)) {
                wrd <- txt[i]  # i-th word in input file       
                wl[[wrd]] <- c(wl[[wrd]],i)    
        }    
        return(wl) 
}
txt <- scan("testconcorda.txt","")
txt
class(txt)
findwords("testconcorda.txt")
#######################
#3.9 Formulas
########################
f<-y~x
class(f)
########################
#3.10 Missing Values
########################
x<-c(1,2,NaN,NA,4)
1/0
0/0
is.na(x)
is.nan(x)
########################
#3.11 attach(), detach(), with()
########################
View(mtcars)
View(mpg)
View(mtcars$mpg)
summary(mtcars$mpg)
plot(mtcars$mpg, mtcars$disp)

attach(mtcars)
summary(mpg)
plot(mpg, disp)
detach(mtcars)

with(mtcars, {
        summary(mpg, disp, wt)
        plot(mpg, disp)
        plot(mpg, wt)
        })

with(mtcars, {
        nokeepstats<-summary(mpg)
        keepstats<<-summary(mpg)
        })
nokeepstats
keepstats

####################################################
#4. Data input
########################
#4.1 Input and Output
########################
install.packages("xlsx")
library(xlsx)
auto<-read.xlsx("C:\\Users\\XXXHHF\\Documents\\R\\workfile\\auto.xlsx",
1,header=TRUE,as.data.frame=TRUE)
View(auto)
save(auto,file="C:\\Users\\XXXHHF\\Documents\\R\\workfile\\auto.RData")
remove(auto)
load("C:\\Users\\XXXHHF\\Documents\\R\\workfile\\auto.RData")

########################
#4.2 Data input
########################
#keyboard
#0
x<-1:20
x
#1
mydata<-data.frame(age=numeric(0),gender=character(0),weight=numeric(0))
View(mydata)
mydata<-edit(mydata)
mydata$gender<-NULL
View(mydata)
#2
x<-numeric(10)
data.entry(x)
y<-character(10)
data.entry(y)
y
#3
aa<-scan()

#Scan
read.table(file="rt.txt")
rt<-scan(file="rt.txt")
rt<-scan(file="rt.txt",sep="\n")
rt<-scan(file="rt.txt",sep="\t")
rt<-scan(file=file.choose(),sep="\t")

#a program
rt<-sapply(1:5,function(i)
        as.numeric(na.omit(
                scan("rt.txt",sep="\t",quiet=T)[(4*i-3):(4*i)])))                	
weight<-scan(file="weight.data")
Forest<-scan(file="ForestData.txt",what=double(),skip=10)

#Read.table
Forest<-read.table(file="ForestData.txt")
View(Forest)
Forest<-read.table(file="ForestData.txt",header=TRUE)
View(Forest)
str(Forest)
names(Forest)
Forest<-read.table(file="ForestData.txt",header=TRUE,stringsAsFactors=FALSE)
Forest<-read.table(file="ForestData.txt",header=TRUE,
                   colClass=c("integer","integer","factor",
                              "character","double","integer","double","double","double"))

#You can save time by using read.delim, because then you can omit header=T
Forest<-read.delim(file="ForestData.txt")

#URL
#1
data2<-read.table("http://www.bio.ic.ac.uk/research/mjcraw/therbook/data/cancer.txt",header=T)
urlcancer<-url("http://www.bio.ic.ac.uk/research/mjcraw/therbook/data/cancer.txt")
data2<-read.table(urlcancer,header=T)
head(data2)
class(data2)

#2
urlhandle=url("http://www.math.smith.edu/r/testdata")
ds=readLines(urlhandle)
class(ds)

#3
ds=read.table("http//www.math.smith.edu/r/testdata")

#4
sp500 <- read.csv("http://ichart.finance.yahoo.com/table.csv?s=%5EGSPC&a=03&b=1&c=1999&d=03&e=1&f=2009&g=m&ignore=.csv")
sp500 <- read.csv(paste("http://ichart.finance.yahoo.com/table.csv?",
                        "s=%5EGSPC&a=03&b=1&c=1999&d=03&e=1&f=2009&g=m&ignore=.csv", sep=""))
paste("Hello,","I","am","Xu.",sep=" ")

#Excel
install.packages("xlsx")
library("xlsx")
Forest<-read.xlsx("ForestData.xlsx",1,header=TRUE,as.data.frame=TRUE)
Forest$month
str(Forest)
levels(Forest$month)
Forest$month<-factor(Forest$month,order=TRUE,
                     levels=c("jan","feb","mar","apr","may","jun",
                              "jul","aug","sep","oct","nov","dec"))
levels(Forest$month)

#Stata
install.packages("foreign")
library("foreign")
auto<-read.dta("auto.dta")
head(auto)

#clipboard
#TRData<-read.table("clipboard",header= T,sep="\t")
country countryisocode year ppp
Afghanistan AFG 2010 23.81
Albania ALB 2010 58.26
Algeria DZA 2010 47.36
TRData<-read.table("clipboard",sep=" ",header=T)
TRData<-read.csv("clipboard")

#Built-in data
data(iris)
str(iris)
View(iris)
data(package=.packages(all.available=TRUE))
try(data(package="MASS") )

#Other command
help(CO2)
help(mtcars)
data(mtcars)
help(iris)
########################
#4.3 Data output
########################
#Saving history
history(Inf)
savehistory(file="history.txt")
loadhistory(file="history.txt")

#Saving graphics
plot(mtcars$wt~mtcars$mpg)
pdf("fig1.pdf")
data(mtcars)
plot(mtcars$wt~mtcars$mpg)
dev.off()

#Saving as Ascii file
For.sep<-Forest[Forest$month=="sep", ]
write.table(For.sep,file="For_sep.txt",sep=" ",quote=FALSE,append=FALSE,na="NA")

#Saving as R data file
save(For.sep,file="For_sep.Rdata")
load("For_sep.Rdata")
save(list=ls(all=TRUE),file="all.Rdata")

#Pasting into an Excel spreadsheet
writeClipboard(as.character(factor.name))
writeClipboard(as.character(CO2))
#Checking files from the command line
file.exists("c:\\temp\\Decay.txt")
########################
#4.4 Annotating datasets
########################
names(patientdata2)
names(patientdata2)[2]<-"Age at hospitalization (in years)"
patientdata2[2]

#value labels
gender<-c(1,2,1,2)
patientdata2<-data.frame(patientID, gender, age, diabetes, status)
patientdata2$gender<-factor(patientdata2$gender,levels=c(1,2),labels=c("male","female"))
patientdata2$gender
View(patientdata2)
########################
#Extended Example: Reading PUMS Census Files
extractpums <- function(pf,flds) {
        dtf <- data.frame()  # data frame to be built
        con <- file(pf,"r")  # connection
        # process the input file
        repeat {  
                hrec <- readLines(con,1)  # read Household record
                if (length(hrec) == 0) break  # end of file, leave loop 
                # get household serial number
                serno <- intextract(hrec,c(2,8))  
                # how many Person records?
                npr <- intextract(hrec,c(106,107))  
                if (npr > 0)
                        for (i in 1:npr) {  
                                prec <- readLines(con,1)  # get Person record
                                # make this person's row for the data frame
                                person <- makerow(serno,prec,flds)  
                                # add it to the data frame
                                dtf <- rbind(dtf,person)  
                        }
        }
        return(dtf)
}

# set up this person's row for the data frame
makerow <- function(srn,pr,fl) {
        l <- list()
        l[["serno"]] <- srn
        for (nm in names(fl)) {
                l[[nm]] <- intextract(pr,fl[[nm]])
        }
        return(l)
}

# extracts an integer field in the string s, in character positions
# rng[1] through rng[2]
intextract <- function(s,rng) {  
        fld <- substr(s,rng[1],rng[2])
        return(as.integer(fld))  
}
###############################################
#HOMEWORK
# (a)
college = read.csv("../data/College.csv")
# (b)
fix(college)
rownames(college) = college[,1]
college = college[,-1]
fix(college)
# (c)
# i.
summary(college)
# ii.
pairs(college[,1:10])
# iii.
plot(college$Private, college$Outstate)
# iv.
Elite = rep("No", nrow(college))
Elite[college$Top10perc>50] = "Yes"
Elite = as.factor(Elite)
college = data.frame(college, Elite)
summary(college$Elite)
plot(college$Elite, college$Outstate)
# v.
par(mfrow=c(2,2))
hist(college$Apps)
hist(college$perc.alumni, col=2)
hist(college$S.F.Ratio, col=3, breaks=10)
hist(college$Expend, breaks=100)
# vi.
par(mfrow=c(1,1))
plot(college$Outstate, college$Grad.Rate)
# High tuition correlates to high graduation rate.
plot(college$Accept / college$Apps, college$S.F.Ratio)
# Colleges with low acceptance rate tend to have low S:F ratio.
plot(college$Top10perc, college$Grad.Rate)
# Colleges with the most students from top 10% perc don't necessarily have
# the highest graduation rate. Also, rate > 100 is erroneous!
########################
end
########################








