

#
library("dplyr")

# load Carseats from ISLR
library("ISLR")
data(package = 'ISLR')
carseats = Carseats # 400 x 11
hitters = Hitters   # 320 x 20

# load the rpart package
library("rpart")
library("rpart.plot")

# model Salary ~ years 
reg.tree = rpart(Salary ~ Years + Hits, data = hitters)
rpart.plot(reg.tree, type = 4)

table(hitters$Years)


# extract the decision boundaries from the regression tree
reg.tree = rpart(Salary ~ Hits, data = hitters)
rpart.plot(reg.tree, type = 4)

# numeric matrix describing the splits:
#     row label: is name of split variable (for purposes of mlike-approx, there
#                is only one name that runs down the rows)
#     count:     the # of obs sent left or right by the split
#     index:     the numeric split point (actual data value)

reg.tree$splits

dim(hitters[!is.na(hitters$Salary),]) # 263 x 20

h = hitters[!is.na(hitters$Salary),]
reg.tree = rpart(Salary ~ Hits, data = h)
rpart.plot(reg.tree, type = 4)
reg.tree$splits

# splits output : count repesents the number of observation sent in the 
# direction that has further splits


dim(h %>% filter(Hits >= 123))              # matches row 2 count
dim(h %>% filter(Hits >= 123, Hits < 160))  # matches row 3 count
dim(h %>% filter(Hits < 160, Hits >= 129))  # matches row 4 count






