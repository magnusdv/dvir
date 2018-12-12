library(donlib)
loadPedSuite(github=TRUE)


#generate
#' Find next candidates
#' 
#' New moves are generated, NULL
#' if none are possible, for one and two step
#' 
#' @param pm A list of singeltons 
#' @param am A list of pedigrees 
#' @param vp Character vector with names of victims
#' @param mp Character vector with names of missing persons
#' @param one logic If TRUE, only one-step
#' @return list of matrices of one and two-moves
#' @examples 
#' vp = c("V1", "V2")
#' pm = singletonList(2, ids = vp)
#' m = marker(pm[[1]], alleles = 1:2, "V1" = 1)
#' pm[[1]] = addMarkers(pm[[1]], m)
#' m = marker(pm[[2]], alleles = 1:2, "V2" = 1:2)
#' pm[[2]] = addMarkers(pm[[2]], m)
#' mp = c("MP1", "MP2")
#' am = nuclearPed(2, children = mp, father = "R1", mother = "R2")
#' m = marker(am, "R1" = 1:2, "R2"= 1:2)
#' am = addMarkers(am, m)
#' # The thre possibilities: a move of two, one or none:
#' tab = generate(pm, am, vp, mp)
#' pm[[1]] = swapSex(pm[[1]], "V1")
#' tab = generate(pm, am, vp,mp)
#' pm[[2]] = swapSex(pm[[2]], "V2")
#' tab = generate(pm, am, vp,mp)
#' @export
#' 





#' 
#reduce

tab = tab[!tab[,1] %in% from, ]
tab = tab[!tab[,2] %in% to, ]
n = 5
pm = singletonList(n)
vp = unlist(lapply(pm, function(x) x$ID))
for (i in 1:n){
  m = marker(pm[[i]], name = "L1")
  pm[[i]] = addMarkers(pm[[i]], m)
}
amat = matrix( rep(1,n*2), ncol = 2)
pm = setAlleles(pm, alleles = amat)

mp = paste("MP", 1:n, sep = "")
am = nuclearPed(n, father ="R1", mother ="R2", 
                 children = mp, sex = c(2,2,1,1,1))
m = marker(am, "R1" = 1, "R2" =1:2, alleles = 1:2)
am = addMarkers(am, m)
windows()
plotPedList(list(pm,am), newdev=F,marker = 1)
###
vp = c("V1", "V2", "V3", "V4")
#pm = singletonList(2)
pm = singletonList(4, ids = vp, sex = c(1,1,1,2))
m = marker(pm[[1]], alleles = 1:2, "V1" = 1)
pm[[1]] = addMarkers(pm[[1]], m)
m = marker(pm[[2]], alleles = 1:2, "V2" = 1)
pm[[2]] = addMarkers(pm[[2]], m)
m = marker(pm[[3]], alleles = 1:2, "V3" = 2)
pm[[3]] = addMarkers(pm[[3]], m)
m = marker(pm[[4]], alleles = 1:2, "V4" = 1)
pm[[4]] = addMarkers(pm[[4]], m)

mp = c("MP1", "MP2", "MP3")
am = nuclearPed(3, children = mp, father = "R1", mother = "R2", sex = c(1,1,2))
m = marker(am, "R1" = 1, "R2"= 1, alleles = 1:2)
am = addMarkers(am, m)

#dviForward2(pm, am, vp, mp, LRlimit = -1)
plotPedList(list(pm, am), marker = 1)
likNULL = prod(LR(list(pm, am), 1)$likelihoodsPerSystem)
reduce (pm, am, vp, mp, likNULL = likNULL, limit = 0)


# Possible actions
# tab = NULL terminate
# dim(tab)[2] ==2 Move one, update
# dim(tab)[2] > 2 Move two, find likMaxTwo
# and make two one steps likMax2, determine best and update
# if(!is.null(tab)){
#   # First do a one step
#   
#   d2 = dim(tab)[2]
#   if (d2 > 2) {
#     res1 = apply(tab ,1,function(x, pm, am) 
#       lik1(x[1:2] , x[3:4], pm, am), pm, am)
#     likres1 = unlist(lapply(res1, function(x) x$lik))
#     no = which.max(likres1)
#     likMax = likres1[no]
#   } else {
#     res1 = apply(tab ,1,function(x, pm, am) 
#       lik1(x[1] , x[2], pm, am), pm, am)
#     likres1 = unlist(lapply(res1, function(x) x$lik))
#     no = which.max(likres1)
#     likMax = likres1[no]
#   }
#   
#   res2 = apply(tab ,1,function(x, pm, am) 
#     lik1(x[1] , x[3], pm, am), pm, am)
#   likres1 = unlist(lapply(res1, function(x) x$lik))
#   likres2 = unlist(lapply(res2, function(x) x$lik))
#   pm = res1$pm
#   vp = res1$vp
#   am = res1$am
#   mp = res1$mp
#   keep1[i+1,1] = res1$from
#   keep1[i+1,2] = res1$to
#   keep2[i+1,1] = res1$lik
#   lr = res1$lik/keep2[i,1]
#   keep2[i+1,2] = lr
#   keep2[i+1,3] = res1$prior
#   keep2[i+1,4] = res1$posterior
#   if(is.na(lr)) lr = LRlimit - 1
#   if(lr > LRlimit) i = i + 1
# }

#singletonList
#' Creates list of singeltons
#' @param n integer
#' @param ids character vector for individuals
#' @param famid character vector for families
#' @param sex integer vector
#' @return list of singletons
#' @export
#' @examples singletonList()
##########
vp = c("V1", "V2")
pm = singletonList(2, ids = vp)
m = marker(pm[[1]], alleles = 1:2, "V1" = 1:2)
pm[[1]] = addMarkers(pm[[1]], m)
m = marker(pm[[2]], alleles = 1:2, "V2" = 1:2)
pm[[2]] = addMarkers(pm[[2]], m)
mp = c("MP1", "MP2")
am = nuclearPed(2, children = mp, father = "R1", mother = "R2")
m = marker(am, "R1" = 1:2, "R2"= 1:2)
am = addMarkers(am, m)
checkInput(pm, am, vp, mp) #lik = 0.5^4
##
.seed(123)
p = c(0.2, 0.8)
nf = 4
vp = paste("V", 1:nf, sep = "")
pm = singletonList(nf, sex = sample(1:2, nf, rep = TRUE))
for (i in 1:nf){
  m = marker(pm[[i]], alleles = 1:2, afreq = p, name ="L1")
  pm[[i]] = addMarkers(pm[[i]], m)
}
amatrix = matrix(sample(1:2, 2*nf, replace = TRUE, prob = p), nrow = nf)
pm = setAlleles(pm, alleles = amatrix)

nt = 2
mp = paste("MP", 1:nt, sep = "")
am = nuclearPed(nt, children = mp, sex = sample(1:2, rep = T, nt))
m = marker(am, alleles = 1:2, afreq = p, name ="L1", 
           "1" = sample(1:2), "2" = sample(1:2) )
am = addMarkers(am, m)
plotPedList(list(pm, am))
generate(pm, am,vp, mp)
pm = pm[[1]]
vp = vp[[1]]
generate(pm, am,vp, mp)
res = loopPair(vp, mp, pm, am)
res = as.matrix(res)
plotPedList(list(pm,am), marker = 1)
from = res[1,1:2]
to = res[1, 3:4]
lik1(from, to, pm, am)
##
pm = singletonList(2)
m = marker(pm[[1]], alleles = 1:2, "V1" = 1:2)
pm[[1]] = addMarkers(pm[[1]], m)
m = marker(pm[[2]], alleles = 1:2, "V2" = 1:2)
pm[[2]] = addMarkers(pm[[2]], m)
mp = c("MP1", "MP2")
am = nuclearPed(2, children = mp, father = "R1", mother = "R2")
m = marker(am, "R1" = 1:2, "R2"= 1:2)
am = addMarkers(am, m)
lik1(vp, mp, pm, am)
                  

# matrix(res[1, c(pe1, pe2, pe3)], ncol=9, byrow = T)
# foo=apply(matrix(res[1,], ncol=9),1, function(x, pe1, pe2,pe3) 
#   rbind(x[pe1], x[pe2], x[pe3]), pe1, pe2, pe3)
# matrix(foo, ncol=9, byrow=F)
rbind(res[1:2,pe1], res[1:2,pe2], res[1:2,pe3])
res[1, pe1]
sex = cbind( sex12[, 1:2], sex34[, 1:2])

setAlleles
am1$FAMID = "am1"
am2 = nuclearPed(2, father =mp[3], mother = mp[4], 
                 children = c("R3", "R4"))
am2$FAMID = "am2"

am = list(am1, am2)
plotPedList(list(pm, am), marker = 1)

### Ex
library(arrangements)
pm = nuclearPed(1, father = "V1", mother = "V2", children = "V3")
vp = labels(pm)
m = marker(pm, alleles = 1:2, "V1" = 1, "V2" =2, "V3" = 1:2)
pm = addMarkers(pm, m)
am = nuclearPed(1, father = "MP1", mother = "MP2", children = "MP3")
m = marker(am, alleles  = 1:2)
am = addMarkers(am, m)
mp = labels(am)

dviForward(pm, am, vp, mp, LRlimit = -1)
foo = loopPair(vp,mp,pm, am)
vp = paste("V",1:20, sep ="");pm = nuclearPed(20, children = vp, sex =rep(1:2,10))
mp = paste("MP",1:20, sep ="");am = nuclearPed(20, children = mp, sex =rep(1:2,10))
foo=loopPair(vp,mp,pm, am)

library(pedtools)
v = list(singleton(1),singleton(2))
getSex(v,"1")
# Trivial example 1
als = 1:2
freq = c(0.01, 0.99)
als = 1:2
v1 = singleton("V1")
m = marker(v1, alleles = als, afreq = freq, "V1" = c(1, 1))
v1 = addMarkers(v1, m)
v2 = singleton("V2")
m = marker(v2, alleles = als, afreq = c(0.1,0.9), "V2" = c(1, 1))
v2 = addMarkers(v2, m)
pm = list(v1, v2)
am = nuclearPed(2, father = "fa", mother = "mo", children = c("MP1", "MP2"))
m = marker(am, alleles = als, afreq = c(0.1,0.9), "fa" = 1:2, "mo" = 1:2)
am = addMarkers(am, m)
plotPedList(list(pm, am), marker = 1)
vp = c("V1", "V2")
mp = c("MP1","MP2")
res = dviForward(pm,am, vp, mp)
res[1]
##########
library(dvi)
data(dvi.nfi)
pm = dvi.nfi$pm
am = dvi.nfi$am
plotPedList(list(pm,am), newdev = TRUE, marker = "TPOX",
            skip.empty.genotypes = TRUE, dev.width = 7, dev.height = 3,
            frametitles = c("Post mortem data (first marker)", 
                            "Ante mortem data (first marker)"))
#### simulation
data(dvi.data)
pm = dvi.data$pm
am = dvi.data$am
vp = dvi.data$vict
mp = dvi.data$miss
rp = c("R1", "R2", "R3")
pmam = c(pm, am)
pmam = setAlleles(pmam, alleles = 0)
dvi.data = profileSim(pmam, N = 1, ids = c(vp, rp), 
                      conditions = 1:13)[[1]]

a = lapply(dvi.data, function(x,vp, mp) dviForward(x[1:5], x[6:7], vp, mp, LRlimit=-1), vp, mp)

a[[1]][[1]]

dviForward(pm, am, vp, mp, LRlimit=-1)


am = dvi.data[6:7]
pm = dvi.data[1:5]
foo=dviForward(pm, am, vp, mp)

plotPedList(list(pm,am), newdev = TRUE, marker = 1:2,
            skip.empty.genotypes = TRUE, dev.width = 7, dev.height = 2,
            frametitles = c("Post mortem data (first marker)", 
                            "Ante mortem data (first marker)"))

res = dviForward(pm, am, vp, mp, LRlimit = 1000)
amx = subset(pm[[5]], "V5" )
amx = setAlleles(amx, alleles=0)
amx$ID = "MP4"
am = list(am[[1]], am[[2]], amx)
mp=c(mp, "MP4")

res = dviForward(pm, am, vp, mp,LRlimit = 10)


pm = dvi.data[1:5]
am = dvi.data[6:7]
res = dviForward(pm, am, vp, mp)

plotPedList(list(pm,am), newdev = TRUE, marker = 1,
            skip.empty.genotypes = TRUE, dev.width = 10, dev.height = 6,
            frametitles = c("Post mortem data", "Ante mortem data"))
am2 = vector("list", 3)
am2[[1]] = am[[1]]
am2[[2]] = am[[2]]
am2[[3]] = pm[[1]]
am2[[3]]$ID = "MP4"
vp = dvi.data$vict
mp = c(dvi.data$miss, "MP4")
pm[[4]] = swapSex(pm[[4]], "V4")
am2[[3]]=setAlleles(am2[[3]], alleles =0)
res = dviForward(pm, am2, vp, mp)


library(arrangements)
rm(list = ls())
source("C:/Users/theg/Dropbox/prj/global/globalR/functions.R")
## General approach

load("C:/Users/theg/Dropbox/prj/global/globalR/forward.RData")
dvi.data = list(pm = pm, vict = vict, am = am, miss = miss )
package.skeleton(name = "dvi", list = c("dvi.data", "dviForward",
       "findsex",    "findsex2", "forward1", "lik1", "loop", "loop2"))       

#windows()
#plotPedList(list(pm,am), newdev = F, marker = 1)
res = dviForward(pm, am, vp, mp)
res = dviForward(dvi.data$pm, dvi.data$am, dvi.data$vict, dvi.data$miss)
##
rm(list = ls())

load("C:/Users/theg/Dropbox/prj/global/globalR/klaas3.RData")
v = paste("V", 1:5, sep="")
res = dviForward(pm, am, vp = v, mp = v)
### Section 1.1 One family
source("functions.R")
load("klaas3.RData")
dviForward(pm=pm, am=am, vp = paste("V1", 1:5, sep=""),
           mp = paste("MP", 1:12, sep=""))
# data pm, am and klaas2 (not used loaded)
# Some noteson the generations of data are in klaas-example.R
v = paste("V", 1:5, sep="")
# Read id.labels
foo = read.table("C:/Users/theg/Dropbox/prj/global/globalR/relabel2.txt", 
      header = F)
windows()
plot(am, marker =1, skip.empty.genotypes = TRUE, 
     deceased = c("NN1", "NN2", "NN3"), cex  = 0.8,
     starred = v, id.labels =as.character(foo[,1]))
plotPedList(pm, marker = 1,newdev = T, 
            skip.empty.genotypes = T, frames =F)
# Perform simulations for correct victimas and unrelated victims
nsim = 1
set.seed(123)
corr = simlik2(pm,am, correct = TRUE, nsim = nsim)
false = simlik2(pm,am, correct = FALSE, nsim = nsim)
# Lokka at simulations results
e1 = false$check
n1 = unlist(lapply(e1, length))
n1 =(1:100)[n1>0]
corr$tab[,sort(c(3*n1-2,3*n1-1))]

save(corr,false, file ="corr-false.RData")
# Example used in global.docx
ex1 = forward(pm,am)

### Section 1.2 Several families
### Generates table combinations in themanuscrip
M = 10; V = 10
nc = matrix(ncol = M, nrow = V)
for (i in 1:V)
  for (j in 1:M){
    gen = generate(M=i,V=j, onlyCount = TRUE, ret = FALSE)
    nc[i,j] = gen$No
  }
dimnames(nc) = list(1:V, 1:M)
###
d = generate(M = 3, V = 5)$d
id =  c("V4", "V5", "MP1", "MP2", "R1", "R2", "R3", "MP3")
load("dvi.null.RData") 
dvi.null = setAlleles(dvi.null, alleles = 0)
nsim = 1
post = simlik(dvi.null, N = nsim, id = id, d, seed = 123)
# post$post now contains the posterior
### Remove symmetries
d = generate(M = 3, V = 5, permute = FALSE)$d
d = rbind(d[,1:3],d[,c(1,3,2)], d[,3:1])
d = unique(d)
mp = paste("MP", 1:3, sep="")

## Calculations on the data loaded
lik = apply(d,1,onelik, dvi.null, mp)
post = lik/sum(lik)
foo = data.frame(d,lik =lik, post =post)

# simulations
nsim = 5
post = simlik(dvi.null, N = nsim, id = id, d, seed = 100)
index = apply(post$post,1, function(x) max(x)>0)
e = data.frame(d[index,], post = post$post[index,])
write.table(e[,1:8], file="e.txt", quote = FALSE)
LRs = apply(e[,-c(1:3)],2,  function(x) x[1]/max(x[-1]))
### Export to Familias

#Mutations
load("dvi.null.mut.RData")
post.mut = simlik(dvi.null.mut, N = nsim, id = id, d, seed = 123)
p1 = post$post
p2 = post.mut$post
###
load("dvi.null.RData") 
plotPedList(dvi, frames = list(c(1:5),6:7), skip.empty.genotypes = TRUE,
            dev.width = 10, dev.height = 3,
            newdev  = TRUE, cex = 1.4, marker = c(1:2),
            frametitles = c("AM data (Two first markers)","PM data"))
g = getAlleles(dvi.null, ids = c("V1", "V2","V3"))
dvi = transferMarkers(dvi.null, dvi.null, ids =id)
mp = c("MP1", "MP2", "MP3")
rownames(g) = mp
dvi = setAlleles(dvi,  ids = mp, alleles = g)

dvi[[6]]= relabel(dvi[[6]], c("FA1", "MO1", "MP1=V1","MP2=V2","R3"))
dvi[[7]] = relabel(dvi[[7]], c("R2", "R3", "MP3=V3"))
plotPedList(dvi, frames = list(c(1:5),c(6:7)), skip.empty.genotypes = TRUE,
            dev.width = 10, dev.height = 3,
            newdev  = TRUE, cex = 1.4, marker = c(1:2),
            frametitles = c("AM data (Two first markers)","PM data"))
### mut data
source("dvi.null.mut.R")
mut = Familias2ped(pedigrees, datamatrix, loci)
V1 = mut[[1]][[1]]
V2 = mut[[1]][[2]]
V3 = mut[[1]][[3]]
V4 = mut[[1]][[4]]
V5 = mut[[1]][[5]]
fam1 = mut[[1]][[6]]
fam2 = mut[[1]][[7]]
dvi.null.mut = list(V1, V2, V3, V4, V5, fam1, fam2)
save(dvi.null.mut, file = "dvi.null.mut.RData")
### OBSOLETE 

load("dvi.null.RData") 
g = getAlleles(dvi.null)
dvi = setAlleles(dvi.null, alleles = 0)
id =  c("V4", "V5", "MP1", "MP2", "R1", "R2","R3", "MP3")
dvi.true = profileSim(dvi, N=1, ids = id, conditions = 1:13,
                      seed = 123)[[1]]
mp = c("MP1", "MP2","MP3")
g = getAlleles(dvi.true,ids = mp)
v = c("V1","V2","V3")
rownames(g) = v
dvi = setAlleles(dvi.true, ids = v, alleles = g)
dvi.null = setAlleles(dvi,ids = mp, alleles = 0)
#save(dvi.true,file="dvi.true.RData")
#save(dvi.null,file="dvi.null.RData")
#########################

lapply(cc[1], onelik, dvi.null, mp)
d= lapply(cc, function(x, dvi.null,mp){
  lapply(x[[1]], onelik, dvi.null, mp)}, dvi.null, mp)
g = getAlleles(dvi.null, id)




permlik = function(id, dvi.null, mp){
  if(all(id==0)){
    lik = LR(dvi.null,1)$likelihoodsPerSystem
    liks = prod(apply(lik, 1, prod))
  } else{
    g = getAlleles(dvi.null, id)
    x = setAlleles(dvi.null, ids = id, alleles = 0)
    mps = unique(permutations(id))
    d1 = dim(mps)[1]
    liks = rep(NA, d1)
    M = length(mp)
    for (i in 1:d1){
      mp0 = mps[i,]
      to = mp[(1:M)[mp0!=0]]
      to = to[order(mp0[mp0!=0])]
      rownames(g) = to
      x2 = setAlleles(x, ids = to, alleles = g)
      lik = LR(x2,1)$likelihoodsPerSystem
      liks[i] = prod(apply(lik, 1, prod))
    }  
  }
  liks
}


b = apply(a,1,  permutations,layout="list")
b = lapply(b, unique)
b = lapply(b, function(x,mp) {lapply(x, function(z,mp) {names(z) = mp; z},mp)}, mp)
sum(unlist(lapply(b,length))) #136
onelik(b[[1]][[3]], dvi.null)

f = function(x, dvi.null) 
  lapply(x, function(z,dvi.null) onelik(z, dvi.null),dvi.null)
f(b[[1]],dvi.null)

lapply(b,f,dvi.null)

id = c(0,0,0)
id = c("0", "0", "0")
id = a[1,]
permlik(id,dvi.null,mp)
apply(a,1,permlik,dvi.null,mp)


####################
mp = paste("MP", 1:3, sep = "") #Missing persons
M = length(mp)
rp = paste("R", 1:3, sep="") # Reference persons
vp = paste("V", 1:5, sep="") # Victim persons
V = length(vp)
N = min(M, V)

#### Make a list of all combinations
vp0 = c(0, paste("V",1:5, sep=""))
a = expand.grid(vp0, vp0, vp0)
a = unique(a)

out = NULL
for(i in 1:2)
  for(j in (i+1):3){
    i1 = (1:dim(a)[1])[(a[,i] == a[,j]) &a[,i]!=0]
    out = c(out,i1)
  }
a = a[-out,]
###
nc = dim(a)[1]
rownames(a) = 1:nc
colnames(a) = mp
liks = rep(NA,nc)
lik = LR(dvi.null,1)$likelihoodsPerSystem
liks[1] = prod(apply(lik, 1, prod))
j=1
for ( i in 2:nc){
  vs = vp[a[i,][a[i,]!=0]]
  g = getAlleles(dvi.null, vs)
  ms = mp[(1:M)[a[i,]>0]]
  rownames(g) = ms[order(a[i,][a[i,]>0])]
  x = setAlleles(dvi.null, ids = ms, alleles = g)
  x = setAlleles(x, ids = vs, alleles = 0)
  j = j+1
  lik = LR(x,1)$likelihoodsPerSystem
  liks[j] = prod(apply(lik, 1, prod))
}





liks = rep(NA, nc)
### Number of combinations checked



lik = LR(dvi.null,1)$likelihoodsPerSystem
liks[1] = prod(apply(lik, 1, prod))
j = 1
foo = rbind(c(0,0,0))
for (k in 1:N){
  idv = combs(1:V, k)
  d1 = dim(idv)[1]
  idm = combs(1:M, k)
  d2 = dim(idm)[1]*length(permn(k))
  j = j+d1*d2
}



  ps = combs(1:k,k)
  rownames(g) = mp
  x = setAlleles(dvi.null, ids = mp, alleles = g)
  x = setAlleles(x, ids = one, alleles = 0)
  j = j+1
  lik = LR(x,1)$likelihoodsPerSystem
  liks[j] = prod(apply(lik, 1, prod))
}
plotPedList(x, frames = list(c(1:5),c(6:7)), skip.empty.genotypes = TRUE,
            dev.width = 10, dev.height = 3,
            newdev  = FALSE, cex = 1.4, marker = 1:2,
            frametitles = c("AM data (2 first markers)","PM data"))
lik = LR(dvi2,1)$likelihoodsPerSystem
lik.true = prod(apply(lik, 1, prod))
g = getAlleles(dvi2, mp)
rownames(g) = vp[1:3]
dvi.null = setAlleles(dvi2, ids = vp[1:3], alleles = g)
dvi.null = setAlleles(dvi.null, ids = mp, alleles = 0)
lik = LR(dvi.null,1)$likelihoodsPerSystem
lik.null = prod(apply(lik, 1, prod))
lik.true/lik.null

s1 = singleton("S1")
m1 = marker(s1,  name = "locus1")
s1 = addMarkers(s1, m1)
s2 = singleton("S2")
m2 = marker(s2,  name ="locus1")
s2 = addMarkers(s2, m2)
peds = list(s1, s2)

peds = setAlleles(peds, ids = "S1", alleles = c(1,2))
g = getAlleles(peds, ids = c("S1"))
setAlleles(peds, ids = "S1", alleles = g)

rm(list = ls())
load("dvi1.RData")
### The following has been loaded:
# id: names of the relevant persons
## MM, i.e., victims V1, V2, V3, V4,V5
## PM reference persons R1, R2, R3 and 
## missing persons MP1,MP2, MP3
## markers names of 13 CODIS markers
## pedlist:se below figure
genotyped = c(paste("V", 1:5, sep= ""), paste("R", 1:3, sep=""))
plotPedList(pedlist, frames = list(c(1:5),c(6:7)), skip.empty.genotypes = TRUE,
            shaded = genotyped, dev.width = 9.5, dev.height = 3,
            newdev  = TRUE, cex = 1.4,
            frametitles = c("AM data","PM data"))

# Simulate data for V4, V5, references and missing persons
genotyped = c("V4", "V5", paste("R", 1:3, sep=""), paste("MP", 1:3, sep = ""))
s = profileSim(pedlist, N = 1, ids = genotyped, cond = markers, seed = 123)
# Find likelihood for the solution; the simulated data
lik = LR(s,1)$likelihoodsPerSystem
lik.true = prod(apply(lik, 1, prod))

# Convert to the situation where no one are identified
# Convert first s to a a pedlist allowing for transfer
save(dvi2, file= "dvi2.RData")
dvi2 = list(V1 = s[[1]][[1]], V2 = s[[1]][[2]],V3 = s[[1]][[3]],V4 = s[[1]][[4]],
                 V5 = s[[1]][[5]],  fam1 = s[[1]][[6]], fam2 = s[[1]][[7]])

# Extract data so that V1 = MP1, V2= MP2 and V3 = MP3
#  V4 and V5 are simulated as  unrelated

V1 = subset(s$fam1[[1]], "MP1")
V1 = relabel(V1, "V1")
V2 = subset(s$fam1[[1]], "MP2")
V2 = relabel(V1, "V2")
V3 = subset(s$fam2[[1]], "MP3")
V3 = relabel(V3, "V3")

s2 = transferMarkers(list(V1, V2,V3), s.pedlist, ids = c("V1","V2", "V3"), 
                     erase = FALSE)
# remove MP1, MP2,MP3
id = c("V1","V2", "V3","V4","V5", "R1",  "R2",  "R3")
s3 = transferMarkers(s2, s2, ids = id, erase = TRUE)

lik0 = LR(s3,1)$likelihoodsPerSystem

###########
rm(list=ls())
#source("example2.R")
#example2 = Familias2ped(pedigrees, datamatrix, loci)
#save(example2, file = "example2.RData")
load("example2.RData")
V1 = example2[[1]][[1]]
m1 = V1$markerdata
navn = lapply(V1$markerdata, function(x) attr(x,"name"))
navn = unlist(navn)
V2 = example2[[1]][[2]]
V3 = example2[[1]][[3]]
V4 = example2[[1]][[4]]
V5 = example2[[1]][[5]]
ped1 = example2[[1]][[6]]
ped2 = example2[[2]][[11]]
ped2 = relabel(ped2, c("R2", "R3", "MP3"))
x = list(V1, V2, V3, V4, V5,  ped1, ped2)
id = c("V1", "V2", "V3", "V4", "V5", "MP1", "MP2", "R1", "R2", "R3", "MP3")
s = profileSim(x, N = 1, ids = id, cond = navn, seed = 123)
lik0 = LR(s,1)$likelihoodsPerSystem

#
fam1 = ped1
fam2 = ped2
pedlist = list(V1, V2, V3, V4, V5,  fam1, fam2)
names(pedlist) = c("V1", "V2", "V3", "V4", "V5",  "fam1", "fam2")
markers = navn
save(pedlist, markers, id, file ="dvi1.RData")
rm(list=ls())
### 

### plot
plotPedList(x[1:5],frames= F, skip.empty.genotypes = TRUE,
            shaded = paste("V",1:5, sep=""), fmar=0.49)
plot(ped1, shaded = "R1", 
     id.labels = c("MP1", "MP2", "R1") )

plot(ped2, shaded = c("R2","R3"), 
     id.labels = c("R2", "R3", "MP3") )



# Example 1
V1 = singleton("V1")
V2 = singleton("V2")
V3 = singleton("V3")
V4 = singleton("V4")
V5 = singleton("V5")
x1 = nuclearPed(3, children=c("MP1", "MP2", "R1"))
x2 = nuclearPed(1, father = "R2",mother = "R3", children= "MP3")


rm(list = ls())
load("glob.RData")
H1 = globMut[[1]][[1]]
m1 = list(H1$markerdata[[1]], H1$markerdata[[2]])
m1 = H1$markerdata
navn = lapply(m1, function(x) attr(x,"name"))
sim1 = profileSim(H1, N = 1, ids =c("S1", "V1", "V2"), 
                  cond = m1, seed = 17)


H2 = globMut[[2]][[7]]
m2 = list(H2$markerdata[[1]], H2$markerdata[[2]])
m2 = H2$markerdata
sim2 = profileSim(H2,N = 1, ids =c("S3", "V3", "V4"), 
                  cond = m2, seed = 177)
for (i in 1:length(navn)){
  attr(sim1[[1]]$markerdata[[i]],"name") = as.character(navn[i])
  attr(sim2[[1]]$markerdata[[i]],"name") = as.character(navn[i])
}

PM = vector("list", 4)
PM[[1]] = subset(sim1[[1]], "V1")
PM[[2]] = subset(sim1[[1]], "V2")
PM[[3]] = subset(sim2[[1]], "V3")
PM[[4]] = subset(sim2[[1]], "V4")

lik1  = rep(NA,6)
k = 0
for (i in 1:3){
  for(j in (i+1):4){
    k = k + 1
    V1 = PM[[i]]
    V1 = relabel(V1, "V1")
    V2 = PM[[j]]
    V2 = relabel(V2, "V2")
    s = setdiff(1:4, c(i,j))
    V3 = PM[[s[1]]]
    V3 = relabel(V3, "V3")
    V4 = PM[[s[2]]]
    V4 = relabel(V4, "V4")
    fra = list(V1, V2, V3, V4)
    til = list(sim1[[1]],sim2[[1]])
    alt = transferMarkers(fra, til, erase = FALSE, matchNames = TRUE)
    foo = apply(LR(alt,1)$likelihoodsPerSystem,1,prod)
    lik1[k] = prod(foo)
  }
}


###

library(donlib)
loadPedSuite( github = TRUE)
# 
library(forrel,quietly = TRUE)
fa = singleton("F")
m = marker(fa, name = "L1")
fa = setMarkers(fa, m)
ch = singleton("C")
m = marker(ch,  name = "L1")
ch = setMarkers(ch, m)
nsim = 5
set.seed(123)
sim = profileSim(list(fa, ch), N = nsim, ids = c("F", "C"), cond = "L1")
simN = vector("list", nsim)
for (i in 1:nsim)
  simN[[i]] = lapply(sim, function(z) z[[i]])
LR(simN, 1)$likelihoodsPerSystem
x = nuclearPed(1, father = "F", children = "C")
m = marker(x,  name = "L1")
x = setMarkers(x,m)
for (i in 1:nsim)
  simN[[i]] = transferMarkers(simN[[i]], x, id = c("F","C"))
LR(simN, 1)$likelihoodsPerSystem





library(forrel)

t1  = nuclearPed(father = "TF", children = "CH")
m = marker(t1, alleles = 1:2, name = "locus1", CH = 1)
t1 = setMarkers(t1, m)

af = singleton("AF")
m = marker(af, alleles = 1:2, name = "locus1")
af = setMarkers(af, m)
mo = singleton("MO", 2)
m = marker(mo, alleles = 1:2, name = "locus1")
mo = setMarkers(mo, m)
ped1 = list(t1, af, mo)
sim = profileSim(ped1, N = 2, ids = c("AF", "MO"), cond = 1)

fra = lapply(sim, function(x) x[[1]])
claim = nuclearPed(father = "AF", mother = "MO", children = "CH")
s = transferMarkers(fra, claim, ids = c("AF", "CH"))
