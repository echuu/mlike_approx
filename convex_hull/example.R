
install.packages("uniformly")
install.packages("rgl")
library('uniformly')
library('rgl')

# vertices of a 3d icosahedron
vs = t(rgl::icosahedron3d()$vb[1:3,])
vs %>% head

install.packages("geometry") 
# calculates a triangulation of convex hull of a set of points
library('geometry')

# convex hull of the points such that no N-spere defined by the N-triangles
# contains any other points from the set
tetrahedra = delaunayn(vs, options = "Qz") # (12 x 3)
head(tetrahedra)

# compute the volumes of eachof the tetrahedra with volume_tetrahedron
volumes <- 
    apply(tetrahedra, 1, 
          function(t){
              volume_tetrahedron(vs[t[1],], vs[t[2],], vs[t[3],], vs[t[4],])
          })

probs = volumes / sum(volumes)

# uniformly sample a point in the icosahedron
#   (0) select a tetrahedron at random, with probability given by the norm vol
#   (1) uniformly sample a point in the picked tetrahedron

i = sample.int(nrow(tetrahedra), 1, prob = probs) # pick the tetra
th = tetrahedra[i,]
# sample from the tetrahedron
runif_in_tetrahedron(1, vs[th[1],], vs[th[2],], vs[th[3],], vs[th[4],])


# repeat for 100 random points
nsims <- 100
sims <- matrix(NA_real_, nrow=nsims, ncol=3)
for(k in 1:nsims){
    th <- tetrahedra[sample.int(nrow(tetrahedra), 1, prob=probs),]
    sims[k,] <- runif_in_tetrahedron(1, vs[th[1],], vs[th[2],], vs[th[3],], vs[th[4],])
}


install.packages("rgl")
library('rgl')

open3d(windowRect=c(100,100,600,600))
#> wgl 
#>   1
shade3d(icosahedron3d(), color="red", alpha=0.3)
points3d(sims)
rglwidget()


