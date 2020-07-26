#==============================================================================================================
#                       Modelo espacio-temporal capturas incidental lobo marino
#==============================================================================================================
rm(list = ls())

# Packages
library(ggplot2)
library(INLA)
library(reshape2)
library(RANN)
library(RandomFields)
library(tidyverse)
library(gridExtra)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(abind)
library(maps)
library(maptools)
library(mapdata)
library(spdep)
library(rgdal)
library(rgeos)
library(ggmap)
library(dplyr)
library(ggpubr)

setwd("....")

dat_fit = read.csv("dat_fit3.csv", header=T)
dim(dat_fit)
head(dat_fit)
glimpse(dat_fit)

options(scipen=999)

# VARIABLES RESPUESTAS
dat_fit$N_EJEMPLAR = as.numeric(dat_fit$N_EJEMPLAR)

dat_fit$ANO = as.factor(dat_fit$ANO)
dat_fit$COD_PESQUERIA = as.factor(dat_fit$COD_PESQUERIA)
dat_fit$TIPO_DE_RED = as.factor(dat_fit$TIPO_DE_RED)
dat_fit$CLASE_LANCE = as.factor(dat_fit$CLASE_LANCE)
dat_fit$MODELO_EXCLUSION = as.factor(dat_fit$MODELO_EXCLUSION)
dat_fit$TRIM = as.factor(dat_fit$TRIM)
dat_fit$DIA = as.numeric(dat_fit$DIA)
dat_fit$distancia_Km = as.numeric(dat_fit$distancia_Km)

# leer archivo shape
new_sp <- readOGR(dsn = "C:/Users/Usuario/Desktop/Projects/IFOP_lobo_2020/mapas/temp", layer = "new_sp2")
shape5 <- readOGR(dsn = "C:/Users/Usuario/Desktop/Projects/IFOP_lobo_2020/mapas/temp", layer = "shape5")
plot(new_sp)

# Plot
points(cbind(dat_fit$LON2, dat_fit$LAT2), col = "red")



#==============================================================================================================
#                             Modelo espacio-temporal zero-inflated Poisson
#                               usando el metodo SPDE para aproximar el GF 
#==============================================================================================================

# MESH
coords <- as.matrix(dat_fit[,21:20])

mesh_poly = inla.mesh.2d(boundary = shape5, max.edge = c(0.5, 2), cutoff = c(0.1), min.angle = c(0.1, 20), offset = 0.5)
plot(mesh_poly, main = '')
mesh_poly$n

#=====================================================================================
# THE NEXT LINES OF CODES PRODUCES A VERY FINE MESH (More exppensive computation)

mesh = inla.mesh.2d(loc=mesh_poly$loc, max.edge = c(0.5, 2)*2, offset=0.5)
plot(mesh, main="")
mesh$n

source('functions-barriers-dt-models-march2017.R')
source('functions-barrier-DT.R')

mesh = dt.mesh.addon.posTri(mesh)


points = SpatialPoints(mesh$posTri)
points@proj4string = shape5@proj4string
barrier3 = over(shape5, points, returnList=T)
barrier3 = unlist(barrier3)
Omega3 = dt.Omega(list(barrier3, 1:mesh$t), mesh)
Omega.SP3 = dt.polygon.omega(mesh, Omega3)


## Visually check correctness
plot(mesh, main="")
plot(Omega.SP3[[1]], add=T, col='darkgrey')
plot(Omega.SP3[[2]], add=T, col='lightblue')
plot(mesh, add=T)
points(coords, col='red', cex = 0.5)

#=====================================================
# Metodo SPDE con hyperparametros para sigma y rango
#=====================================================
spde = inla.spde2.pcmatern(mesh, prior.range = c(500, .5), prior.sigma = c(2, 0.01))


#Matrz de proyeccion A
A = inla.spde.make.A(mesh=mesh, loc=coords)
dim(A)

# Definimos las variables respuestas
z = (dat_fit$N_EJEMPLAR>0) + 0
y = ifelse(dat_fit$N_EJEMPLAR>0, dat_fit$N_EJEMPLAR, NA)




#=========================================================
#     Modelo seleccionado zeroinflatedbinomialn0
#=========================================================

# Hacemos análisis respecto a este nuevo modelo e introducciendo las variables:
# i)    Trimestre
# ii)   día (como una función)
# iii)  Distancia lobera

# Controles en INLA
cres = list(return.marginals.predictor = TRUE, return.marginals.random = TRUE)
cinla <- list(strategy = 'adaptive', int.strategy = 'eb')  # Estrategias de estimación

# Valores iniciales de los hyperparámetros
ini.zb <- c(1.834,-6.141,-1.079,-0.336)


#=======================
# STACK for estimation
#=======================
# Stack data for z
stk.z <- inla.stack(tag='est.z', data=list(z=z, y=cbind(z, NA)), A=list(A, 1),
                    effects=list(list(i.z =1:spde$n.spde),
                                 list(z.b0=rep(1,length(z)),
                                      ano   = dat_fit$ANO,
                                      dia   = dat_fit$DIA,
                                      trim  = dat_fit$TRIM,
                                      pesq  = dat_fit$COD_PESQUERIA,
                                      red   = dat_fit$TIPO_DE_RED,
                                      lance = dat_fit$CLASE_LANCE,
                                      espe  = dat_fit$ESPECIE_OBJETIVO_LANCE,
                                      dist  = dat_fit$distancia_Km,
                                      mode  = dat_fit$MODELO_EXCLUSION)))

# Stack data for y
stk.y <- inla.stack(tag='est.y', data=list(r=y, y=cbind(NA, y)), A=list(A, 1),
                    effects=list(list(i.y =1:spde$n.spde),
                                 list(y.b0=rep(1,length(y)),
                                      ano   = dat_fit$ANO,
                                      dia   = dat_fit$DIA,
                                      trim  = dat_fit$TRIM,
                                      pesq  = dat_fit$COD_PESQUERIA,
                                      red   = dat_fit$TIPO_DE_RED,
                                      lance = dat_fit$CLASE_LANCE,
                                      espe  = dat_fit$ESPECIE_OBJETIVO_LANCE,
                                      dist  = dat_fit$distancia_Km,
                                      mode  = dat_fit$MODELO_EXCLUSION)))

#=======================
# STACK for prediction
#=======================

# Stack data predict for z
stk.zp <- inla.stack(tag='pred.z', data=list(z=NA, y=cbind(z, NA)),
                     A=list(A, 1),
                     effects=list(list(i.z =1:spde$n.spde),
                                  list(z.b0=rep(1,length(z)),
                                       ano   = dat_fit$ANO,
                                       dia   = dat_fit$DIA,
                                       trim  = dat_fit$TRIM,
                                       pesq  = dat_fit$COD_PESQUERIA,
                                       red   = dat_fit$TIPO_DE_RED,
                                       lance = dat_fit$CLASE_LANCE,
                                       espe  = dat_fit$ESPECIE_OBJETIVO_LANCE,
                                       dist  = dat_fit$distancia_Km,
                                       mode  = dat_fit$MODELO_EXCLUSION)))

# Stack data predict for y
stk.yp <- inla.stack(tag='pred.y', data=list(r=NA, y=cbind(NA, y)),
                     A=list(A, 1),
                     effects=list(list(i.y =1:spde$n.spde),
                                  list(y.b0=rep(1,length(y)),
                                       ano   = dat_fit$ANO,
                                       dia   = dat_fit$DIA,
                                       trim  = dat_fit$TRIM,
                                       pesq  = dat_fit$COD_PESQUERIA,
                                       red   = dat_fit$TIPO_DE_RED,
                                       lance = dat_fit$CLASE_LANCE,
                                       espe  = dat_fit$ESPECIE_OBJETIVO_LANCE,
                                       dist  = dat_fit$distancia_Km,
                                       mode  = dat_fit$MODELO_EXCLUSION)))

stk.zy <- inla.stack(stk.z, stk.y, stk.zp, stk.yp)



#======================================
#    Formula modelo elegido con GLM
#======================================
formula = y ~ 0  + z.b0 + y.b0 + ano + pesq + red + mode +
                                 f(i.z, model=spde) + f(i.y, copy='i.z')

m0  = inla(formula, 
           family=c('binomial', 'zeroinflatednbinomial0'),
           data=inla.stack.data(stk.zy), 
           control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
           control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE, link = 1),
           control.results = cres, control.inla = cinla,
           control.mode = list(theta = ini.zb, restart = TRUE),
           verbose=TRUE)


#======================
# Formula base  + trim
#======================
formula1 = y ~ 0  + z.b0 + y.b0 + ano + pesq + red + mode + trim +
                                  f(i.z, model=spde) + f(i.y, copy='i.z')
m1 = inla(formula1, 
          family=c('binomial', 'zeroinflatednbinomial0'),
          data=inla.stack.data(stk.zy), 
          control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
          control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE, link = 1),
          control.results = cres, control.inla = cinla,
          control.mode = list(theta = ini.zb, restart = TRUE),
          verbose=TRUE)




#========================
# Formula base + día ('rw1')
#========================
formula2 = y ~ 0  + z.b0 + y.b0 + ano + pesq + red + mode +
  f(dia, model='rw1') +                      
  f(i.z, model=spde) + f(i.y, copy='i.z')

m2 = inla(formula2, 
          family=c('binomial', 'zeroinflatednbinomial0'),
          data=inla.stack.data(stk.zy), 
          control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
          control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE, link = 1),
          control.mode = list(theta = ini.zb, restart = TRUE),
          control.results = cres, control.inla = cinla,
          verbose=TRUE)


#============================
# Formula base + especie objetivo
#============================
formula3 = y ~ 0  + z.b0 + y.b0 + ano + pesq + red + mode + espe + 
  f(i.z, model=spde) + f(i.y, copy='i.z')

m3 = inla(formula3, 
          family=c('binomial', 'zeroinflatednbinomial0'),
          data=inla.stack.data(stk.zy), 
          control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
          control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE, link = 1),
          control.mode = list(theta = ini.zb, restart = TRUE),
          control.results = cres, control.inla = cinla,
          verbose=TRUE)

#============================
# Formula base + distancia lobera
#============================
formula4 = y ~ 0  + z.b0 + y.b0  + ano +  pesq + red + mode + dist + 
  f(i.z, model=spde) + f(i.y, copy='i.z')

m4 = inla(formula4, 
          family=c('binomial', 'zeroinflatednbinomial0'),
          data=inla.stack.data(stk.zy), 
          control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
          control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE, link = 1),
          control.mode = list(theta = ini.zb, restart = TRUE),
          control.results = cres, control.inla = cinla,
          verbose=TRUE)


m0$waic$waic
m1$waic$waic
m2$waic$waic
m3$waic$waic
m4$waic$waic



#=============================================
#       Vizualización de los resultados
#=============================================

# CPO para comparar modelos
slcpo <- function(m, na.rm = TRUE) {
  - sum(log(m$cpo$cpo), na.rm = na.rm)
}


c("Model_base" = slcpo(m0), 
  "Model_base + trim" = slcpo(m1), 
  "Model_base + dia" = slcpo(m2), 
  "Model_base + espe" = slcpo(m3), 
  "Model_base + dist" = slcpo(m4),
  "Model_base + dist" = slcpo(m5))

#===================================
m1$waic$waic # modelo base + trim

# Resumen mejor modelo
summary(m1)

# Desviacion estandar de sigma y rango
round(m1$summary.hyperpar[-1, ], 3)
round(m1$summary.random[-1, ], 3)



# Plot interceptos
z.b0 =  m1$marginals.fixed[[1]]
plot_z.b0 = ggplot(data.frame(inla.smarginal(z.b0)), aes(x, y)) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Distribución posterior para los interceptos", y = "z.b0", x = "") + 
  theme_bw()

y.b0 =  m1$marginals.fixed[[2]]
plot_y.b0 = ggplot(data.frame(inla.smarginal(y.b0)), aes(x, y)) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  labs(title = "", y = "y.b0", x = "") + 
  theme_bw() 
grid.arrange(plot_z.b0, plot_y.b0, ncol = 2)


# Estimacion de variabza marginal
marg.variance = inla.tmarginal(function(x) 1/x, m1$marginals.hyperpar$`Precision for ano`)

# Visualización de resultados
index <- inla.stack.index(stk.zy, tag = "pred.z")$data

pred_mean = exp(m1$summary.fitted.values[index, "mean"])
pred_ll = exp(m1$summary.fitted.values[index, "0.025quant"])
pred_ul = exp(m1$summary.fitted.values[index, "0.975quant"])

#coords = cbind(dat_fit$LON2, dat_fit$LAT2)

dpm <- rbind(data.frame(Longitud = coords[, 1], Latitud = coords[, 2], value = pred_mean, variable = "media_pred"),
             data.frame(Longitud = coords[, 1], Latitud = coords[, 2], value = pred_ll, variable = "sup_media_pred"),
             data.frame(Longitud = coords[, 1], Latitud = coords[, 2], value = pred_ul, variable = "inf_media_pred"))
dpm$variable <- as.factor(dpm$variable)

dpm$value

x11()
ggplot(dpm) + geom_tile(aes(Longitud, Latitud, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "N° Capturas incidentales",
    low = "blue", high = "orange"
  ) +
  theme_bw()

summary(pred_mean)
summary(pred_ll)
summary(pred_ul)

# Media posterior para varianza y el rango
r.f <- inla.spde2.result(m1, 'i.z', spde, do.transf=TRUE)   # "i.z" es el indice espacial declarado en data.stack

# Media marginal para la varianza
inla.emarginal(function(x) x, r.f$marginals.log.variance.nominal[[1]])

# Media marginal para el rango
inla.emarginal(function(x) x, r.f$marginals.range.nominal[[1]])


# Grafico para las marginales de los parametros (PRIMERA FORMA)
par(mfrow=c(1,2), mar=c(6,4,2,2), oma=c(1, 1, 1, 1))
plot.default(r.f$marginals.variance.nominal[[1]], type='l',xlab=expression(sigma[i.z]^2), ylab='Density', cex.lab=1, cex.axis=1)
plot.default(r.f$marginals.range.nominal[[1]], type='l',xlab='Practical range', ylab='Density', cex.lab=1, cex.axis=1)


library(brinla)
bri.hyperpar.summary(m1)
#bri.hyperpar.summary(inla.hyperpar(m1)) # Costoso computacionalmente!

# Para los efectos fijos (incluidos el intercepto)
bri.fixed.plot(m1)



## Testeo de la significancia del efecto espacial (solo corremos el efecto espacial dentro de la formula)
formula_spatial = y ~ 0  + z.b0 + y.b0 +  f(i.z, model=spde) + f(i.y, copy='i.z')

r0.sead <- inla(formula_spatial,
                family=c('binomial', 'zeroinflatedpoisson0'),
                data=inla.stack.data(stk.zy), 
                control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
                control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE),
                control.results = cres, control.inla = cinla,
                verbose=TRUE)

# Comparamos ambos DIC (modelos 6 y modelo solo con efecto espacial)
c(r0.sead$waic$waic, m1$waic$waic)

#================================
#   Proyecion del rando field 
#================================
library(dplyr)
library(raster)
lin <- as(new_sp, "SpatialLinesDataFrame")  
pts <- as.data.frame(as(lin, "SpatialPointsDataFrame"))
border = cbind(pts$coords.x1, pts$coords.x2)
min(border)

border.ll <- SpatialPolygons(list(Polygons(list(Polygon(border)), '0')), 
                             proj4string = CRS("+proj=longlat +datum=WGS84"))
border <- spTransform(border.ll,  
                      CRS("+proj=longlat +units=km +zone=18 +south"))
bbox(border)

wh <- apply(bbox(border), 1, diff)
nxy <- round(500 * wh / wh[1])
projgrid <- inla.mesh.projector(mesh, xlim = bbox(border)[1, ],
                                ylim = bbox(border)[2, ],
                                dims = c(1000, 1000))

# Media y desviación estandar del GMRF
xmean    <- inla.mesh.project(projgrid, m1$summary.random$i.y$mean)  # media GMRF
xsd      <- inla.mesh.project(projgrid, m1$summary.random$i.y$sd)    # desviación estandar GMRF


x11()
par(mfrow=c(1,2), mar=c(4,4,3,5))
image.plot(x=projgrid$x, y=projgrid$y, z=xmean, asp=1,xlab='Longitud', ylab='Latitud', axes=T, cex.axis=0.9, axis.args = list(cex.axis = 0.9))
plot(shape5, add=T)
title(main="Media posterior del GMRF")
image.plot(x=projgrid$x, y=projgrid$y, z=xsd, asp=1,xlab='Longitud', ylab='', axes=T, cex.axis=0.9, axis.args = list(cex.axis = 0.9))
plot(shape5)
title(main="Desviación estandar posterior del GMRF")
dev.off()

