#===========================================================
#                   Modelos preliminares
#===========================================================


## Histograma de la variable respuesta
library(ggplot2)
library(pscl)
library(boot)
library(MASS)
library(dplyr)
library(mgcv)
library(mgcViz)


# Determinamos las variables de importancia en función del menor AIC 
# (función step)

setwd("C:/Users/Usuario/Desktop/Projects/UV_lobo_2020")

dat_fit = read.csv("dat_fit.csv", header=T)
dim(dat_fit)
head(dat_fit)
glimpse(dat_fit)

options(scipen=999)

# VARIABLES RESPUESTA
dat_fit$loboc = as.numeric(dat_fit$loboc)

dat_fit$buq = as.factor(dat_fit$buq)
dat_fit$ano = as.factor(dat_fit$ano)
dat_fit$trim = as.factor(dat_fit$trim)
dat_fit$cod_pesq = as.factor(dat_fit$cod_pesq)
dat_fit$nom_pesq = as.factor(dat_fit$nom_pesq)
dat_fit$repro = as.factor(dat_fit$repro)
dat_fit$estado_mar = as.factor(dat_fit$estado_mar)

# Nuevo data set que contiene sólo las variables de ineteres
dat_fit = dat_fit[, c(3:12, 14, 17, 20, 21, 67, 70, 71, 73:75)]
dim(dat_fit)

glimpse(dat_fit)
head(dat_fit)

#================================================================================
#                          JUREL INDUSTRIAL CENTRO SUR
#================================================================================


# Vamos a filtrar por JUREL INDUSTRIAL CENTRO SUR (jics)
dat_jics <- dat_fit %>% filter(nom_pesq %in% "JUREL INDUSTRIAL CENTRO SUR")
summary(dat_jics) # cod_pesq tiene 5 niveles pero solo 1 para hacer el contraste
                  # nom_pesq tiene 7 niveles pero solo 1 para hacer el contraste
                  # esp_obj es sólo 1 
                  # mes será reemplazado por trim
                  

dat_jics %>% filter(cod_pesq == "112") # No hay información y es lo mismo para los demas niveles 
                                       # 115, 118, 120
dat_jics = droplevels(dat_jics) 
summary(dat_jics$cod_pesq)



# Sacamos el nombre de la pesquería para poder realizar los constrastes
str(dat_jics)
dat_jics = dat_jics[, c(-1, -2, -3, -6)] # Sacamos del data.frame 'cod_pesq', 'nom_pesq', esp_obj'
head(dat_jics)
table(dat_jics$loboc)
length(dat_jics$num_lan)




full = glm(loboc~., data = dat_jics, family = poisson)
step = stepAIC(full, trace = FALSE)
step$anova


backward  = stepAIC(glm(loboc~.,data = dat_jics, family = poisson),direction="backward")
forward   = stepAIC(glm(loboc~.,data = dat_jics, family = poisson),direction="forward")
both      = stepAIC(glm(loboc~.,data = dat_jics, family = poisson),direction="both")


backward$anova
forward$anova
both$anova

# Comparar devianzas
step$deviance
backward$deviance
forward$deviance
both$deviance   # Tienen igual devianza así que podemos ocupar cualquier

# Menos AIC es forward
formula(step)
formula(backward)
formula(forward)
formula(both)

#===============================================================================================================
#                                         MODELOS PRELIMINARES
#===============================================================================================================


p = ggplot(dat_jics, aes(x=loboc))+ geom_histogram(color="black", fill="grey")
  p  + theme_bw()

p + scale_color_brewer(palette="black") +
  theme_classic()+theme(legend.position="top") +
  labs(title = "Jurel Industrial Centro Sur", y = "Frecuencia", x = "N° ejemplares capturados") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=11, face = "bold"),
        axis.title=element_text(size=11,face="bold")) +
  theme(plot.title = element_text(size = 14, face = "bold"))



fr <- table(dat_jics$loboc) %>% data.frame
names(fr) <- c('Capturas', 'freq')
fr$Capturas <- as.numeric(as.character(fr$Capturas)) #convert factor to numeric
ggplot(fr, aes(x = Capturas, y = freq)) +
  labs(title = "Jurel Industrial Centro Sur", x = "Capturas", y = "N°") +
  geom_col() +
  theme_bw() +
  lims(y = c(0, 1700)) + 
  geom_line() + 
  geom_text(aes(x = Capturas, y = freq, label = freq, vjust = -1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(angle = 0)) 



## Poisson GLM
# model_1 = glm(loboc ~ num_lan + ano + repro + hora_l + estado_mar + tipo_agreg + 
#                 tsm + captura.to + dist_lob + prof_med + trim + total_aves, 
#                 family = 'poisson', data = dat_jics, na.action = na.omit)
# summary(model_1)
# model_1$aic
# 
# 
# ## Verificiamos si existe baja/sobredispersión en el modelo
# E2 = resid(model_1, type = "pearson")
# N  = nrow(dat_jics)
# p  = length(coef(model_1))   
# sum(E2^2) / (N - p)


## Parece ser que sí existe soberdispersión (varianza no es igual a la media) 
## por tanto probaremos con una regresión binomial negativa

model_2 = glm.nb(loboc ~ ano + repro + dist_lob + prof_med + trim + dis_cost, 
                 data = dat_jics, na.action = na.omit, maxit = 2000)
summary(model_2)



# Comprobamos dispersión nuevamente
E2_2 = resid(model_2, type = "pearson")
N  = nrow(dat_jics)
p  <- length(coef(model_2)) + 1  # '+1' is for variance parameter in NB
sum(E2_2^2) / (N - p)

m2 <- update(model_2, . ~ . - num_lan)  # Hacer esto sólo si le agrego la variable ano
anova(model_2, m2)

# Luce mejor que el modelo de conteo Poisson pero esto igualmente puede generar
# problemas en la inferencia de los parámetros

# Modelo zero-inflated Poisson
model_3 <- zeroinfl(loboc ~ num_lan + ano + repro + dist_lob + prof_med + trim | ## Predictor for the Poisson process
                      num_lan + ano + repro + dist_lob + prof_med + trim,  ## Predictor for the Bernoulli process;
                      dist = 'poisson',
                      data = dat_jics)

summary(model_3) # error en la estimación de mat var-cov
model_3


# Modelo zero-inflated Poisson 2
model_3_2 <- zeroinfl(loboc ~ num_lan + ano + repro + dist_lob + prof_med + trim | ## Predictor for the Poisson process
                        num_lan + ano + repro + dist_lob + prof_med + trim,  ## Predictor for the Bernoulli process;
                        dist = 'poisson',
                        data = dat_jics)

summary(model_3_2) # 
model_3_2$loglik


# Estadístico de dispersión 2
E2_32 <- resid(model_3_2, type = "pearson")
N  <- nrow(dat_jics)
p  <- length(coef(model_3_2))  
sum(E2_32^2) / (N - p)



#=============================================================
#               Zero inflated negative binomial
#=============================================================

# 
model_4 <- zeroinfl(loboc ~ ano + repro + dist_lob + prof_med + trim,
                    dist = 'negbin',
                    data = dat_jics)
summary(model_4)


# 
model_4_2 <- zeroinfl(loboc ~ ano + repro + dist_lob + prof_med + trim,
                      dist = 'negbin',
                      data = dat_jics)
summary(model_4_2)


# Estadístico de dispersión
E2_4 <- resid(model_4, type = "pearson")
N  <- nrow(dat_jics)
p  <- length(coef(model_4)) + 1 # 
test1= sum(E2_4^2) / (N - p)



E2_42 <- resid(model_4_2, type = "pearson")
N  <- nrow(dat_jics)
p  <- length(coef(model_4_2)) + 1 # 
test2= sum(E2_42^2) / (N - p)


# Comparación modelos bajo/sobredispersión
c(test1, test2)

#===================================================================
# Modelado con GAM para incluir latitud y longitud en la estimación
#===================================================================
library(mgcv)
m0 <- gam(loboc ~ ano + repro + dist_lob + prof_med + trim + 
                             s(lon,lat, k=200), 
                             data = dat_jics,
                             family = nb, method = "REML")


m1 <- gam(loboc ~ ano + repro + dist_lob + prof_med + trim + 
                            s(lon,lat, bs="gp",k=200,m=c(5,5)), 
                            data = dat_jics,
                            family = nb, method = "ML")



m2 <- gam(loboc ~ ano + repro + dist_lob + prof_med + trim + 
                            s(lon,lat, bs="ds", m=c(5,.5),k=200), 
                            data = dat_jics,
                            family = nb, method = "ML")


m3 <- gam(loboc ~ ano + repro + dist_lob + prof_med + trim + 
                            s(lon,lat, bs="tp", k=200), 
                            data = dat_jics,
                            family = nb, method = "ML")


AIC(m0, m1, m2, m3)
anova(m0, m1, m2, m3, test="Chisq")
gam.check(m2, k.rep=1000)

summary(m3)


AIC(model_2, m3)
summary(m3)


library(mgcViz)
m3 <- getViz(m3)
print(plot(m3) + labs(title = NULL), pages = 1)

print(plot(m3, select = 2:7), pages = 1)


print(plot(m3), ask = FALSE)
print(plot(m3, allTerms = TRUE), pages = 1)


vis.gam(x = m3, # GAM object
        view = c("lon","lat"), # variables
        plot.type = "persp") # kind of plot

qq(m3, method = "simul1", a.qqpoi = list("shape" = 1), a.ablin = list("linetype" = 2))
qq(m3, rep = 20, showReps = T, CI = "none", a.qqpoi = list("shape" = 19), a.replin = list("alpha" = 0.2))


o <- qq(m3, rep = 10, method = "simul1", CI = "normal", showReps = TRUE,
        ngr = 1e2, a.replin = list(alpha = 0.1), a.qqpoi = list(shape = 19))
o 

gridPrint(o, zoom(o, xlim = c(2, 2.5), ylim = c(2, 2.5)), ncol = 2)

vis.gam(x = m3, # GAM object
        view = c("lon","lat"), # variables
        plot.type = "contour", main='Efecto espacial',
        color = "heat") # kind of plot

plotRGL(sm(m2, 1), fix = c("z" = 0), residuals = TRUE)








#================================================================================
#                      ANCHOVETA ARTESANAL ZONA NORTE
#================================================================================

# Vamos a filtrar por ANCHOVETA ARTESANAL ZONA NORTE
dat_aazn <- dat_fit %>% filter(nom_pesq %in% "ANCHOVETA ARTESANAL ZONA NORTE")
glimpse(dat_aazn) 

#dat_aazn = droplevels(dat_aazn) 

# Sacamos el nombre de la pesquería para poder realizar los constrastes
str(dat_aazn)
dat_aazn = dat_aazn[, c(-1, -2, -3, -6)] # Sacamos del data.frame 'cod_pesq', 'nom_pesq', esp_obj'
head(dat_aazn)
table(dat_aazn$loboc)
length(dat_aazn$num_lan)



p = ggplot(dat_aazn, aes(x=loboc))+ geom_histogram(color="black", fill="grey")
p  + theme_bw()

p + scale_color_brewer(palette="black") +
  theme_classic()+theme(legend.position="top") +
  labs(title = "Anchoveta Artesanal Zona Norte", y = "Frecuencia", x = "N° ejemplares capturados") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=11, face = "bold"),
        axis.title=element_text(size=11,face="bold")) +
  theme(plot.title = element_text(size = 14, face = "bold"))



#==================================
# Modelos lineales generalizados
# Poisson GLM
#==================================
model_1 = glm(loboc ~ num_lan + ano + repro + hora_l + estado_mar + tipo_agreg + 
              tsm + captura.to + dist_lob + prof_med + trim + total_aves, 
              family = 'poisson', data = dat_aazn, na.action = na.omit)
summary(model_1)



## Verificiamos si existe baja/sobredispersión en el modelo
E2 = resid(model_1, type = "pearson")
N  = nrow(dat_aazn)
p  = length(coef(model_1))   
sum(E2^2) / (N - p)

#==================================
# Modelos lineales generalizados
# glm Binomial negativo
#==================================
model_2 = glm.nb(loboc ~ ano + repro + estado_mar + lat + lon + 
                 tsm + captura.to + dist_lob, 
                 data = dat_aazn, na.action = na.omit, maxit = 1000)

# Comprobamos dispersión nuevamente
E2_2 = resid(model_2, type = "pearson")
N  = nrow(dat_aazn)
p  <- length(coef(model_2)) + 1  # '+1' is for variance parameter in NB
sum(E2_2^2) / (N - p)

# Modelo zero-inflated Poisson
model_3 <- zeroinfl(loboc ~ ano + repro + estado_mar + lat + lon + 
                    tsm + captura.to + dist_lob | ## Predictor for the Poisson process
                    ano + repro + estado_mar + lat + lon + 
                    tsm + captura.to + dist_lob,  ## Predictor for the Bernoulli process;
                    dist = 'poisson',
                    data = dat_aazn)

summary(model_3) # error en la estimación de mat var-cov


# Modelo zero-inflated Poisson 2
model_3_2 <- zeroinfl(loboc ~ num_lan + ano + repro  + 
                      tsm + dist_lob | ## Predictor for the Poisson process
                      num_lan + ano + repro  + 
                      tsm + dist_lob,  ## Predictor for the Bernoulli process;
                      dist = 'poisson',
                      data = dat_aazn)

# Estadístico de dispersión 2
E2_32 <- resid(model_3_2, type = "pearson")
N  <- nrow(dat_aazn)
p  <- length(coef(model_3_2))  
sum(E2_32^2) / (N - p)



#=============================================================
#               Zero inflated negative binomial
#=============================================================

# Incluimos sin variable 'estado_mar'
model_4 <- zeroinfl(loboc ~ ano + repro + estado_mar + lat + lon + 
                    tsm + captura.to + dist_lob,
                    dist = 'negbin',
                    data = dat_aazn)
summary(model_4) # Problema de singularidad (probableemnte en estimación de desvest)


# Incluimos estado_mar
model_4_2 <- zeroinfl(loboc ~ ano + repro + estado_mar + tsm + captura.to + dist_lob,
                      dist = 'negbin',
                      data = dat_aazn)
summary(model_4_2)


# Estadístico de dispersión
E2_4 <- resid(model_4, type = "pearson")
N  <- nrow(dat_aazn)
p  <- length(coef(model_4)) + 1 # 
test1= sum(E2_4^2) / (N - p)



E2_42 <- resid(model_4_2, type = "pearson")
N  <- nrow(dat_aazn)
p  <- length(coef(model_4_2)) + 1 # 
test2= sum(E2_42^2) / (N - p)


# Comparación modelos bajo/sobredispersión
c(test1, test2)



#===================================================================
# Modelado con GAM para incluir latitud y longitud en la estimación
#===================================================================
m0 <- gam(loboc ~ ano + repro + estado_mar + tsm + captura.to + dist_lob +
                  s(lon,lat, k=200), 
                  data = dat_aazn,
                  family = nb, method = "ML")


m1 <- gam(loboc ~ ano + repro + estado_mar + tsm + captura.to + dist_lob +
                  s(lon,lat, bs="gp",k=200,m=c(5,5)), 
                  data = dat_aazn,
                  family = nb, method = "ML")



m2 <- gam(loboc ~ ano + repro + estado_mar + tsm + captura.to + dist_lob +
                  s(lon,lat, bs="ds", m=c(5, 0.5),k=200), 
                  data = dat_aazn,
                  family = nb, method = "ML")

m3 <- gam(loboc ~ ano + repro + estado_mar + tsm + captura.to + dist_lob +
                  s(lon,lat, bs="tp", k=200), 
                  data = dat_aazn,
                  family = nb, method = "ML")


AIC(m0, m1, m2, m3)
anova(m0, m1, m2, m3, test="Chisq")
gam.check(m2, k.rep=1000)

AIC(model_2, m2)



summary(m2)
plot(m2)


library(mgcViz)
m2 <- getViz(m2)
print(plot(m2) + labs(title = NULL), pages = 1)

print(plot(m2, select = 2:7), pages = 1)


print(plot(m2), ask = FALSE)
print(plot(m2, allTerms = TRUE), pages = 1)


vis.gam(x = m2, # GAM object
          view = c("lon","lat"), # variables
          plot.type = "persp") # kind of plot

qq(m2, method = "simul1", a.qqpoi = list("shape" = 1), a.ablin = list("linetype" = 2))
qq(m2, rep = 20, showReps = T, CI = "none", a.qqpoi = list("shape" = 19), a.replin = list("alpha" = 0.2))


o <- qq(m2, rep = 10, method = "simul1", CI = "normal", showReps = TRUE,
        ngr = 1e2, a.replin = list(alpha = 0.1), a.qqpoi = list(shape = 19))
o 

gridPrint(o, zoom(o, xlim = c(2, 2.5), ylim = c(2, 2.5)), ncol = 2)

vis.gam(x = m2, # GAM object
        view = c("lon","lat"), # variables
        plot.type = "contour", main='Efecto espacial',
        color = "heat") # kind of plot

plotRGL(sm(m2, 1), fix = c("z" = 0), residuals = TRUE)

# Mejor modelo m2


#================================================================================
#                     ANCHOVETA JUREL ARTESANAL CENTRO NORTE
#================================================================================

# Vamos a filtrar por JUREL INDUSTRIAL CENTRO SUR (jics)
dat_ajacn <- dat_fit %>% filter(nom_pesq %in% "ANCHOVETA JUREL ARTESANAL CENTRO NORTE")
glimpse(dat_ajacn) 

#dat_ajacn = droplevels(dat_ajacn) 

# Sacamos el nombre de la pesquería para poder realizar los constrastes
dat_ajacn = dat_ajacn[, c(-1, -2, -3, -6)] # Sacamos del data.frame 'cod_pesq', 'nom_pesq', esp_obj'
head(dat_ajacn)
table(dat_ajacn$loboc)
length(dat_ajacn$loboc)


p = ggplot(dat_ajacn, aes(x=loboc))+ geom_histogram(color="black", fill="grey")
p  + theme_bw()

p + scale_color_brewer(palette="black") +
  theme_classic()+theme(legend.position="top") +
  labs(title = "Anchoveta Jurel Artesanal Centro Norte", y = "Frecuencia", x = "N° ejemplares capturados") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=11, face = "bold"),
        axis.title=element_text(size=11,face="bold")) +
  theme(plot.title = element_text(size = 14, face = "bold"))

table(dat_ajacn$num_lan)

#==================================
# Modelos lineales generalizados
# Poisson GLM
#==================================
model_1 = glm(loboc ~ num_lan + ano + repro  + lat + lon + 
              prof_med, 
              family = 'poisson', data = dat_ajacn, na.action = na.omit)
summary(model_1)

## Verificiamos si existe baja/sobredispersión en el modelo
E2 = resid(model_1, type = "pearson")
N  = nrow(dat_ajacn)
p  = length(coef(model_1))   
sum(E2^2) / (N - p)

#==================================
# Modelos lineales generalizados
# glm Binomial negativo
#==================================
model_2 = glm.nb(loboc ~ num_lan + ano + repro  + lat + lon, 
                 data = dat_ajacn, na.action = na.omit, maxit = 2000)

summary(model_2)

# Comprobamos dispersión nuevamente
E2_2 = resid(model_2, type = "pearson")
N  = nrow(dat_ajacn)
p  <- length(coef(model_2)) + 1  # '+1' is for variance parameter in NB
sum(E2_2^2) / (N - p)



# Modelo zero-inflated Poisson
model_3 <- zeroinfl(loboc ~ num_lan + ano + repro  + lat + lon| ## Predictor for the Poisson process
                    num_lan + ano + repro  + lat + lon + 
                    prof_med,  ## Predictor for the Bernoulli process;
                    dist = 'poisson',
                    data = dat_ajacn)

# Estadístico de dispersión 2
E2_3 <- resid(model_3, type = "pearson")
N  <- nrow(dat_ajacn)
p  <- length(coef(model_3))  
sum(E2_3^2) / (N - p)





#=============================================================
#               Zero inflated negative binomial
#=============================================================

# Incluimos sin variable 'estado_mar'
model_4 <- zeroinfl(loboc ~ num_lan + repro  + lat + lon,
                    dist = 'negbin',
                    data = dat_ajacn)
summary(model_4) 


# Estadístico de dispersión
E2_4 <- resid(model_4, type = "pearson")
N  <- nrow(dat_ajacn)
p  <- length(coef(model_4)) + 1 # 
sum(E2_4^2) / (N - p)




# Comparación modelos bajo/sobredispersión
c(test1, test2)



#===================================================================
# Modelado con GAM para incluir latitud y longitud en la estimación
#===================================================================
m0 <- gam(loboc ~ num_lan + ano + repro +
          s(lon,lat, k=100), 
          data = dat_ajacn,
          family = nb, method = "ML")


m1 <- gam(loboc ~ num_lan + ano + repro +
            s(lon,lat, bs="gp",k=100,m=c(1,1)), 
          data = dat_ajacn,
          family = nb, method = "ML")



m2 <- gam(loboc ~ num_lan + ano + repro +
          s(lon,lat, bs="ds", m=c(1,.5),k=100), 
          data = dat_ajacn,
          family = nb, method = "ML")


m3 <- gam(loboc ~ num_lan + ano + repro +
          s(lon,lat, bs="tp", k=100), 
          data = dat_ajacn,
          family = nb, method = "ML")

  
AIC(m0, m1, m2, m3)
anova(m0, m1, m2, m3, test="Chisq")
gam.check(m2, k.rep=1000)

# Analizar modelo elegido entre model_2 y m3
summary(m3)

AIC(model_2, m3) # NOs quedamos con m3 !!!!!!!!!!

library(mgcViz)
m3 <- getViz(m3)

print(plot(m3, select = 2:7), pages = 1)


print(plot(m3, allTerms = TRUE), pages = 1)


vis.gam(x = m3, # GAM object
        view = c("lon","lat"), # variables
        plot.type = "persp") # kind of plot

qq(m3, method = "simul1", a.qqpoi = list("shape" = 1), a.ablin = list("linetype" = 2))
qq(m3, rep = 20, showReps = T, CI = "none", a.qqpoi = list("shape" = 19), a.replin = list("alpha" = 0.2))

o <- qq(m3, rep = 10, method = "simul1", CI = "normal", showReps = TRUE,
        ngr = 1e2, a.replin = list(alpha = 0.1), a.qqpoi = list(shape = 19))
o 

gridPrint(o, zoom(o, xlim = c(2, 2.5), ylim = c(2, 2.5)), ncol = 2)

vis.gam(x = m3, # GAM object
        view = c("lon","lat"), # variables
        plot.type = "contour", main='Efecto espacial',
        color = "heat") # kind of plot


#================================================================================
#                      SARDINA-ANCHOVETA ARTESANAL CENTRO SUR
#================================================================================

# Vamos a filtrar por JUREL INDUSTRIAL CENTRO SUR (jics)
dat_saacs <- dat_fit %>% filter(nom_pesq %in% "SARDINA-ANCHOVETA ARTESANAL CENTRO SUR")
glimpse(dat_saacs) 

#dat_saacs = droplevels(dat_saacs) 

# Sacamos el nombre de la pesquería para poder realizar los constrastes
str(dat_saacs)
dat_saacs = dat_saacs[, c(-1, -2, -3, -6)] # Sacamos del data.frame 'cod_pesq', 'nom_pesq', esp_obj'
head(dat_saacs)
table(dat_saacs$loboc)
length(dat_saacs$loboc)

mean(dat_saacs$loboc)
4.052138^0 * exp(-4.052138) / factorial(0)

p = ggplot(dat_saacs, aes(x=loboc))+ geom_histogram(color="black", fill="grey")
p  + theme_bw()

p + scale_color_brewer(palette="black") +
  theme_classic()+theme(legend.position="top") +
  labs(title = "Sardina Anchoveta Artesanal Centro Sur", y = "Frecuencia", x = "N° ejemplares capturados") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=11, face = "bold"),
        axis.title=element_text(size=11,face="bold")) +
  theme(plot.title = element_text(size = 14, face = "bold"))


#==================================
# Modelos lineales generalizados
# Poisson GLM
#==================================
model_1 = glm(loboc ~ num_lan + ano + repro + estado_mar + tipo_agreg + 
              tsm + captura.to + dist_lob + prof_med + trim + total_aves, 
              family = 'poisson', data = dat_saacs, na.action = na.omit)
summary(model_1)

## Verificiamos si existe baja/sobredispersión en el modelo
E2 = resid(model_1, type = "pearson")
N  = nrow(dat_saacs)
p  = length(coef(model_1))   
sum(E2^2) / (N - p)

#==================================
# Modelos lineales generalizados
# glm Binomial negativo
#==================================
model_2 = glm.nb(loboc ~ repro + dis_cost + dist_lob + prof_med, 
                 data = dat_saacs, na.action = na.omit, maxit = 2000)

summary(model_2)

# Comprobamos dispersión nuevamente
E2_2 = resid(model_2, type = "pearson")
N  = nrow(dat_saacs)
p  <- length(coef(model_2)) + 1  # '+1' is for variance parameter in NB
sum(E2_2^2) / (N - p)

m2 <- update(model_2, . ~ . + ano)  # Hacer esto sólo si le agrego la variable ano
anova(model_2, m2)

(est <- cbind(Estimate = coef(model_2), confint(model_2)))
exp(est)

# Modelo zero-inflated Poisson
model_3 <- zeroinfl(loboc ~ repro + dis_cost + dist_lob + prof_med | ## Predictor for the Poisson process
                      repro + dis_cost + dist_lob + prof_med,  ## Predictor for the Bernoulli process;
                      dist = 'poisson',
                      data = dat_saacs)

# Estadístico de dispersión 2
E2_3 <- resid(model_3, type = "pearson")
N  <- nrow(dat_saacs)
p  <- length(coef(model_3))  
sum(E2_3^2) / (N - p)


library(jtools)
summ(model_2, exp = T) #could exp the coefficients but this is

# Para el modelo zero-inflado Poisson
tmp <- summary(model_3) #doing this manually
tmp$coefficients$count[-7, 1] %>% exp()



#=============================================================
#               Zero inflated negative binomial
#=============================================================

# Incluimos sin variable 'estado_mar'
model_4 <- zeroinfl(loboc ~ repro + dis_cost + dist_lob + prof_med,
                    dist = 'negbin',
                    data = dat_saacs)
summary(model_4) 


# Estadístico de dispersión
E2_4 <- resid(model_4, type = "pearson")
N  <- nrow(dat_saacs)
p  <- length(coef(model_4)) + 1 # 
sum(E2_4^2) / (N - p)



# Incluimos sin variable 'estado_mar'
model_4_2 <- zeroinfl(loboc ~ repro + dist_lob + prof_med,
                      dist = 'negbin',
                      data = dat_saacs)
summary(model_4_2) 


# Estadístico de dispersión
E2_4_2 <- resid(model_4_2, type = "pearson")
N  <- nrow(dat_saacs)
p  <- length(coef(model_4_2)) + 1 # 
sum(E2_4_2^2) / (N - p)


AIC(model_4, model_4_2)

# Para el modelo zero-inflado Poisson
tmp <- summary(model_4) #doing this manually
tmp$coefficients$count[-7, 1] %>% exp()

tmp <- summary(model_4_2) #doing this manually
tmp$coefficients$count[-7, 1] %>% exp()

(est <- cbind(Estimate = coef(model_2), confint(model_2)))
exp(est)

AIC(model_2, model_4, model_4_2)


# ELEGIDO MODEL_2


m2 <- gam(loboc ~ repro + dis_cost + dist_lob + prof_med +
            s(lon,lat, bs="ds", m=c(1,.5),k=100), 
            data = dat_saacs,
            family = nb, method = "ML")

AIC(model_2, model_4, m2)

# Best model
summary(m2)


library(mgcViz)
m2 <- getViz(m2)

print(plot(m2, select = 2:7), pages = 1)


print(plot(m2, allTerms = TRUE), pages = 1)


vis.gam(x = m2, # GAM object
        view = c("lon","lat"), # variables
        plot.type = "persp") # kind of plot

qq(m2, method = "simul1", a.qqpoi = list("shape" = 1), a.ablin = list("linetype" = 2))
qq(m2, rep = 20, showReps = T, CI = "none", a.qqpoi = list("shape" = 19), a.replin = list("alpha" = 0.2))


o <- qq(m2, rep = 10, method = "simul1", CI = "normal", showReps = TRUE,
        ngr = 1e2, a.replin = list(alpha = 0.1), a.qqpoi = list(shape = 19))
o 

gridPrint(o, zoom(o, xlim = c(2, 2.5), ylim = c(2, 2.5)), ncol = 2)

vis.gam(x = m2, # GAM object
        view = c("lon","lat"), # variables
        plot.type = "contour", main='Efecto espacial',
        color = "heat") # kind of plot








#================================================================================
#                      SARDINA-ANCHOVETA INDUSTRIAL CENTRO SUR
#================================================================================

# Vamos a filtrar por JUREL INDUSTRIAL CENTRO SUR (jics)
dat_saics <- dat_fit %>% filter(nom_pesq %in% "SARDINA-ANCHOVETA INDUSTRIAL CENTRO SUR")
glimpse(dat_saics) 

#dat_saics = droplevels(dat_saics) 

# Sacamos el nombre de la pesquería para poder realizar los constrastes
str(dat_saics)
dat_saics = dat_saics[, c(-1, -2, -3, -6)] # Sacamos del data.frame 'cod_pesq', 'nom_pesq', esp_obj'
table(dat_saics$loboc)
length(dat_saics$loboc)

mean(dat_saics$loboc)
5.448052^0 * exp(-5.448052) / factorial(0)


p = ggplot(dat_saics, aes(x=loboc))+ geom_histogram(color="black", fill="grey")
p  + theme_bw()

p + scale_color_brewer(palette="black") +
  theme_classic()+theme(legend.position="top") +
  labs(title = "Sardina Anchoveta Industrial Centro Sur", y = "Frecuencia", x = "N° ejemplares capturados") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=11, face = "bold"),
        axis.title=element_text(size=11,face="bold")) +
  theme(plot.title = element_text(size = 14, face = "bold"))




#==================================
# Modelos lineales generalizados
#==================================

#__________________________________
#           Modelo Poisson
#==================================
model_1 = glm(loboc ~ num_lan + ano + repro + estado_mar + tipo_agreg + dis_cost + 
              tsm + captura.to + dist_lob + prof_med + trim + total_aves, 
              family = 'poisson', data = dat_saics, na.action = na.omit)

## Verificiamos si existe baja/sobredispersión en el modelo
E2 = resid(model_1, type = "pearson")
N  = nrow(dat_saics)
p  = length(coef(model_1))   
sum(E2^2) / (N - p)


#__________________________________
#       glm Binomial negativo
#==================================
model_2 = glm.nb(loboc ~ ano + repro + estado_mar + dist_lob  + total_aves + trim, 
                 data = dat_saics, na.action = na.omit, maxit = 2000)

# Comprobamos dispersión nuevamente
E2_2 = resid(model_2, type = "pearson")
N  = nrow(dat_saics)
p  <- length(coef(model_2)) + 1  # '+1' is for variance parameter in NB
sum(E2_2^2) / (N - p)


m2 <- update(model_2, . ~ . - ano)  # Hacer esto sólo si le agrego la variable ano
anova(model_2, m2)

(est <- cbind(Estimate = coef(model_2), confint(model_2)))
exp(est)

# Modelo zero-inflated Poisson
model_3 <- zeroinfl(loboc ~ ano + repro + tipo_agreg + dist_lob + prof_med | ## Predictor for the Poisson process
                    ano + repro + tipo_agreg + dist_lob + prof_med,  ## Predictor for the Bernoulli process;
                    dist = 'poisson',
                    data = dat_saics)

# Estadístico de dispersión 2
E2_3 <- resid(model_3, type = "pearson")
N  <- nrow(dat_saics)
p  <- length(coef(model_3))  
sum(E2_3^2) / (N - p)


summ(model_2, exp = T) #could exp the coefficients but this is

# Para el modelo zero-inflado Poisson
tmp <- summary(model_3) #doing this manually
tmp$coefficients$count[-7, 1] %>% exp()



#=============================================================
#               Zero inflated negative binomial
#=============================================================

# Incluimos sin variable 'estado_mar'
model_4 <- zeroinfl(loboc ~ num_lan + ano + trim + dist_lob,
                    dist = 'negbin',
                    data = dat_saics)
summary(model_4) 


# Estadístico de dispersión
E2_4 <- resid(model_4, type = "pearson")
N  <- nrow(dat_saics)
p  <- length(coef(model_4)) + 1 # 
sum(E2_4^2) / (N - p)


model_4_2 <- zeroinfl(loboc ~ ano + trim + dist_lob + dis_cost,
                    dist = 'negbin',
                    data = dat_saics)
summary(model_4_2) 


# Estadístico de dispersión
E2_4_2 <- resid(model_4_2, type = "pearson")
N  <- nrow(dat_saics)
p  <- length(coef(model_4_2)) + 1 # 
sum(E2_4_2^2) / (N - p)



# Para el modelo zero-inflado Poisson
tmp <- summary(model_4) #doing this manually
tmp$coefficients$count[-7, 1] %>% exp()

(est <- cbind(Estimate = coef(model_2), confint(model_2)))
exp(est)

AIC(model_2, model_4, model_4_2)

summary(model_2)


# ELEGIDO MODEL_2 por la interpreatibilidad de los resultados

#===================================================================
# Modelado con GAM para incluir latitud y longitud en la estimación
#===================================================================
m0 <- gam(loboc ~ num_lan + ano + trim + dist_lob +
          s(lon, lat, k = 100), 
          data = dat_saics,
          family = nb, method = "ML")


m1 <- gam(loboc ~ num_lan + ano + trim + dist_lob +
          s(lon,lat, bs="gp",k=100, m=c(5,5)), 
          data = dat_saics,
          family = nb, method = "ML")



m2 <- gam(loboc ~ num_lan + ano + trim + dist_lob +
            s(lon,lat, bs="ds", m=c(5,.5),k=100), 
          data = dat_saics,
          family = nb, method = "ML")


m3 <- gam(loboc ~ num_lan + ano + trim + dist_lob +
            s(lon,lat, bs="tp", k=100), 
          data = dat_saics,
          family = nb, method = "ML")

AIC(m0, m1, m2, m3)
anova(m0, m1, m2, m3, test="Chisq")
gam.check(m2, k.rep=1000)

AIC(model_2, model_4, m2)
summary(model_2)


# Modelo elegido model_2









#================================================================================
#                     ANCHOVETA INDUSTRIAL ZONA NORTE
#================================================================================

# Vamos a filtrar por JUREL INDUSTRIAL CENTRO SUR (jics)
dat_aizn <- dat_fit %>% filter(nom_pesq %in% "ANCHOVETA INDUSTRIAL ZONA NORTE")
glimpse(dat_aizn) 

#dat_aizn = droplevels(dat_aizn) 

# Sacamos el nombre de la pesquería para poder realizar los constrastes
str(dat_aizn)
dat_aizn = dat_aizn[, c(-1, -2, -3, -6)] # Sacamos del data.frame 'cod_pesq', 'nom_pesq', esp_obj'
table(dat_aizn$loboc)
length(dat_aizn$loboc)

mean(dat_aizn$loboc)
1.521852^0 * exp(-1.521852) / factorial(0)


p = ggplot(dat_aizn, aes(x=loboc))+ geom_histogram(color="black", fill="grey")
p  + theme_bw()

p + scale_color_brewer(palette="black") +
  theme_classic()+theme(legend.position="top") +
  labs(title = "Anchoveta Industrial Zona Norte", y = "Frecuencia", x = "N° ejemplares capturados") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=11, face = "bold"),
        axis.title=element_text(size=11,face="bold")) +
  theme(plot.title = element_text(size = 14, face = "bold"))



#==================================
# Modelos lineales generalizados
#==================================

#__________________________________
#           Modelo Poisson
#==================================
model_1 = glm(loboc ~ num_lan + ano + repro + estado_mar + tipo_agreg + 
              tsm + captura.to + dist_lob + prof_med + trim + total_aves, 
              family = 'poisson', data = dat_aizn, na.action = na.omit)

## Verificiamos si existe baja/sobredispersión en el modelo
E2 = resid(model_1, type = "pearson")
N  = nrow(dat_aizn)
p  = length(coef(model_1))   
sum(E2^2) / (N - p)


#__________________________________
#       glm Binomial negativo
#==================================
model_2 = glm.nb(loboc ~ ano + dis_cost + tsm + prof_med + dist_lob, 
                 data = dat_aizn, na.action = na.omit, maxit = 2000)

# Comprobamos dispersión nuevamente
E2_2 = resid(model_2, type = "pearson")
N  = nrow(dat_aizn)
p  <- length(coef(model_2)) + 1  # '+1' is for variance parameter in NB
sum(E2_2^2) / (N - p)


m2 <- update(model_2, . ~ . - ano)  # Hacer esto sólo si le agrego la variable ano
anova(model_2, m2)

(est <- cbind(Estimate = coef(model_2), confint(model_2)))
exp(est)

# Modelo zero-inflated Poisson
model_3 <- zeroinfl(loboc ~ ano + dis_cost + tsm + prof_med + dist_lob | ## Predictor for the Poisson process
                    ano + dis_cost + tsm + prof_med + dist_lob,  ## Predictor for the Bernoulli process;
                    dist = 'poisson',
                    data = dat_aizn)

# Estadístico de dispersión 2
E2_3 <- resid(model_3, type = "pearson")
N  <- nrow(dat_aizn)
p  <- length(coef(model_3))  
sum(E2_3^2) / (N - p)


summ(model_2, exp = T) #could exp the coefficients but this is

# Para el modelo zero-inflado Poisson
tmp <- summary(model_3) #doing this manually
tmp$coefficients$count[-7, 1] %>% exp()



#=============================================================
#               Zero inflated negative binomial
#=============================================================

# Incluimos sin variable 'estado_mar'
model_4 <- zeroinfl(loboc ~ ano + dis_cost + tsm + prof_med + dist_lob,
                    dist = 'negbin',
                    data = dat_aizn)
summary(model_4) 


# Estadístico de dispersión
E2_4 <- resid(model_4, type = "pearson")
N  <- nrow(dat_aizn)
p  <- length(coef(model_4)) + 1 # 
sum(E2_4^2) / (N - p)


# Para el modelo zero-inflado Poisson
tmp <- summary(model_4) #doing this manually
tmp$coefficients$count[-7, 1] %>% exp()

(est <- cbind(Estimate = coef(model_2), confint(model_2)))
exp(est)

AIC(model_2, model_4)
summary(model_2)

# ELEGIDO MODEL_4

#===================================================================
# Modelado con GAM para incluir latitud y longitud en la estimación
#===================================================================
m0 <- gam(loboc ~ ano + dis_cost + tsm + prof_med + dist_lob + 
            s(lon, lat, k = 100), 
          data = dat_aizn,
          family = nb, method = "ML")


m1 <- gam(loboc ~ ano + dis_cost + tsm + prof_med + dist_lob +
            s(lon,lat, bs="gp",k=100, m=c(5,5)), 
          data = dat_aizn,
          family = nb, method = "ML")



m2 <- gam(loboc ~ ano + dis_cost + tsm + prof_med + dist_lob +
            s(lon,lat, bs="ds", m=c(5,.5),k=100), 
          data = dat_aizn,
          family = nb, method = "ML")


m3 <- gam(loboc ~ ano + dis_cost + tsm + prof_med + dist_lob +
            s(lon,lat, bs="tp", k=100), 
          data = dat_aizn,
          family = nb, method = "ML")

AIC(m0, m1, m2, m3)
anova(m0, m1, m2, m3, test="Chisq")
gam.check(m2, k.rep=1000)

AIC(model_2, m1)

# Modelo elegido m2
summary(m1)


library(mgcViz)
m1 <- getViz(m1)

print(plot(m1, select = 2:7), pages = 1)


print(plot(m1, allTerms = TRUE), pages = 1)


vis.gam(x = m1, # GAM object
        view = c("lon","lat"), # variables
        plot.type = "persp") # kind of plot

qq(m1, method = "simul1", a.qqpoi = list("shape" = 1), a.ablin = list("linetype" = 2))
qq(m1, rep = 20, showReps = T, CI = "none", a.qqpoi = list("shape" = 19), a.replin = list("alpha" = 0.2))

o <- qq(m1, rep = 10, method = "simul1", CI = "normal", showReps = TRUE,
        ngr = 1e2, a.replin = list(alpha = 0.1), a.qqpoi = list(shape = 19))
o 

gridPrint(o, zoom(o, xlim = c(2, 2.5), ylim = c(2, 2.5)), ncol = 2)

vis.gam(x = m1, # GAM object
        view = c("lon","lat"), # variables
        plot.type = "contour", main='Efecto espacial',
        color = "heat") # kind of plot










#================================================================================
#                      SARDINA AUSTRAL ARTESANAL CENTRO SUR
#================================================================================

# Vamos a filtrar por JUREL INDUSTRIAL CENTRO SUR (jics)
dat_sausacs <- dat_fit %>% filter(nom_pesq %in% "SARDINA AUSTRAL ARTESANAL CENTRO SUR")
glimpse(dat_sausacs) 

#dat_sausacs = droplevels(dat_sausacs) 

# Sacamos el nombre de la pesquería para poder realizar los constrastes
str(dat_sausacs)
dat_sausacs = dat_sausacs[, c(-1, -2, -3, -6)] # Sacamos del data.frame 'cod_pesq', 'nom_pesq', esp_obj'
# 'mes' (vamos a probar con trim)
table(dat_sausacs$lobo)
length(dat_sausacs$lobo)

mean(dat_sausacs$loboc)
2.815217^0 * exp(-2.815217) / factorial(0)

p = ggplot(dat_sausacs, aes(x=loboc))+ geom_histogram(color="black", fill="grey")
p  + theme_bw()

p + scale_color_brewer(palette="black") +
  theme_classic()+theme(legend.position="top") +
  labs(title = "Sardina Austral Artesanal Centro Sur", y = "Frecuencia", x = "N° ejemplares capturados") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=11, face = "bold"),
        axis.title=element_text(size=11,face="bold")) +
  theme(plot.title = element_text(size = 14, face = "bold"))



#==================================
# Modelos lineales generalizados
#==================================

#__________________________________
#           Modelo Poisson
#==================================
model_1 = glm(loboc ~ num_lan + ano + repro + estado_mar + tipo_agreg + 
                tsm + captura.to + dist_lob + prof_med + trim + total_aves, 
              family = 'poisson', data = dat_sausacs, na.action = na.omit)

## Verificiamos si existe baja/sobredispersión en el modelo
E2 = resid(model_1, type = "pearson")
N  = nrow(dat_saacs)
p  = length(coef(model_1))   
sum(E2^2) / (N - p)


#__________________________________
#       glm Binomial negativo
#==================================
model_2 = glm.nb(loboc ~ ano + dis_cost + prof_med + dist_lob + trim, 
                 data = dat_sausacs, na.action = na.omit, maxit = 2000)

# Comprobamos dispersión nuevamente
E2_2 = resid(model_2, type = "pearson")
N  = nrow(dat_sausacs)
p  <- length(coef(model_2)) + 1  # '+1' is for variance parameter in NB
sum(E2_2^2) / (N - p)


m2 <- update(model_2, . ~ . - num_lanc)  # Hacer esto sólo si le agrego la variable ano
anova(model_2, m2)

(est <- cbind(Estimate = coef(model_2), confint(model_2)))
exp(est)

# Modelo zero-inflated Poisson
model_3 <- zeroinfl(loboc ~ ano + dis_cost + prof_med + dist_lob + trim | ## Predictor for the Poisson process
                    ano + dis_cost + prof_med + dist_lob + trim,  ## Predictor for the Bernoulli process;
                    dist = 'poisson',
                    data = dat_sausacs)

# Estadístico de dispersión 2
E2_3 <- resid(model_3, type = "pearson")
N  <- nrow(dat_sausacs)
p  <- length(coef(model_3))  
sum(E2_3^2) / (N - p)


summ(model_2, exp = T) #could exp the coefficients but this is

# Para el modelo zero-inflado Poisson
tmp <- summary(model_3) #doing this manually
tmp$coefficients$count[-7, 1] %>% exp()



#=============================================================
#               Zero inflated negative binomial
#=============================================================

# Incluimos sin variable 'estado_mar'
model_4 <- zeroinfl(loboc ~ ano + dis_cost + prof_med + dist_lob + trim,
                    dist = 'negbin',
                    data = dat_sausacs)
summary(model_4) 


# Estadístico de dispersión
E2_4 <- resid(model_4, type = "pearson")
N  <- nrow(dat_sausacs)
p  <- length(coef(model_4)) + 1 # 
sum(E2_4^2) / (N - p)


# Para el modelo zero-inflado Poisson
tmp <- summary(model_4) #doing this manually
tmp$coefficients$count[-7, 1] %>% exp()

(est <- cbind(Estimate = coef(model_2), confint(model_2)))
exp(est)

AIC(model_2, model_4)


# ELEGIDO MODEL_2

#===================================================================
# Modelado con GAM para incluir latitud y longitud en la estimación
#===================================================================
m0 <- gam(loboc ~ ano + dis_cost + prof_med + dist_lob + trim +
            s(lon, lat, k = 100), 
          data = dat_sausacs,
          family = nb, method = "ML")


m1 <- gam(loboc ~ ano + dis_cost + prof_med + dist_lob + trim +
            s(lon,lat, bs="gp",k=100, m=c(5,1)), 
          data = dat_sausacs,
          family = nb, method = "ML")



m2 <- gam(loboc ~ ano + dis_cost + prof_med + dist_lob + trim +
            s(lon,lat, bs="ds", m=c(5,.5),k=100), 
          data = dat_sausacs,
          family = nb, method = "ML")


m3 <- gam(loboc ~ ano + dis_cost + prof_med + dist_lob + trim +
            s(lon,lat, bs="tp", k=100), 
          data = dat_sausacs,
          family = nb, method = "ML")

AIC(m0, m1, m2, m3)
anova(m0, m1, m2, m3, test="Chisq")
gam.check(m2, k.rep=1000)

AIC(model_2, m1)

summary(m1)

# Modelo elegido m1

library(mgcViz)
m1 <- getViz(m1)

print(plot(m1, select = 2:7), pages = 1)


print(plot(m1, allTerms = TRUE), pages = 1)


vis.gam(x = m1, # GAM object
        view = c("lon","lat"), # variables
        plot.type = "persp") # kind of plot

qq(m1, method = "simul1", a.qqpoi = list("shape" = 1), a.ablin = list("linetype" = 2))
qq(m1, rep = 20, showReps = T, CI = "none", a.qqpoi = list("shape" = 19), a.replin = list("alpha" = 0.2))

o <- qq(m1, rep = 10, method = "simul1", CI = "normal", showReps = TRUE,
        ngr = 1e2, a.replin = list(alpha = 0.1), a.qqpoi = list(shape = 19))
o 

gridPrint(o, zoom(o, xlim = c(2, 2.5), ylim = c(2, 2.5)), ncol = 2)

vis.gam(x = m1, # GAM object
        view = c("lon","lat"), # variables
        plot.type = "contour", main='Efecto espacial',
        color = "heat") # kind of plot

gam.check(m1)



