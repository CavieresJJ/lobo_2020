#===========================================================
#                   Modelos preliminares
#===========================================================


## Histograma de la variable respuesta
library(ggplot2)
library(pscl)
library(boot)
library(MASS)


# Determinamos las variables de importancia en función del menor AIC 
# (función step)

setwd(".....")

dat_fit = read.csv("dat_fit.csv", header=T)
dim(dat_fit)
head(dat_fit)
glimpse(dat_fit)

options(scipen=999)

# VARIABLES RESPUESTA
dat_fit$N_EJEMPLAR = as.numeric(dat_fit$N_EJEMPLAR)

dat_fit$ANO = as.factor(dat_fit$ANO)
dat_fit$MES = as.factor(dat_fit$MES)
dat_fit$TRIM = as.factor(dat_fit$TRIM)
dat_fit$COD_PESQUERIA = as.factor(dat_fit$COD_PESQUERIA)
dat_fit$TIPO_DE_RED = as.factor(dat_fit$TIPO_DE_RED)
dat_fit$CLASE_LANCE = as.factor(dat_fit$CLASE_LANCE)
dat_fit$MODELO_EXCLUSION = as.factor(dat_fit$MODELO_EXCLUSION)
dat_fit$INTENSIDAD_VIENTO_AR = as.factor(dat_fit$INTENSIDAD_VIENTO_AR)
dat_fit$ESPECIE_OBJETIVO_LANCE = as.factor(dat_fit$ESPECIE_OBJETIVO_LANCE)

dat_fit = dat_fit[, c(-1, -14, -15)]

glimpse(dat_fit)
head(dat_fit)


full = glm(N_EJEMPLAR~.,data = dat_fit, family = poisson)
step = stepAIC(full, trace = FALSE)
step$anova


backward  = stepAIC(glm(N_EJEMPLAR~.,data = dat_fit, family = poisson),direction="backward")
forward   = stepAIC(glm(N_EJEMPLAR~.,data = dat_fit, family = poisson),direction="forward")
both      = stepAIC(glm(N_EJEMPLAR~.,data = dat_fit, family = poisson),direction="both")


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

## Poisson GLM
model_1 = glm(N_EJEMPLAR ~ ANO + MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                PROFUNDIDAD_PROM, family = 'poisson', data = dat_fit, na.action = na.omit)
summary(model_1)
model_1$aic


## Verificiamos si existe baja/sobredispersión en el modelo
E2 = resid(model_1, type = "pearson")
N  = nrow(dat_fit)
p  = length(coef(model_1))   
sum(E2^2) / (N - p)


## Parece ser que sí existe soberdispersión (varianza no es igual a la media) 
## por tanto probaremos con una regresión binomial negativa

model_2 = glm.nb(N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                   MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                   PROFUNDIDAD_PROM, data = dat_fit, na.action = na.omit, maxit = 1000)
summary(model_2)



# Comprobamos dispersión nuevamente
E2_2 = resid(model_2, type = "pearson")
N  = nrow(dat_fit)
p  <- length(coef(model_2)) + 1  # '+1' is for variance parameter in NB
sum(E2_2^2) / (N - p)

# Luce mejor que el modelo de conteo Poisson pero esto igualmente puede generar
# problemas en la inferencia de los parámetros

# Modelo zero-inflated Poisson
model_3 <- zeroinfl(N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                      MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                      PROFUNDIDAD_PROM | ## Predictor for the Poisson process
                      ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                      MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                      PROFUNDIDAD_PROM,  ## Predictor for the Bernoulli process;
                      dist = 'poisson',
                      data = dat_fit)

summary(model_3) # error en la estimación de mat var-cov
model_3


# Modelo zero-inflated Poisson 2
model_3_2 <- zeroinfl(N_EJEMPLAR ~ ANO + MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                        MODELO_EXCLUSION + PROFUNDIDAD_PROM  | ## Predictor for the Poisson process
                        ANO + MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                        MODELO_EXCLUSION + PROFUNDIDAD_PROM,  ## Predictor for the Bernoulli process;
                        dist = 'poisson',
                        data = dat_fit)

summary(model_3_2) # 
model_3_2$loglik


# Hay dos variables que presentan problemas para ajustar un modelozero-inflado Poisson
# INTENSIDAD_VIENTO_AR y ESPECIE_OBJETIVO_LANCE



# Estadístico de dispersión 2
E2_32 <- resid(model_3_2, type = "pearson")
N  <- nrow(dat_fit)
p  <- length(coef(model_3_2))  
sum(E2_32^2) / (N - p)



#=============================================================
#               Zero inflated negative binomial
#=============================================================

# Incluimos ANO y TRIM
model_4 <- zeroinfl(N_EJEMPLAR ~ ANO + MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                    MODELO_EXCLUSION +PROFUNDIDAD_PROM ,
                    dist = 'negbin',
                    data = dat_fit)
summary(model_4)


# Incluimos sólo 'ANO'
model_4_2 <- zeroinfl(N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                      MODELO_EXCLUSION + PROFUNDIDAD_PROM,
                      dist = 'negbin',
                      data = dat_fit)
summary(model_4_2)


# Estadístico de dispersión
E2_4 <- resid(model_4, type = "pearson")
N  <- nrow(dat_fit)
p  <- length(coef(model_4)) + 1 # 
test1= sum(E2_4^2) / (N - p)



E2_42 <- resid(model_4_2, type = "pearson")
N  <- nrow(dat_fit)
p  <- length(coef(model_4_2)) + 1 # 
test2= sum(E2_42^2) / (N - p)



library(lmtest)
lrtest(model_4, model_4_2)  

# Comparación modelos bajo/sobredispersión
c(test1, test2)
#====================================================
# GLM BAYESIANO
# Controles en INLA
cres = list(return.marginals.predictor = TRUE, return.marginals.random = TRUE)
cinla <- list(strategy = 'adaptive', int.strategy = 'eb')  # Estrategias de estimación

# Initial values of parameters
ini.zb <- c(1.834,-6.141)

formula = N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + MODELO_EXCLUSION
  
inla = inla(formula,
            family = 'zeroinflatednbinomial0',
            data = dat_fit, 
            control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
            control.predictor=list(compute=TRUE, link = 1),
            control.results = cres, control.inla = cinla,
            control.mode = list(theta = ini.zb, restart = TRUE),
            verbose=TRUE)

summary(inla)

library(brinla)
bri.fixed.plot(inla)


z.b0 =  inla$marginals.fixed[[1]]
plot_z.b0 = ggplot(data.frame(inla.smarginal(z.b0)), aes(x, y)) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Distribución posterior para el intercepto", y = "Intercepto", x = "") + 
  theme_bw()

