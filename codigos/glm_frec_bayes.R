#===========================================================
#                   Modelos preliminares
#===========================================================


## Histograma de la variable respuesta
library(ggplot2)
library(pscl)
library(boot)
library(MASS)
library(INLA)


# Determinamos las variables de importancia en funci?n del menor AIC 
# (funci?n step)

setwd("C:/Users/Usuario/Desktop/Projects/Ejecución/IFOP_lobo_2020/data")

dat_fit = read.csv("dat_fit5.csv", header=T)
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
#dat_fit$ESPECIE_OBJETIVO_LANCE2 = as.factor(dat_fit$ESPECIE_OBJETIVO_LANCE2)
dim(dat_fit)

dat_fit = dat_fit[, c(-1, -3, -4, -21, -22)]

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
both$deviance   

# Comparar AIC
step$aic
backward$aic
forward$aic
both$aic

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
model_1$aic


model_1_1 = glm(N_EJEMPLAR ~ ANO + MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                MODELO_EXCLUSION + PROFUNDIDAD_PROM, 
                family = 'poisson', data = dat_fit, na.action = na.omit)
model_1_1$aic





## Verificiamos si existe baja/sobredispersi?n en el modelo
E2 = resid(model_1, type = "pearson")
N  = nrow(dat_fit)
p  = length(coef(model_1))   
sum(E2^2) / (N - p)



## Verificiamos si existe baja/sobredispersi?n en el modelo
E2 = resid(model_1_1, type = "pearson")
N  = nrow(dat_fit)
p  = length(coef(model_1_1))   
sum(E2^2) / (N - p)






#=======================================================================================
#                              Predictor 1 
#                           Binomial Negativa
#=======================================================================================
model_2 = glm.nb(N_EJEMPLAR ~ ANO + MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                   MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                   PROFUNDIDAD_PROM, data = dat_fit, na.action = na.omit, maxit = 1000)
anova(model_2)
#termplot(model_2)

# Comprobamos dispersi?n nuevamente
E2_2 = resid(model_2, type = "pearson")
N  = nrow(dat_fit)
p  <- length(coef(model_2)) + 1  # '+1' is for variance parameter in NB
sum(E2_2^2) / (N - p)



#=======================================================================================
#                              Predictor 2 
#                           Binomial Negativa
#=======================================================================================
model_2_1 = glm.nb(N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                   MODELO_EXCLUSION + PROFUNDIDAD_PROM, data = dat_fit, na.action = na.omit, maxit = 1000)
anova(model_2_1)

# Comprobamos dispersi?n nuevamente
E2_2 = resid(model_2_1, type = "pearson")
N  = nrow(dat_fit)
p  <- length(coef(model_2_1)) + 1  # '+1' is for variance parameter in NB
pred2_bn = sum(E2_2^2) / (N - p)


AIC(model_2, model_2_1)
anova(model_2, model_2_1)

# termplot(model_2)
# termplot(model_2_1)


# Luce mejor que el modelo de conteo Poisson pero esto igualmente puede generar
# problemas en la inferencia de los par?metros


#=================================================================================
#                          Modelo zero-inflated Poisson 
#                                  Predictor 1
#=================================================================================

model_3 <- zeroinfl(N_EJEMPLAR ~ ANO + MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                      MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                      PROFUNDIDAD_PROM | ## Predictor for the Poisson process
                      ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                      MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                      PROFUNDIDAD_PROM,  ## Predictor for the Bernoulli process;
                      dist = 'poisson',
                      data = dat_fit)

summary(model_3) # error en la estimaciin de mat var-cov
model_3$loglik


#=================================================================================
#                          Modelo zero-inflated Poisson 
#                                  Predictor 2
#=================================================================================
model_3_2 <- zeroinfl(N_EJEMPLAR ~ ANO + MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                        MODELO_EXCLUSION + PROFUNDIDAD_PROM  | ## Predictor for the Poisson process
                        ANO + MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                        MODELO_EXCLUSION + PROFUNDIDAD_PROM,  ## Predictor for the Bernoulli process;
                        dist = 'poisson',
                        data = dat_fit)

model_3_2$loglik


# Hay dos variables que presentan problemas para ajustar un modelozero-inflado Poisson
# INTENSIDAD_VIENTO_AR y ESPECIE_OBJETIVO_LANCE



# Estadistico dispersion para predictor 2
E2_32 <- resid(model_3_2, type = "pearson")
N  <- nrow(dat_fit)
p  <- length(coef(model_3_2))  
sum(E2_32^2) / (N - p)



#=================================================================================
#                         Zero inflated negative binomial
#                                  Predictor 1
#=================================================================================
model_4 <- zeroinfl(N_EJEMPLAR ~ ANO + MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                    MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                    PROFUNDIDAD_PROM,
                    dist = 'negbin',
                    data = dat_fit)

summary(model_4) #==> Problemas de estimacion


# Estadistico de dispersiOn
E2_4 <- resid(model_4, type = "pearson")
N  <- nrow(dat_fit)
p  <- length(coef(model_4)) + 1 # 
pred1 = sum(E2_4^2) / (N - p)



#=================================================================================
#                         Zero inflated negative binomial
#                                  Predictor 2
#=================================================================================

# Con variable 'MES'
model_4_1 <- zeroinfl(N_EJEMPLAR ~ ANO + MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                        MODELO_EXCLUSION + PROFUNDIDAD_PROM,
                      dist = 'negbin',
                      data = dat_fit)
summary(model_4_1) #==> Problemas de estimacion


E2_4_1 <- resid(model_4_1, type = "pearson")
N  <- nrow(dat_fit)
p  <- length(coef(model_4_1)) + 1 # 
pred2 = sum(E2_4_1^2) / (N - p)




# Sin variable 'MES'
model_4_2 <- zeroinfl(N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                    MODELO_EXCLUSION + PROFUNDIDAD_PROM,
                    dist = 'negbin',
                    data = dat_fit)
summary(model_4_2)


E2_4_2 <- resid(model_4_2, type = "pearson")
N  <- nrow(dat_fit)
p  <- length(coef(model_4_2)) + 1 # 
pred2_no_mes = sum(E2_4_2^2) / (N - p)




# Con variable 'TRIM'
model_4_3 <- zeroinfl(N_EJEMPLAR ~ ANO + TRIM + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                        MODELO_EXCLUSION + PROFUNDIDAD_PROM,
                      dist = 'negbin',
                      data = dat_fit)
summary(model_4_3)


E2_4_3 <- resid(model_4_3, type = "pearson")
N  <- nrow(dat_fit)
p  <- length(coef(model_4_3)) + 1 # 
pred2_trim = sum(E2_4_3^2) / (N - p)




# Comparacion loglike
library(lmtest)
lrtest(model_4, model_4_1, model_4_2, model_4_3)  

# Comparacion AIC
AIC(model_2, model_2_1, model_4, model_4_1, model_4_2, model_4_3)


# Comparacion AIC (Tabla informe)
AIC(model_2, model_2_1, model_4, model_4_2)



# Comparacion baja/dispersion
c(pred2_bn, pred1, pred2, pred2_no_mes, pred2_trim)


#===========================================
#         Mejor predictor lineal
#===========================================

# N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
# MODELO_EXCLUSION + PROFUNDIDAD_PROM




# # Muestra todos los coeficientes estimados (sin intercepto)
# coef(zeroinfl(N_EJEMPLAR ~ -1 + ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
#                         MODELO_EXCLUSION + PROFUNDIDAD_PROM,
#                       dist = 'negbin',
#                       data = dat_fit))
# 
# # No muestra todos los coeficientes (con intercepto)
# coef(zeroinfl(N_EJEMPLAR ~ 1 + ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
#                 MODELO_EXCLUSION + PROFUNDIDAD_PROM,
#               dist = 'negbin',
#               data = dat_fit))
# Muestra todos los coeficientes estimados


# plot(density(resid(model_2, type='response')), main = "Modelo Binomail Negativo")
# plot(density(resid(model_4_2, type='response')), main = "Modelo Zero-Inflated Binomail Negativo", col='red')
# 
# lines(density(resid(model_4_2, type='response')), main = "Modelo Zero-Inflated Binomail Negativo", col='red')
# 
# 
# plot(density(resid(model_2, type='pearson')))
# lines(density(resid(model_4_2, type='pearson')), col='red')
# 
# plot(density(resid(model_2, type='deviance')))
# lines(density(resid(model_4_2, type='deviance')), col='red')
# 
# 
# 
# 
# par(mfrow=c(1,2))
# plot(density(dat_fit$N_EJEMPLAR), xlim=c(0, 50), ylim=c(0, 5), main='M_0 y_hat')
# lines(density(predict(model_2, type='response')), col='red')
# 
# plot(density(dat_fit$N_EJEMPLAR), xlim=c(0, 50), ylim=c(0, 5), main='M_1 y_hat')
# lines(density(predict(m1, type='response')), col='red')
# 
# 
# par(mfrow=c(1,2))
# scatter.smooth(1:10555, rstandard(model_2, type='deviance'), col='gray')
# scatter.smooth(1:10555, rstandard(model_4_2, type='deviance'), col='gray')
# 
# 
# par(mfrow=c(1,2))
# qqnorm(statmod::qresid(model_2)); qqline(statmod::qresid(model_2))
# qqnorm(statmod::qresid(model_4_2)); qqline(statmod::qresid(model_4_2))


#===========================================================================
#                              GLM BAYESIANO
#===========================================================================
# Controles en INLA
library(INLA)
library(INLAutils)
cres = list(return.marginals.predictor = TRUE, return.marginals.random = TRUE)
cinla <- list(strategy = 'adaptive', int.strategy = 'eb')  # Estrategias de estimaci?n

# Initial values of parameters
ini.zb <- c(1.834)
ini.zb2 <- c(1.834, -1)


#=========================================================================
#                     Predictor 1 (sin variable mes)
#=========================================================================
formula1 = N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
           MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
           PROFUNDIDAD_PROM


inla1 = inla(formula1,
             family = 'nbinomial',
             data = dat_fit, 
             control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
             control.predictor=list(compute=TRUE, link = 1),
             control.results = cres, control.inla = cinla,
             control.mode = list(theta = ini.zb, restart = TRUE),
             verbose=TRUE)

summary(inla1)


inla1$waic$waic


#=========================================================================
#                   Predictor 2 (sin variable mes)
#=========================================================================
formula2 = N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + MODELO_EXCLUSION + PROFUNDIDAD_PROM
  
inla2 = inla(formula2,
            family = 'nbinomial',
            data = dat_fit, 
            control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
            control.predictor=list(compute=TRUE, link = 1),
            control.results = cres, control.inla = cinla,
            control.mode = list(theta = ini.zb, restart = TRUE),
            verbose=TRUE)

summary(inla2)

inla2$waic$waic



#=========================================================================
#           Predictor 1 => zeroinflatednbinomial0 (sin variable mes)
#=========================================================================

ini.zb2 <- c(1.834, -6.02)

formula3 = N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
           MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
           PROFUNDIDAD_PROM

inla3 = inla(formula3,
             family = 'zeroinflatednbinomial0',
             data = dat_fit, 
             control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
             control.predictor=list(compute=TRUE, link = 1),
             control.results = cres, control.inla = cinla,
             control.mode = list(theta = ini.zb2, restart = TRUE),
             verbose=TRUE)

summary(inla3)





#=========================================================================
#           Predictor 2 => zeroinflatednbinomial0 (sin variable mes)
#=========================================================================
formula4 = N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + MODELO_EXCLUSION + 
           PROFUNDIDAD_PROM

inla4 = inla(formula4,
             family = 'zeroinflatednbinomial0',
             data = dat_fit, 
             control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
             control.predictor=list(compute=TRUE, link = 1),
             control.results = cres, control.inla = cinla,
             control.mode = list(theta = ini.zb2, restart = TRUE),
             verbose=TRUE)


#==========================
# summary para cada modelo
#==========================

summary(inla1)
summary(inla2)
summary(inla3) # Problemas en la estimacion del WAIC y número efectivo de parametros
summary(inla4)


# Comparacion WAIC
inla1$waic$waic
inla2$waic$waic
inla3$waic$waic
inla4$waic$waic

inla1$dic$dic
inla2$dic$dic
inla3$dic$dic
inla4$dic$dic


slcpo <- function(m, na.rm = TRUE) {
  - sum(log(m$cpo$cpo), na.rm = na.rm)
}


c("Model_1" = slcpo(inla1), 
  "Model_2" = slcpo(inla2), 
  "Model_3" = slcpo(inla3), 
  "Model_4" = slcpo(inla4))


#=====================================================================================================
# Mejor modelo ==> Predictor 2 con familia 'nbinomial'
# N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + MODELO_EXCLUSION + PROFUNDIDAD_PROM
#=====================================================================================================


library(brinla)
bri.fixed.plot(inla2)


z.b0 =  inla1$marginals.fixed[[1]]
plot_z.b0 = ggplot(data.frame(inla.smarginal(z.b0)), aes(x, y)) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Distribución posterior para el intercepto", y = "Intercepto", x = "") + 
  theme_bw()


source("Efxplot_2.R")
x11()
Efxplot_2(list(m5))

summary(model_2)


# compare models with cross validation (lowest is best)
cval <- function(result=inla1){
  test <- sum(result$cpo$failure)
  if(test != 0) print("some cpo estimates failed the quality test. do not use.")
  return(-2 * sum(log(result$cpo$cpo), na.rm = T))
}
cval(inla1); cval(inla2); cval(inla3); cval(inla4)

# pick a model and get credible intervals for fixed effects
inla1$summary.fixed[,c(1,3,5)]

# get posterior probability that a parameter estimate is positive
1 - inla.pmarginal(0, inla1$marginals.fixed$COD_PESQUERIA9)

# get variation explained by random effects
brinla::bri.hyperpar.summary(inla1)[, c(1,3,5)]

# look at model fit to data
hist(inla2$cpo$pit) # should be uniformish
plot(inla2$summary.fitted.values$mean, dat_fit$N_EJEMPLAR); abline(0,1)
# x values are posterior mean fitted values that include fixed and random effects