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

# Vamos a seleccionar las variables que creemos importantes para el análisis
dat_fit = dat2[, c("N_EJEMPLAR", "MES", "ANO", "COD_PESQUERIA", "TIPO_DE_RED", "CLASE_LANCE", "MODELO_EXCLUSION", "INTENSIDAD_VIENTO_AR", 
                   "ESPECIE_OBJETIVO_LANCE", "PROFUNDIDAD_PROM", "LAT2", "LON2")]


glimpse(dat_fit)

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
step$deviance; backward$deviance; forward$deviance; both$deviance   # Tienen igual devianza así que podemos ocupar cualquier

# Menos AIC es forward
formula(step)


#===============================================================================================================
#                                         MODELOS PRELIMINARES
#===============================================================================================================

## Poisson GLM
model_1 = glm(N_EJEMPLAR ~ MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                PROFUNDIDAD_PROM, family = 'poisson', data = dat_fit, na.action = na.omit)
summary(model_1)


## Verificiamos si existe baja/sobredispersión en el modelo
E2 = resid(model_1, type = "pearson")
N  = nrow(dat2)
p  = length(coef(model_1))   
sum(E2^2) / (N - p)


## Parece ser que sí existe soberdispersión (varianza no es igual a la media) 
## por tanto probaremos con una regresión binomial negativa

model_2 = glm.nb(N_EJEMPLAR ~ MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                   MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                   PROFUNDIDAD_PROM, data = dat_fit, na.action = na.omit, maxit = 1000)
summary(model_2)


# Comprobamos dispersión nuevamente
E2_2 = resid(model_2, type = "pearson")
N  = nrow(dat2)
p  <- length(coef(model_2)) + 1  # '+1' is for variance parameter in NB
sum(E2_2^2) / (N - p)

# Luce mejor que el modelo de conteo Poisson pero esto igualmente puede generar
# problemas en la inferencia de los parámetros

# Modelo zero-inflated Poisson
model_3 <- zeroinfl(N_EJEMPLAR ~ MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                      MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                      PROFUNDIDAD_PROM | ## Predictor for the Poisson process
                      MES + ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                      MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                      PROFUNDIDAD_PROM,  ## Predictor for the Bernoulli process;
                      dist = 'poisson',
                      data = dat_fit)

summary(model_3) # error en la estimación de mat var-cov


#=================================================================================
#                     Modelo zero-inflated Poisson
#=================================================================================
model_3 <- zeroinfl(N_EJEMPLAR ~ MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                      MODELO_EXCLUSION + INTENSIDAD_VIENTO_AR + ESPECIE_OBJETIVO_LANCE + 
                      PROFUNDIDAD_PROM,
                      dist = 'poisson',
                      data = dat_fit, na.action = na.omit)

summary(model_3) # error en la estimación de mat var-cov




# Modelo zero-inflated Poisson 2
model_3_2 <- zeroinfl(N_EJEMPLAR ~ MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                        MODELO_EXCLUSION + PROFUNDIDAD_PROM  | ## Predictor for the Poisson process
                        MES + ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + MODELO_EXCLUSION + 
                        PROFUNDIDAD_PROM,  ## Predictor for the Bernoulli process;
                        dist = 'poisson',
                        data = dat_fit, na.action = na.omit)

summary(model_3_2) # 



# Hay dos variables que presentan problemas para ajustar un modelozero-inflado Poisson
# INTENSIDAD_VIENTO_AR y ESPECIE_OBJETIVO_LANCE



# Estadístico de dispersión 2
E2_32 <- resid(model_3_2, type = "pearson")
N  <- nrow(dat2)
p  <- length(coef(model_3_2))  
sum(E2_32^2) / (N - p)










#=============================================================
#               Zero inflated negative binomial
#=============================================================

# Incluimos el 'MES'
model_4 <- zeroinfl(N_EJEMPLAR ~ MES + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                    MODELO_EXCLUSION + PROFUNDIDAD_PROM,
                    dist = 'negbin',
                    data = dat_fit, na.action = na.omit)
summary(model_4)


# Incluimos el 'ANO'
model_4_2 <- zeroinfl(N_EJEMPLAR ~ ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                      MODELO_EXCLUSION + PROFUNDIDAD_PROM,
                      dist = 'negbin',
                      data = dat_fit, na.action = na.omit)
summary(model_4_2)



# Incluimos el 'ANO' y 'MES' 
model_4_3 <- zeroinfl(N_EJEMPLAR ~ MES + ANO + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                        MODELO_EXCLUSION + PROFUNDIDAD_PROM,
                        dist = 'negbin',
                        data = dat_fit, na.action = na.omit)

summary(model_4_3) ##===> error de estimación mat var-cov



# Estadístico de dispersión
E2_4 <- resid(model_4, type = "pearson")
N  <- nrow(dat2)
p  <- length(coef(model_4)) + 1 # '+1' is due to theta
test1= sum(E2_4^2) / (N - p)



E2_42 <- resid(model_4_2, type = "pearson")
N  <- nrow(dat2)
p  <- length(coef(model_4_2)) + 1 # '+1' is due to theta
test2= sum(E2_42^2) / (N - p)

library(lmtest)
lrtest(model_3_2, model_4, model_4_2)  # Mejor modelo con variable 'MES'

# Comparación modelos bajo/sobredispersión
c(test1, test2)

# Me quedo con el modelo que contiene a la variable 'AÑO' ya que su test está
# cercano a 1 y puede explicar la variación anual de las capturas.


mnull <- update(model_4_2, . ~ 1)
pchisq(2 * (logLik(model_4_2) - logLik(mnull)), df = 6, lower.tail = FALSE)





# Vamos a probar categorizando las variables "MES" para que tenga menos niveles

levels(dat_fit$MES) # Tiene 13 niveles y vamos a dejar en 4 (TRIMESTRES)
# Niveles de 1 a 3   =  (1)
# Niveles de 4 A 6   =  (2)
# Niveles de 7 A 9   =  (3)
# Niveles de 10 A 12 =  (4)



library(varhandle)
dat_fit$MES2  = dat_fit$MES

dat_fit$MES2 = unfactor(dat_fit$MES)


dat_fit$TRIM = with(dat_fit$MES2, ifelse(dat_fit$MES2 < 4, "1", ifelse(dat_fit$MES2 < 7, "2", ifelse(dat_fit$MES2 < 10, "3", "4"))))
class(dat2$INTENSIDAD_VIENTO_AR3)
dat2$INTENSIDAD_VIENTO_AR3 = as.factor(dat2$INTENSIDAD_VIENTO_AR3)
levels(dat2$INTENSIDAD_VIENTO_AR3)


dat_fit$TRIM = ifelse(dat_fit$MES2 < 4 ,1,
                      ifelse(dat_fit$MES2 < 7, 2,
                      ifelse(dat_fit$MES2 < 10, 3, 4)))

dat_fit$TRIM = as.factor(dat_fit$TRIM)

model_4_4 <- zeroinfl(N_EJEMPLAR ~ ANO + TRIM + COD_PESQUERIA + TIPO_DE_RED + CLASE_LANCE + 
                      MODELO_EXCLUSION + PROFUNDIDAD_PROM,
                      dist = 'negbin',
                      data = dat_fit)

summary(model_4_4) #----> Se produce error en estimación de las desviaciones estandard asociadas a 
                   # la desviación estandar



