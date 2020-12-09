rm(list=ls(all=TRUE))

#====================================================================================================================================
#                                         Proyecto IFOP 2020 - lobos marinos
#====================================================================================================================================

setwd("C:/Users/Usuario/Desktop/Projects/Ejecución/IFOP_lobo_2020/data")

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


# load data 
dat = read.csv("data_total2.csv", header=T)
dim(dat)
head(dat)
str(dat)

options(scipen=99)

#================================================================================================
# OBSERVACIONES:
# 101_1: Lobo marino común estado vivo
# 101_2: Lobo marino común estado muerto  
# 971_1: Lobo fino austral estado vivo
# 971_2: Lobo fino austral estado muerto
# 101: Total de ejemplares de lobo marino común capturados incidentalmente (vivos + muertos)
# 971: Total de ejemplares de lobo fino austral capturados incidentalmente  (vivos + muertos)
#================================================================================================


# Creamos una nueva variable para el total de lobos capturados incidentalmente
dat$X101_1 = as.numeric(dat$X101_1)
dat$X101_2 = as.numeric(dat$X101_2)
dat$X971_1 = as.numeric(dat$X971_1)
dat$X971_2 = as.numeric(dat$X971_2)

dat$X101 = as.numeric(dat$X101)
dat$X971 = as.numeric(dat$X971)


dat$N_EJEMPLAR = dat$X101 + dat$X971  # Variable conteo 'N_EJEMPLAR'
table(dat$N_EJEMPLAR)


dat$ESTADO_C1 = ifelse(dat$X101_1 == 0, 0, 1)  # Lobo común muerto = 0, vivo = 1
dat$ESTADO_F1 = ifelse(dat$X971_1 == 0, 0, 1)  # Lobo fino  muerto = 0, vivo = 1



# Creamos otra variable para ver la proporción de vivos y/muertos
dat$N_EJEMPLAR2 = ifelse(dat$N_EJEMPLAR > 0, "Con captura incidental", "Sin captura incidental")
table(dat$N_EJEMPLAR2)


dat$COD_BARCO = as.factor(dat$COD_BARCO)
dat$COD_PESQUERIA = as.factor(dat$COD_PESQUERIA)
dat$NOMBRE_PESQUERIA = as.factor(dat$NOMBRE_PESQUERIA)
dat$CLASE_LANCE = as.factor(dat$CLASE_LANCE)
dat$ESPECIE_OBJETIVO = as.factor(dat$ESPECIE_OBJETIVO)
dat$ESPECIE_OBJETIVO_LANCE = as.factor(dat$ESPECIE_OBJETIVO_LANCE)
dat$TIPO_DE_RED = as.factor(dat$TIPO_DE_RED)
dat$INTENSIDAD_VIENTO_AR = as.factor(dat$INTENSIDAD_VIENTO_AR)
dat$MODELO_EXCLUSION = as.factor(dat$MODELO_EXCLUSION)


# Creamos variables 'ANO', 'MES' y 'DIA'
dat$FECHA <-  as.Date(dat$FECHA_LANCE,'%m/%d/%Y')
dat$ANO <- as.numeric(format(dat$FECHA,'%Y'))
dat$MES <- as.numeric(format(dat$FECHA,'%m'))
dat$DIA <- as.numeric(format(dat$FECHA,'%d'))


# Vamos a categorizar la variable "MES" en "TRIM"
library(varhandle)

dat$TRIM = ifelse(dat$MES < 4 ,1,ifelse(dat$MES < 7, 2, ifelse(dat$MES < 10, 3, 4)))

# Filtramos las observaciones que efectivamente fueron registradas por el observador
# basada en la información de la variable 'OBS_CIAMT'
dat2 = dat %>% filter(OBS_CIAMT == "1")
length(dat2$COD_BARCO) # Tenemos 10560 observaciones
summary(dat2$LATITUD)  # Celdas con NA
summary(dat2$LONGITUD)  # Celdas con NA

dat2 = dat2 %>% filter(!is.na(dat2$LONGITUD))
dim(dat2)
summary(dat2$LATITUD) # Ya no hay celdas +

LAT = as.numeric(substring(dat2$LATITUD,1,2))+as.numeric(substring(dat2$LATITUD,3,4))/60
as.numeric(substring(dat2$LATITUD,5,6))/3600

LON <- as.numeric(substring(dat2$LONGITUD,1,2))+as.numeric(substring(dat2$LONGITUD,3,4))/60
as.numeric(substring(dat$LONGITUD,5,6))/3600


londec <- (-1*(LON))
latdec <- (-1*(LAT))

dat2$LON2 = londec
dat2$LAT2 = latdec

summary(dat2$LON2)
summary(dat2$LAT2)

# Revisamos si existen NA en las variables relevantes
glimpse(dat2)


#=======================================
# Variable 'ESPECIE_OBJETIVO_LANCE'
#=======================================
# De acuerdo a lo conversado con Marcelo, se debería utilizar la variable de mayor
# resolución, por tanto consideramos la variable 'ESPECIE_OBJETIVO_LANCE'
table(dat2$ESPECIE_OBJETIVO_LANCE, exclude = NULL)

#=======================================
# Variable 'COD_EXCLUSION'
#=======================================
table(dat2$COD_EXCLUSION, exclude = NULL)

#=======================================
# Variable 'CLASE_LANCE'
#=======================================
table(dat2$CLASE_LANCE, exclude = NULL)


#=======================================
# Variable 'TIPO_DE_RED'
#=======================================
table(dat2$TIPO_DE_RED, exclude = NULL)



#===================================================================================================
#                              Creación o transformación de variables
#===================================================================================================
table(dat2$ESPECIE_OBJETIVO_LANCE, exclude = NULL)

table(dat2$TIPO_DE_RED, dat2$ESPECIE_OBJETIVO_LANCE, exclude = NULL)

summary(dat2$ESPECIE_OBJETIVO_LANCE)
dat2$ESPECIE_OBJETIVO_LANCE <- factor(dat2$ESPECIE_OBJETIVO_LANCE, exclude = NULL, 
                           levels = c("1", "2", "3", "4", "5", "6", "7", "27", "29", "45", "96", "956", NA), 
                           labels = c("1", "2", "3", "4", "5", "6", "7", "27", "29", "45", "96", "956", "4"))
summary(dat2$ESPECIE_OBJETIVO_LANCE)



dat2 %>% filter(ESPECIE_OBJETIVO_LANCE == "956") # No hay información, la vamos a eliminar
dat2 = subset(dat2, ESPECIE_OBJETIVO_LANCE != "956")

dat2 = droplevels(dat2) 





# 1	Arrastre de Fondo
# 2	Media Agua
# 3	Multiproposito

# El 'ESPECIE_OBJETIVO_LANCE' tiene 12 NA asociado con 'TIPO_DE_RED' = 1
#                             tiene 2  NA asociado con 'TIPO_DE_RED' = 2

# 'TIPO_DE_RED' tiene 38 NA asociado a 'ESPECIE_OBJETIVO_LANCE' = 4


# Transformación NA en varilable 'CLASE_LANCE' 
# NA (celdas vacias) son pasados a codigo "1"
summary(dat2$CLASE_LANCE)
dat2$CLASE_LANCE <- factor(dat2$CLASE_LANCE, exclude = NULL, 
                           levels = c("1", "2", NA), 
                           labels = c("1", "2", "1"))
summary(dat2$CLASE_LANCE)



#==================================================================================
# En primera instancia ibamos a considerar la variable 'ESTADO_DEL_MAR'
table(dat2$ESTADO_DEL_MAR, exclude = NULL)
# Como se aprecia tenemos muchos NA que pueden causar confusión, y como
# la variable 'ESTADO_DEL_MAR' esta en función de la variable 'INTENSIDAD_DEL_VIENTO'
# consideramos esta variable ya que tiene pocos NA (25)

table(dat2$INTENSIDAD_VIENTO_AR, exclude = NULL)

dat2$INTENSIDAD_VIENTO_AR = factor(dat2$INTENSIDAD_VIENTO_AR, exclude = NULL, 
                           levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", NA), 
                           labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"))
table(dat2$INTENSIDAD_VIENTO_AR) # donde "13" es un nuevo nivel "sin información"


dat2 %>% filter(INTENSIDAD_VIENTO_AR == "12") # No hay información, la vamos a eliminar
dat2 = subset(dat2, INTENSIDAD_VIENTO_AR != "12")

dat2 = droplevels(dat2) 
summary(dat2$INTENSIDAD_VIENTO_AR)


#=================================
# Variable 'MODELO_DE_EXCLUSION'
#=================================

#   1 =  Con Dispositivo
#   2 =  Sin Dispositivo
#  	3 =  Con dispositivo modif.
# 	9 =  No observado
table(dat2$MODELO_EXCLUSION, exclude = NULL)


#=================================
# Variable 'COD_EXCLUSION'
#=================================

#  1 Ventana Cuadrada
#  2 Tamaño Malla
#  3 Rejilla Metálica
#  4 Vertical cuadrado y Tamaño malla grande
#  5 Vertical cuadrado y Rejilla Metálica
#  6 Tamaño malla grande y Rejilla Metálica
#  7 Abertura Superior de Mallas
#  99 No observado
table(dat2$COD_EXCLUSION, exclude = NULL)


table(dat2$MODELO_EXCLUSION, dat2$COD_EXCLUSION, exclude = NULL)

#=======================================
# NUEVA VARIABLE QUE UNIFICA ESTAS DOS
#=======================================
# dat2$MOD_EX = ifelse(dat2$COD_EXCLUSION == 99, "9", dat2$MODELO_EXCLUSION)
# table(dat2$MOD_EX, exclude = NULL)


# Tengo ciertas dudas respecto a esta nueva variable ya que si comparamos
# esta nueva variable con 'MODELO_EXCLUSION', los observados = 1 disminuyen y 
# 2 = no observados tambien disminuyen, por tanto se creara una nueva categoria
# en 'MODELO_EXCLUSION' -> 4 que será 'S/I'.

dat2$MODELO_EXCLUSION = factor(dat2$MODELO_EXCLUSION, exclude = NULL, 
                                   levels = c("1", "2", "3", "9", NA), 
                                   labels = c("1", "2", "3", "9", "4"))
table(dat2$MODELO_EXCLUSION) # donde "4" es un nuevo nivel "sin información"

# dat2$MODELO_EXCLUSION = ordered(dat2$MODELO_EXCLUSION, levels = c("1", "2", "3", "4", "9"))
# table(dat2$MODELO_EXCLUSION)



#=============================================================================
# Variable 'TIPO_DE_RED'
#=============================================================================
table(dat2$TIPO_DE_RED, exclude = NULL) # Tenemos 38 NA

# Vemos a que especie objetivo corresponden esos lances
table(dat2$TIPO_DE_RED, dat2$ESPECIE_OBJETIVO_LANCE, exclude = NULL)
# corresponden a la especie 4 y asignamos el tipo de red capturado generalmente
# a esa especie (TIPO_RED = 2)

dat2$TIPO_DE_RED = factor(dat2$TIPO_DE_RED, exclude = NULL, 
                               levels = c("1", "2", "3", NA), 
                               labels = c("1", "2", "3", "2"))
table(dat2$TIPO_DE_RED, exclude = NULL) # donde "4" es un nuevo nivel "sin información"


# Imputación de datos para la PROFUNDIAD_PROM
dat2$PROFUNDIDAD_PROM[is.na(dat2$PROFUNDIDAD_PROM)] = mean(dat2$PROFUNDIDAD_PROM, na.rm=TRUE)
summary(dat2$PROFUNDIDAD_PROM)

table(dat2$PROFUNDIDAD_PROM)
hist(dat2$PROFUNDIDAD_PROM)



#=============================
# Nueva especie objetivo
#=============================

dat2$ESPECIE_OBJETIVO_LANCE2 = factor(dat2$ESPECIE_OBJETIVO_LANCE, exclude = NULL, 
                          levels = c("1", "2",  "3", "4", "5",  "6", "7", "27", "29", "45", "96"), 
                          labels = c("1", "2",  "3", "4", "1",  "1", "1", "27", "1", "1", "1"))
table(dat2$ESPECIE_OBJETIVO_LANCE2, exclude = NULL) # donde "4" es un nuevo nivel "sin información"




#================================
#        COD_PESQUERIA
#================================
library(plyr)
count(dat2, "COD_PESQUERIA")
count(dat2, c("ANO", "COD_PESQUERIA", "N_EJEMPLAR"))   # captura por buque

dat2 %>%
  group_by(COD_PESQUERIA) %>%
  summarise_at(vars(N_EJEMPLAR), funs(mean(., na.rm=TRUE)))


aggregate(N_EJEMPLAR ~ COD_PESQUERIA + ANO, data = dat2, FUN= "length" ) # N observaciones por año para cada pesquería
aggregate(N_EJEMPLAR ~ COD_PESQUERIA + ANO, data = dat2, FUN= "mean" ) # media de lobos capturados por año y pesquería
aggregate(N_EJEMPLAR ~ COD_PESQUERIA + ANO, data = dat2, FUN= "sum" ) # Suma de lobos capturados por año y pesquería




ddply(dat2, c("ANO", "COD_PESQUERIA"), summarise,
      media_n_ejemplar=mean(N_EJEMPLAR),
      max_n_ejemplar=max(N_EJEMPLAR),
      n.obs=length(N_EJEMPLAR))


ddply(dat2, c("ANO", "CLASE_LANCE"), summarise,
      media_n_ejemplar=mean(N_EJEMPLAR),
      max_n_ejemplar=max(N_EJEMPLAR),
      n.obs=length(N_EJEMPLAR))



#====================
#  ESPECIE OBJETIVO
#====================
count(dat2, "ESPECIE_OBJETIVO_LANCE")
count(dat2, c("ANO", "ESPECIE_OBJETIVO_LANCE", "N_EJEMPLAR"))   # captura por buque

aggregate(N_EJEMPLAR ~ ESPECIE_OBJETIVO_LANCE + ANO, data = dat2, FUN= "length" ) # N observaciones por año para cada pesquería
aggregate(N_EJEMPLAR ~ ESPECIE_OBJETIVO_LANCE + ANO, data = dat2, FUN= "mean" ) # media de lobos capturados por año y pesquería
aggregate(N_EJEMPLAR ~ ESPECIE_OBJETIVO_LANCE + ANO, data = dat2, FUN= "sum" ) # Suma de lobos capturados por año y pesquería


levels(dat2$COD_PESQUERIA) <- c("Fábrica", "Hielera")

ddply(dat2, c("ANO", "COD_PESQUERIA"), summarise,
      max_n_ejemplar=max(N_EJEMPLAR),
      n.obs=length(N_EJEMPLAR))

ddply(dat2, c("ANO", "CLASE_LANCE"), summarise,
      max_n_ejemplar=max(N_EJEMPLAR),
      n.obs=length(N_EJEMPLAR))


library(data.table)
setDT(dat2)[, list(Sum_N_ejemplar = sum(N_EJEMPLAR)), keyby = list(ANO, COD_PESQUERIA)][,Porcentaje := paste0(round(Sum_N_ejemplar/sum(Sum_N_ejemplar), 2)*100, "%"), by = ANO][]



# Calcula el porcentaje de lobos capturados para un año y por pesquería
setDT(dat2)[, list(Sum_N_ejemplar = sum(N_EJEMPLAR)), keyby = list(ANO, ESPECIE_OBJETIVO_LANCE)][,Porcentaje := paste0(round(Sum_N_ejemplar/sum(Sum_N_ejemplar), 2)*100, "%"), by = ANO][]


# Calcula el porcentaje de lobos capturados para un año y codigo de exclusión (REVISAR)
setDT(dat2)[, list(Sum_N_ejemplar = sum(N_EJEMPLAR)), keyby = list(ANO, MODELO_EXCLUSION)][,Porcentaje := paste0(round(Sum_N_ejemplar/sum(Sum_N_ejemplar), 2)*100, "%"), by = ANO][]




aggregate(N_EJEMPLAR ~ COD_PESQUERIA, data = dat2, FUN= "sum" )

x11()
library(scales)
ggplot(data = dat2, aes(x = N_EJEMPLAR2, fill = N_EJEMPLAR2)) +
  geom_bar(alpha = 1/2, aes(y = (..count..)/sum(..count..))) +
  scale_fill_manual(values = c("blue", "orangered2"))+
  geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust = -0.25) +
  scale_y_continuous(labels = percent) +
  labs(title = "Proporción de lances con captura incidental", y = "Porcentaje(%)", x = "") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=14, face = "bold"),
        axis.title=element_text(size=14,face="bold")) + 
  theme(plot.title = element_text(size = 16))

# Grafico FABRICA
dat_fabrica <- dat2 %>% filter(NOMBRE_PESQUERIA %in% "fabrica")

plot_fabrica = ggplot(data = dat_fabrica, aes(x = N_EJEMPLAR2, fill = N_EJEMPLAR2)) +
  geom_bar(alpha = 1/2, aes(y = (..count..)/sum(..count..))) +
  scale_fill_manual(values = c("blue", "orangered2"))+
  geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust = -0.25) +
  scale_y_continuous(labels = percent) +
  labs(title = "Flota FÁBRICA", y = "Porcentaje(%)", x = "") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=14, face = "bold"),
        axis.title=element_text(size=14,face="bold")) + 
  theme(plot.title = element_text(size = 17))


# Grafico HIELERA
dat_hielera <- dat2 %>% filter(NOMBRE_PESQUERIA %in% "hielera")

plot_hielera = ggplot(data = dat_hielera, aes(x = N_EJEMPLAR2, fill = N_EJEMPLAR2)) +
  geom_bar(alpha = 1/2, aes(y = (..count..)/sum(..count..))) +
  scale_fill_manual(values = c("blue", "orangered2"))+
  geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust = -0.25) +
  scale_y_continuous(labels = percent) +
  labs(title = "Flota HIELERA", y = "Porcentaje(%)", x = "") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=14, face = "bold"),
        axis.title=element_text(size=14,face="bold")) + 
  theme(plot.title = element_text(size = 17))


# Gráfico Fábrica y Hielera
grid.arrange(plot_fabrica, plot_hielera, ncol = 2)




# Gráfico Fábrica por año
ggplot(dat_fabrica, aes(x= N_EJEMPLAR2,  group=ANO)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.1, size = 3) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(title = "Proporción captura incidental FABRICA", y = "Porcentaje(%)", x = "")  +
  facet_wrap(~ANO) +
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=16,face="bold")) + 
        theme(strip.text = element_text(face="bold", size=12),
        strip.background = element_rect(fill="white", colour="grey",size=0.2))

# Gráfico Hielera por año
ggplot(dat_hielera, aes(x= N_EJEMPLAR2,  group=ANO)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.1, size = 3) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(title = "Proporción captura incidental HIELERA", y = "Porcentaje(%)", x = "")  +
  facet_wrap(~ANO) +
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=16,face="bold")) + 
  theme(strip.text = element_text(face="bold", size=12),
        strip.background = element_rect(fill="white", colour="grey",size=0.2))




#==================================
#     Estadísticos de resumen
#==================================
library(summarytools)
descr(dat2)


table(dat2$COD_PESQUERIA, dat2$CLASE_LANCE, dat2$N_EJEMPLAR2)





#====================================================================================================
# Grafico CLASE_LANCE
dat_LANCE1 <- dat2 %>% filter(CLASE_LANCE %in% "1")


# Gráfico 'CLASE_LANCE1'
ggplot(dat_LANCE1, aes(x= N_EJEMPLAR2,  group=ANO)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.1, size = 3) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(title = "Proporción captura incidental CLASE DE LANCE 1", y = "Porcentaje(%)", x = "")  +
  facet_wrap(~ANO) +
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=16,face="bold")) + 
  theme(strip.text = element_text(face="bold", size=12),
        strip.background = element_rect(fill="white", colour="grey",size=0.2))



# Gráfico 'CLASE_LANCE2'
dat_LANCE2 <- dat2 %>% filter(CLASE_LANCE %in% "2")

ggplot(dat_LANCE2, aes(x= N_EJEMPLAR2,  group=ANO)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.1, size = 3) +
  labs(title = "Proporción captura incidental CLASE DE LANCE 2", y = "Porcentaje(%)", x = "")  +
  facet_wrap(~ANO) +
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=12, face = "bold"),
        axis.title=element_text(size=16,face="bold")) + 
  theme(strip.text = element_text(face="bold", size=12),
        strip.background = element_rect(fill="white", colour="grey",size=0.2))
#====================================================================================================







# Test para ver correlación entre variables categoricas
cv.test = function(x,y) {
  CV = sqrt(chisq.test(x, y, correct=FALSE)$statistic /
              (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  print.noquote("Cramér V / Phi:")
  return(as.numeric(CV))
}

with(dat2, cv.test(COD_PESQUERIA, COD_BARCO))
with(dat2, cv.test(CLASE_LANCE, TIPO_DE_RED))
with(dat2, cv.test(CLASE_LANCE, MODELO_EXCLUSION))
with(dat2, cv.test(CLASE_LANCE, INTENSIDAD_VIENTO_AR))



# Cramer desde librería
library('rcompanion')
cramerV(dat2$CLASE_LANCE, dat2$INTENSIDAD_VIENTO_AR, bias.correct = FALSE)
cramerV(dat2$ESPECIE_OBJETIVO_LANCE, dat2$TIPO_DE_RED, bias.correct = FALSE)
cramerV(dat2$ESPECIE_OBJETIVO_LANCE, dat2$COD_BARCO, bias.correct = FALSE)
cramerV(dat2$TIPO_DE_RED, dat2$CLASE_LANCE, bias.correct = FALSE)
cramerV(dat2$TIPO_DE_RED, dat2$MODELO_EXCLUSION, bias.correct = FALSE)
cramerV(dat2$ESPECIE_OBJETIVO_LANCE, dat2$CLASE_LANCE, bias.correct = FALSE)

library(gmodels)
library(summarytools)
#CrossTable(dat2$NOMBRE_PESQUERIA, dat2$N_EJEMPLAR2, expected = TRUE, format = "SPSS")
summarytools::freq(dat2$NOMBRE_PESQUERIA, order = "freq")
summarytools::freq(dat2$TIPO_DE_RED, order = "freq")
summarytools::freq(dat2$CLASE_LANCE, report.nas = FALSE, headings = FALSE)
summarytools::freq(dat2$INTENSIDAD_VIENTO_AR, report.nas = FALSE, headings = FALSE)

# Rmarkdown
summarytools::freq(c(dat2$COD_PESQUERIA, dat2$N_EJEMPLAR), report.nas = FALSE, totals = FALSE,
                   cumul = FALSE, style = "rmarkdown", headings = FALSE)


summarytools::freq(dat2$N_EJEMPLAR, report.nas = FALSE, totals = FALSE,
                   cumul = FALSE, headings = FALSE)


library(tidyverse)
library(iml)
library(skimr)
library(DataExplorer)
library(ggpubr)
library(univariateML)
library(recipes)
library(kableExtra)
library(magrittr)



# Vamos a seleccionar las variables que creemos importantes para el análisis
dat2 = dat2[, c("N_EJEMPLAR", "ESTADO_C1", "ESTADO_F1", "ANO", "MES", "DIA", "TRIM", "COD_BARCO", "COD_PESQUERIA", 
                "TIPO_DE_RED", "CLASE_LANCE", "MODELO_EXCLUSION", "INTENSIDAD_VIENTO_AR", 
                "distancia_Km", "Nlobos_lobera", "lobera_cercana",
                "ESPECIE_OBJETIVO_LANCE2", "ESPECIE_OBJETIVO_LANCE", "PROFUNDIDAD_PROM", "LAT2", "LON2")]

kable(dat2, "latex", booktabs = T)


skim(dat2)

skim(dat2) %>%
  yank("numeric")

iris_list <- skimr::skim_to_list(dat2)

iris_list %>% .$factor %>% 
  knitr::kable( caption = 'Factor Variables', format = 'latex', booktabs = T ) %>% 
  kable_styling(latex_options = 'HOLD_position', font_size = 9)

iris_list  %>% .$numeric %>% 
  knitr::kable( caption = 'Numeric Variables', booktabs = T ) %>% 
  kable_styling(latex_options = 'HOLD_position', font_size = 9)


# Número de datos ausentes por variable
dat2 %>% map_dbl(.f = function(x){sum(is.na(x))})


plot_missing(
  data    = dat2, 
  #title   = "Porcentaje de valores ausentes",
  ggtheme = theme_bw(),
  theme_config = list(xlab = ",", axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.text=element_text(size=6),legend.position = "none") , 
)





plot_bar(
  dat2,
  ncol    = 4,
  title   = "Número de observaciones por grupo",
  ggtheme = theme_bw(),
  theme_config = list(
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(colour = "black", size = 12, face = 2),
    legend.position = "none"
  )
)




levels(dat2$COD_PESQUERIA) <- c("Fábrica", "Hielera")

p = ggplot(dat2, aes(x=N_EJEMPLAR))+ geom_histogram(color="black", fill="grey")+
  facet_grid(COD_PESQUERIA ~ .)
p  + theme_bw() 

p + scale_color_brewer(palette="black") + 
  theme_classic()+theme(legend.position="top") + 
  labs(title = "", y = "Frecuencia", x = "N° ejemplares capturados") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=9, face = "bold"),
        axis.title=element_text(size=9,face="bold")) + 
  theme(plot.title = element_text(size = 10))

#==========================================================
# Escribimos los datos que vamos a utilizar para modelar
write.csv(dat2 ,"dat_fit5.csv")
