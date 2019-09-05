
# Paquetes ----------------------------------------------------------------


require(MortalityLaws)
require(lifecontingencies)
library(ggplot2)
library(xtable)
library(gamlss)


# Lectura de datos --------------------------------------------------------

D = read.table(file = "Tabla_vida_hombre.txt",header=T)
D$qx.1000 <- D$qx.1000/1000
colnames(D) <- c("Edad", "qx", "lx", "dx", "ex")
plot(D$Edad,D$qx, xlab = "Edad", ylab = expression(q[x]))
head(D)

mx = D$qx/(1-D$qx/2)


# Función generadora de vidas aleatorias ----------------------------------

una.muerte <- function(Edad){
  e <- c(Edad:100)
  lx <- D$lx[Edad-19]
  tqx = 1-(D$lx[(Edad-19):length(D$lx)]/lx)
  Fx <- tqx/tqx[length(tqx)]
  u <- runif(1)
  i <- 1
  while (u > Fx[i]) {
    i <- i + 1
  }
  return(e[i])
}


# Revisión de la función --------------------------------------------------

muerte <- function(n, Edad){
  replicate(n = n, una.muerte(Edad = Edad))
}
muerte(10, c(30, 25, 40))
eval <- muerte(10000, 30)
hist(eval)

# Evaluar valor presente del pago -----------------------------------------

vpSeguro <- function(x,k,i){
  vp <- 1
  for (j in 1:(k+1-x)) {
    vp <- vp/(1+i[j])
  }
  return(vp)
}


# Calculo prima neta  -----------------------------------------------------

prima.neta <- function(x, tpx.ley, i, parametros){
  pn <- 0
  v <- c()
  temp <- 1
  for (j in 1:(81-x)) {
    temp <- temp/(1+i[j])
    v <- c(v, temp)
  }
  for (k in 1:(81-x)) {
    pn <- pn + v[k]*tpx.ley((k-1),x,parametros)*(1-tpx.ley(1, (k-1)+x, parametros))
  }
  return(pn)
}

# Tasa de interés ---------------------------------------------------------

historico <- read.table("Tint2.txt")
porcent <- historico$V1/100
fit <- fitDist(porcent, type = "real0to1")
fit$fits
pari <- c(fit$mu, fit$sigma, fit$nu, fit$tau)
simi <- rGB1(80, mu = pari[1], sigma = pari [2], nu = pari[3], tau = pari[4])
hist(simi)
nombres.par <- c("mu", "sigma", "nu", "tau")
tabla.i <- data.frame(Parametro = nombres.par, Valor_estimado = pari)
xtable(tabla.i)

# FunciónECM --------------------------------------------------------------

one.ECM <- function(x, tpx, parametros, tipo.descuento = "Constante"){
  
  # Vector de tasas de descuento
  ifelse(tipo.descuento=="Constante",{
    i <- rep(0.0373287, (81))
  },
  {
    i <- rGB1((81),  mu = pari[1], sigma = pari [2], nu = pari[3], tau = pari[4])
  })
  
  # Generación de un tiempo de vida aleatorio
  k <- una.muerte(Edad = x)
  
  # Hallar Prima para este tiempo de vida
  vp <- vpSeguro(x = x, k = k, i = i)
  
  # Hallar la prima neta
  A <- prima.neta(x, tpx, i, parametros)
 
  # Error cuadrado
  E <- (vp-A)^2
  
  return(E)
}


# Ajuste Leyes de Mortalidad ----------------------------------------------


# Ajuste de Ley de Gompertz -----------------------------------------------


Mgompertz <- MortalityLaw(x   =  D$Edad, 
                       mx  = mx, 
                       law = 'gompertz', 
                       opt.method = 'LF2')


parametrosg =  Mgompertz$coefficients

# tpx Gompertz

tpx.gompertz <- function(t, x, parametros){
  fn <- function(s){gompertz(s+x, parametros)$hx}
  a = ifelse(t == 0, 0, integrate(fn, 0, t)$value)
  v = exp(-a)
  return(v)
}


# Ajuste de Ley de Makeham ------------------------------------------------


Mmakeham <- MortalityLaw(x = D$Edad, mx = mx, law = "makeham", opt.method = "LF2")
parametrosm <- Mmakeham$coefficients

# tpx Makeham

tpx.makeham <- function(t, x, parametros){
  fn <- function(s){makeham(s+x, parametros)$hx}
  a = ifelse(t == 0, 0, integrate(fn, 0, t)$value)
  v = exp(-a)
  return(v)
}


# Ajuste Ley de Siller ----------------------------------------------------


Msiler <- MortalityLaw(x   =  D$Edad, 
                       mx  = mx, 
                       law = 'siler', 
                       opt.method = 'LF2')


parametross =  Msiler$coefficients

#tpx con Siler

tpx.siler = function(t,x,parametros){
  fn = function(s){siler(s+x,parametross)$hx}
  a = ifelse(t==0,0,integrate(fn,0,t)$value)
  v = exp(-a)
  return(v)
}


# Ajuste de Ley de Heligman and Pollard ----------------------------------

MHP <- MortalityLaw(x = D$Edad, mx = mx, law = "HP", opt.method = "LF2")
parametroshp <- MHP$coefficients


#tpx con HP

tpx.HP = function(t,x,parametros){
  fn = function(s){HP(s+x,parametros)$hx}
  a = ifelse(t==0,0,integrate(fn,0,t)$value)
  v = exp(-a)
  return(v)
}


# Revisión de las funciones -----------------------------------------------

tpx.gompertz(50, 30, parametrosg)
tpx.makeham(50,30, parametrosm)
tpx.siler(50, 30, parametros)
tpx.HP(50, 30, parametroshp)


# Función Para generar los datos simulados -------------------------------

simulacion <- function(Nsim, x, Ley.Mortalidad, Metodo.optimizacion, tipo.descuento = "Constante"){
  
  # Ajustando la ley de mortalidad
  Ley <- MortalityLaw(x = D$Edad, mx = mx, law = Ley.Mortalidad, opt.method = Metodo.optimizacion)
  
  # Obteniendo el vector de parámetros
  parametros <- Ley$coefficients
  
  # Asignando función para integrar
  ifelse(Ley.Mortalidad == "gompertz",{
    FUN <- function(y, pars) gompertz(y, pars)
  }, ifelse(Ley.Mortalidad == "makeham", {
    FUN <- function(y, pars) makeham(y, pars)
  }, ifelse(Ley.Mortalidad == "siler",{
    FUN <- function(y, pars) siler(y, pars)
  },{
    FUN <- function(y, pars) HP(y, pars)
  })))
  
  # Creando la función tpx
  tpx = function(t,x,parametros){
    fn = function(s){FUN(s+x,parametros)$hx}
    a = ifelse(t==0,0,integrate(fn,0,t)$value)
    v = exp(-a)
    return(v)
  }
  
  # Simular los ECM
  ECM <- replicate(n = Nsim, one.ECM(x = x, tpx = tpx, parametros = parametros, tipo.descuento = tipo.descuento))
  resultados <- cbind(mean(ECM),sd(ECM)/sqrt(Nsim), x, Ley.Mortalidad, Metodo.optimizacion, tipo.descuento)
  write(x = t(resultados), file = "resultados.txt", ncolumns = 6, append = T)
  
}


# Asignando valores de los parametros -------------------------------------

Nsim <- 10000
xs <- c(30:49)
Leyes <- c("gompertz", "makeham", "siler")
Metodo <- c("poissonL", "LF1", "LF2","LF3", "LF4", "LF5", "LF6")
tasa <- c("Constante", "Variable")

grid1 <- expand.grid(xs=xs, Leyes=Leyes, Metodos=Metodo, Tasa=tasa)
grid2 <- expand.grid(xs=xs, Leyes="HP", Metodos="LF2", Tasa=tasa)
params <- rbind(grid1, grid2)


# Simulación --------------------------------------------------------------

set.seed(12345)

lapply(1:NROW(params), function(i){
  cat(i," ")
  simulacion(Nsim=Nsim, x = params[i,1], Ley.Mortalidad = as.character(params[i,2]), 
             Metodo.optimizacion = as.character(params[i,3]), tipo.descuento = as.character(params[i,4]))
})


# Lectura de datos --------------------------------------------------------

datos <- read.table(file = "resultados.txt", header = F)
colnames(datos) <- c("ECM", "Desv.est", "Edad", "Ley.Mortalidad", "Metodo.optimizacion", "Tipo.interes")

max(datos$Desv.est)
max(datos$ECM)
datos


# Separacion por valores de parametros ------------------------------------

sep.tasa <- split(datos, datos$Tipo.interes)
t.constante <- sep.tasa$Constante
t.variable <- sep.tasa$Variable

sep.ley1 <- split(t.constante, t.constante$Ley.Mortalidad)
sep.ley2 <- split(t.variable, t.variable$Ley.Mortalidad)

gompertz.c <- sep.ley1$gompertz
makeham.c <- sep.ley1$makeham
siler.c <- sep.ley1$siler
hp.c <- sep.ley1$HP

gompertz.v <- sep.ley2$gompertz
makeham.v <- sep.ley2$makeham
siler.v <- sep.ley2$siler
hp.v <- sep.ley2$HP



# Función para gráficar todos los métodos ---------------------------------

graf <- function(datos){
  par(mfrow=c(3,3))
  layout(matrix(c(1,2,3,4,5,6,0,7,0), ncol = 3, byrow = T))
  lapply(split(datos, datos$Metodo.optimizacion),
         function(x){
           with(data = x, {
             plot(x = Edad, y = ECM, type = "l", lty = 2, col = "blue", 
                  main = paste(Metodo.optimizacion[1], "con tasa", Tipo.interes[1]), 
                  ylab = "ECM", 
                  xlab = "Edad",
                  ylim = c(0, 0.15), las=1)
             media <- mean(ECM)
             abline(h = media, col = "black")
             text((max(Edad)-2),media,paste("Prom = ",round(media,digits=4)),pos=3,font=2,cex=0.9)
             lines(x = Edad, y = Desv.est, col = "red", lty = 2)
           })
         })
  par(mfrow=c(1,1))
}


# Gráficos Gompertz --------------------------------------------------------

graf(gompertz.c)
graf(gompertz.v)


# Gráficos makeham --------------------------------------------------------

graf(makeham.c)
graf(makeham.v)


# Gráficos siler ----------------------------------------------------------

graf(siler.c)
graf(siler.v)

# Gráficos para ley HP ----------------------------------------------------

par(mfrow=c(1,2))

# Tasa constante
with(data = hp.c,{
  plot(x = Edad, y = ECM, type = "l", lty = 2, col = "blue", las = 1, 
       main = "Tasa constante", 
       ylab = "ECM", 
       xlab = "Edad",
       ylim = c(0, 0.1))
  media <- mean(ECM)
  abline(h = media, col = "black")
  text((max(Edad)-2),media+0.01,paste("Prom = ",round(media,digits=4)),pos=3,font=2,cex=0.9)
  lines(x = Edad, y = Desv.est, col = "red", lty = 2)
})

# Tasa variable
with(data = hp.v,{
  plot(x = Edad, y = ECM, type = "l", lty = 2, col = "blue", las = 1,
       main = "Tasa variable", 
       ylab = "ECM", 
       xlab = "Edad",
       ylim = c(0, 0.1))
  media <- mean(ECM)
  abline(h = media, col = "black")
  text((max(Edad)-2),media+0.01,paste("Prom = ",round(media,digits=4)),pos=3,font=2,cex=0.9)
  lines(x = Edad, y = Desv.est, col = "red", lty = 2)
})

par(mfrow=c(1,1))



# Mejores ajustes ---------------------------------------------------------

# Ley de Gompertz

mgc <- split(gompertz.c, gompertz.c$Metodo.optimizacion)$LF4
mgv <- split(gompertz.v, gompertz.v$Metodo.optimizacion)$LF4


# Ley de Makeham

mmc <- split(makeham.c, makeham.c$Metodo.optimizacion)$LF6
mmv <- split(makeham.v, makeham.v$Metodo.optimizacion)$LF6


# Ley de Makeham

msc <- split(siler.c, siler.c$Metodo.optimizacion)$LF2
msv <- split(siler.v, siler.v$Metodo.optimizacion)$LF2


# Graficando para tasa constante ------------------------------------------


# ECM

par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)

plot(x = hp.c$Edad, y = hp.c$ECM, type = "l", lty = 2, col = "blue",
       ylab = "ECM", 
       xlab = "Edad",
       ylim = c(0, 0.08), las = 1)
points(x = hp.c$Edad, y = hp.c$ECM, col = "blue", pch = 5, cex = 0.5)
lines(x = mgc$Edad, y = mgc$ECM, lty = 3, col = "red")
points(x = mgc$Edad, y = mgc$ECM, col = "red", pch = 4, cex = 0.5)
lines(x = mmc$Edad, y = mmc$ECM, lty = 4, col = "darkgreen")
points(x = mmc$Edad, y = mmc$ECM, lty = 4, col = "darkgreen", pch = 3, cex = 0.5)
lines(x = msc$Edad, y = msc$ECM, lty = 5, col = "black")
points(x = msc$Edad, y = msc$ECM, lty = 5, col = "black", pch = 6, cex = 0.5)
legend("topright", inset = c(-0.16,0), 
       legend = c("ley de \nGompertz", "ley de \nMakeham", "ley de \nSiler", "ley de \nHeligman- \nPollard"), 
       col = c("red", "darkgreen", "black", "blue"), 
       lty = c(3,4,5,2),
       pch = c(5, 4, 3, 6),
       ncol = 1, cex = 0.7, text.width = 1, text.font = 4, 
       bg = "transparent", box.lwd = 1, box.lty = 1)
par(mfrow= c(1,1))

#Desv.estandar

par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
plot(x = hp.c$Edad, y = hp.c$Desv.est, type = "l", lty = 2, col = "blue", 
     ylab = "Desviación Estandar del estimador", 
     xlab = "Edad",
     ylim = c(0, 0.005), las = 1)
points(x = hp.c$Edad, y = hp.c$Desv.est, col = "blue", pch = 5, cex = 0.5)
lines(x = mgc$Edad, y = mgc$Desv.est, lty = 3, col = "red")
points(x = mgc$Edad, y = mgc$Desv.est, col = "red", pch = 4, cex = 0.5)
lines(x = mmc$Edad, y = mmc$Desv.est, lty = 4, col = "darkgreen")
points(x = mmc$Edad, y = mmc$Desv.est, lty = 4, col = "darkgreen", pch = 3, cex = 0.5)
lines(x = msc$Edad, y = msc$Desv.est, lty = 5, col = "black")
points(x = msc$Edad, y = msc$Desv.est, lty = 5, col = "black", pch = 6, cex = 0.5)
legend("topright", inset = c(-0.16,0), 
       legend = c("ley de \nGompertz", "ley de \nMakeham", "ley de \nSiler", "ley de \nHeligman- \nPollard"), 
       col = c("red", "darkgreen", "black", "blue"), 
       lty = c(3,4,5,2),
       pch = c(5, 4, 3, 6),
       ncol = 1, cex = 0.7, text.width = 1, text.font = 4, 
       bg = "transparent", box.lwd = 1, box.lty = 1)


# Graficando para tasa variable ------------------------------------------


# ECM

plot(x = hp.v$Edad, y = hp.v$ECM, type = "l", lty = 2, col = "blue", 
     ylab = "ECM", 
     xlab = "Edad",
     ylim = c(0, 0.08), las = 1)
points(x = hp.v$Edad, y = hp.v$ECM, col = "blue", pch = 5, cex = 0.5)
lines(x = mgv$Edad, y = mgv$ECM, lty = 3, col = "red")
points(x = mgv$Edad, y = mgv$ECM, col = "red", pch = 4, cex = 0.5)
lines(x = mmv$Edad, y = mmv$ECM, lty = 4, col = "darkgreen")
points(x = mmv$Edad, y = mmv$ECM, lty = 4, col = "darkgreen", pch = 3, cex = 0.5)
lines(x = msv$Edad, y = msv$ECM, lty = 5, col = "black")
points(x = msv$Edad, y = msv$ECM, lty = 5, col = "black", pch = 6, cex = 0.5)
legend("topright", inset = c(-0.16,0), 
       legend = c("ley de \nGompertz", "ley de \nMakeham", "ley de \nSiler", "ley de \nHeligman- \nPollard"), 
       col = c("red", "darkgreen", "black", "blue"), 
       lty = c(3,4,5,2),
       pch = c(5, 4, 3, 6),
       ncol = 1, cex = 0.7, text.width = 1, text.font = 4, 
       bg = "transparent", box.lwd = 1, box.lty = 1)


#Desv.estandar


plot(x = hp.v$Edad, y = hp.v$Desv.est, type = "l", lty = 2, col = "blue", 
     ylab = "Desviación Estandar del estimador", 
     xlab = "Edad",
     ylim = c(0, 0.005), las = 1)
points(x = hp.v$Edad, y = hp.v$Desv.est, col = "blue", pch = 5, cex = 0.5)
lines(x = mgv$Edad, y = mgv$Desv.est, lty = 3, col = "red")
points(x = mgv$Edad, y = mgv$Desv.est, col = "red", pch = 4, cex = 0.5)
lines(x = mmv$Edad, y = mmv$Desv.est, lty = 4, col = "darkgreen")
points(x = mmv$Edad, y = mmv$Desv.est, lty = 4, col = "darkgreen", pch = 3, cex = 0.5)
lines(x = msv$Edad, y = msv$Desv.est, lty = 5, col = "black")
points(x = msv$Edad, y = msv$Desv.est, lty = 5, col = "black", pch = 6, cex = 0.5)
legend("topright", inset = c(-0.16,0), 
       legend = c("ley de \nGompertz", "ley de \nMakeham", "ley de \nSiler", "ley de \nHeligman- \nPollard"), 
       col = c("red", "darkgreen", "black", "blue"), 
       lty = c(3,4,5,2),
       pch = c(5, 4, 3, 6),
       ncol = 1, cex = 0.7, text.width = 1, text.font = 4, 
       bg = "transparent", box.lwd = 1, box.lty = 1)
