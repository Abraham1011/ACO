library(ggplot2)


tou_inicio <- function(n_ciudades,evaporation_rate){
  tau <- (1/(1 - evaporation_rate)) * 1/1e6
  M <- matrix(tau,
              ncol = n_ciudades,nrow = n_ciudades)
  return(M)
}

prob_ij <- function(feromona, distancia, alpha, beta){
  prob <- feromona^alpha * (1/distancia)^beta
  prob[prob==Inf] <- 1e-40
  prob <- prob/sum(prob)
  return(prob)
}

objfun <- function(ruta, d) {
  n_ciudades <- length(ruta)
  fitness <- sum(d[cbind(ruta[-n_ciudades], ruta[-1])])
  fitness <- fitness + d[ruta[n_ciudades], ruta[1]]
  
  return(fitness)
}

ACO <- function(n_ants, alpha, beta, Q, evaporation_rate,
                d, max_iteraciones, pbest){
  n_ciudades <- dim(d)[1]
  tau <- tou_inicio(n_ciudades, evaporation_rate)
  best_rutas <- matrix(nrow = max_iteraciones, 
                       ncol = n_ciudades)
  pdec <- pbest^(1/n_ciudades)
  avg <- n_ciudades/2
  for (nit in 1:max_iteraciones) {
    rutas <- matrix(nrow = n_ants, ncol = n_ciudades)
    for(n_ant in 1:n_ants){
      ruta <- c(1)
      for(i in 2:n_ciudades){
        pij <- prob_ij(feromona = tau[ruta[i-1],],
                       distancia = d[ruta[i-1],],
                       alpha = alpha, beta = beta)
        pij[ruta[1:(i-1)]] <- 0
        ruta[i] <- sample(1:n_ciudades, 1, replace = FALSE, prob = pij)
      }
      rutas[n_ant,] <- ruta
    }
    
    tau <- tau * (1 - evaporation_rate)
    
    fitness <- apply(rutas, 1, function(r){objfun(ruta = r,d = d)})
    best_ant <- which.min(fitness)
    best_rutas[nit,] <- rutas[best_ant,]
    
    ruta <- rutas[best_ant,]
    for(i in 1:(length(ruta)-1)){
      Lk <- objfun(ruta = ruta,d = d)
      tau[ruta[i], ruta[i+1]] <- tau[ruta[i], ruta[i+1]] + 1 / Lk
      tau[ruta[i+1], ruta[i]] <- tau[ruta[i], ruta[i+1]]  
    }
    tmax <- 1/(1 - evaporation_rate) * 1/min(fitness)
    tmin <- (tmax * (1 - pdec))/((avg - 1)*pdec)
    if(tmin > tmax){
      tmin = tmax
    }
    tau[tau < tmin] <- tmin
    tau[tau > tmax] <- tmax
  }
  fitness <- apply(best_rutas, 1, function(r){objfun(ruta = r,d = d)})
  best_ruta <- best_rutas[which.min(fitness),]
  return(list(best = best_ruta, rutas = best_rutas,
              fitness = min(fitness)))
}


#Ejemplo
set.seed(1011)
cord <- data.frame(lon = rnorm(100,0,1),lat = rnorm(100,0,1))
d <- as.matrix(dist(cord))
plot(cord)


res <- ACO(n_ants = 100,alpha = 2, beta = 5, Q = 1,
           evaporation_rate = 0.9,d = d, 
           max_iteraciones = 200, pbest = 0.01)

df <- data.frame(
  x = cord$lon[res$best],
  y = cord$lat[res$best]
)
ggplot(data = df) +
  geom_path(aes(x = x, y = y), color = "#25ced1",
            size = 1) +
  geom_point(aes(x = x, y = y), color = "#ff8a5b", size = 3) +
  labs(title = "Mejor ruta encontrada",
       x = "", y = "") + 
  theme_void( 
    base_size = 10, 
    base_line_size =2, 
    base_rect_size =2)+ 
  theme(plot.title=element_text(hjust=0.5,size=20)) 


res$fitness

#4.651719
#8.006542

rutas_hist <- res$rutas
df <- data.frame()
for(i in 1:dim(rutas_hist)[1]){
  df_aux <- data.frame(ciclo = rep(i,dim(rutas_hist)[2]),
                       lon = cord$lon[rutas_hist[i,]],
                       lat = cord$lat[rutas_hist[i,]])
  df <- rbind(df,df_aux)
}

df$ciclo <- factor(df$ciclo)
library(gganimate)
p <- ggplot(data = df) +
  geom_path(aes(x = lon, y = lat, group = ciclo),
            color = "#25ced1", size = 1) +  
  geom_point(aes(x = lon, y = lat), 
             color = "#ff8a5b", size = 3) +  
  theme_void( 
    base_size = 10, 
    base_line_size =2, 
    base_rect_size =2)+ 
  theme(plot.title = element_text(hjust=0.5,size=20)) +  
  labs(title = "Mejor ruta encontrada",
       x = "", y = "") + 
  transition_manual(ciclo) +  
  labs(title = "Mejor ruta en el ciclo: {current_frame}") 

animate(p,fps = 100/20,renderer = magick_renderer(),
        height = 7, width = 7, units = "in",
        res = 100)



# Guardar el archivo como gif
#anim_save("ruta_animada.gif")


