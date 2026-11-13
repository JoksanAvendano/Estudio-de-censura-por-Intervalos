# ====================================================================
# EFECTO CENSURA POR INTERVALOS (CODA + BAYESIANO)
# ====================================================================

library(coda)
library(doParallel)
library(foreach)

set.seed(1234)
M <- 5000          # Número de simulaciones
n <- 20            
theta_true <- 20.4 
N_mcmc <- 10000    # Muestras MCMC para estabilizar coda

# Corrección asintótica de Díaz-Francés para n=20
c_95 <- 0.1465 - (1 / (32.0361 + 37.0825 * n))
corte_log <- log(c_95)
deltas <- c("Semanal" = 0.25, "Quincenal" = 0.5, "Mensual" = 1)

# ====================================================================
# 1. DEFINICIÓN DE FUNCIONES (Deben cargarse antes de paralelizar)
# ====================================================================

obtener_ic_exacto <- function(tiempos_obs, theta_emv) {
  log_lik <- function(theta) {
    if (is.na(theta) || theta <= 0) return(-Inf)
    sum(-log(theta) - tiempos_obs/theta)
  }
  max_ll <- log_lik(theta_emv)
  if(is.infinite(max_ll) || is.na(max_ll)) return(c(NA, NA))
  
  f_raiz <- function(theta) log_lik(theta) - (max_ll + corte_log)
  inf <- tryCatch(uniroot(f_raiz, lower = 0.1, upper = theta_emv)$root, error = function(e) NA)
  sup <- tryCatch(uniroot(f_raiz, lower = theta_emv, upper = 500)$root, error = function(e) NA)
  return(c(inf, sup))
}

obtener_ic_discreto <- function(delta, tiempos_obs) {
  ll_func <- function(theta) {
    if (is.na(theta) || theta <= 0) return(-Inf)
    k_ev <- ceiling(tiempos_obs / delta)
    L_i <- (k_ev - 1) * delta
    R_i <- k_ev * delta
    prob_ev <- exp(-L_i/theta) - exp(-R_i/theta)
    prob_ev[prob_ev <= 0] <- 1e-10
    return(sum(log(prob_ev)))
  }
  
  T_medio <- (ceiling(tiempos_obs / delta) - 0.5) * delta
  opt <- tryCatch(
    optim(par = mean(T_medio), fn = ll_func, method = "Brent", lower=0.1, upper=500, control=list(fnscale=-1)),
    error = function(e) NULL
  )
  if (is.null(opt)) return(c(NA, NA))
  
  max_ll <- opt$value
  f_raiz <- function(theta) ll_func(theta) - (max_ll + corte_log)
  inf <- tryCatch(uniroot(f_raiz, lower = 0.1, upper = opt$par)$root, error = function(e) NA)
  sup <- tryCatch(uniroot(f_raiz, lower = opt$par, upper = 500)$root, error = function(e) NA)
  return(c(inf, sup))
}

obtener_hdi_exacto_coda <- function(tiempos_obs, prior="jeffreys") {
  n_obs <- length(tiempos_obs)
  S <- sum(tiempos_obs)
  alpha <- ifelse(prior == "jeffreys", n_obs, n_obs + 2.1)
  beta_p <- ifelse(prior == "jeffreys", S, S + 22)
  
  muestras_post <- 1 / rgamma(N_mcmc, shape = alpha, rate = beta_p)
  return(as.numeric(HPDinterval(as.mcmc(muestras_post), prob = 0.95)))
}

obtener_hdi_discreto_ambos_coda <- function(delta, tiempos_obs) {
  theta_grid <- seq(1, 150, length.out = 4000) 
  k_ev <- ceiling(tiempos_obs / delta)
  prob_ev <- outer((k_ev - 1) * delta, theta_grid, function(l, th) exp(-l/th)) - 
    outer(k_ev * delta, theta_grid, function(r, th) exp(-r/th))
  prob_ev[prob_ev <= 0] <- 1e-10
  log_lik <- colSums(log(prob_ev))
  
  extraer_coda <- function(log_prior) {
    log_post <- log_lik + log_prior
    post <- exp(log_post - max(log_post))
    post <- post / sum(post)
    muestras_post <- sample(theta_grid, size = N_mcmc, replace = TRUE, prob = post)
    return(as.numeric(HPDinterval(as.mcmc(muestras_post), prob = 0.95)))
  }
  
  return(list(
    jeffreys = extraer_coda( -log(theta_grid) ),
    ig = extraer_coda( -3.1 * log(theta_grid) - 22/theta_grid )
  ))
}

# ====================================================================
# 2. CONFIGURACIÓN DE PARALELISMO
# ====================================================================
# Detectar cuántos núcleos tiene tu CPU y dejar 1 libre para que la PC no se congele
num_cores <- parallel::detectCores() - 1 
cl <- makeCluster(num_cores)
registerDoParallel(cl)

cat("Iniciando Monte Carlo en", num_cores, "núcleos. Por favor espera...\n")

# ====================================================================
# 3. BUCLE MONTE CARLO PARALELIZADO
# ====================================================================
# Utilizamos foreach en lugar del bucle FOR clásico. 
# Los resultados de cada iteración se "apilarán" (rbind) en un data.frame masivo.

resultados_crudos <- foreach(i = 1:M, .combine = rbind, .packages = c("coda")) %dopar% {
  
  # Generar muestra para esta iteración
  tiempos_obs <- rexp(n, rate = 1/theta_true)
  theta_emv <- mean(tiempos_obs)
  
  # Función auxiliar rápida para dar formato a los resultados de esta iteración
  evaluar <- function(nombre, ic) {
    if (any(is.na(ic))) return(NULL)
    data.frame(
      Metodo = nombre,
      Len = ic[2] - ic[1],
      Low = ifelse(theta_true < ic[1], 1, 0),
      Up  = ifelse(theta_true > ic[2], 1, 0),
      Dentro = ifelse(theta_true >= ic[1] && theta_true <= ic[2], 1, 0),
      IC_inf = ic[1],
      IC_sup = ic[2]
    )
  }
  
  res_iter <- list() # Guardar resultados de la iteración actual
  
  # A. CASO EXACTO
  res_iter[[1]] <- evaluar("Frec_Exacto", obtener_ic_exacto(tiempos_obs, theta_emv))
  res_iter[[2]] <- evaluar("BayesJeff_Exacto", obtener_hdi_exacto_coda(tiempos_obs, "jeffreys"))
  res_iter[[3]] <- evaluar("BayesIG_Exacto", obtener_hdi_exacto_coda(tiempos_obs, "ig"))
  
  # B. CASO POR INTERVALOS
  idx <- 4
  for (d_name in names(deltas)) {
    delta <- deltas[[d_name]]
    res_iter[[idx]] <- evaluar(paste0("Frec_", d_name), obtener_ic_discreto(delta, tiempos_obs))
    
    hpd_disc <- obtener_hdi_discreto_ambos_coda(delta, tiempos_obs)
    res_iter[[idx+1]] <- evaluar(paste0("BayesJeff_", d_name), hpd_disc$jeffreys)
    res_iter[[idx+2]] <- evaluar(paste0("BayesIG_", d_name), hpd_disc$ig)
    idx <- idx + 3
  }
  
  # Devolver todos los métodos evaluados en esta simulación
  do.call(rbind, res_iter)
}

# Apagar el cluster de núcleos al terminar
stopCluster(cl)

# ====================================================================
# 4. AGREGACIÓN Y REPORTE FINAL
# ====================================================================
# Ahora resumimos el data.frame masivo sacando promedios por Método

resultados_finales <- aggregate(
  cbind(Len, Low, Dentro, Up, IC_inf, IC_sup) ~ Metodo, 
  data = resultados_crudos, 
  FUN = mean
)

# Dar formato de presentación (Porcentajes e Intervalos legibles)
resultados_finales$Abajo <- resultados_finales$Low * 100
resultados_finales$Dentro <- resultados_finales$Dentro * 100
resultados_finales$Arriba <- resultados_finales$Up * 100
resultados_finales$Intervalo_Promedio <- paste0("[", round(resultados_finales$IC_inf, 2), ", ", round(resultados_finales$IC_sup, 2), "]")

# Limpiar columnas no necesarias para la impresión
resultados_finales <- resultados_finales[, c("Metodo", "Len", "Abajo", "Dentro", "Arriba", "Intervalo_Promedio")]
colnames(resultados_finales)[2] <- "Longitud"

print(resultados_finales)
