# ====================================================================
# Simulación de cesnura por la derecha con K censurados
# ====================================================================

set.seed(1234)
# Valores para la simulacion
M <- 10000         
n <- 20            
theta_true <- 20.4 
k_valores <- c(0, 1, 2, 3, 4, 5, 6, 7) # Ajustado a la realidad de la distribución

c_95 <- 0.1465 - (1 / (32.0361 + 37.0825 * n))
corte_log <- log(c_95)

obtener_ic_perfil <- function(status, tiempos_obs, theta_emv) {
  log_lik <- function(theta) {
    if (is.na(theta) || theta <= 0) return(-Inf) # Protección contra NAs
    sum(status * (-log(theta) - tiempos_obs/theta) + (1 - status) * (-tiempos_obs/theta))
  }
  max_ll <- log_lik(theta_emv)
  
  # Falla controlada si la verosimilitud se rompe
  if(is.infinite(max_ll) || is.na(max_ll)) return(c(NA, NA))
  
  f_raiz <- function(theta) log_lik(theta) - (max_ll + corte_log)
  
  inf <- tryCatch(uniroot(f_raiz, lower = 0.1, upper = theta_emv)$root, error = function(e) NA)
  sup <- tryCatch(uniroot(f_raiz, lower = theta_emv, upper = 500)$root, error = function(e) NA)
  return(c(inf, sup))
}

# --------------------------------------------------------------------
#  Cálculo Analítico Exacto del HPD para Inversa-Gamma
# --------------------------------------------------------------------
obtener_hpd_ig <- function(alpha, beta_p) {
  # Igualamos las alturas de la densidad en el límite Inferior (L) y Superior (U)
  f_obj <- function(L) {
    prob_L <- pgamma(1/L, shape = alpha, rate = beta_p, lower.tail = FALSE)
    prob_U <- prob_L + 0.95
    if(prob_U >= 0.999999) return(Inf)
    
    U <- 1 / qgamma(prob_U, shape = alpha, rate = beta_p, lower.tail = FALSE)
    
    log_f_L <- -(alpha + 1)*log(L) - beta_p/L
    log_f_U <- -(alpha + 1)*log(U) - beta_p/U
    
    return(log_f_L - log_f_U)
  }
  
  # Límite máximo donde podríamos empezar a buscar L
  L_max <- 1 / qgamma(0.05, shape = alpha, rate = beta_p, lower.tail = FALSE)
  
  L_opt <- tryCatch(uniroot(f_obj, lower = 1e-6, upper = L_max * 0.999)$root, error = function(e) NA)
  if(is.na(L_opt)) return(c(NA, NA))
  
  prob_L <- pgamma(1/L_opt, shape = alpha, rate = beta_p, lower.tail = FALSE)
  U_opt <- 1 / qgamma(prob_L + 0.95, shape = alpha, rate = beta_p, lower.tail = FALSE)
  
  return(c(L_opt, U_opt))
}

resultados_finales <- data.frame()

for (k in k_valores) {
  stats <- list(f    = list(dentro=0, low=0, up=0, len=0, ic=c(0,0)),
                b_J  = list(dentro=0, low=0, up=0, len=0, ic=c(0,0)),
                b_IG = list(dentro=0, low=0, up=0, len=0, ic=c(0,0)))
  
  for (i in 1:M) {
    T_r <- rexp(n, rate = 1/theta_true)
    status <- rep(1, n)
    tiempos_obs <- T_r 
    
    candidatos <- which(T_r > theta_true)
    k_efectivo <- min(k, length(candidatos))
    
    if (k_efectivo > 0) {
      # --- CORRECCIÓN DEL BUG DE R EN LA FUNCIÓN SAMPLE ---
      if (length(candidatos) == 1) {
        idx_cens <- candidatos
      } else {
        idx_cens <- sample(candidatos, k_efectivo)
      }
      
      status[idx_cens] <- 0
      # Simulación uniforme rigurosa por la derecha
      tiempos_obs[idx_cens] <- runif(k_efectivo, min = theta_true, max = T_r[idx_cens])
    }
    
    r <- sum(status)
    S_t <- sum(tiempos_obs) 
    
    if (r == 0) next # Evitar división por cero si ocurre un evento extremo
    
    theta_emv <- S_t/r
    
    # 1. Frecuentista
    ic_p <- obtener_ic_perfil(status, tiempos_obs, theta_emv)
    if (!any(is.na(ic_p))) {
      stats$f$len <- stats$f$len + (ic_p[2] - ic_p[1]); stats$f$ic <- stats$f$ic + ic_p
      if (theta_true < ic_p[1]) stats$f$low <- stats$f$low + 1 
      else if (theta_true > ic_p[2]) stats$f$up <- stats$f$up + 1 
      else stats$f$dentro <- stats$f$dentro + 1
    }
    
    # 2. Bayesiano Jeffreys ( usando HPD exacto)
    # Posterior: Inversa-Gamma(alpha = r, beta = S_t)
    hpd_J <- obtener_hpd_ig(alpha = r, beta_p = S_t)
    if (!any(is.na(hpd_J))) {
      stats$b_J$len <- stats$b_J$len + (hpd_J[2] - hpd_J[1]); stats$b_J$ic <- stats$b_J$ic + hpd_J
      if (theta_true < hpd_J[1]) stats$b_J$low <- stats$b_J$low + 1 
      else if (theta_true > hpd_J[2]) stats$b_J$up <- stats$b_J$up + 1 
      else stats$b_J$dentro <- stats$b_J$dentro + 1
    }
    
    # 3. Bayesiano IG (usando HPD exacto)
    # Posterior: Inversa-Gamma(alpha = 2.1 + r, beta = 22 + S_t)
    hpd_IG <- obtener_hpd_ig(alpha = 2.1 + r, beta_p = 22 + S_t)
    if (!any(is.na(hpd_IG))) {
      stats$b_IG$len <- stats$b_IG$len + (hpd_IG[2] - hpd_IG[1]); stats$b_IG$ic <- stats$b_IG$ic + hpd_IG
      if (theta_true < hpd_IG[1]) stats$b_IG$low <- stats$b_IG$low + 1 
      else if (theta_true > hpd_IG[2]) stats$b_IG$up <- stats$b_IG$up + 1 
      else stats$b_IG$dentro <- stats$b_IG$dentro + 1
    }
  }
  
  for (met in names(stats)) {
    resultados_finales <- rbind(resultados_finales, data.frame(
      Censura_k = k, Metodo = met,
      Longitud = stats[[met]]$len/M,
      Abajo = (stats[[met]]$low/M)*100,
      Dentro = (stats[[met]]$dentro/M)*100,
      Arriba = (stats[[met]]$up/M)*100,
      Intervalo_Promedio = paste0("[", round(stats[[met]]$ic[1]/M, 2), ", ", round(stats[[met]]$ic[2]/M, 2), "]")
    ))
  }
}

print(resultados_finales)