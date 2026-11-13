# ====================================================================
# SCRIPT MAESTRO: CENSURA SELECTIVA (k) + CENSURA POR INTERVALOS
# COMPARATIVA COMPLETA: FRECUENTISTA VS BAYESIANO 
# ====================================================================

set.seed(1234)
M <- 5000          # Iteraciones
n <- 20            
theta_true <- 20.4 
k_valores <- c(0, 1, 3, 5, 7) 

# Corrección asintótica de Díaz-Francés
c_95 <- 0.1465 - (1 / (32.0361 + 37.0825 * n))
corte_log <- log(c_95)

# --------------------------------------------------------------------
# 1. FUNCIÓN PERFILADA: DATOS EXACTOS
# --------------------------------------------------------------------
obtener_ic_exacto <- function(status, tiempos_obs, theta_emv) {
  log_lik <- function(theta) {
    if (is.na(theta) || theta <= 0) return(-Inf)
    sum(status * (-log(theta) - tiempos_obs/theta) + (1 - status) * (-tiempos_obs/theta))
  }
  max_ll <- log_lik(theta_emv)
  if(is.infinite(max_ll) || is.na(max_ll)) return(c(NA, NA))
  
  f_raiz <- function(theta) log_lik(theta) - (max_ll + corte_log)
  inf <- tryCatch(uniroot(f_raiz, lower = 0.1, upper = theta_emv)$root, error = function(e) NA)
  sup <- tryCatch(uniroot(f_raiz, lower = theta_emv, upper = 500)$root, error = function(e) NA)
  return(c(inf, sup))
}

# --------------------------------------------------------------------
# 2. FUNCIÓN PERFILADA: DATOS DISCRETOS (Frecuentista)
# --------------------------------------------------------------------
obtener_ic_discreto <- function(delta, status, tiempos_obs) {
  ll_func <- function(theta) {
    if (is.na(theta) || theta <= 0) return(-Inf)
    
    # Fallas: P( T en (L, R] )
    k_ev <- ceiling(tiempos_obs[status == 1] / delta)
    L_i <- (k_ev - 1) * delta
    R_i <- k_ev * delta
    prob_ev <- exp(-L_i/theta) - exp(-R_i/theta)
    prob_ev[prob_ev <= 0] <- 1e-10 
    ll_ev <- sum(log(prob_ev))
    
    # Censuras: Sobrevive a la inspección C_disc
    C_disc <- floor(tiempos_obs[status == 0] / delta) * delta
    ll_cens <- sum(-C_disc / theta)
    
    return(ll_ev + ll_cens)
  }
  
  emv_inicial <- sum(tiempos_obs)/max(1, sum(status))
  opt <- tryCatch(
    optim(par = emv_inicial, fn = ll_func, method = "Brent", lower = 0.1, upper = 500, control = list(fnscale = -1)),
    error = function(e) NULL
  )
  
  if (is.null(opt)) return(c(NA, NA))
  theta_emv_disc <- opt$par
  max_ll <- opt$value
  
  f_raiz <- function(theta) ll_func(theta) - (max_ll + corte_log)
  inf <- tryCatch(uniroot(f_raiz, lower = 0.1, upper = theta_emv_disc)$root, error = function(e) NA)
  sup <- tryCatch(uniroot(f_raiz, lower = theta_emv_disc, upper = 500)$root, error = function(e) NA)
  return(c(inf, sup))
}

# --------------------------------------------------------------------
# 3. FUNCIONES BAYESIANAS: HPD CONTINUO Y HPD DISCRETO EXACTO
# --------------------------------------------------------------------

# Para el caso de Datos Exactos (Continuos): La conjugación sí se mantiene perfectamente.
obtener_hpd_ig <- function(alpha, beta_p) {
  f_obj <- function(L) {
    prob_L <- pgamma(1/L, shape = alpha, rate = beta_p, lower.tail = FALSE)
    prob_U <- prob_L + 0.95
    if(prob_U >= 0.999999) return(Inf)
    U <- 1 / qgamma(prob_U, shape = alpha, rate = beta_p, lower.tail = FALSE)
    log_f_L <- -(alpha + 1)*log(L) - beta_p/L
    log_f_U <- -(alpha + 1)*log(U) - beta_p/U
    return(log_f_L - log_f_U)
  }
  L_max <- 1 / qgamma(0.05, shape = alpha, rate = beta_p, lower.tail = FALSE)
  L_opt <- tryCatch(uniroot(f_obj, lower = 1e-6, upper = L_max * 0.999)$root, error = function(e) NA)
  if(is.na(L_opt)) return(c(NA, NA))
  
  prob_L <- pgamma(1/L_opt, shape = alpha, rate = beta_p, lower.tail = FALSE)
  U_opt <- 1 / qgamma(prob_L + 0.95, shape = alpha, rate = beta_p, lower.tail = FALSE)
  return(c(L_opt, U_opt))
}

# Para el caso Discreto: Verosimilitud exacta evaluada matricialmente sobre grilla
obtener_hpd_discreto_exacto <- function(delta, status, tiempos_obs) {
  # Definir una grilla adaptativa para soportar casos extremos de censura (k=7)
  emv_aprox <- sum(tiempos_obs) / max(1, sum(status))
  theta_grid <- seq(1, max(200, emv_aprox * 5), length.out = 5000)
  
  # 1. Log-Verosimilitud de las Fallas (Intervalos)
  idx_fail <- which(status == 1)
  if (length(idx_fail) > 0) {
    k_ev <- ceiling(tiempos_obs[idx_fail] / delta)
    L_i <- (k_ev - 1) * delta
    R_i <- k_ev * delta
    P_fail <- outer(L_i, theta_grid, function(l, th) exp(-l/th)) - 
      outer(R_i, theta_grid, function(r, th) exp(-r/th))
    P_fail[P_fail <= 0] <- 1e-10
    ll_fail <- colSums(log(P_fail))
  } else {
    ll_fail <- rep(0, length(theta_grid))
  }
  
  # 2. Log-Verosimilitud de las Censuras (Por la derecha en tiempo discreto)
  idx_cens <- which(status == 0)
  if (length(idx_cens) > 0) {
    C_disc <- floor(tiempos_obs[idx_cens] / delta) * delta
    ll_cens <- colSums(outer(C_disc, theta_grid, function(c, th) -c/th))
  } else {
    ll_cens <- rep(0, length(theta_grid))
  }
  
  log_lik <- ll_fail + ll_cens
  
  # Función interna para calcular HPD dado el prior
  extraer_hpd_grilla <- function(log_prior) {
    log_post <- log_lik + log_prior
    post <- exp(log_post - max(log_post)) # Estabilizador numérico
    post <- post / sum(post)
    
    ord <- order(post, decreasing = TRUE)
    cum_post <- cumsum(post[ord])
    idx_95 <- which(cum_post >= 0.95)[1]
    intervalo <- theta_grid[ord[1:idx_95]]
    
    return(c(min(intervalo), max(intervalo)))
  }
  
  return(list(
    jeffreys = extraer_hpd_grilla( -log(theta_grid) ),
    ig       = extraer_hpd_grilla( -3.1 * log(theta_grid) - 22/theta_grid )
  ))
}


# --------------------------------------------------------------------
# 4. EJECUCIÓN MONTE CARLO
# --------------------------------------------------------------------
resultados_finales <- data.frame()
deltas <- c("Semanal" = 0.25, "Quincenal" = 0.5, "Mensual" = 1)

for (k in k_valores) {
  stats <- list(f_Exacto = list(dentro=0, low=0, up=0, len=0, ic=c(0,0)),
                b_Jeff_Exacto = list(dentro=0, low=0, up=0, len=0, ic=c(0,0)),
                b_IG_Exacto = list(dentro=0, low=0, up=0, len=0, ic=c(0,0)))
  
  for (d_name in names(deltas)) {
    stats[[paste0("f_", d_name)]] <- list(dentro=0, low=0, up=0, len=0, ic=c(0,0))
    stats[[paste0("b_Jeff_", d_name)]] <- list(dentro=0, low=0, up=0, len=0, ic=c(0,0))
    stats[[paste0("b_IG_", d_name)]] <- list(dentro=0, low=0, up=0, len=0, ic=c(0,0))
  }
  
  for (i in 1:M) {
    T_r <- rexp(n, rate = 1/theta_true)
    status <- rep(1, n)
    tiempos_obs <- T_r 
    
    # Aplicación de Censura por la Derecha Selectiva
    candidatos <- which(T_r > theta_true)
    k_efectivo <- min(k, length(candidatos))
    
    if (k_efectivo > 0) {
      idx_cens <- if(length(candidatos) == 1) candidatos else sample(candidatos, k_efectivo)
      status[idx_cens] <- 0
      tiempos_obs[idx_cens] <- runif(k_efectivo, min = theta_true, max = T_r[idx_cens])
    }
    
    r <- sum(status); S_t <- sum(tiempos_obs) 
    if (r == 0) next 
    theta_emv <- S_t/r
    
    # --- 1. EVALUACIONES EXACTAS ---
    # Frecuentista
    ic_ex <- obtener_ic_exacto(status, tiempos_obs, theta_emv)
    if (!any(is.na(ic_ex))) {
      stats$f_Exacto$len <- stats$f_Exacto$len + (ic_ex[2] - ic_ex[1]); stats$f_Exacto$ic <- stats$f_Exacto$ic + ic_ex
      if (theta_true < ic_ex[1]) stats$f_Exacto$low <- stats$f_Exacto$low + 1 else if (theta_true > ic_ex[2]) stats$f_Exacto$up <- stats$f_Exacto$up + 1 else stats$f_Exacto$dentro <- stats$f_Exacto$dentro + 1
    }
    
    # Bayesiano Exacto (HPD Analítico)
    hpd_J_ex <- obtener_hpd_ig(alpha = r, beta_p = S_t)
    if (!any(is.na(hpd_J_ex))) {
      stats$b_Jeff_Exacto$len <- stats$b_Jeff_Exacto$len + (hpd_J_ex[2] - hpd_J_ex[1]); stats$b_Jeff_Exacto$ic <- stats$b_Jeff_Exacto$ic + hpd_J_ex
      if (theta_true < hpd_J_ex[1]) stats$b_Jeff_Exacto$low <- stats$b_Jeff_Exacto$low + 1 else if (theta_true > hpd_J_ex[2]) stats$b_Jeff_Exacto$up <- stats$b_Jeff_Exacto$up + 1 else stats$b_Jeff_Exacto$dentro <- stats$b_Jeff_Exacto$dentro + 1
    }
    
    hpd_IG_ex <- obtener_hpd_ig(alpha = 2.1+r, beta_p = 22+S_t)
    if (!any(is.na(hpd_IG_ex))) {
      stats$b_IG_Exacto$len <- stats$b_IG_Exacto$len + (hpd_IG_ex[2] - hpd_IG_ex[1]); stats$b_IG_Exacto$ic <- stats$b_IG_Exacto$ic + hpd_IG_ex
      if (theta_true < hpd_IG_ex[1]) stats$b_IG_Exacto$low <- stats$b_IG_Exacto$low + 1 else if (theta_true > hpd_IG_ex[2]) stats$b_IG_Exacto$up <- stats$b_IG_Exacto$up + 1 else stats$b_IG_Exacto$dentro <- stats$b_IG_Exacto$dentro + 1
    }
    
    # --- 2. EVALUACIONES POR INTERVALOS (Bucle sobre Deltas) ---
    for (d_name in names(deltas)) {
      delta <- deltas[[d_name]]
      
      # Frecuentista Discreto (Perfilado Exacto vía optim)
      ic_d <- obtener_ic_discreto(delta, status, tiempos_obs)
      if (!any(is.na(ic_d))) {
        fname <- paste0("f_", d_name)
        stats[[fname]]$len <- stats[[fname]]$len + (ic_d[2] - ic_d[1]); stats[[fname]]$ic <- stats[[fname]]$ic + ic_d
        if (theta_true < ic_d[1]) stats[[fname]]$low <- stats[[fname]]$low + 1 else if (theta_true > ic_d[2]) stats[[fname]]$up <- stats[[fname]]$up + 1 else stats[[fname]]$dentro <- stats[[fname]]$dentro + 1
      }
      
      # Bayesiano Discreto (HPD de Verosimilitud Discreta Exacta)
      hpd_d <- obtener_hpd_discreto_exacto(delta, status, tiempos_obs)
      
      # Jeffreys
      if (!any(is.na(hpd_d$jeffreys))) {
        bjname <- paste0("b_Jeff_", d_name)
        stats[[bjname]]$len <- stats[[bjname]]$len + (hpd_d$jeffreys[2] - hpd_d$jeffreys[1]); stats[[bjname]]$ic <- stats[[bjname]]$ic + hpd_d$jeffreys
        if (theta_true < hpd_d$jeffreys[1]) stats[[bjname]]$low <- stats[[bjname]]$low + 1 else if (theta_true > hpd_d$jeffreys[2]) stats[[bjname]]$up <- stats[[bjname]]$up + 1 else stats[[bjname]]$dentro <- stats[[bjname]]$dentro + 1
      }
      
      # Inversa-Gamma
      if (!any(is.na(hpd_d$ig))) {
        bigname <- paste0("b_IG_", d_name)
        stats[[bigname]]$len <- stats[[bigname]]$len + (hpd_d$ig[2] - hpd_d$ig[1]); stats[[bigname]]$ic <- stats[[bigname]]$ic + hpd_d$ig
        if (theta_true < hpd_d$ig[1]) stats[[bigname]]$low <- stats[[bigname]]$low + 1 else if (theta_true > hpd_d$ig[2]) stats[[bigname]]$up <- stats[[bigname]]$up + 1 else stats[[bigname]]$dentro <- stats[[bigname]]$dentro + 1
      }
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
  cat(sprintf("Simulación completada para k = %d\n", k))
}

print(resultados_finales)
