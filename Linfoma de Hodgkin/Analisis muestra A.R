# ====================================================================
# 1. CARGA DE DATOS (Grupo A - Hodgkin's Disease)
# ====================================================================

# Cargar la librería 'coda' para el análisis MCMC
if(!require(coda)) {
  install.packages("coda")
  library(coda)
}
tiempos <- c(1.25, 1.41, 4.98, 5.25, 5.38, 6.92, 8.89, 10.98, 11.18, 13.11, 
             13.21, 16.33, 19.77, 21.08, 21.84, 22.07, 31.38, 32.62, 37.18, 42.92)
status <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1)


r <- sum(status)       # Número de fallas observadas (Tamaño de muestra efectivo)
sum_t <- sum(tiempos)  # Suma total de tiempos en riesgo
n <-length(tiempos)

# ====================================================================
# 2. DEFINICIÓN DE FUNCIONES DE LOG-VEROSIMILITUD (Escala: theta = media)
# ====================================================================

# Log-verosimilitud para el Caso Exacto
log_lik_exacta <- function(theta) {
  if (theta <= 0) return(-Inf)
  return(-r * log(theta) - (sum_t / theta))
}

# Log-verosimilitud para Casos Discretos (Censura por Intervalo)
log_lik_discreta <- function(theta, delta) {
  if (theta <= 0) return(-Inf)
  
  L <- floor(tiempos / delta) * delta
  R <- L + delta
  
  # Fallas
  prob_falla <- exp(-L[status == 1] / theta) - exp(-R[status == 1] / theta)
  # Censuras
  prob_censura <- exp(-L[status == 0] / theta)
  # No usar 0 
  prob_falla[prob_falla <= 0] <- 1e-10
  prob_censura[prob_censura <= 0] <- 1e-10
  
  return(sum(log(prob_falla)) + sum(log(prob_censura)))
}

# ====================================================================
# 3. ESTIMACIÓN UNIFICADA Y CÁLCULO DE INTERVALOS DE CONFIANZA
# ====================================================================

# Nivel de corte 'c' para el 95% de confianza
# (Basado en la aproximación chi-cuadrada con corrección por tamaño de muestra)
c_95 <- 0.1465 - (1 / (32.0361 + 37.0825 * n))
corte_log <- log(c_95) 

# Función genérica: Calcula los límites del IC mediante la búsqueda de raíces
# de la log-verosimilitud perfilada.
obtener_ic_unificado <- function(fun_loglik, theta_emv, lower = 5, upper = 60) {
  max_ll <- fun_loglik(theta_emv)
  
  # La función a la que le buscaremos el cero
  f_raiz <- function(theta) fun_loglik(theta) - (max_ll + corte_log)
  
  # Buscar raíces a los lados del EMV
  limite_inf <- uniroot(f_raiz, lower = lower, upper = theta_emv)$root
  limite_sup <- uniroot(f_raiz, lower = theta_emv, upper = upper)$root
  
  return(c(limite_inf, limite_sup))
}

# --- CASO 1: EXACTO (Con censura) ---
opt_ex <- optim(par = 20, fn = log_lik_exacta, method = "Brent", 
                lower = 5, upper = 50, control = list(fnscale = -1))
theta_ex <- opt_ex$par
ic_ex <- obtener_ic_unificado(log_lik_exacta, theta_ex)

# --- CASO 2: SEMANAL ---
# Usamos delta_sem calculado previamente
opt_sema <- optim(par = 20, fn = function(t) log_lik_discreta(t, delta_sem), 
                  method = "Brent", lower = 5, upper = 50, control = list(fnscale = -1))
theta_sema <- opt_sema$par
ic_sema <- obtener_ic_unificado(function(t) log_lik_discreta(t, delta_sem), theta_sema)

# --- CASO 3: QUINCENAL ---
opt_quin <- optim(par = 20, fn = function(t) log_lik_discreta(t, 0.5), 
                  method = "Brent", lower = 5, upper = 50, control = list(fnscale = -1))
theta_quin <- opt_quin$par
ic_quin <- obtener_ic_unificado(function(t) log_lik_discreta(t, 0.5), theta_quin)

# --- CASO 4: MENSUAL ---
opt_men <- optim(par = 20, fn = function(t) log_lik_discreta(t, 1), 
                 method = "Brent", lower = 5, upper = 50, control = list(fnscale = -1))
theta_men <- opt_men$par
ic_men <- obtener_ic_unificado(function(t) log_lik_discreta(t, 1), theta_men)

# ====================================================================
# 4. IMPRESIÓN DE RESULTADOS
# ====================================================================
cat("========================================================\n")
cat("   ESTIMACIÓN DE THETA (TIEMPO MEDIO DE VIDA EN MESES)  \n")
cat("========================================================\n")
cat(sprintf("%-12s | %-10s | %-20s | %-10s\n", "Escala", "EMV", "IC 95%", "Distancia"))
cat("--------------------------------------------------------\n")
cat(sprintf("%-12s | %-10.4f | [%.4f, %.4f] | %-10.4f\n", "Exacto", theta_ex, ic_ex[1], ic_ex[2], ic_ex[2]-ic_ex[1]))
cat(sprintf("%-12s | %-10.4f | [%.4f, %.4f] | %-10.4f\n", "Semanal", theta_sema, ic_sema[1], ic_sema[2], ic_sema[2]-ic_sema[1]))
cat(sprintf("%-12s | %-10.4f | [%.4f, %.4f] | %-10.4f\n", "Quincenal", theta_quin, ic_quin[1], ic_quin[2], ic_quin[2]-ic_quin[1]))
cat(sprintf("%-12s | %-10.4f | [%.4f, %.4f] | %-10.4f\n", "Mensual", theta_men, ic_men[1], ic_men[2], ic_men[2]-ic_men[1]))
cat("========================================================\n")
# ====================================================================
# 5. GRÁFICA COMPARATIVA DE LOG-VEROSIMILITUD
# ====================================================================
theta_grid <- seq(10, 40, length.out = 500)

y_ex <- sapply(theta_grid, log_lik_exacta)
y_qu <- sapply(theta_grid, function(x) log_lik_discreta(x, 0.5))
y_se <- sapply(theta_grid, function(x) log_lik_discreta(x, delta_sem))

y_ex_norm <- y_ex - max(y_ex)
y_qu_norm <- y_qu - max(y_qu)
y_se_norm <- y_se - max(y_se)

plot(theta_grid, y_ex_norm, type = "l", lwd = 2, col = "black",
     main = expression(paste("Efecto de la Discretización ", theta)),
     xlab = expression(theta ~ " (Tiempo Medio en Meses)"), 
     ylab = "Log-Verosimilitud Relativa",
     ylim = c(-5, 0)) 

lines(theta_grid, y_se_norm, col = "blue", lwd = 2, lty = 2)
lines(theta_grid, y_qu_norm, col = "red", lwd = 2, lty = 3)

abline(v = theta_ex, col = "black", lty = 4)
abline(v = theta_sema, col = "blue", lty = 4)
abline(v = theta_quin, col = "red", lty = 4)

legend("topright", 
       legend = c("Exacta (Continuo)", "Semanal (0.23 meses)", "Quincenal (0.5 meses)"),
       col = c("black", "blue", "red"), lty = c(1, 2, 3), lwd = 2, cex = 0.7)

# ====================================================================
# 5. B GRÁFICA COMPARATIVA DE VEROSIMILITUD (L en lugar de Log-L)
# ====================================================================
# Crear un grid de valores de theta (de 5 a 40 meses)
theta_grid <- seq(5, 45, length.out = 500)
delta_sem <- 7 / 30.41

# 1. Evaluar las funciones LOG-verosimilitud en el grid
# Asegúrate de tener cargada la función log_lik_hibrida
log_y_ex <- sapply(theta_grid, log_lik_exacta)
log_y_se <- sapply(theta_grid, function(x) log_lik_discreta(x, delta_sem))
log_y_qu <- sapply(theta_grid, function(x) log_lik_discreta(x, 0.5))
log_y_me <- sapply(theta_grid, function(x) log_lik_discreta(x, 1))

# 2. Normalizar la LOG-verosimilitud (Restar el máximo)
# El pico de estas nuevas variables será exactamente 0
log_y_ex_norm <- log_y_ex - max(log_y_ex)
log_y_se_norm <- log_y_se - max(log_y_se)
log_y_qu_norm <- log_y_qu - max(log_y_qu)
log_y_me_norm <- log_y_me - max(log_y_me)

# 3. Transformar a VEROSIMILITUD RELATIVA (Aplicando exponencial)
# L(theta) = exp( log_L(theta) )
L_ex <- exp(log_y_ex_norm)
L_se <- exp(log_y_se_norm)
L_qu <- exp(log_y_qu_norm)
L_me <- exp(log_y_me_norm)

# 4. Graficar la Verosimilitud
plot(theta_grid, L_ex, type = "l", lwd = 2, col = "black",
     main = expression(paste("Función de Verosimilitud Relativa para ", theta)),
     xlab = expression(theta ~ " (Tiempo  en Meses)"), 
     ylab = expression("Verosimilitud Relativa  " * L(theta)),
     ylim = c(0, 1)) # El eje Y ahora va de 0 a 1

# Añadir las curvas discretas (Híbridas)
lines(theta_grid, L_se, col = "blue", lwd = 2, lty = 2)
lines(theta_grid, L_qu, col = "red", lwd = 2, lty = 3)
lines(theta_grid, L_me, col = "purple", lwd = 2, lty = 3)

# Opcional: Encontrar las modas en el grid para trazar las líneas verticales
theta_hat_ex <- theta_grid[which.max(L_ex)]
theta_hat_se <- theta_grid[which.max(L_se)]
theta_hat_qu <- theta_grid[which.max(L_qu)]

abline(v = theta_hat_ex, col = "black", lty = 4)
abline(v = theta_hat_se, col = "blue", lty = 4)
abline(v = theta_hat_qu, col = "red", lty = 4)

# Añadir leyenda
legend("topright", 
       legend = c("Exacta (Continuo)", "Semanal ", "Quincenal ", "Mensual"),
       col = c("black", "blue", "red", "purple"), lty = c(1, 2, 3), lwd = 2, cex = 0.7)


# ====================================================================
# 6. DEFINICIÓN DE PREVIAS Y ALGORITMO MCMC
# ====================================================================

# Datos a priori segundo caso
r2 <- 16
sum_t2 <-  sum(c(1.25, 1.41, 4.98, 5.25, 5.38, 6.92, 8.89, 10.98, 11.18, 13.11, 
                 13.21, 16.33, 19.77, 21.08, 21.84, 22.07, 31.38, 32.62, 37.18, 42.92))

# Función unificada de la Previa. Permite devolver la densidad (para graficar) o el logaritmo (para MCMC)
calc_prior <- function(theta, tipo_previa = "invgamma_debil", log_val = TRUE) {
  if (theta <= 0) return(ifelse(log_val, -Inf, 0))
  
  if (tipo_previa == "invgamma_debil") {
    alpha <- 2.1
    beta <- 22
    # Se añade la constante de normalización real de la Inv-Gamma para poder graficarla
    log_dens <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(theta) - (beta / theta)
    
  } else if (tipo_previa == "jeffreys") {
    # Al ser impropia, no integra a 1. Su escala gráfica no será comparable con densidades propias.
    log_dens <- -log(theta)
    
  } else if (tipo_previa == "invgamma_muestra") {
    alpha <- 2.1 + r2
    beta <- 22 + sum_t2
    log_dens <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(theta) - (beta / theta)
    
  } else {
    stop("Tipo de previa no reconocido")
  }
  
  return(ifelse(log_val, log_dens, exp(log_dens)))
}

# Vectorizamos la función para que curve() pueda evaluarla correctamente en la gráfica
calc_prior_vec <- Vectorize(calc_prior, vectorize.args = "theta")

# El MCMC ahora recibe el tipo de previa como argumento
run_mcmc <- function(iter, delta, tipo_verosimilitud = "exacto", tipo_previa = "invgamma_debil") {
  theta_chain <- numeric(iter)
  theta_chain[1] <- 20  
  sd_prop <- 2.0        
  
  for (i in 2:iter) {
    theta_curr <- theta_chain[i-1]
    theta_prop <- rnorm(1, mean = theta_curr, sd = sd_prop)
    
    if (theta_prop <= 0) {
      log_alpha <- -Inf
    } else {
      if (tipo_verosimilitud == "exacto") {
        log_lik_curr <- log_lik_exacta(theta_curr)
        log_lik_prop <- log_lik_exacta(theta_prop)
      } else {
        log_lik_curr <- log_lik_discreta(theta_curr, delta)
        log_lik_prop <- log_lik_discreta(theta_prop, delta)
      }
      
      # Log-Posterior = Log-Likelihood + Log-Prior según la selección
      log_post_curr <- log_lik_curr + calc_prior(theta_curr, tipo_previa, log_val = TRUE)
      log_post_prop <- log_lik_prop + calc_prior(theta_prop, tipo_previa, log_val = TRUE)
      
      log_alpha <- log_post_prop - log_post_curr
    }
    
    if (log(runif(1)) < log_alpha) theta_chain[i] <- theta_prop
    else theta_chain[i] <- theta_curr
  }
  return(theta_chain)
}

# ====================================================================
# 7. EJECUCIÓN DEL ANÁLISIS
# ====================================================================
set.seed(123) 
iteraciones <- 50000
burn_in <- 10000 

# SELECCIONA AQUÍ LA PREVIA QUE DESEAS USAR: 
# Opciones: "invgamma_debil", "jeffreys", "invgamma_muestra"
previa_elegida <- "jeffreys" 

cadena_ex <- run_mcmc(iteraciones, delta = 0, tipo_verosimilitud = "exacto", tipo_previa = previa_elegida)
cadena_se <- run_mcmc(iteraciones, delta = 7/30.41, tipo_verosimilitud = "discreto", tipo_previa = previa_elegida)
cadena_qu <- run_mcmc(iteraciones, delta = 0.5, tipo_verosimilitud = "discreto", tipo_previa = previa_elegida)
cadena_me <- run_mcmc(iteraciones, delta = 1, tipo_verosimilitud = "discreto", tipo_previa = previa_elegida)

post_ex <- cadena_ex[(burn_in + 1):iteraciones]
post_se <- cadena_se[(burn_in + 1):iteraciones]
post_qu <- cadena_qu[(burn_in + 1):iteraciones]
post_me <- cadena_me[(burn_in + 1):iteraciones]

# Función para obtener la Moda (MAP) basada en la densidad de Kernel
get_mode <- function(chain) {
  d <- density(chain)
  return(d$x[which.max(d$y)])
}

hpd_ex <- HPDinterval(as.mcmc(post_ex), prob = 0.95)
hpd_se <- HPDinterval(as.mcmc(post_se), prob = 0.95)
hpd_qu <- HPDinterval(as.mcmc(post_qu), prob = 0.95)
hpd_me <- HPDinterval(as.mcmc(post_me), prob = 0.95)

# ====================================================================
# 8. RESULTADOS: MEDIA, MODA (MAP) E INTERVALOS HPD AL 95%
# ====================================================================
res_bayes <- data.frame(
  Escenario = c("Exacto", "Semanal", "Quincenal","Mensual"),
  Media = c(mean(post_ex), mean(post_se), mean(post_qu), mean(post_me)),
  Moda_MAP = c(get_mode(post_ex), get_mode(post_se), get_mode(post_qu), get_mode(post_me)),
  HPD_Inf = c(hpd_ex[1], hpd_se[1], hpd_qu[1], hpd_me[1]),
  HPD_Sup = c(hpd_ex[2], hpd_se[2], hpd_qu[2],hpd_me[2]),
  Distancia = c(hpd_ex[2]-hpd_ex[1], hpd_se[2]-hpd_se[1], hpd_qu[2]-hpd_qu[1], hpd_me[2]-hpd_me[1])
)

cat("\n========================================================================\n")
cat("          ESTIMADORES BAYESIANOS (TIEMPO MEDIO EN MESES)                \n")
cat("========================================================================\n")
cat(sprintf("%-12s | %-12s | %-12s | %-20s | %-12s\n", "Escenario", "Media", "Moda (MAP)", "HPD 95%", "Distancia"))
cat("------------------------------------------------------------------------\n")
for(i in 1:4) {
  cat(sprintf("%-12s | %-12.4f | %-12.4f | [%.4f, %.4f] | %-12.4f\n", 
              res_bayes$Escenario[i], res_bayes$Media[i], res_bayes$Moda_MAP[i], 
              res_bayes$HPD_Inf[i], res_bayes$HPD_Sup[i], res_bayes$Distancia[i]))
}
cat("========================================================================\n")

# ====================================================================
# 9. GRÁFICA COMPARATIVA DE DENSIDADES POSTERIORES Y PREVIA
# ====================================================================
# Ajustamos ylim dinámicamente para asegurar que las densidades o previas más picudas no se recorten
max_y <- max(density(post_ex)$y, density(post_se)$y, density(post_qu)$y, density(post_me)$y)

if(previa_elegida != "jeffreys") {
  max_prior <- max(calc_prior_vec(seq(1, 45, length.out=500), tipo_previa = previa_elegida, log_val = FALSE))
  max_y <- max(max_y, max_prior)
}

plot(density(post_ex), lwd=2, col="black", 
     main=expression(paste("Densidades Posteriores y Previa de ", theta)),
     xlab=expression(theta~" (meses)"), ylab="Densidad",
     xlim=c(0, 45), ylim=c(0, max_y * 1.1))

lines(density(post_se), lwd=2, col="blue", lty=2)
lines(density(post_qu), lwd=2, col="red", lty=3)
lines(density(post_me), lwd=2, col="purple", lty=4)

# Lógica para graficar la previa elegida
if (previa_elegida != "jeffreys") {
  curve(calc_prior_vec(x, tipo_previa = previa_elegida, log_val = FALSE), 
        add = TRUE, col = "darkgreen", lwd = 2, lty = 1)
  
  legend("topright", 
         legend=c("Exacta", "Semanal", "Quincenal", "Mensual", "Previa"),
         col=c("black", "blue", "red", "purple", "darkgreen"), 
         lty=c(1, 2, 3, 4, 1), lwd=2, cex = 0.7)
} else {
  # Para la previa de Jeffreys (1/x), usamos un factor de escala arbitrario 
  # para que su FORMA sea visible en la misma ventana de las densidades propias.
  # Multiplicamos 1/x por un porcentaje del valor máximo de Y.
  factor_escala <- max_y * 1 
  
  curve((1/x) * factor_escala, add = TRUE, col = "darkgreen", lwd = 2, lty = 1)
  
  legend("topright", 
         legend=c("Exacta", "Semanal", "Quincenal", "Mensual", "Previa Jeffreys"),
         col=c("black", "blue", "red", "purple", "darkgreen"), 
         lty=c(1, 2, 3, 4, 1), lwd=2, cex = 0.7)
}

