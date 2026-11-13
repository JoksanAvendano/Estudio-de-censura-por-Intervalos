#Criterio de Información de Akaike (AIC) para ver el mejor modelo.


# Cargar la librería estándar de supervivencia en R
if(!require(survival)) {install.packages("survival"); library(survival)}

# Datos de la muestra de Hodgkin
tiempos <- c(1.25, 1.41, 4.98, 5.25, 5.38, 6.92, 8.89, 10.98, 11.18, 13.11, 
             13.21, 16.33, 19.77, 21.08, 21.84, 22.07, 31.38, 32.62, 37.18, 42.92)
status <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1)
# Crear el objeto de supervivencia
surv_obj <- Surv(tiempos, status)


# ---------------------------------------------------------
#  PRUEBA FRECUENTISTA: Comparación de Modelos vía AIC
# ---------------------------------------------------------

# Ajustar modelo Exponencial
ajuste_exp <- survreg(
  surv_obj ~ 1,
  dist = "exponential"
)

# Ajustar modelo Weibull
ajuste_weibull <- survreg(
  surv_obj ~ 1,
  dist = "weibull"
)

# Ajustar modelo Log-Normal
ajuste_lognormal <- survreg(
  surv_obj ~ 1,
  dist = "lognormal"
)

# Ajustar modelo Log-Logístico
ajuste_loglogistic <- survreg(
  surv_obj ~ 1,
  dist = "loglogistic"
)

# ---------------------------------------------------------
# CÁLCULO DEL AIC
# Fórmula: AIC = -2 * LogLik + 2 * k
# ---------------------------------------------------------

# Número de parámetros:
# Exponencial  -> 1
# Weibull      -> 2
# Lognormal    -> 2
# Loglogistic  -> 2

aic_exp <- -2 * ajuste_exp$loglik[1] + 2 * 1

aic_weibull <- -2 * ajuste_weibull$loglik[1] + 2 * 2

aic_lognormal <- -2 * ajuste_lognormal$loglik[1] + 2 * 2

aic_loglogistic <- -2 * ajuste_loglogistic$loglik[1] + 2 * 2

# ---------------------------------------------------------
# TABLA COMPARATIVA
# ---------------------------------------------------------

tabla_aic <- data.frame(
  Modelo = c(
    "Exponencial",
    "Weibull",
    "Log-Normal",
    "Log-Logístico"
  ),
  AIC = c(
    aic_exp,
    aic_weibull,
    aic_lognormal,
    aic_loglogistic
  )
)

# Ordenar de mejor a peor (menor AIC = mejor)
tabla_aic <- tabla_aic[order(tabla_aic$AIC), ]

# ---------------------------------------------------------
# RESULTADOS
# ---------------------------------------------------------

cat("\n=======================================================\n")
cat("          COMPARACIÓN DE MODELOS (AIC)                 \n")
cat("=======================================================\n")

print(tabla_aic)

cat("-------------------------------------------------------\n")

# Mejor modelo
mejor_modelo <- tabla_aic$Modelo[1]

cat(sprintf(
  "El modelo con menor AIC es: %s\n",
  mejor_modelo
))

# Diferencia respecto al segundo mejor modelo
delta_aic <- tabla_aic$AIC[2] - tabla_aic$AIC[1]

cat(sprintf(
  "Diferencia de AIC respecto al segundo mejor: %.2f\n",
  delta_aic
))

# Interpretación
if (delta_aic > 2) {
  
  cat(
    "CONCLUSIÓN: Existe evidencia significativa a favor del mejor modelo.\n"
  )
  
} else {
  
  cat(
    "CONCLUSIÓN: La diferencia entre los mejores modelos no es contundente.\n"
  )
  
  cat(
    "Por parsimonia, podría preferirse el modelo más simple.\n"
  )
}

cat("=======================================================\n")