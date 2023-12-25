# Importation des bibliothèques
library(ggplot2)
library(forecast)

########################## Description des données ############################

# Lecture des données
data <- read.table("donnees.csv", header = TRUE, sep = ",", skip = 0)

# Conversion des températures en °C
data$Value <- (data$Value - 32)/1.8

# Ajout de la colonne des indices temporelles
data["t"] <- 1:278

# Estimation de la température en fonction du temps par régression linéaire
# afin de dégager la tendance
model <- lm(Value ~ t, data = data)
summary(model) 
# On trouve X = a + b*t avec a = 17.167564 et b = 0.006426
# On constate une tendance croissance

# Affichage de la série et de la tendance
plot(Value ~ t, data = data, type = "l", main = "Série et tendance", 
     xlab = "Temps", ylab = "Températures (°C)", col.lab = 'blue', 
     col.main = 'blue')
abline(model, col = "red")

# Création de la série temporelle
serie <- ts(data$Value, start = c(2000, 1), end = c(2023, 2), frequency = 12)
print(serie)
plot(serie, main = "Séries des températures moyennes mensuelles à San Diego", 
     xlab = "Temps", ylab = "Température (°C)", col.main = 'blue')

# Affichage des séries saisonnières
monthplot(serie)

# Affichage par années
colors <- rainbow(23)
plot(serie[1:12], type = "l", ylim = c(min(serie), max(serie)), 
     xlab = "Mois", ylab = "Température (°C)")
for (i in 1:23) lines(serie[(1+12*i):(12*(i+1))], col = colors[i])
legend("topright", legend = paste("Année", 2000:2022), col = colors, lty = 1)

# Autocorrélation
auto_corr <- acf(serie, lag.max = 60, type = 'correlation', plot = T)
# On voit une saisonnalité de période 12
# Et confirmation d'une tendance (fonction d'auto-corrélation qui décroit)

# Nuages de points (X_t, X_(t+h)) pour h = 1, h = 3, h = 6 et h = 12
par(mfrow=c(2,2)) 
matrice_1 <- matrix(0, 277, 2)
cpt <- 1
for (i in 1:277){
  matrice_1[cpt, 1] <- data$Value[(i-1)*12 + 1]
  matrice_1[cpt, 2] <- data$Value[i*12 + 1]
  cpt <- cpt + 1
}
matrice_1 <- data.frame(matrice_1)
plot(matrice_1$X1, matrice_1$X2, main = 'Nuage de points pour h = 1', 
     xlab = 'X_t', ylab = 'X_(t+1)', col.main = 'blue', col.lab = 'blue')
matrice_3 <- matrix(0, 273, 2)
cpt <- 1
for (t in 1:3) {
  for (i in 1:91){
    matrice_3[cpt, 1] <- data$Value[(i-1)*3 + t]
    matrice_3[cpt, 2] <- data$Value[i*3 + t]
    cpt <- cpt + 1
  }
}
matrice_3 <- data.frame(matrice_3)
plot(matrice_3$X1, matrice_3$X2, main = 'Nuage de points pour h = 3', 
     xlab = 'X_t', ylab = 'X_(t+3)', col.main = 'blue', col.lab = 'blue')
matrice_6 <- matrix(0, 270, 2)
cpt <- 1
for (t in 1:6) {
  for (i in 1:45){
    matrice_6[cpt, 1] <- data$Value[(i-1)*3 + t]
    matrice_6[cpt, 2] <- data$Value[i*3 + t]
    cpt <- cpt + 1
  }
}
matrice_6 <- data.frame(matrice_6)
plot(matrice_6$X1, matrice_6$X2, main = 'Nuage de points pour h = 6', 
     xlab = 'X_t', ylab = 'X_(t+6)', col.main = 'blue', col.lab = 'blue')
matrice_12 <- matrix(0, 264, 2)
cpt <- 1
for (t in 1:12) {
  for (i in 1:22){
    matrice_12[cpt, 1] <- data$Value[(i-1)*12 + t]
    matrice_12[cpt, 2] <- data$Value[i*12 + t]
    cpt <- cpt + 1
  }
}
matrice_12 <- data.frame(matrice_12)
plot(matrice_12$X1, matrice_12$X2, main = 'Nuage de points pour h = 12', 
     xlab = 'X_t', ylab = 'X_(t+12)', col.main = 'blue', col.lab = 'blue')

# Séparation des données en données d'entrainement et de validation
train <- window(serie, start = c(2000, 1), end = c(2017, 12))
test <- window(serie, start = c(2018, 1), end = c(2023, 2))


################### Désaisonnalisation par moyenne mobile ###################

# Filtre de moyenne mobile d'ordre 12
filtre <- rep(1/12, 11)
filtre <- c(1/24, filtre, 1/24)
M12 <- filter(train, filter = filtre, sides = 2)
ts.plot(train, M12) # Superposition

# Calcul des coefficients saisonniers
S <- train - M12
print(S)
c_temp <- tapply(S, cycle(S), mean, na.rm = TRUE)
print(c_temp) # On constate qu'on a bien les mêmes coefs que sur excel
c_tmp_mean <- mean(c_temp)
print(c_tmp_mean)
c = c_temp - c_tmp_mean # On centre les coefficients saisonniers
print(c)
print(mean(c)) # On voit que la moyenne des coefficients est nulle
print(sum(c)) # De plus, leur somme est nulle

# Voyons les coeffients saisonniers donnés par la fonction décompose
decomp_train <- decompose(train, type = "additive")
coefs_decomp <- decomp_train$seasonal[1:12]
print(coefs_decomp)
print(c - coefs_decomp) # La différence fait 0: ce sont les mêmes

# Calcul de la série corrigées des variations saisonnières (sur train)
Xcvs <- train - decomp_train$seasonal

# Ajustement d'une tendance sur la série corrigée des variations saisonnières
df_train <- data.frame(Xcvs)
df_train$t = 1:216
model1 <- lm(Xcvs ~ t, data = df_train)
summary(model1)
a0 <- model1$coefficients[1]
a1 <- model1$coefficients[2]
# On trouve Xcvs_t_pred = 16.929 + 0.0096*t comme sur Excel

# Prévision sur l'ensemble de test
df_test <- data.frame(test)
df_test$t <- 217:278
df_test$coef_saisonniersMM <- c(rep(c, 5), c[1], c[2])
df_test$predMM <- a0 + a1*df_test$t + df_test$coef_saisonniersMM

# Affichage des prédictions sur l'ensemble de test
ts_test_pred <- ts(df_test, start = c(2018, 1), end = c(2023, 2), 
                   frequency = 12)
plot(df_test$t, df_test$test, type = "l", col = "blue", xlab = "Temps", 
     ylab = "Températures (°C)", col.lab = 'blue', 
     main = "Série de test: janvier 2018 à février 2023", col.main = 'blue')
lines(df_test$t, df_test$predMM, col = "red", lwd = "4")
legend("topright", 
       legend = c("Série initiale", "Série prédite par moyenne mobile"), 
       col = c("blue", "red"), lty = 1)


############## Désaisonnalisation par régression linéaire #####################
##################### avec variables indicatrices #############################

# Création des variables indicatrices
D <- data[1:216, c('t', 'Value')] # Récupération des données d'entrainement
tmp <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
D$C1 <- rep(tmp, 18)
D$C2 <- c(0, rep(tmp, 17), 1, rep(0, 10))
D$C3 <- c(0, 0, rep(tmp, 17), 1, rep(0, 9))
D$C4 <- c(rep(0, 3), rep(tmp, 17), 1, rep(0, 8))
D$C5 <- c(rep(0, 4), rep(tmp, 17), 1, rep(0, 7))
D$C6 <- c(rep(0, 5), rep(tmp, 17), 1, rep(0, 6))
D$C7 <- c(rep(0, 6), rep(tmp, 17), 1, rep(0, 5))
D$C8 <- c(rep(0, 7), rep(tmp, 17), 1, rep(0, 4))
D$C9 <- c(rep(0, 8), rep(tmp, 17), 1, rep(0, 3))
D$C10 <- c(rep(0, 9), rep(tmp, 17), 1, rep(0, 2))
D$C11 <- c(rep(0, 10), rep(tmp, 17), 1, 0)
D$C12 <- c(rep(0, 11), rep(tmp, 17), 1)

# Régréssion linéaire sans constante
model2 <- lm(formula = Value ~ t+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11+C12+0, 
             data = D)
summary(model2)
print(model2$coefficients)
# Mêmes coefficients obetenus sur Excel
beta2 <- model2$coefficients[1] # coefficient directeur
coefs_prov <- model2$coefficients[2:13] # coefficients saisonniers provisoires

# Calcul de la constante beta1 = moyenne des c_i trouvés
beta1 <- mean(coefs_prov)
print(beta1)

# On centre les coefficients
coefs_sais <- coefs_prov - beta1
print(coefs_sais) # Mêmes coeffcients que sur Excel
print(sum(coefs_sais)) # La somme des coefficients saisonniers est bien nulle

# Prédictions sur la période 2018-2023
df_test$coef_saisonniersRVI <- c(rep(coefs_sais, 5), coefs_sais[1], 
                                 coefs_sais[2])
df_test$predRVI <- a0 + a1*df_test$t + df_test$coef_saisonniersRVI

# Affichage des prédictions sur la période 2018-2023
plot(df_test$t, df_test$test, type = "l", col = "blue", xlab = "Temps", 
     ylab = "Températures (°C)", col.lab = 'blue', 
     main = "Série de test: janvier 2018 à février 2023", col.main = 'blue')
lines(df_test$t, df_test$predRVI, col = "yellow", lwd = "4")
legend("topright", legend = c("Série initiale", "Série prédite par RVI"), 
       col = c("blue", "yellow"), lty = 1)

# Comparaison coefficients saisonniers MM et RVI
print(df_test$coef_saisonniersMM[1:12])
print(df_test$coef_saisonniersRVI[1:12])
# Les coefficients saisonniers trouvés sont quasiment les mêmes avec des 
# variations de l'ordre de 10^(-2)
print(df_test$coef_saisonniersMM[1:12] - df_test$coef_saisonniersRVI[1:12])

# Affichage des prédictions avec les deux méthodes
plot(df_test$t, df_test$test, type = "l", col = "blue", xlab = "Temps", 
     ylab = "Températures (°C)", col.lab = 'blue', 
     main = "Série de test: janvier 2018 à février 2023", col.main = 'blue')
lines(df_test$t, df_test$predMM, col = "red", lwd = "4")
lines(df_test$t, df_test$predRVI, col = "yellow", lwd = "4")
legend("topright", legend = c("Série initiale", "Série prédite par MM", 
                  "Série prédite par RVI"), 
       col = c("blue", "red", "yellow"), lty = 1)
# On remarque que les séries prédites avec les deux méthodes sont presques
# identiques


######################## Lissage exponentiel simple ##########################

X_cvs <- a0 + a1*(df_train$t)

# Ajustement d'un LES (lissage exponentiel simple) 
LES <- ets(X_cvs, model = "ANN", additive.only = TRUE) # LES pour modèle additif
summary(LES) 
# R a trouvé alpha = 0.9999 comme constante de lissage la plus adapté
# Cela signifie qu'on ne donne presque aucune importance aux observations
# anciennes.
# Par défaut, alpha est ici estimé par maximum de vraissemblance

# Avec une valeur de alpha
LES <- ets(X_cvs, model = "ANN", additive.only = TRUE, alpha = 0.1)
summary(LES)
# Nous souhaitons obtenir de meilleurs prévisions en donnant plus d'importance 
# aux valeurs passées donc en réduisant la valeur d'alpha.
# De plus on remarque que la prévision est la meilleure en terme d'erreur
# absolue moyenne pour alpha = 0.1

# Prévisions sur l'ensemble de test (sans saisonnalités)
test_pred <- forecast(LES, h = 62)
print(test_pred)
plot(test_pred, main = 'Lissage exponentiel simple', col.main = "blue",
     xlab = "Temps", ylab = "Températures (°C)", col.lab = 'blue', lwd = 2)
# Dans le cas du LES, la prévision est constante, 18.91 ici.
# On peut voir les intervalles de confiance à 80% et 95%
# On obtient des résultats plus proches de la réalité en ajoutant les 
# coefficients saisonniers trouvés précédemment à la prédiction

# Ajout des coefficients saisonniers aux predictions
plot(x = df_test$t, y = test_pred$mean + df_test$coef_saisonniersRVI,
     type = 'l', main = 'Prévisions', col = 'red', ylim = c(13,26),
     col.main = 'blue', xlab = 'Temps', ylab = 'Températures (°C)', lwd = 2,
     col.lab = 'blue')
points(x = df_test$t, y = df_test$test, type = 'l')
points(x = df_test$t, y = rep(test_pred$mean[1], 62), type = 'l', col = 'blue')
legend("topright", legend = c("Prédictions", "Série réelle", 'LES sans saison'), 
       col = c("red", "black", 'blue'), lty = 1)

# Calcul de l'erreur absolue moyenne
eps_abs <- rep(0, 62)
for (i in 1:62) eps_abs[i] = abs(df_test$test[i] - (test_pred$mean[1] + 
                                   df_test$coef_saisonniersRVI[i]))
MAE_LES <- mean(eps_abs)
print(MAE_LES)
# L'ecart moyen entre les valeurs réelles et les valeurs prédite est d'environ
# 0.84


########################## Lissage Holt Winters #############################

# Création du modèle de Holt Winters
HW_model <- HoltWinters(train, alpha = NULL, beta = NULL, gamma = NULL,
                        seasonal = "add")
print(HW_model$coefficients)
# Les coefficients ressemblent à ceux trouvés avec les autres méthodes
print(HW_model)
# R trouve alpha = 0.37, beta = 0.006 et gamma = 0.30

# Prédictions sur l'ensemble de test
HW_pred <- predict(HW_model, n.ahead = 62)

# Affichage des prédictions
plot(x = df_test$t, y = HW_pred, type = 'l', main = 'Prévisions Holt Winters'
     , col = 'red', ylim = c(13,26), col.main = 'blue', xlab = 'Temps', 
     ylab = 'Températures (°C)', lwd = 2, col.lab = 'blue')
points(x = df_test$t, y = df_test$test, type = 'l')
legend("topright", legend = c("Prédictions HW", "Série réelle"), 
       col = c("red", "black"), lty = 1)

# Calcul de l'erreur absolue moyenne
eps_abs_HW <- rep(0, 62)
for (i in 1:62) eps_abs_HW[i] = abs(df_test$test[i] - HW_pred[i])
MAE_HW <- mean(eps_abs_HW)
print(MAE_HW)
# L'ecart moyen entre les valeurs réelles et les valeurs prédites est d'environ
# 0.92.

# Affichage des deux prédictions: par LES et par HW
plot(x = df_test$t, y = HW_pred, type = 'l', main = 'Prévisions avec LES et HW'
     , col = 'red', ylim = c(13,26), col.main = 'blue', xlab = 'Temps', 
     ylab = 'Températures (°C)', lwd = 2, col.lab = 'blue')
points(x = df_test$t, y = test_pred$mean + df_test$coef_saisonniersRVI, 
       type = 'l', col = 'blue', lwd = 2)
points(x = df_test$t, y = df_test$test, type = 'l', col = 'black')
legend("topright", 
       legend = c("Prédictions HW", "Prédictions LES", "Série réelle"), 
       col = c("red", "blue", "black"), lty = 1)
# A vue d'oeil, les deux prédictions sont assez similaires et on ne discerne pas
# de différence quant à leur efficacité

# Comparaison des erreurs absolues moyennes
print(MAE_LES)
print(MAE_HW)
# Elle est légèrement supérieure pour le LES que pour la méthode de HW


