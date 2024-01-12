#Import des librairies utiles
library(mvtnorm)
library(glmnet)
library(xtable)
library(ggplot2)
library(gridExtra)
library(knitr)
library(dplyr)
library(kableExtra)
library(corrplot)
library(pls)
library(stargazer)

#Lecture des données

data=read.table("transportmod3.txt",h=T)
View(data)
summary(data)

########################### Statistiques descriptives ##########################

- "Valeurs descriptives"

# Sélection des dolonnes à décrire
colonnes_a_analyser <- c("CO2", "Incid", "CA")
donnees_selectionnees <- data[colonnes_a_analyser]

# Affichage des statistiques descriptives
res_summary <- summary(donnees_selectionnees)
table_latex <- xtable(res_summary)

# Génération du code LaTeX
print(table_latex)


- "Graphe des coorélations"

# Définition un seuil pour les fortes corrélations 
threshold <- 0.7
correlation_matrix <- cor(data)

# Récupérer les noms de lignes et de colonnes de corrélations fortes
strong_correlations <- as.data.frame(which(abs(correlation_matrix) > threshold, arr.ind = TRUE))

# Afficher la carte thermique pour les corrélations fortes
ggplot(strong_correlations, aes(Var2, Var1, fill = correlation_matrix[abs(correlation_matrix) > threshold])) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "Correlation Heatmap for Strong Correlations")

# Créer un graphique en violon
ggplot(data, aes(x = CA)) +
  geom_bar(fill = "#69b3a2", color = "#e9ecef") +
  labs(title = "Countplot of Variable Binaire (CA)",
       x = "CA Category",
       y = "Count")



# Création de l'histogramme
hist(data$Incid, main = "Histogramme de ma_variable", xlab = "Valeurs de ma_variable", col = "skyblue", border = "black")

# Création de l'histogramme avec ggplot2
ggplot(data, aes(x = CA3)) +
  geom_histogram(fill = "skyblue", color = "black", bins = 3, alpha = 0.7) +
  labs(title = "Histogramme de CA3",
       x = "Valeurs de CA3",
       y = "Fréquence") +
  theme_minimal()

############################# Conditionnenment de XX ###########################

X=as.matrix(data[,5:64])

# Calculer la matrice X'X
XTX <- t(X) %*% X

# Calculer les valeurs propres de X'X
eigen_XTX <- eigen(XTX)
cat("Valeurs propres de X'X:", eigen_XTX$values, "\n")

# Calculer l'indice de conditionnement
condition_index <- max(eigen_XTX$values) / min(eigen_XTX$values)

# Afficher l'indice de conditionnement
cat("Indice de conditionnement:", condition_index, "\n")

########################### Standardisation de la donnée #######################

#Centrer et réduire les variables explicatives
Y=data$CO2
X=as.matrix(data[,5:64])
XS=scale(X)

########################### Modélisation du CO2 ################################


" Modèle RIDGE " 

# On a un modèle linéaire. On va prendre family = 'gaussian'
rescv=cv.glmnet(XS,Y,family='gaussian',alpha=0, lambda=seq(0.1, 100, 0.1))
plot(rescv)

lambda=rescv$lambda.min
lambda

# Stabilisation du lambda avec une boucle
lambda=0
for (j in 1:10)
{
  rescv=cv.glmnet(XS,Y,family='gaussian',alpha=0, lambda=seq(0.1, 100, 0.1))
  print(rescv$lambda.min)
  lambda=lambda+rescv$lambda.min
}
#Redéfinition du seuil
seuil=lambda/10
rescv=glmnet(XS,Y,family='gaussian',alpha=0, lambda=seuil) 
#On retire le coefficient -1
coef=coefficients(rescv)[-1] 

# Choix des variables en focntion du plot
plot(sort(abs(coef)))
which(abs(coef)>0.6) 

# Choix automatique des variables qui ressortent à l'aide d'un boxplot 
resbox=boxplot(coef)
which(coef>resbox$stats[5,])

selecridge=order(coef,decreasing=TRUE)[1:5]
colnames(X[,selecridge])

'Les ports retenus sont: 3  5  8 34 58'


################################### Regression LASSO 
# Modèle de LASSO

lambda=0
for (j in 1:10)
{
  reslass=cv.glmnet(XS,Y,family='gaussian',alpha=1, lambda=seq(0.1, 100, 0.1))
  print(reslass$lambda.1se)
  lambda=lambda+reslass$lambda.1se
}
seuil=lambda/10 # Définition du sueillage par calcul de la moyenne
# On va alors redefinir un meilleur seuil et relancer notre regression de ridge
reslass=glmnet(XS,Y,family='gaussian',alpha=1, lambda=seuil)

coef=coefficients(reslass)[-1]
plot(sort(abs(coef)))
which(coef!=0)
"On retient les meme ports sauf le 7 qui vient s'ajouter"
plot(res)


################################### ELASTINET 
"On procède de la meme façon en stabilisant notre modèle sur un intervalle de test de valeur"
lambda=0
for (j in 1:10)
{
  resenet=cv.glmnet(XS,Y,family='gaussian',alpha=0.5,lambda=seq(0.1, 100, 0.1))
  print(resenet$lambda.1se)
  lambda=lambda+resenet$lambda.1se
}
seuil=lambda/10
resenet=glmnet(XS,Y,family='gaussian',alpha=0.5,lambda=seuil)

coef=coefficients(resenet)[-1]
plot(sort(abs(coef)))
which(coef!=0)
coef
# regardons le boxplot des coef non nuls 
coefnonnul=coef[coef!=0]
boxplot(coefnonnul)
selecnet=order(coef,decreasing=TRUE)[1:8]


################################### OLS 
'
Modele OLS sur les variables retenues: nous utilisons les coeff retenus dans la regression elastinet pour ce cas
Il nous faudrait penser à justifier les résultats et pourquoi un tel choix dans notre modélisation
Faire attention à standardiser une nouvelle fois nos variables avant de les incluires dans notre regression ols 
'
reslm=lm(Y~scale(X[,selecnet]))
summary(reslm)

################################### Modèle PCA 
library(pls)
library(caret)

# Effectuez le PCR
pcr_model <- pcr(Y ~ XS, scale = TRUE)

# Résumé du modèle PCR
summary(pcr_model)

# Initialisation du vecteur pour stocker les RMSEP
rmsep_values <- numeric(36)

# Boucle sur le nombre de composantes
for (i in 1:36) {
  pcr_model <- pcr(Y ~ XS, ncomp = i, scale = TRUE)
  
  # Prédiction sur l'ensemble de données
  predicted_values <- predict(pcr_model, as.data.frame(X))
  
  # Calcul du RMSEP
  rmsep_values[i] <- sqrt(mean((Y - predicted_values)^2))
}

# Tracer le graphique du RMSEP
plot(1:36, rmsep_values, type = "b", xlab = "Nombre de composantes", ylab = "RMSEP")
# Trouver le nombre optimal de composantes
optimal_components <- which.min(rmsep_values)

# Seuillage pour les coefficients
threshold <- 0.3 # Choisissez votre seuil
significant_coef <- abs (coef(pcr_model)) > threshold

# Afficher les coefficients significatifs
significant_coef

par(mfrow=c(3,1))

plot(respca,ncomp=2)
abline(0,1)

plot(respca,ncomp=4)
abline(0,1)

plot(respca,ncomp=8)
abline(0,1)

plot(respca,ncomp=12)
abline(0,1)

plot(respca,ncomp=16)
abline(0,1)

plot(respca,ncomp=20)
abline(0,1)


respca2=plsr(Y ~XS, ncomp=16)
summary(respca2)

coefpca=coef(respca2)
par(mfrow=c(1,1))
plot(coefpca) # on pourrait faire du seuillage ici
plot(sort(abs(coefpca)))
selecpca=order(abs(coefpca),decreasing=TRUE)[1:9] 

colnames(X[,selecpca])
################################### Modèle partial least square

# Réaliser la régression PLS
pls_model <- plsr(Y ~ XS, ncomp = 16)  # Vous pouvez ajuster ncomp selon votre analyse

# Afficher un résumé du modèle PLS
summary(pls_model)

# Visualiser les coefficients PLS
coef_pls <- coef(pls_model, ncomp = 16)  # Vous pouvez ajuster ncomp ici aussi
print(coef_pls)


############################################################################## modélisation du nombre d’incidents ############################################################################## 
'On va redefinir nos variable car on change de modèle
Notre Y sera ici le nombre incidents
les X restent les ports
'
Y=data$CO2
X=as.matrix(data[,5:64])
XS=scale(X)# Standardisation de la donées
################################### RIDGE DE POISSON 
# Modèle de Ridge de poisson
# Sélectionner le meilleur modèle en utilisant la validation croisée
lambda=0
for (j in 1:10)
{
  rescv_poisson <- cv.glmnet(XS, Y, family = "poisson", alpha=0, lambda=seq(0.1, 100, 0.1))
  # Trouver le meilleur lambda (paramètre de régularisation)
  print(rescv_poisson$lambda.min)
  lambda=lambda+rescv_poisson$lambda.min
}
seuil=lambda/10

# on refait notre modèle avec le nouveau seuil trouvé
rescv_poisson <- glmnet(XS, Y, family = "poisson", alpha=0, lambda=seuil)

# Obtenir le vecteur des coefficients du meilleur modèle
best_coeffs <- coef(rescv_poisson)
# Sélectionner les variables non nulles du meilleur modèle
selected_vars <- which(best_coeffs != 0)
# Obtenez le nom des variables sélectionnées (en supposant que vous avez des noms de variables)
selected_var_names <- colnames(XS)[selected_vars]

# Deuxieme méthode
coef=coefficients(rescv_poisson)[-1] # on retire le coefficient 1
'la plupart de nos coeff sont différents de zéro donc une autre technique serait optimal
notamment un choix basé sur le seuillage ??
'
plot(sort(abs(coef)))
# Mise en place du seuillage
which(abs(coef)>0.02)

################################### LASSO DE POISSON 
# Modèle de Lasso de poisson
lambda=0
for (j in 1:10)
{
  Lasso_poisson <- cv.glmnet(XS, Y, family = "poisson", alpha=1, lambda=seq(0.1, 100, 0.1))
  # Trouver le meilleur lambda (paramètre de régularisation)
  print(Lasso_poisson$lambda.1se)
  lambda=lambda+Lasso_poisson$lambda.1se
}
seuil=lambda/10

Lasso_poisson=glmnet(XS,Y,family='gaussian',alpha=1, lambda=seuil)

coef=coefficients(Lasso_poisson)[-1]
plot(sort(abs(coef)))
which(coef!=0)
"On retient les meme ports sauf le 7 qui vient s'ajouter"
plot(res)

# Extraire les coefficients du modèle Lasso (s = 0)
lasso_coefficients <- coef(Lasso_poisson, s = 0)
# Obtenir les indices des variables retenues
selected_indices <- which(lasso_coefficients != 0)
# Obtenir les noms des variables retenues
selected_variables <- colnames(lasso_coefficients)[selected_indices]
# Afficher les noms des variables retenues
print(selected_variables)

################################### ELASTICNET DE POISSON 
# Modèle de ELASTICNET de poisson

lambda=0
for (j in 1:10)
{
  Elastinet_poisson <- cv.glmnet(XS, Y, family = "poisson", alpha=0.5, lambda=seq(0.1, 100, 0.1))
  # Trouver le meilleur lambda (paramètre de régularisation)
  print(Elastinet_poisson$lambda.1se)
  lambda=lambda+Elastinet_poisson$lambda.1se
}
seuil=lambda/10

Lasso_poisson=glmnet(XS,Y,family='gaussian',alpha=0.5, lambda=seuil)

################################### REGRESSION DE POISSON 
# Modèle de Regression de poisson
# Exemple de régression de Poisson
# Supposons que votre dataframe s'appelle "df" et que vous voulez sélectionner les colonnes "col1", "col2" et "col3"
ports_selectionnees <- X[, c("Port.Port.3","Port.Port.5", "Port.Port.7", "Port.Port.58")]

model_poisson <- glm(Y ~ ports_selectionnees, family = poisson)


############################################################################## modélisation de CA ############################################################################## 
################################ RIDGE
Y=data$CA
X=as.matrix(data[,5:64])
XS=scale(X)# Standardisation de la donées

# Supposons que y est initialement une variable binaire (0 ou 1)
Y <- as.factor(Y)

# Définir "1" comme le niveau par défaut
levels(Y) <- c("0", "1")

# Réalisez une régression logistique Ridge avec validation croisée
cv.ridge <- cv.glmnet(XS, Y, alpha = 0, family = "binomial",lambda=seq(0.1, 20, 0.2))
best_lambda <- cv.ridge$lambda.min  # ou cv.ridge$lambda.1se selon votre choix
plot(cv.ridge)

# Entraînez le modèle final avec le meilleur lambda
ridge_model <- glmnet(XS, Y, alpha = 0, lambda = best_lambda, family = "binomial")

# Sélection des variables
coefridge=coefficients(ridge_model)[-1]
plot(sort(abs(coefridge))) 
boxplot(coefridge)
selec3=order(abs(coefridge),decreasing=TRUE)[1:3]
selec3

selec12=order(abs(coefridge),decreasing=TRUE)[1:6]
selec12

colnames(X[,selec3])
colnames(X[,selec12])

################################ LASSO

rescv=cv.glmnet(XS,Y,family='binomial',alpha=1,lambda=seq(0.1, 20, 0.2))
plot(rescv)
seuil=rescv$lambda.1se # s?lection => lambda.1SE 
# le lambda.1se varie car il d?pend de la partition utilis?e pour la CV 
# l'id?al 
reslasso=glmnet(X,Y,family='binomial',alpha=1,lambda=seuil)
reslasso
coeflasso=coefficients(reslasso)[-1]
plot(sort(abs(coeflasso))) 
selec6=order(abs(coeflasso),decreasing=TRUE)[1:6]
colnames(X[,selec6])

################################ AUC & Interprétation
# Supposons que 'final_lasso_model' soit votre modèle LASSO final

# Obtenez les probabilités prédites
probas_lasso <- predict(reslasso, newx = as.matrix(XS), s = seuil, type = "response")

# Supposons que vous avez déjà les probabilités prédites dans une variable 'probas'
# et que 'Y' est votre variable de réponse binaire

library(pROC)

roc_curve <- roc(Y, probas_lasso)
roc_curve
auc_value <- auc(roc_curve)

# Supposons que 'final_lasso_model' soit votre modèle LASSO final

# Obtenez les coefficients
coefficients_lasso <- coef(reslasso, s = seuil)

# Affichez les coefficients
print(coefficients_lasso)


# Supposons que vous avez déjà les coefficients dans une variable 'coefficients'
variable_significative <- names(coeflasso)[which.max(abs(coeflasso))]
cat("La variable la plus significative est :", variable_significative)

############################################################################## modélisation de CA3 ############################################################################## 
# Supposons que X soit votre matrice de variables explicatives et Y votre variable réponse CA3

# Convertir la variable réponse en facteur
Y=data$CA3
X=as.matrix(data[,5:64])
XS=scale(X)# Standardisation de la donées
Y <- as.factor(Y)

# Appliquer la régression multinomiale pénalisée avec la validation croisée
cv_fit_multinom <- cv.glmnet(x = XS, y = Y, family = "multinomial", alpha = 1, lambda=seq(0.1, 20, 0.2))

# Obtenir le meilleur lambda à partir de la validation croisée
best_lambda_multinom <- cv_fit_multinom$lambda.min

# Appliquer la régression multinomiale pénalisée avec le meilleur lambda
fit_multinom_best_lambda <- glmnet(x = XS, y = Y, family = "multinomial", alpha = 1, lambda = best_lambda_multinom)
fit_multinom_best_lambda

# Sélectionner les variables
selected_vars <- coef(fit_multinom_best_lambda)
selected_vars <- rownames(selected_vars[selected_vars != 0, ])

# Afficher les variables sélectionnées
print(selected_vars)

########################## Regression polytomique ordonnée

library(ordinal)

# Supposons que Y_ordered est votre variable réponse ordonnée et X_Select est votre matrice d'variables explicatives
# Convertir Y_ordered en une variable ordonnée si ce n'est pas déjà fait
Y_ordered <- as.ordered(Y)

# Créer une matrice X_Select (remplacez cela par votre propre matrice)
X_Select <- data[, c("Port.Port.5", "Port.Port.24", "Port.Port.29",
                     "Port.Port.40","Port.Port.43","Port.Port.54")]

# Appliquer la régression polytomique ordonnée
fit_ord <- clm(Y_ordered ~ ., data = cbind(X_Select, Y_ordered))

# Afficher les résultats
summary(fit_ord)

stargazer(fit_ord, title = "Régression polytomique ordonnée", out = "table.tex")



##############################################  HISTOGRAMME ##############################################  
library(palmerpenguins)
histo_mass <- ggplot(data)+
  geom_histogram(aes(x=Port.Port.60), fill="darkorchid4", color="darkgray", bins=50)+
  labs(title = "Incid", subtitle = "Histogram")+
  ylab("Count")+theme_light()
histo_mass

a <- ggplot(data, aes(x = Port.Port.25))
a
# Histogram with density plot
a + geom_histogram(aes(y = ..density..), 
                   colour="black", fill="white") +
  geom_density(alpha = 0.2, fill = "#FF6666")

a +geom_histogram(aes(y = after_stat(density)), bins = 30, colour = "black", fill = "white") +
  geom_density(alpha = 0.2, fill = "#FF6666")





################################### Modèle PLS ###################################
respls=plsr(Y~XS,ncomp=20)
summary(respls)

par(mfrow=c(3,1))

plot(respls,ncomp=2)
abline(0,1)

plot(respls,ncomp=4)
abline(0,1)

plot(respls,ncomp=8)
abline(0,1)

plot(respls,ncomp=12)
abline(0,1)

plot(respls,ncomp=16)
abline(0,1)

plot(respls,ncomp=20)
abline(0,1)

plot(respls,ncomp=24)
abline(0,1)

respls2=plsr(Y ~XS, ncomp=12)
summary(respls2)

coefpls=coef(respls2)
par(mfrow=c(1,1))
plot(coefpls) # on pourrait faire du seuillage ici
plot(sort(abs(coefpls)))
selecpls=order(abs(coefpls),decreasing=TRUE)[1:9] 

colnames(X[,selecpls])






respls2=plsr(Y ~X, ncomp=10)
coefpls=coef(respls2)
