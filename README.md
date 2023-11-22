# M2-Baboon-Project
Estimating Chacma baboon population characteristics in the Garden Route National Park (South Africa) from opportunistic camera-trap data

This project aims to utilize camera traps for biodiversity monitoring, focusing on the Chacma Baboon in the northern region of Knysna, South Africa. The primary goal is to develop a model to assess the characteristics of an unmarked population, using images to identify individuals and analyze their social interactions and distribution.
The complete report and related materials are available in the project's repository.

## Background
Datas used in this project consist of photographs from camera traps, provided by SANParks (South African National Parks). I personally analyzed these photos to create a usable dataset, aimed at addressing our research question.
Since completing my Master's thesis until now, additional data has been analyzed and slight changes have been made to the code. Although the overall project remains the same, it is important to note that some results obtained in the R and Python scripts may differ from those in the report. The project was initially coded in R (with RStudio) and was transcribed into Python by myself for training.

## Important Note on Running the Bayesian Model

To successfully run the Bayesian model included in our code, it is essential to have JAGS (Just Another Gibbs Sampler) installed on your system. JAGS is utilized for Bayesian analysis and is a crucial component for executing the model.

Please be aware that running the Bayesian model can be time-consuming. To save you from the lengthy process of installation and waiting for the model to run, I have provided the results directly in the Bayesian_results.csv file. This file contains the output of the Bayesian analysis and can be consulted for a quick overview of the model's results.

To download JAGS, please visit the file page of the mcmc-jags project on SourceForge: https://mcmc-jags.sourceforge.io/

## Repository Structure

- `/R`: Original R scripts and notebooks.
- `/Python`: Transcribed Python scripts.
- `/Data`: Dataset.

# M2-Baboon-Project
Estimating Chacma baboon population characteristics in the Garden Route National Park (South Africa) from opportunistic camera-trap data

Ce projet vise à exploiter les pièges photographiques pour surveiller la biodiversité, en se concentrant sur le Babouin Chacma dans la région nord de Knysna, en Afrique du Sud. L'objectif principal est de développer un modèle pour évaluer les caractéristiques d'une population non marquée, en utilisant des images pour identifier les individus et analyser leurs interactions sociales et leur répartition. 
Le rapport complet et les documents associés sont disponibles dans le dépôt du projet. Le projet a été initialement codé en R (avec RStudio) et a par la suite été transcrit par moi-même en Python dans le but de m'entrainer.


## Contexte

Les données utilisées dans ce projet sont des photographies issues de pièges photographiques, fournies par les SANParks (Parcs Nationaux Sud-Africains). J'ai personnellement analysé ces photos pour créer un jeu de données exploitable, dans le but de répondre à notre problématique de recherche.
Depuis la finalisation de mon rapport de Master 2 jusqu'à aujourd'hui, de nouvelles données ont été analysées et des modifications mineures ont été apportées au code. Bien que l'essence globale du projet reste inchangée, il est important de noter que certains résultats obtenus dans les scripts R et Python pourraient différer de ceux présentés dans le rapport. 

## Note Importante sur l'Exécution du Modèle Bayésien

Pour exécuter avec succès le modèle bayésien inclus dans le code, il est essentiel d'avoir installé JAGS (Just Another Gibbs Sampler) sur votre système. JAGS est utilisé pour l'analyse bayésienne et est un composant crucial pour l'exécution du modèle.

Veuillez noter que l'exécution du modèle bayésien peut être longue. Pour vous éviter le processus d'installation et d'attente, j'ai fourni directement les résultats dans le fichier Bayesian_results.csv. Ce fichier contient la sortie de l'analyse bayésienne et peut être consulté pour un aperçu rapide des résultats du modèle sans avoir besoin de l'exécuter.

Pour télécharger JAGS, veuillez visiter la page de fichiers du projet mcmc-jags sur SourceForge : https://mcmc-jags.sourceforge.io/

## Structure du Répertoire

- `/R`: Scripts et notebooks R originels.
- `/Python`: Scripts Python transcrits.
- `/Data`: Jeu de données.
