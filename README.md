# AnalysisofIBMs

Scripts for the paper entitled "Analysis of an Individual-Based Influenza Epidemic Model using Random Forest Metamodels and Adaptive Sequential Sampling" by Edali &amp; Yucel (2020)

This repository contains four files (two R scripts and two data files):

* FluTE_DataGeneration.r: This R script generates a training and a test for experimentation.
* FluTE_Analysis.r: This R script performs the three-step analysis procedure proposed in the paper by using a training and a test set.
* flute_tr.csv: For those who do not want to generate a training set, this file contains a ready-to-use training set with 30 instances.
* flute_test.csv: For those who do not want to generate a test set, this file contains a ready-to-use test set with 1000 instances.

In order to run the FluTE model, the user is referred to https://github.com/dlchao/FluTE
