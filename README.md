# AnalysisofIBMs

Scripts for the paper entitled ["Analysis of an Individual-Based Influenza Epidemic Model using Random Forest Metamodels and Adaptive Sequential Sampling"](https://onlinelibrary.wiley.com/doi/full/10.1002/sres.2763) by Edali &amp; Yücel (2020).

This repository contains four files (two R scripts and two data files):

* FluTE_DataGeneration.r: This R script generates a training and a test for experimentation.
* FluTE_Analysis.r: This R script performs the three-step analysis procedure proposed in the paper by using a training and a test set.
* flute_tr.csv: For those who do not want to generate a training set, this file contains a ready-to-use training set with 30 instances.
* flute_test.csv: For those who do not want to generate a test set, this file contains a ready-to-use test set with 1000 instances.

As a result, the analyst has two options:

1. The analyst may generate a training and a test set by running FluTE_DataGeneration.r script. Then, s/he must feed these dataset to FluTE_Analysis.r script in order to perform the analysis procedure proposed in the manuscript.
2. The analyst may directly use the provided training and test set by feeding them to FluTE_Analysis.r script in order to perform the analysis procedure proposed in the manuscript.

## Reminders and Warnings

* Both script files contain comments to ease the application process.
* In order to run the FluTE model, the user is referred to https://github.com/dlchao/FluTE
* The rule extraction phase requires [Gurobi solver](https://www.gurobi.com/). Although it is commercial, Gurobi provides one-year free licence to academics/students through emails with .edu extension if they make their registration in-campus.

