# Quick Explanation
Main codes related to the manuscript "Predictive Non-linear Dynamics via Neural Networks and Recurrence Plots", authored by Luiza Lober¹, Matheus Palmero² and Francisco A. Rodrigues.

Codes in C produce the Recurrence Plots (RPs) related to the non-linear dynamical system of choice. The Makefile should compile and run all necessary dependencies (just run `make`), and the RPs will be placed on the `results_path` (in rec_plots_chaotic_systems.c).

Notebooks .ipynb are responsible for the Convolutional Neural Network (CNN) part. It will provide the main results on the predicted control parameters of the models.  

making_gif.py produces the GIFs uploaded to the supplemental material.

Correspondence: ¹luiza.lober@usp.br; ²palmero@usp.br
