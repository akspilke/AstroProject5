# AstroProject5

This repository contains all programs, figures and films from the Astro version of Project 5 in Computational Physics at UiO, simulating
the collapse of an open galactic cluster. We encourage the reader to have a look at the films, as illustrate whats happening in this project.
In order to run the C++ program, select a number of initial particles N and a size of your "cluster", an initial radius R0. The program produces
three different textfiles (one for positions, one for energies and one for a histogram of densities), comment out the desired writeToFile commands.
There is a function for calculating the energy of bound objects or the energy of all objects, see comments in the code. 

The Python program only contains plotting commands, change the ones you want to look at to "True". The variable name for the if statement before
each plot should give an indication for what is being plotted; "plot_galaxy" plots the trajectories of star particles, "energy_plot" plots the energies
of all particles (smoothed and unsmoothed depending on text file), "energy_bouded" plots the energies of only the bound particles (smooth/unsmooth),
"virial" plots the test of the virial theorem, "Density_Hist" plots histograms of the density, and "Density_Hist_fixmass_profile"  plots density 
distributions for different N and compares to the density function described in the report, and finally "ploteps" compares different values for the
smoothing parameter epsilon.

Figures are in the final_figs and films are in the Movies folder. The text files are in the text file folder, but we were not able to upload all text files due to the limitation of 100MB on github, let us know if you want to have a look at these! 

Thanks for all the help this semester we have learnt a lot!
<3 C++ <3

Tiffany & Andri
