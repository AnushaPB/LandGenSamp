# knit File S4 and S5
# note 1: this creates the figure PDFs which are placed in the "figures" directory
Rscript -e "rmarkdown::render(here::here('p4_analysis', 'FileS6_GEA.Rmd')"
Rscript -e "rmarkdown::render(here::here('p4_analysis', 'FileS7_IBDIBE.Rmd')"
Rscript -e "rmarkdown::render(here::here('p4_analysis', 'FileS5_simulation_summary_statistics.Rmd')"

# note 2: supplementary plots/other data visualizations were created from running notebooks interactively
# analysis_megatable.Rmd - used to create a giant table of results from the mixed effects models
# example_simulations.Rmd - used to create plots of examples from simulations (e.g. Figures 1, S1 and S2)
