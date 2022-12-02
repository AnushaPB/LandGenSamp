#job script to create some plots to look at sampling schemes

#run script
Rscript sampling_viz_test.R

#export png to bdrive
rclone copy sampling_plots.pdf bdrive:plots
rclone copy site_sampling_plots.pdf bdrive:plots
