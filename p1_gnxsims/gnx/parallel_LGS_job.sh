#run this script in LandGenSamp/p1_gnxsims/gnx

#create genomic architecture file 
Rscript create_genomic_architecture

#create nnloci dir if it does not exist
mkdir -p nnloci

#run simulation script
python3 parallel_LGS.py > parallel_LGS.pyout

#run data cleanup script
bash data_cleanup.sh

#backup data to bdrive
rclone copy LGS_data bdrive:new_LGS_data

date "+%H:%M:%S   %d/%m/%y"
