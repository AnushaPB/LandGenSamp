#run this script in LandGenSamp/p1_gnxsims/parallel

python3 parallel_LGS.py > parallel_LGS.pyout

cd ..

rclone copy parallel/ bdrive:new_data

date "+%H:%M:%S   %d/%m/%y"
