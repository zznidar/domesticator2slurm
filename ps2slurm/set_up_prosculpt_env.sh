#!/bin/bash
echo SETTING UP PROSCULPT RUNNING AS ${USER}@${HOSTNAME}

#activate AF2 enviorment
source /home/folivieri/miniforge3/etc/profile.d/conda.sh
conda activate /home/folivieri/.conda/envs/prosculpt

if [ $HOSTNAME != hpc.ki.si ]
then
    echo SETTING UP PROXY
    #Enable internet access on compute nodes
    source /home/aljubetic/bin/setup_proxy_settings.sh
fi

echo RUNNING IN:
pwd

export NVIDIA_VISIBLE_DEVICES=$GPU_DEVICE_ORDINAL
export CUDA_VISIBLE_DEVICES=$GPU_DEVICE_ORDINAL
#export CUDA_VISIBLE_DEVICES=$NVIDIA_VISIBLE_DEVICES
#export PATH:$PATH:/usr/local/cuda-11.3
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib/:$LD_LIBRARY_PATH
echo Free device is $CUDA_VISIBLE_DEVICES