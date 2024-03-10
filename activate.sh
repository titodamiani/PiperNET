#!/bin/bash

#check if conda environment pipernet already exists
if conda env list | grep -q 'pipernet'; then

    #if env exists, print message and activate environment
    echo "Conda environment 'pipernet' already exists. Activating existing environment..."
    conda activate pipernet

else
    #if env doesn't exist, create it and install packages in requirements.txt
    echo "Creating conda environment pipernet..."
    conda create -y --name pipernet
    echo "Installing packages..."
    conda install --file requirements.txt -y
    conda activate pipernet
fi

#export cwd to PYTHONPATH
echo "Exporting current working directory to PYTHONPATH..."
export PYTHONPATH=$(pwd):$PYTHONPATH
echo "PiperNET ready!"