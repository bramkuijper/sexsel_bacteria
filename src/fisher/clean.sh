#!/usr/bin/env bash

echo "Do you want to delete all batch files? (y/N)"
read -p 'Option: ' delete_batch_var

if [[ "$delete_batch_var" == "y" ]]; then
    echo "deleting batch files..."
    rm -f batch_file_*
    rm -f batch_record_*
fi


echo "Do you want to delete all graphs? (y/N)"
read -p 'Option: ' delete_graph_var

if [[ "$delete_graph_var" == "y" ]]; then
    echo "deleting graphs..."
    rm -f graph_sim_bact*
fi

echo "Do you want to delete all simulation output files? (y/N)"
read -p 'Option: ' delete_sim

if [[ "$delete_sim" == "y" ]]; then
    echo "deleting simulation output files..."
    rm -f sim_bact_sexsel_*
fi
