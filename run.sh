#!/bin/sh

echo "-database /home/rojan/software/rosetta/rosetta32/rosetta32ss/rosetta_database" > flags
# echo "-database /data/rojan/public/software/rosetta/rosetta32_no_mpi/rosetta_database" > flags

./NEFILIM -r @flags \
          --path_to_scores ./inputs/unit_apmds.lst \
          --path_to_models ./inputs/unit_model.lst \
          --path_to_native ./inputs/1opd_atom.pdb \
          --window_size 9 \
          --cutoff 0.10 \
          --proportion 1 \
          --number_of_fragment 10 \
          --fragment_per_cluster 2 \
          --number_of_template 3 \
          --cluster_radius 0.60
