#!/bin/bash
#OAR --project avjc
#OAR -n PE_1000_50_2
#OAR -l /core=1,walltime=200:00:00
source ~/activate_R.sh
# Lancement...
R CMD BATCH /home/vbrault/MixSegTime/Plan_t_fixe_dgp/PE_n_1000_d_50_alpha_2.R

