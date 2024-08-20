#!/bin/bash
#OAR --project avjc
#OAR -n PE_100_50_2
#OAR -l /core=1,walltime=200:00:00
source ~/activate_R.sh
# Lancement...
R CMD BATCH /home/vbrault/MixSegTime/Plan_t_fixe_dgp/PE_n_100_d_50_alpha_2.R

