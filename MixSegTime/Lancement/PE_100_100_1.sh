#!/bin/bash
#OAR --project avjc
#OAR -n PE_100_100_1
#OAR -l /core=1,walltime=200:00:00
source ~/activate_R.sh
# Lancement...
R CMD BATCH /home/vbrault/MixSegTime/Plan_t_fixe_dgp/PE_n_100_d_100_alpha_1.R

