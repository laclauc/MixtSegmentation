#!/bin/bash
#OAR --project avjc
#OAR -n Comp_1000_100_10
#OAR -l /core=1,walltime=200:00:00
source ~/activate_R.sh
# Lancement...
R CMD BATCH /home/vbrault/MixSegTime/Plan_Comp/PE_n_1000_d_100_alpha_10.R

