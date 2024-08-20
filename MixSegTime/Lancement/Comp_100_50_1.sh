#!/bin/bash
#OAR --project avjc
#OAR -n Comp_100_50_1
#OAR -l /core=1,walltime=200:00:00
source ~/activate_R.sh
# Lancement...
R CMD BATCH /home/vbrault/MixSegTime/Plan_Comp/PE_n_100_d_50_alpha_1.R

