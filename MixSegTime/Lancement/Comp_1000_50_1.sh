#!/bin/bash
#OAR --project avjc
#OAR -n Comp_1000_50_1
#OAR -l /core=1,walltime=200:00:00
source ~/activate_R.sh
# Lancement...
R CMD BATCH /home/vbrault/MixSegTime/Plan_Comp/PE_n_1000_d_50_alpha_1.R

