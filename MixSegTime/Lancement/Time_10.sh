#!/bin/bash
#OAR --project avjc
#OAR -n Time_10
#OAR -l /core=1,walltime=200:00:00
source ~/activate_R.sh
# Lancement...
R CMD BATCH /home/vbrault/MixSegTime/Plan_t_fixe_time/PE_alpha_10.R

