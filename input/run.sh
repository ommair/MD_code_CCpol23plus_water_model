#!/bin/bash
#
#SBATCH --job-name=MY_MD_mpi
#SBATCH --output=myjob-%A_%a.out
#SBATCH --error=myjob-%A_%a.err
#SBATCH --partition=kolos2
#
###SBATCH --constraint="infiniband"
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
##SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:0
#SBATCH --mem-per-cpu=1024MB

echo "***********************************************"
hostname
ulimit -a

##SBATCH --ntasks-per-core=1
##SBATCH --exclusive


# this job requests 32 cores. Cores can be selected from various nodes.

df -h 
ls 
echo $JOBDIR

#module load openmpi/3.1.3-intel 
module load openmpi/3.0.0-intel

EXE=/data/ommair/My_MD_MC_Codes/MD_FCC_VV_NBind_ES_3B_MPI_NVE_NVT_NPT/Src_NVE_NVT_NPT_090722_test12_etot/execute  #ommair
# /data/szalewicgrp/ommair/dl_poly_classic_test/dl_class_1.10/execute
#WRKDIR=/scratch/local/prayerz/$SLURM_JOB_ID

##WRKDIR=/scratch/local/ommair/$SLURM_JOB_ID
##ILEDIR=$PWD

#OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_NUM_THREADS

#export OMP_NUM_THREADS=32

##mkdir -p $WRKDIR


##cp $ILEDIR/* $WRKDIR/
#cp $ILEDIR/STATIS $WRKDIR/

##cd $WRKDIR/
##echo $HOSTNAME
##hostname
##which mpirun 

#mpirun -np ${OMP_NUM_THREADS} $EXE  
#mpirun -np 32 $EXE
mpirun $EXE
##rm slurm*.out
##cp OUTPUT $ILEDIR/
##cp REVCON $ILEDIR/
##cp HISTORY.xyz $ILEDIR/
##cp fort.401 $ILEDIR/
