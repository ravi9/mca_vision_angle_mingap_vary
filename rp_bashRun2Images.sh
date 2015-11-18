#!/bin/bash
#$ -q devel
#$ -N matlab_MCAVision
#$ -l h_rt=00:58:59,pcpus=1
#$ -j y
#$ -o /work/r/ravi1/EMT/rp_anglepaper_mingap/gridlogs/output.$JOB_ID.log

module purge
module load apps/matlab/r2012b
cd /home/r/ravi1
matlab -nodisplay -r "cd /work/r/ravi1/EMT/rp_anglepaper_mingap; rp_gridMain_Run2Img('$1', '$2', $3)"
