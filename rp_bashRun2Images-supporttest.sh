#!/bin/bash
#$ -q devel
#$ -N matlab_MCAVision
#$ -l h_rt=00:58:59,pcpus=1
#$ -j y
#$ -o /work/r/ravi1/EMT_TEST/ravi_mca-vision-angle-paper-100img/gridlogs/output.$JOB_ID.log

module purge
module load apps/matlab/r2012b
cd /home/r/ravi1
matlab -nodisplay -r "cd /work/r/ravi1/EMT_TEST/ravi_mca-vision-angle-paper-100img; rp_gridMain_Run2Img('007.pgm.atts', '008.pgm.atts', 6)"
