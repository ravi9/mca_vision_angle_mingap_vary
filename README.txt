README for rp_anglepaper_mingap_orig.

Code Folder: rp_anglepaper_mingap_orig
Github URL: https://github.com/ravi9/mca_vision_angle_mingap_vary

Input Data Folder: Data

Different execution runs by varying the minimum gap between the nanomagnets is performed in this code. The code execution is the nanomagnetic simulation experiments of optimization by varying the magnetization angles.

Prerequisites:
1.	Access to CIRCE grid. To request account on CIRCE, see here: http://www.rc.usf.edu/
2.	Knowledge of submitting jobs to CIRCE computing grid. Beginner guide to submit jobs on CIRCE computing grid are available here: https://cwa.rc.usf.edu/projects/research-computing/wiki/Guide_to_Slurm


Steps to run the comparison:
1.  First choose the mingap parameter and to launch experiments run rp_qsubCmds_<mingap>.sh. For example to run the experiments with mingap 35, run rp_qsubCmds_35.sh

Execution Workflow:
The dataset (Data) contains 101 images as attribute files (.atts). The submission script (rp_qsubCmds_35.sh) launches 51 jobs. Each job has 2 images as input along with the sparsity value. Each job invokes script rp_bashRun2Images.sh.
The rp_bashRun2Images.sh script is used to submit the job for SLURM(CIRCE grid).
The rp_bashRun2Images.sh script loads Matlab, and launches the MATLAB script rp_gridMain_Run2Img.m

The Matlab script rp_gridMain_Run2Img.m launches the script: main.m, which runs the nanomagnetic simulation.

To summarize, the order of execution workflow for mingap 35:
1.	rp_ qsubCmds _35.sh
2.	rp_bashRun2Images.sh
3.	rp_gridMain_Run2Img.m
4.	main.m


Description of files:

1.	Results_mingap_<mingap parameter>: Result folder for a specific mingap parameter run
2.	Data: Contains the dataset. 101 Images in .atts format
3.	rp_qsubCmds__<mingap parameter>.sh: Shell script to launch the experiments using qsub command
4.	rp_bashRun2Images.sh: Shell script which loads Matlab, and launches the MATLAB script rp_gridMain_Run2Img.m
5.	rp_gridMain_Run2Img.m: matlab script which launches the script: main.m, which runs the nanomagnetic simulation.
6.	main.m: Matlab script which start execution of nanomagnetic simulation.
7.	Main_driver.m: Matlab script to run main.m for testing purposes on laptop.
8.	gridlogs: This folder contains the logs generated while running the jobs.

