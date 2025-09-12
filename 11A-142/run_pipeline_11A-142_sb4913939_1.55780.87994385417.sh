
track_name=11A-142_sb4913939_1.55780.87994385417
echo Starting $track_name at $(date)


project_folder="/reduction/erickoch/M51/VLA/11A-142/11A-142/"

cd $project_folder

this_casa="/py3opt/casa-6.5.4-9-pipeline-2023.1.0.124/bin/casa"

conda activate pyuvdata_ewk

###############
ms_name=${track_name}.ms

cd $track_name

tar -xf ${ms_name}.tgz
tar -xf weblog.tgz
timestamp=$(date +%Y-%m-%d_%H-%M-%S)
mkdir run_${timestamp}
mv ${ms_name} run_${timestamp}
cd run_${timestamp}

$this_casa -c ~/M51_HI/11A-142/11A-142_sb4913939_1.55780.87994385417_manual_lband_script.py ${ms_name}

# Make additonal QA plots:
mkdir -p products

$this_casa -c ~/M51_HI/11A-142//continuum_qaproducts.py ${ms_name}

cd products
# Copy the initial pipeline MS import from the NRAO job
cp -r $project_folder/$track_name/pipeline-2024* weblog

ipython -c "import qaplotter; qaplotter.make_all_plots(msname='${ms_name}')"

cd $project_folder

echo Finished $track_name at $(date)
