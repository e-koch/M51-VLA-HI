
project_folder="/reduction/erickoch/M51/VLA/11A-142/11A-142/"

cd $project_folder

this_casa="/py3opt/casa-6.5.4-9-pipeline-2023.1.0.124/bin/casa"

conda activate pyuvdata_ewk

###############
# Track 1

# NOTE: needs custom script because the SPW setup is different than the other tracks.

# track_name=11A-142_sb4913939_1.55780.87994385417
# echo $track_name

# ms_name=${track_name}.ms

# cd $track_name

# tar -xf ${ms_name}.tgz
# timestamp=$(date +%Y-%m-%d_%H-%M-%S)
# mkdir run_${timestamp}
# mv ${ms_name} run_${timestamp}
# cd run_${timestamp}

# $this_casa -c ~/M51_HI/11A-142/11A-142_manual_lband_script.py ${ms_name}
# mkdir -p products
# # Make additonal QA plots:
# $this_casa -c ~/M51_HI/11A-142//continuum_qaproducts.py ${ms_name}

# cd products
# ipython -c "import qaplotter; qaplotter.make_all_plots(msname='${ms_name}')"


# cd $project_folder

###############
# Track 2-12

cat ~/M51_HI/11A-142/track_names.txt |  parallel -n 4 ~/M51_HI/11A-142/run_pipeline_11A-142.sh
