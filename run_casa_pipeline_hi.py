
import os
import sys

from pathlib import Path
import numpy as np

data_path = Path("/scratch/public/sao/erickoch/M51_VLA/")

this_idx = int(sys.argv[-1])

all_lines = ['hi21cm_2p5kms', 'hi21cm_5kms', 'hi21cm_10kms', 'oh1665', 'oh1667']

this_line_product = all_lines[this_idx-1]

# Spectral coverage and resolution changed. Don't bother with some configs
# for certain lines
if all_lines[this_idx-1] in ['hi21cm_2p5kms', 'oh1665', 'oh1667']:
    this_config = ['A']
else:
    this_config = ['A+B+C+D', 'B+C+D', 'C+D', 'D', 'A']
    # this_config = ['B+C+D']

scripts_dir = '/home/erickoch/M51_HI/'

master_key = os.path.join(scripts_dir, 'keys/master_key.txt')

if not os.path.exists(master_key):
    raise ValueError("master_key does not exist at {}".format(master_key))

# Set the logging
from phangsPipeline import phangsLogger as pl
# reload(pl)
pl.setup_logger(level='DEBUG', logfile=None)
# Imports

#sys.path.insert(1, )
from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerVis as uvh
from phangsPipeline import handlerImaging as imh
from phangsPipeline import handlerPostprocess as pph

# Initialize key handler

this_kh = kh.KeyHandler(master_key = master_key)
this_uvh = uvh.VisHandler(key_handler = this_kh)
this_imh = imh.ImagingHandler(key_handler = this_kh)
this_pph = pph.PostProcessHandler(key_handler= this_kh)

this_uvh.set_dry_run(False)
this_uvh.set_targets(only=['m51'])


this_uvh.set_interf_configs(only=this_config)
this_imh.set_interf_configs(only=this_config)
this_pph.set_interf_configs(only=this_config)


these_line_products = [this_line_product]


this_uvh.set_line_products(only=these_line_products)
this_imh.set_line_products(only=these_line_products)
this_pph.set_line_products(only=these_line_products)

this_uvh.set_no_cont_products(True)
this_imh.set_no_cont_products(True)
this_pph.set_no_cont_products(True)


do_uvprocessing = False
do_imaging = True
do_postproc = False


for this_target, this_proj, this_array, this_obsnum in this_kh.loop_over_input_ms():
    print(this_kh.get_file_for_input_ms(target=this_target,
                                        project=this_proj,
                                        array_tag=this_array,
                                        obsnum=this_obsnum))

target_list = this_uvh.get_targets()
product_list = this_uvh.get_all_products()
config_list = this_uvh.get_interf_configs()

print(target_list)
print(product_list)
print(config_list)

just_projects = None
strict_config = False

# Our first loop goes over the individual measurement sets,
# splits, and continuum subtracts the data. At this stage we
# have no knowledge of configs except that selection may
# reduce the number of input measurement sets.

for this_target, this_project, this_array_tag, this_obsnum in \
        this_uvh._kh.loop_over_input_ms(
            target=target_list,
            config=config_list,
            project=just_projects,
            strict_config=strict_config):

    print(this_target, this_project, this_array_tag, this_obsnum)

    for this_product in product_list:
        print(this_product)

##############################################################################
# Stage the uv data
##############################################################################

# This step loses some archival data due to channel width

if do_uvprocessing:
    this_uvh.loop_stage_uvdata(do_copy=True, do_contsub=True,
                               do_extract_line=False, do_extract_cont=False,
                               do_remove_staging=False, overwrite=True,
                               statwt_line=True, statwt_cont=True,
                               strict_config=False)

    this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False,
                               do_extract_line=True, do_extract_cont=False,
                               do_remove_staging=False, overwrite=True,
                               statwt_line=True, statwt_cont=True,
                               strict_config=False)

    this_uvh.loop_stage_uvdata(do_copy=False, do_contsub=False,
                               do_extract_line=False, do_extract_cont=False,
                               do_remove_staging=True, overwrite=True,
                               statwt_line=True, statwt_cont=True,
                               strict_config=False)

# Imaging
if do_imaging:

    # Enable splitting to avoid memory failured with casa subcube
    # operations in tclean
    nchunks = 9

    these_parts = [f'm51_{ii}' for ii in range(1, nchunks+1)]
    for config_name in this_config:
        this_imh.set_interf_configs(only=[config_name])

        # Split the concatenated MS file:
        this_ms_file = data_path / f"imaging/m51/m51_{config_name}_{this_line_product}.ms"
        if not this_ms_file.exists():
            raise ValueError(f"Cannot find parent MS {this_ms_file}")

        msmd.open(str(this_ms_file))
        nchan = msmd.nchan(0)
        msmd.close()

        chan_per_chunk = int(np.ceil(nchan / nchunks))
        print(f"Splitting into approx chunks of {chan_per_chunk}")

        # We'll actually do the splitting with explicity numbering for each. This ensures
        # we don't accidentally miss a channel!
        chans = np.arange(nchan)
        chan_chunks = np.array_split(chans, nchunks)


        for ii, (this_part, these_chans) in enumerate(zip(these_parts,
                                                          chan_chunks)):
            print(f"On part {this_part}")

            if this_part != these_parts[-1]:
                these_chans = np.append(these_chans, these_chans.max()+1)

            these_chans_str = ";".join([str(num) for num in these_chans])
            print(these_chans_str)

            this_ms_file_part = data_path / f"imaging/m51/{this_part}_{config_name}_{this_line_product}.ms"
            if not this_ms_file_part.exists():
                split(vis=str(this_ms_file),
                      outputvis=str(this_ms_file_part),
                      spw=f"*:{these_chans_str}",
                      datacolumn='DATA')

            this_imh.set_targets(only=[this_part])

            # Don't both cleaning the OH:
            if "oh" in all_lines[this_idx-1]:
                this_imh.loop_imaging(do_all=False,
                                    recipe='phangsalma',
                                    do_export_to_fits=True,
                                    do_dirty_image=True,)
            else:
                this_imh.loop_imaging(do_all=True,
                                    recipe='phangsalma',)


# Postprocessing (in this case, just downsampling)
if do_postproc:

    raise ValueError("Need to add mosaicking along spectral dimension first.")

    this_pph.loop_postprocess(do_prep=True,
                            do_feather=False,
                            do_mosaic=False,
                            do_cleanup=True,
                            )
