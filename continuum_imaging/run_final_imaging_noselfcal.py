
from pathlib import Path
import os
import sys

target = sys.argv[-3]

orig_data_path = Path(sys.argv[-2])

data_path = Path(sys.argv[-1])

search_str = "11A-142*target"

my_vis = []

cont_spws = '3~7'

for this_vis in orig_data_path.glob(search_str):
       output_vis_name = str(data_path / f"{this_vis.name}.cont")
       split(vis=str(this_vis), outputvis=output_vis_name,
             spw=cont_spws, datacolumn='DATA')

       my_vis.append(output_vis_name)

if len(my_vis) == 0:
     raise ValueError(f"No MSs found with {search_str}")

# my_weighting = "uniform"

my_weighting = "briggs"
my_robust = 0.5

weight_str = f"{my_weighting}"
if "briggs" in weight_str:
    weight_str += f"_{my_robust}".replace(".", "p")

if my_weighting == "briggs" or my_weighting == "natural":
       my_imsize = 16384
       my_cell_val = 0.2
else:
       my_imsize = 16384
       my_cell_val = 0.15

my_cell = f'{my_cell_val}arcsec'

my_imagename = f"{target}_A_L_{weight_str}"

# Stage 1: Clean bright sources. Single scale
nsigma_0 = 20.

tclean(vis=my_vis,
       field='', spw='',
       datacolumn='corrected',
       imagename=my_imagename,
       imsize=my_imsize,
       cell=my_cell,
       specmode='mfs',
       reffreq='1.5GHz',
       gridder='wproject',
       wprojplanes=256,
       pblimit=-0.1,
       normtype='flatnoise',
       deconvolver='mtmfs',
       scales=[0],
       nterms=2,
       smallscalebias=0.6,
       restoration=True,
       pbcor=False,
       weighting=my_weighting,
       robust=my_robust,
       cycleniter=2000,
       cyclefactor=5.0,
       niter=10000,
       gain=0.1,
       nsigma=nsigma_0,
       nmajor=-1,
       usemask='auto-multithresh',
       mask='',
       pbmask=0.0,
       sidelobethreshold=4.0,
       noisethreshold=15.0,
       lownoisethreshold=5.0,
       negativethreshold=0.0,
       smoothfactor=1.0,
       minbeamfrac=0.8,
       cutthreshold=0.01,
       growiterations=10,
       dogrowprune=True,
       minpercentchange=1.0,
       verbose=True,
       fastnoise=True,
       restart=True,
       savemodel='none',
       calcres=True,
       calcpsf=True,
       psfcutoff=0.35,
       parallel=False)


# Stage 2: Deep MS clean of the central source.

# Need to remove the previous mask
os.system(f"mv {my_imagename}.mask {my_imagename}.mask_stage1")

# Backup stage 1
for this_image in Path(".").glob(f"{my_imagename}.*"):
    os.system(f"cp -r {this_image} {this_image}_stage1")


nsigma_1 = 3.


# Beams are ~1".
nbeams = [0, 1, 2]
pix_per_beam = 5
scales_1 = [int(beam * pix_per_beam * my_cell_val) for beam in nbeams]

# Define the central box for the deeper MS clean
y = my_imsize // 2
x = my_imsize // 2

# 8'
ang_width_arcmin = 8
ang_width_pix = int((ang_width_arcmin * 60) / my_cell_val)

region_filename = f'centerbox[[{x}pix, {y}pix], [{ang_width_pix}pix, {ang_width_pix}pix]]'


tclean(vis=my_vis,
       field='', spw='',
       datacolumn='corrected',
       imagename=my_imagename,
       imsize=my_imsize,
       cell=my_cell,
       specmode='mfs',
       reffreq='1.5GHz',
       gridder='wproject',
       wprojplanes=256,
       pblimit=-0.1,
       normtype='flatnoise',
       deconvolver='mtmfs',
       scales=scales_1,
       nterms=2,
       smallscalebias=0.6,
       restoration=True,
       pbcor=False,
       weighting=my_weighting,
       robust=my_robust,
       cycleniter=2000,
       cyclefactor=3.0,
       niter=10000,
       gain=0.1,
       nsigma=nsigma_1,
       nmajor=-1,
       usemask='user',
       mask=region_filename,
       pbmask=0.0,
       verbose=True,
       fastnoise=True,
       restart=True,
       savemodel='none',
       calcres=False,
       calcpsf=False,
       psfcutoff=0.35,
       parallel=False)


# Stage 3: Deep single-scale clean of the central PB.

# Need to remove the previous mask
os.system(f"mv {my_imagename}.mask {my_imagename}.mask_stage2")

# Backup stage 1
for this_image in Path(".").glob(f"{my_imagename}.*"):
    if "stage1" in str(this_image):
        continue
    os.system(f"cp -r {this_image} {this_image}_stage2")


nsigma_2 = 3.
scales_2 = [0]


tclean(vis=my_vis,
       field='', spw='',
       datacolumn='corrected',
       imagename=my_imagename,
       imsize=my_imsize,
       cell=my_cell,
       specmode='mfs',
       reffreq='1.5GHz',
       gridder='wproject',
       wprojplanes=256,
       pblimit=-0.1,
       normtype='flatnoise',
       deconvolver='mtmfs',
       scales=scales_2,
       nterms=2,
       smallscalebias=0.6,
       restoration=True,
       pbcor=False,
       weighting=my_weighting,
       robust=my_robust,
       cycleniter=2000,
       cyclefactor=4.0,
       niter=10000,
       gain=0.1,
       nsigma=nsigma_2,
       nmajor=-1,
       usemask='pb',
       mask='',
       pbmask=0.5,
       verbose=True,
       fastnoise=True,
       restart=True,
       savemodel='none',
       calcres=False,
       calcpsf=False,
       psfcutoff=0.35,
       parallel=False)


widebandpbcor(my_vis[0],
              imagename=my_imagename,
              nterms=2,
              threshold='5e-5Jy',
              action='pbcor',
              reffreq='1.5GHz',
              pbmin=0.2,
              field='',
              spwlist=[3, 4, 5, 6, 7],
              chanlist=[64] * 5,
              weightlist=[1.0] * 5)
