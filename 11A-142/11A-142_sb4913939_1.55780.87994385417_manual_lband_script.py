
'''
Separate script to handle the different SPW setup for the first track.
'''

import os
import sys

import numpy as np

from casatasks import (listobs, gaincal, applycal,
                       flagdata, flagmanager, gencal, setjy,
                       bandpass, fluxscale, rmtables,
                       hanningsmooth, split)

from casaplotms import plotms


myvis = sys.argv[-1]

obs_dict = listobs(vis=myvis)

bpcal = '1331+305=3C286'
phasecal = 'J1327+4326'
target = "M51a"

all_cals = f"{bpcal},{phasecal}"


# plotants(vis=myvis)

myrefant = 'ea12,ea05,ea13,ea21'
# ea28,

# Priorcals
gencal(vis=myvis,
       caltable=f'{myvis}.antpos',
       caltype='antpos')

priorcals = []
if os.path.exists(f'{myvis}.antpos'):
    priorcals.append(f'{myvis}.antpos')

# Gain curve
# gencal(vis=myvis,
#        caltable=f'{myvis}.gc',
#        caltype='gc')
# priorcals.append(f'{myvis}.gc')

# Requantizer gains
# gencal(vis=myvis,
#        caltable=f'{myvis}.rq',
#        caltype='rq')
# priorcals.append(f'{myvis}.rq')

# Switched power
# gencal(vis=myvis,
#        caltable=f'{myvis}.swpow',
#        caltype='swpow')
# priorcals.append(f'{myvis}.swpow')

# Apply Hanning smoothing.
hanningsmooth(vis=myvis,
              outputvis=f"temphanning.ms",
              datacolumn='data')

os.system(f"rm -rf {myvis}")
os.system(f"mv temphanning.ms {myvis}")


# Flag spw edges:
flagmanager(vis=myvis, mode='save', versionname='inital_flags')

# Flag unused C-band SPWs
flagdata(vis=myvis, mode='manual', spw='0~1', flagbackup=False)
flagmanager(vis=myvis, mode='save', versionname='flag_cband', comment='Flag Cband')


# 128 channels with 3.906 kHz width in spws 2~33
all_spws = "2~33"
flagdata(vis=myvis, mode='manual', spw=f'{all_spws}:0~10', flagbackup=False)
flagdata(vis=myvis, mode='manual', spw=f'{all_spws}:117~127', flagbackup=False)

flagmanager(vis=myvis, mode='save', versionname='flag_edges', comment='Flag edges')

# Duck says
flagdata(vis=myvis, mode='quack', quackinterval=5.0, quackmode='beg', flagbackup=False)
flagmanager(vis=myvis, mode='save', versionname='quack', comment='Quack 5s beg')

# Manual flags
try:
    flag_file = f"/home/erickoch/M51_HI/11A-142/manual_flags/{myvis.split('.ms')[0]}_flags.txt"
    if not os.path.exists(flag_file):
        raise ValueError(f"Flag file {flag_file} not found")
    flagdata(vis=myvis, mode='list', flagbackup=False,
            inpfile=flag_file)
    flagmanager(vis=myvis, mode='save', versionname='manual_flags',
            comment=f'{myvis.split(".ms")[0]}_flags.txt')
except ValueError as exc:
      print(exc)


# setjy
setjy(vis=myvis, field=bpcal)


# Initial phase solution:
rmtables(f'{myvis}.finalBPinitialgain.tbl')
gaincal(vis=myvis,
        caltable=f'{myvis}.finalBPinitialgain.tbl',
        field=all_cals,
        refant=myrefant,
        spw='',
        gaintype='G',
        calmode='p',
        solint='int',
        minsnr=5,
        gaintable=priorcals)

os.system(f"rm phasecal_init*.png*")
plotms(vis=f'{myvis}.finalBPinitialgain.tbl',
       xaxis='time',
       yaxis='phase',
       coloraxis='corr',
       iteraxis='antenna',
       plotrange=[-1,-1,-180,180],
       dpi=200, width=2048, height=2048,
       plotfile='phasecal_init.png')

# Delay
rmtables(f'{myvis}.finaldelay.tbl')
gaincal(vis=myvis,
        caltable=f'{myvis}.finaldelay.tbl',
        field=bpcal,
        scan='3~7',
        refant=myrefant,
        gaintype='K',
        solint='inf',
        combine='scan',
        minsnr=5,
        gaintable=priorcals + [f'{myvis}.finalBPinitialgain.tbl'])

os.system(f"rm delaycal*.png*")
plotms(vis=f'{myvis}.finaldelay.tbl',
       xaxis='antenna1',
       yaxis='delay',
       coloraxis='baseline',
       plotfile='delaycal.png', dpi=200, width=2048, height=2048)




# Bandpass
# NOTE: using the phase cal for bandpass as it's unresolved and gives decent solutions.

# Poor BP amps
# flagdata(vis=myvis, mode='manual', antenna='ea06', spw='0~8', flagbackup=False)


bandpass(vis=myvis,caltable=f'{myvis}.finalBPcal.tbl',
         field=bpcal,
         spw='',
         scan='3~7',
         refant=myrefant,
         combine='scan',
         solint='inf',
         bandtype='B',
         gaintable=priorcals+[f'{myvis}.finalBPinitialgain.tbl',
                              f'{myvis}.finaldelay.tbl'])

os.system(f"rm bpcal_amp*.png* bpcal_phase*.png*")

plotms(vis=f"{myvis}.finalBPcal.tbl",
       xaxis='chan',
       yaxis='amp',
       coloraxis='corr',
       iteraxis='spw',
       gridrows=2,
       gridcols=2,
       plotfile='bpcal_amp.png', dpi=200, width=2048, height=2048)

plotms(vis=f"{myvis}.finalBPcal.tbl",
       xaxis='chan',
       yaxis='phase',
       coloraxis='corr',
       plotrange=[-1,-1,-180,180],
       iteraxis='spw',
       gridrows=2,
       gridcols=2,
       plotfile='bpcal_phase.png', dpi=200, width=2048, height=2048)


# Gain cal

# Short int
gaincal(vis=myvis,
        caltable=f'{myvis}.phaseshortgaincal.tbl',
        field=all_cals,
        spw='',
        solint='15s',
        refant=myrefant,
        gaintype='G',
        calmode='p',
        solnorm=False,
        gaintable=[f'{myvis}.finaldelay.tbl',
                   f'{myvis}.finalBPcal.tbl'] + priorcals,
        interp=['', 'nearest'] + ['']*len(priorcals))

os.system(f"rm phasecal_short*.png*")
plotms(vis=f'{myvis}.phaseshortgaincal.tbl',
       xaxis='time',
       yaxis='phase',
       gridrows=1,gridcols=2,
       iteraxis='antenna',
       coloraxis='spw',
       plotrange=[-1,-1,-180,180],
       plotfile='phasecal_short.png', dpi=200, width=2048, height=2048)

# Long int
gaincal(vis=myvis,
        caltable=f'{myvis}.finalphasegaincal.tbl',
        field=all_cals,
        spw='',
        solint='120s',
        refant=myrefant,
        gaintype='G',
        calmode='p',
        solnorm=False,
        gaintable=[f'{myvis}.finaldelay.tbl',
                   f'{myvis}.finalBPcal.tbl'] + priorcals,
        interp=['', 'nearest'] + ['']*len(priorcals))

os.system(f"rm phasecal_long*.png*")

plotms(vis=f'{myvis}.finalphasegaincal.tbl',
       xaxis='time',
       yaxis='phase',
       gridrows=1,gridcols=2,
       iteraxis='antenna',
       coloraxis='spw',
       plotrange=[-1,-1,-180,180],
       plotfile='phasecal_long.png', dpi=200, width=2048, height=2048)

# Amp
gaincal(vis=myvis,
        caltable=f'{myvis}.finalampgaincal.tbl',
        field=all_cals,
        spw='',
        solint='300s',
        refant=myrefant,
        gaintype='G',
        calmode='a',
        solnorm=False,
        gaintable=[f'{myvis}.finaldelay.tbl',
                   f'{myvis}.finalBPcal.tbl',
                   f'{myvis}.phaseshortgaincal.tbl'] + priorcals,
        interp=['', 'nearest', 'linear'] + ['']*len(priorcals))


os.system(f"rm gaincal_amp*.png*")

plotms(vis=f'{myvis}.finalampgaincal.tbl',
       xaxis='time',
       yaxis='amp',
       gridrows=1,gridcols=2,
       iteraxis='spw',
       coloraxis='antenna1',
       plotfile='gaincal_amp.png', dpi=200, width=2048, height=2048)


# Fluxscale:
fluxresults = fluxscale(vis=myvis,
                        caltable=f'{myvis}.finalampgaincal.tbl',
                        refspwmap=[-1],
                        transfer=phasecal,
                        fluxtable=f'{myvis}.finalfluxgaincal.tbl',
                        reference=bpcal,
                        fitorder=1,
                        incremental=False)

fluxdict_numpyfile = "{0}.fluxresults.npy".format(myvis)
if os.path.exists(fluxdict_numpyfile):
    os.remove(fluxdict_numpyfile)

np.save(fluxdict_numpyfile,
        fluxresults, allow_pickle=True)

# In CASA
applycal(vis=myvis,
         field=bpcal,
         gaintable=[f'{myvis}.finalfluxgaincal.tbl',
                    f'{myvis}.phaseshortgaincal.tbl',
                    f'{myvis}.finaldelay.tbl',
                    f'{myvis}.finalBPcal.tbl',] + priorcals,
         gainfield=[bpcal, bpcal, '', ''] + ['']*len(priorcals),
         interp=['linear', 'nearest','',''] + ['']*len(priorcals),
         calwt=False)

applycal(vis=myvis,
         field=phasecal,
         gaintable=[f'{myvis}.finalfluxgaincal.tbl',
                    f'{myvis}.phaseshortgaincal.tbl',
                    f'{myvis}.finaldelay.tbl',
                    f'{myvis}.finalBPcal.tbl',] + priorcals,
         gainfield=[phasecal, phasecal,'',''] + ['']*len(priorcals),
         interp=['linear', 'nearest','',''] + ['']*len(priorcals),
         calwt=False)

applycal(vis=myvis,
         field=target,
         gaintable=[f'{myvis}.finalfluxgaincal.tbl',
                    f'{myvis}.finalphasegaincal.tbl',
                    f'{myvis}.finaldelay.tbl',
                    f'{myvis}.finalBPcal.tbl'] + priorcals,
         gainfield=[phasecal, phasecal, '', ''] + ['']*len(priorcals),
         interp=['linear', 'linear','',''] + ['']*len(priorcals),
         calwt=False)

# statwt(vis=myvis, datacolumn='data')

# Split out the science target
split(vis=myvis, outputvis=f'{myvis}.target',
      field=target,
      datacolumn='corrected', keepflags=False)

# This CASA pipescript is meant for use with CASA 6.5.4 and pipeline 2023.1.0.124
# context = h_init()
# context.set_state('ProjectSummary', 'observatory', 'Karl G. Jansky Very Large Array')
# context.set_state('ProjectSummary', 'telescope', 'EVLA')
# try:
#     hifv_importdata(vis=[myvis],
#                     datacolumns={'data': 'raw',
#                                  'corrected': 'regcal_contline_all'})
#     hif_makeimlist(intent='PHASE,BANDPASS', specmode='cont')
#     hif_makeimages(hm_masking='centralregion')
#     hifv_exportdata()
#     hifv_flagtargetsdata()
#     hif_mstransform()
#     hif_checkproductsize(maximsize=16384)
#     hif_makeimlist(specmode='cont', datatype='regcal')
#     hif_makeimages(hm_cyclefactor=3.0, hm_nsigma=5.0)
#     # hif_selfcal()
#     # hif_makeimlist(specmode='cont', datatype='selfcal')
#     # hif_makeimages(hm_cyclefactor=3.0)
#     # hifv_pbcor()
#     #hifv_exportdata(imaging_products_only=True)
# finally:
#     h_save()
