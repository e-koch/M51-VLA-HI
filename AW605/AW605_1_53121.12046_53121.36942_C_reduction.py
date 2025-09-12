
importvla("AW605_1_53121.12046_53121.36942.exp",
          vis='AW605_1_53121.12046_53121.36942.ms')

myvis='AW605_1_53121.12046_53121.36942.ms'



caltab_plot_path = 'calibration_plots'
if not os.path.exists(caltab_plot_path):
    os.mkdir(caltab_plot_path)

out_root = myvis[:-3]

flagmanager(vis=myvis, mode='save', versionname='imported',
            comment='flagging after import', merge='replace')

# Systematic flagging:

flagdata(vis=myvis, mode='quack', quackinterval=30, flagbackup=False)
flagdata(vis=myvis, mode='shadow', flagbackup=False)
flagmanager(vis=myvis, mode='save', versionname='quack_shadow')


# Find refant
plotants(vis=myvis)

# For M51, we only want spw 2, 3.
# Use spw 0,1 for M101.

myvis_m51 = 'AW605_1_53121.12046_53121.36942_m51.ms'

split(vis=myvis, outputvis=myvis_m51, spw='2,3')

myvis = myvis_m51

ref_ant = "VA11,VA15,VA15,VA04,VA12,VA25"

# Set fieldnames
fluxcal = '1331+305'
bpcal = '1331+305'
phasecal = '1252+565'
targets = 'NGC5194'

allcals = ",".join([fluxcal, phasecal])

setjy(vis=myvis, field=fluxcal)


# Phase self-cal on BP calibrator
phaseshortgaincal_table = out_root+".intphase_bpcal.gcal"
os.system("rm -rf "+phaseshortgaincal_table)

gaincal(vis=myvis, field=bpcal, caltable=phaseshortgaincal_table,
        append=False, calmode='p',
        solint='int', minsnr=2.0, refant=ref_ant)

plotms(vis=phaseshortgaincal_table,
       xaxis='time',
       yaxis='phase',
       coloraxis='spw',
       iteraxis='antenna',
       ydatacolumn='data',
       plotfile="{0}/{1}.png".format(caltab_plot_path,
                                     phaseshortgaincal_table),
       gridrows=2, gridcols=2,
       yselfscale=True,
       xconnector='line', timeconnector=True,
       showgui=False, overwrite=True, dpi=400,
       exprange='all')

bandpass_table = out_root+".bpcal"
os.system("rm -rf "+bandpass_table)
bandpass(vis=myvis, field=bpcal, selectdata=False,
         caltable=bandpass_table,
         gaintable=phaseshortgaincal_table,
         solint='inf', solnorm=True,
         refant=ref_ant, bandtype='B')

# Plot bandpass phase and amplitude solutions.
plotms(vis=bandpass_table, xaxis='freq',
       yaxis='phase',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
    #    plotfile=f"{caltab_plot_path}/{bandpass_table}.phase.png",
       plotfile="{0}/{1}.phase.png".format(caltab_plot_path, bandpass_table),
       gridrows=2, gridcols=2,
       yselfscale=True,
       showgui=False, overwrite=True, dpi=400, exprange='all')
plotms(vis=bandpass_table, xaxis='freq',
       yaxis='amp',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
    #    plotfile=f"{caltab_plot_path}/{bandpass_table}.amp.png",
       plotfile="{0}/{1}.amp.png".format(caltab_plot_path, bandpass_table),
       gridrows=2, gridcols=2,
       yselfscale=True,
       showgui=False, overwrite=True, dpi=400, exprange='all')


gain_phase_scan_table = out_root+".scanphase.gcal"
os.system("rm -rf "+gain_phase_scan_table)

gaincal(vis=myvis, field=allcals, selectdata=True,
        caltable=gain_phase_scan_table, append=False,
        gaintable=bandpass_table, gainfield=bpcal,
        calmode='p',solint='inf',
        minsnr=2.0, minblperant=2, refant=ref_ant,
        spwmap=[])

plotms(vis=gain_phase_scan_table,xaxis='time',
       yaxis='phase',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       xconnector='line', timeconnector=True,
       plotfile="{0}/{1}.png".format(caltab_plot_path,
                                     gain_phase_scan_table),
       gridrows=4, gridcols=2,
       yselfscale=True, showgui=False, overwrite=True, dpi=400, exprange='all')


gain_phase_int_table = out_root+".intphase.gcal"
os.system("rm -rf "+gain_phase_int_table)

gaincal(vis=myvis, field=allcals, selectdata=True,
        caltable=gain_phase_int_table, append=False,
        gaintable=bandpass_table, gainfield=bpcal,
        calmode='p',solint='int',
        minsnr=2.0, minblperant=2, refant=ref_ant,
        spwmap=[])

plotms(vis=gain_phase_int_table,xaxis='time',
       yaxis='phase',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       xconnector='line', timeconnector=True,
       plotfile="{0}/{1}.png".format(caltab_plot_path,
                                     gain_phase_int_table),
       gridrows=4, gridcols=2,
       yselfscale=True, showgui=False, overwrite=True, dpi=400, exprange='all')


gain_amp_scan_table = out_root+".amp.gcal"
os.system("rm -rf "+gain_amp_scan_table)
gaincal(vis=myvis, field=allcals, selectdata=True,
        caltable=gain_amp_scan_table, append=False,
        gaintable=[bandpass_table, gain_phase_int_table],
        gainfield=[bpcal, fluxcal+','+phasecal],
        solint='inf', minsnr=2.0, minblperant=2,
        calmode='a',refant=ref_ant,
        spwmap = [])

plotms(vis=gain_amp_scan_table,xaxis='time',
       yaxis='amp',
       coloraxis='spw',iteraxis='antenna',ydatacolumn='data',
       xconnector='line', timeconnector=True,
    #    plotfile=f"{caltab_plot_path}/{gain_amp_scan_table}.png",
       plotfile="{0}/{1}.png".format(caltab_plot_path, gain_amp_scan_table),
       gridrows=4, gridcols=2,
       yselfscale=True, showgui=False, overwrite=True, dpi=400, exprange='all')

fluxboot_table = out_root+".fcal"
os.system("rm -rf "+fluxboot_table)
# fluxscale(vis=myvis, caltable=gain_amp_scan_table, fluxtable=fluxboot_table,
#           transfer=phasecal, reference=fluxcal,
#           refspwmap=[])
fluxresults = fluxscale(vis=myvis, caltable=gain_amp_scan_table,
                        refspwmap=[-1], transfer=phasecal,
                        fluxtable=fluxboot_table, reference=fluxcal,
                        fitorder=1)


flagmanager(vis=myvis,mode='save',versionname='beforeapplycal')

# Now apply the solutions:
interpmode = ['linear', 'linear', 'nearest']

# ... the source
applycal(vis=myvis, field=targets,
            gaintable=[fluxboot_table, gain_phase_scan_table, bandpass_table],
            interp=interpmode, gainfield=[phasecal,phasecal,bpcal],
            spwmap=[])

# ... the phase calibrator
applycal(vis=myvis, field=phasecal,
            gaintable=[fluxboot_table, gain_phase_scan_table, bandpass_table],
            interp=interpmode, gainfield=[phasecal,phasecal,bpcal],
            spwmap=[])

# ... the flux calibrator
applycal(vis=myvis, field=fluxcal,
            gaintable=[fluxboot_table, gain_phase_int_table, bandpass_table],
            interp=interpmode, gainfield=[fluxcal,fluxcal,bpcal],
            spwmap=[[],[],[]])


# Additional manual flagging

# Flag SPW edges
flagdata(vis=myvis, spw="*:0~4", flagbackup=False)
flagdata(vis=myvis, spw="*:58~62", flagbackup=False)

# Poor antenna in 2nd science scan
flagdata(vis=myvis, field=targets, scan='8', antenna='VA26', flagbackup=False)
flagdata(vis=myvis, field=targets, scan='8', antenna='VA16&&VA23', flagbackup=False)
flagdata(vis=myvis, field=targets, scan='8', antenna='VA20', flagbackup=False)

flagdata(vis=myvis, scan='8',
         timerange='04:41:20~04:42:30',
         flagbackup=False)

flagdata(vis=myvis, scan='12',
         timerange='05:37:30~05:38:30',
         flagbackup=False)

# Consistently high amps on targets
flagdata(vis=myvis, field=targets, antenna='VA02,VA04', flagbackup=False)

flagmanager(vis=myvis, versionname='manual_flagging',
            mode='save')


# Diagnostic visibility plots

plot_path = "calibrator_plots"
if not os.path.exists(plot_path):
    os.mkdir(plot_path)

for myfield in allcals.split(","):

    # Amp/Phase vs freq

    plotms(vis=myvis,xaxis='freq',
            yaxis='amp',field=myfield,avgtime='1e8',avgscan=True,
            coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
            gridrows=4, gridcols=2, showgui=False, overwrite=True, dpi=400, exprange='all',
            yselfscale=True,
            plotfile="{0}/{1}.{2}.aftercal.amp_freq.png".format(plot_path,
                                                                myvis,
                                                                myfield.lower()))

    plotms(vis=myvis,xaxis='freq',
            yaxis='phase',field=myfield,avgtime='1e8',avgscan=True,
            coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
            gridrows=4, gridcols=2, showgui=False, overwrite=True, dpi=400, exprange='all',
            yselfscale=True,
            plotfile="{0}/{1}.{2}.aftercal.pha_freq.png".format(plot_path,
                                                                myvis,
                                                                myfield.lower()))

    # Amp/Phase vs time

    plotms(vis=myvis,xaxis='time',
            yaxis='amp', field=myfield, avgchannel='1e8', avgscan=False,
            coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
            gridrows=4, gridcols=2, showgui=False, overwrite=True, dpi=400, exprange='all',
            yselfscale=True,
            plotfile="{0}/{1}.{2}.aftercal.amp_time.png".format(plot_path,
                                                                myvis,
                                                                myfield.lower()))

    plotms(vis=myvis,xaxis='time',
            yaxis='phase',field=myfield, avgchannel='1e8', avgscan=False,
            coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
            gridrows=4, gridcols=2, showgui=False, overwrite=True, dpi=400, exprange='all',
            yselfscale=True,
            plotfile="{0}/{1}.{2}.aftercal.pha_time.png".format(plot_path,
                                                                myvis,
                                                                myfield.lower()))

    # Amp/Phase vs uvdist

    plotms(vis=myvis,xaxis='uvdist',
            yaxis='amp', field=myfield, avgchannel='8', avgscan=False, avgtime='1e8',
            coloraxis='ant1',iteraxis='spw',ydatacolumn='corrected',
            gridrows=4, gridcols=3, showgui=False, overwrite=True, dpi=400, exprange='all',
            yselfscale=True,
            plotfile="{0}/{1}.{2}.aftercal.amp_uvdist.png".format(plot_path,
                                                                  myvis,
                                                                  myfield.lower()))

    plotms(vis=myvis,xaxis='uvdist',
            yaxis='phase',field=myfield, avgchannel='8', avgscan=False, avgtime='1e8',
            coloraxis='ant1',iteraxis='spw',ydatacolumn='corrected',
            gridrows=1, gridcols=3, showgui=False, overwrite=True, dpi=400, exprange='all',
            yselfscale=True,
            plotfile="{0}/{1}.{2}.aftercal.pha_uvdist.png".format(plot_path,
                                                                  myvis,
                                                                  myfield.lower()))

    # Real vs Imag
    plotms(vis=myvis,xaxis='imag',
            yaxis='real',field=myfield, avgchannel='8', avgscan=False, avgtime='1e8',
            coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
            gridrows=4, gridcols=2, showgui=False, overwrite=True, dpi=400, exprange='all',
            plotfile="{0}/{1}.{2}.aftercal.amp_pha.png".format(plot_path,
                                                               myvis,
                                                               myfield.lower()))

    # Amp/Phase vs baseline

    plotms(vis=myvis,xaxis='baseline',
            yaxis='amp', field=myfield, avgchannel='1e8', avgscan=False, avgtime='1e8',
            coloraxis='spw',iteraxis=None, ydatacolumn='corrected',
            gridrows=1, gridcols=1, showgui=False, overwrite=True, dpi=400, exprange='all',
            plotfile="{0}/{1}.{2}'.aftercal.amp_baseline.png".format(plot_path,
                                                                     myvis,
                                                                     myfield.lower()))

    plotms(vis=myvis,xaxis='baseline',
            yaxis='phase',field=myfield, avgchannel='1e8', avgscan=False, avgtime='1e8',
            coloraxis='spw',iteraxis=None, ydatacolumn='corrected',
            gridrows=1, gridcols=1, showgui=False, overwrite=True, dpi=400, exprange='all',
            plotfile="{0}/{1}.{2}'.aftercal.pha_baseline.png".format(plot_path,
                                                                     myvis,
                                                                     myfield.lower()))



scitarget_plot_path = "target_plots"
if not os.path.exists(scitarget_plot_path):
    os.mkdir(scitarget_plot_path)


for myfield in targets.split(","):

    # Amp vs freq
    plotms(vis=myvis,xaxis='freq',
            yaxis='amp',field=myfield,avgtime='1e8',avgscan=True,
            coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
            gridrows=4, gridcols=2, showgui=False, overwrite=True, dpi=400, exprange='all',
            yselfscale=True,
            # plotfile=f"{scitarget_plot_path}/{myvis}.{myfield.lower()}.aftercal.amp_freq.png")
            plotfile="{0}/{1}.{2}.aftercal.amp_freq.png".format(scitarget_plot_path,
                                                                myvis,
                                                                myfield.lower()))

    # Amp vs time
    plotms(vis=myvis,xaxis='time',
            yaxis='amp', field=myfield, avgchannel='1e8', avgscan=False,
            coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
            gridrows=4, gridcols=2, showgui=False, overwrite=True, dpi=400, exprange='all',
            yselfscale=True,
            # plotfile=f"{scitarget_plot_path}/{myvis}.{myfield.lower()}.aftercal.amp_time.png")
            plotfile="{0}/{1}.{2}.aftercal.amp_time.png".format(scitarget_plot_path,
                                                                myvis,
                                                                myfield.lower()))

    # Amp vs uvdist
    plotms(vis=myvis,xaxis='uvdist',
            yaxis='amp', field=myfield, avgchannel='1e8', avgscan=False, avgtime='1e8',
            coloraxis='spw',iteraxis='antenna',ydatacolumn='corrected',
            gridrows=1, gridcols=3, showgui=False, overwrite=True, dpi=400, exprange='all',
            yselfscale=True,
            # plotfile=f"{scitarget_plot_path}/{myvis}.{myfield.lower()}.aftercal.amp_uvdist.png")
            plotfile="{0}/{1}.{2}.aftercal.amp_uvdist.png".format(scitarget_plot_path,
                                                                  myvis,
                                                                  myfield.lower()))


# Check if more flagging is needed with a test cube
os.system('rm -r test_imaging/{}_HI_dirty.*'.format(targets))

tclean(vis=myvis, spw='*', specmode='cube', field=targets,
       imagename='test_imaging/{}_HI_dirty'.format(targets),
       cell='3arcsec', imsize=512,
       start=1, nchan=-1, width=1,
       niter=0)


# Split out the target and regrid to a single SPW
mstransform(vis=myvis,outputvis=f'{myvis}.target',
            datacolumn='CORRECTED',field=targets,
            spw='*:4~57', combinespws=True, regridms=True,
            nchan=-1, width=1, outframe='LSRK', veltype='radio',)

tclean(vis=f'{myvis}.target', spw='', specmode='cube', field=targets,
       imagename='test_imaging/{}_HI_dirty_regrid'.format(targets),
       cell='3arcsec', imsize=512,
       start=1, nchan=-1, width=1,
       niter=0)


# uvcontsub
out = uvcontsub(vis=f"{myvis}.target", outputvis=f"{myvis}.target.contsub",
                fitspec='0:0~30;90~99')

tclean(vis=f"{myvis}.target.contsub", spw='', specmode='cube', field=targets,
       imagename='test_imaging/{}_HI_dirty_contsub'.format(targets),
       cell='3arcsec', imsize=512,
       start=1, nchan=-1, width=1,
       niter=0)
