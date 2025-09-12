
'''
Add intents to pre-EVLA data sets.
This is needed to run the PHANGS pipeline.
'''

from pathlib import Path

data_path = Path('/reduction/erickoch/M51/VLA/calibrated_ms/')


for myvis in (data_path / "AW605").glob('AW605*'):
    print(myvis.name)

    # Add intents
    defintent(vis=str(myvis),
              outputvis=str(myvis),
              intent='OBSERVE_TARGET',
              mode='set')

    # For C and D configs, combinespws and split out M51
    if "_C.ms" in myvis.name or "_D.ms" in myvis.name:
        mstransform(vis=str(myvis),
                    outputvis=str(data_path / f"{myvis.name}.target"),
                    combinespws=True,
                    field='NGC5194')
