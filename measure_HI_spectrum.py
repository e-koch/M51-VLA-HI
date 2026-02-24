
import re
import gc

import numpy as np
import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import Table
from regions import Regions
from spectral_cube import SpectralCube, OneDSpectrum

from pathlib import Path


alpha_HI = ((1.82e18 * u.cm**-2) * c.m_p).to(u.solMass / u.pc**2) / (u.K * u.km / u.s)


data_path = Path("/Users/ekoch/storage/M51/VLA")

# Distance used for mass conversion
dist = 8.6 * u.Mpc

# THINGS HI mass
hi_mass_db08 = 25.4e8 * u.solMass
dist_db08 = 8.0 * u.Mpc

# FEASTS HI mass (Wang+24)
hi_mass_w24 = 10**9.59 * u.solMass
dist_w24 = 8.0 * u.Mpc

# Read the region file
reg_path = data_path / "m51_hi_avg_spectrum_region.reg"
reg = Regions.read(str(reg_path), format="ds9")[0]

# Output directory for saved spectra
spec_out_dir = data_path / "avg_spectra"
spec_out_dir.mkdir(exist_ok=True)

# Distance-corrected reference masses
hi_mass_db08_corr = hi_mass_db08 * (dist / dist_db08) ** 2
hi_mass_w24_corr = hi_mass_w24 * (dist / dist_w24) ** 2

names = []
hi_masses = []
frac_db08s = []
frac_w24s = []

for folder in sorted(data_path.glob("m51_*")):
    if not folder.is_dir():
        continue

    # Native-resolution cube: exclude smoothed, noise, coverage, moment files
    all_fits = list(folder.glob("*.fits"))
    native_cubes = [
        f for f in all_fits
        if not re.search(
            r'_\d+as\.fits$|_noise\.fits$|_coverage|_feathered|_broad|_strict|mom\d?\.fits$|emom',
            f.name)
    ]

    if len(native_cubes) != 1:
        print(f"Skipping {folder.name}: "
              f"found {len(native_cubes)} native cubes: {[f.name for f in native_cubes]}")
        continue

    cube_file = native_cubes[0]
    print(f"Processing {cube_file.name} ...")

    cube = SpectralCube.read(cube_file)
    cube.allow_huge_operations = True
    cube = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio')

    # Header quantities needed for K → Jy conversion
    jytok = cube.header['JYTOK']               # K per (Jy/beam)
    bmaj_deg = cube.header['BMAJ']
    bmin_deg = cube.header['BMIN']
    cdelt_deg = abs(cube.header['CDELT2'])
    pixel_area_arcsec2 = (cdelt_deg * 3600.0) ** 2
    beam_area_arcsec2 = (np.pi / (4 * np.log(2))
                         * bmaj_deg * 3600.0 * bmin_deg * 3600.0)
    beams_per_pixel = pixel_area_arcsec2 / beam_area_arcsec2
    pixel_area_pc2 = (np.radians(cdelt_deg) * dist.to(u.pc).value) ** 2

    # Apply the broad mask before region selection
    broad_mask_files = [f for f in folder.glob("*_broadmask.fits")
                        if not re.search(r'_\d+as_broadmask\.fits$', f.name)]
    if len(broad_mask_files) != 1:
        print(f"  Warning: found {len(broad_mask_files)} broad mask files, skipping mask")
    else:
        broad_mask = fits.getdata(broad_mask_files[0]).astype(bool)
        cube = cube.with_mask(broad_mask)
        print(f"  Applied broad mask: {broad_mask_files[0].name}")

    # Crop and mask spatially to the region, then compute the sum spectrum.
    # subcube_from_regions crops to the bounding box and masks outside the
    # region, keeping only one cube's worth of data in memory at a time.
    cube_sub = cube.subcube_from_regions([reg])
    sum_k_spec = cube_sub.sum(axis=(1, 2))    # K, summed over spatial pixels

    # Free the cube before loading the next one
    del cube, cube_sub
    gc.collect()

    # Convert to Jy and save.
    # S_total[v] [Jy] = Σ_pix T_B(pix,v) [K] * (pixel_area / beam_area) / JYTOK
    sum_jy_spec = OneDSpectrum(
        sum_k_spec.value * beams_per_pixel / jytok * u.Jy,
        wcs=sum_k_spec.wcs,
    )
    spec_outfile = spec_out_dir / f"{folder.name}_sum_spectrum.fits"
    sum_jy_spec.write(spec_outfile, overwrite=True)
    print(f"  Saved spectrum → {spec_outfile.name}")

    # HI mass: integrate sum spectrum over velocity, scale by pixel area.
    # Σ_v sum_T[v] dv  [K km/s · pix]  ×  pixel_area [pc²/pix]  ×  α_HI  →  M_HI
    chan_width = np.abs(np.diff(sum_k_spec.spectral_axis)[0]).to(u.km / u.s)
    integrated_sum = np.nansum(sum_k_spec.value) * chan_width   # km/s (K implicit)

    hi_mass = (integrated_sum          # km/s (from chan_width; K from sum_k_spec)
               * u.K                   # K is the unit of sum_k_spec values
               * pixel_area_pc2 * u.pc**2
               * alpha_HI).to(u.solMass)

    frac_db08 = (hi_mass / hi_mass_db08_corr).decompose().value
    frac_w24 = (hi_mass / hi_mass_w24_corr).decompose().value

    names.append(folder.name)
    hi_masses.append(hi_mass.value)
    frac_db08s.append(frac_db08)
    frac_w24s.append(frac_w24)

    print(f"  mass={hi_mass:.3e},  frac_dB08={frac_db08:.3f},  frac_W24={frac_w24:.3f}")


table = Table(
    [names, hi_masses, frac_db08s, frac_w24s],
    names=['config', 'hi_mass_Msun', 'frac_of_dB08', 'frac_of_W24'],
)
table['hi_mass_Msun'].unit = u.solMass

out_file = data_path / 'hi_mass_from_spectra.ecsv'
table.write(out_file, overwrite=True)
print(f"\nTable saved to {out_file}")
print(table)
