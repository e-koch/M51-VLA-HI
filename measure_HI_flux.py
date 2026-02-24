
import re

import numpy as np
import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from regions import Regions

from pathlib import Path


alpha_HI = ((1.82e18 * u.cm**-2) * c.m_p).to(u.solMass / u.pc**2) / (u.K * u.km / u.s)


data_path = Path("/Users/ekoch/storage/M51/VLA")


# THINGS HI mass
hi_mass_db08 = 25.4e8 * u.solMass
dist_db08 = 8.0 * u.Mpc

# FEASTS HI mass
# Wang+24
hi_mass_w24 = 10**9.59 * u.solMass
dist_w24 = 8.0 * u.Mpc
hi_flux_w24 = 493.97 * u.Jy * u.km / u.s

# Missed fraction compared to THINGS is 0.35

# Distance used for mass conversion
dist = 8.6 * u.Mpc


# Read the region file:

reg_path = data_path / "m51_hi_avg_spectrum_region.reg"

reg = Regions.read(str(reg_path), format="ds9")[0]


names = []
hi_fluxes_Kkms = []
hi_flux_errs_Kkms = []
hi_fluxes = []
hi_flux_errs = []
hi_masses = []
hi_mass_errs = []
frac_db08s = []
frac_w24s = []

for folder in sorted(data_path.glob("m51_*")):
    if not folder.is_dir():
        continue

    # Find native-resolution broad mom0/emom0 files; exclude smoothed
    # variants like *_30as_broad_mom0.fits, *_60as_broad_mom0.fits, etc.
    all_mom0 = list(folder.glob("*_broad_mom0.fits"))
    all_emom0 = list(folder.glob("*_broad_emom0.fits"))

    native_mom0 = [f for f in all_mom0
                   if not re.search(r'_\d+as_broad_mom0\.fits$', f.name)]
    native_emom0 = [f for f in all_emom0
                    if not re.search(r'_\d+as_broad_emom0\.fits$', f.name)]

    if len(native_mom0) != 1 or len(native_emom0) != 1:
        print(f"Skipping {folder.name}: "
              f"found {len(native_mom0)} mom0, {len(native_emom0)} emom0")
        continue

    mom0_file = native_mom0[0]
    emom0_file = native_emom0[0]

    with fits.open(mom0_file) as hdul:
        hdr = hdul[0].header
        data_mom0 = hdul[0].data.squeeze()

    wcs = WCS(hdr).celestial
    jytok = hdr['JYTOK']
    bmaj_deg = hdr['BMAJ']
    bmin_deg = hdr['BMIN']
    cdelt_deg = abs(hdr['CDELT2'])

    with fits.open(emom0_file) as hdul:
        data_emom0 = hdul[0].data.squeeze()

    # Build pixel mask from sky region
    pix_reg = reg.to_pixel(wcs)
    mask_img = pix_reg.to_mask(mode='center').to_image(data_mom0.shape)
    if mask_img is None:
        print(f"Region does not overlap image for {folder.name}, skipping.")
        continue
    mask = mask_img.astype(bool)

    # Sum within region; propagate emom0 errors in quadrature
    sum_mom0 = np.nansum(data_mom0[mask])            # K km/s * npix
    sum_emom0 = np.sqrt(np.nansum(data_emom0[mask] ** 2))  # K km/s * npix

    # Beam and pixel geometry
    pixel_arcsec = cdelt_deg * 3600.0
    pixel_area_arcsec2 = pixel_arcsec ** 2
    beam_area_arcsec2 = (np.pi / (4 * np.log(2))
                         * bmaj_deg * 3600.0
                         * bmin_deg * 3600.0)
    beams_per_pixel = pixel_area_arcsec2 / beam_area_arcsec2

    # HI flux in Jy km/s: T[K] / JYTOK = S[Jy/beam], then * (pix/beam ratio)
    hi_flux = (sum_mom0 * beams_per_pixel / jytok) * u.Jy * u.km / u.s
    hi_flux_err = (sum_emom0 * beams_per_pixel / jytok) * u.Jy * u.km / u.s

    # Pixel area in pc^2 for mass via alpha_HI
    pixel_area_pc2 = (np.radians(cdelt_deg) * dist.to(u.pc).value) ** 2

    hi_mass = (sum_mom0 * u.K * u.km / u.s
               * pixel_area_pc2 * u.pc**2
               * alpha_HI).to(u.solMass)
    hi_mass_err = (sum_emom0 * u.K * u.km / u.s
                   * pixel_area_pc2 * u.pc**2
                   * alpha_HI).to(u.solMass)

    hi_mass_db08_corr = hi_mass_db08 * (dist / dist_db08) ** 2
    hi_mass_w24_corr = hi_mass_w24 * (dist / dist_w24) ** 2
    frac_db08 = (hi_mass / hi_mass_db08_corr).decompose().value
    frac_w24 = (hi_mass / hi_mass_w24_corr).decompose().value

    names.append(folder.name)
    hi_fluxes_Kkms.append(sum_mom0)
    hi_flux_errs_Kkms.append(sum_emom0)
    hi_fluxes.append(hi_flux.value)
    hi_flux_errs.append(hi_flux_err.value)
    hi_masses.append(hi_mass.value)
    hi_mass_errs.append(hi_mass_err.value)
    frac_db08s.append(frac_db08)
    frac_w24s.append(frac_w24)

    print(f"{folder.name}:  flux={hi_flux:.1f},  mass={hi_mass:.3e},  "
          f"frac_dB08={frac_db08:.3f},  frac_W24={frac_w24:.3f}")


table = Table(
    [names, hi_fluxes_Kkms, hi_flux_errs_Kkms,
     hi_fluxes, hi_flux_errs, hi_masses, hi_mass_errs,
     frac_db08s, frac_w24s],
    names=['config',
           'hi_flux_K_kms', 'hi_flux_err_K_kms',
           'hi_flux_Jy_kms', 'hi_flux_err_Jy_kms',
           'hi_mass_Msun', 'hi_mass_err_Msun',
           'frac_of_dB08', 'frac_of_W24'],
)
table['hi_flux_K_kms'].unit = u.K * u.km / u.s
table['hi_flux_err_K_kms'].unit = u.K * u.km / u.s
table['hi_flux_Jy_kms'].unit = u.Jy * u.km / u.s
table['hi_flux_err_Jy_kms'].unit = u.Jy * u.km / u.s
table['hi_mass_Msun'].unit = u.solMass
table['hi_mass_err_Msun'].unit = u.solMass

out_file = data_path / 'hi_flux_summary.ecsv'
table.write(out_file, overwrite=True)
print(f"\nTable saved to {out_file}")
print(table)
