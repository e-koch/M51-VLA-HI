
import re

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.patches import Circle
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from regions import Regions
from spectral_cube import OneDSpectrum

from pathlib import Path


data_path = Path("/Users/ekoch/storage/M51/VLA")
spec_dir = data_path / "avg_spectra"

reg_path = data_path / "m51_hi_avg_spectrum_region.reg"
reg = Regions.read(str(reg_path), format="ds9")[0]
reg_center = SkyCoord(reg.center.ra, reg.center.dec, frame='icrs')
reg_radius_arcsec = reg.radius.to(u.arcsec).value

# -------------------------------------------------------------------------
# Collect A+B+C+D configurations
# -------------------------------------------------------------------------

configs = []
for folder in sorted(data_path.glob("m51_A+B+C+D*")):
    if not folder.is_dir():
        continue

    all_mom0 = list(folder.glob("*_broad_mom0.fits"))
    native = [f for f in all_mom0
              if not re.search(r'_\d+as_broad_mom0\.fits$', f.name)]
    if len(native) != 1:
        print(f"Skipping {folder.name}: {len(native)} native mom0 files")
        continue

    with fits.open(native[0]) as hdul:
        hdr = hdul[0].header
        data = hdul[0].data.squeeze()

    bmaj_arcsec = hdr['BMAJ'] * 3600.0
    bmin_arcsec = hdr['BMIN'] * 3600.0
    beam_avg = round(0.5 * (bmaj_arcsec + bmin_arcsec), 1)
    cdelt_arcsec = abs(hdr['CDELT2']) * 3600.0

    m = re.search(r'5kms_(\d+)as', folder.name)
    taper = m.group(1) if m else '?'

    # Spectrum file (named by folder, not by cube file)
    spec_file = spec_dir / f"{folder.name}_sum_spectrum.fits"

    configs.append({
        'folder': folder,
        'data': data,
        'hdr': hdr,
        'wcs': WCS(hdr).celestial,
        'cdelt_arcsec': cdelt_arcsec,
        'beam_avg': beam_avg,
        'taper': taper,
        'spec_file': spec_file,
    })

# -------------------------------------------------------------------------
# Build figure: 2×3 grid
# -------------------------------------------------------------------------

# "talk"-level text sizes
plt.rcParams.update({
    'font.size': 19,
    'axes.titlesize': 20,
    'axes.labelsize': 19,
    'xtick.labelsize': 17,
    'ytick.labelsize': 17,
    'legend.fontsize': 17,
})

fig = plt.figure(figsize=(18, 12), layout='constrained')
fig.get_layout_engine().set(h_pad=0.02, w_pad=0.02, hspace=0.0, wspace=0.0)
gs = gridspec.GridSpec(2, 3, figure=fig)

# Field of view half-width: slightly larger than the region
fov_arcsec = reg_radius_arcsec * 1.15

# Common colourscale across all moment-0 panels
# Compute vmax only within the source region to avoid edge/noise spikes
def _reg_p99(cfg):
    pix_reg = reg.to_pixel(cfg['wcs'])
    mask = pix_reg.to_mask(mode='center').to_image(cfg['data'].shape).astype(bool)
    return np.nanpercentile(cfg['data'][mask], 99)

common_vmax = np.max([_reg_p99(cfg) for cfg in configs])
common_norm = mcolors.Normalize(vmin=0, vmax=common_vmax)

# Colours for spectrum overlay
spec_colors = plt.cm.plasma(np.linspace(0.1, 0.80, len(configs)))

# -------------------------------------------------------------------------
# Moment-0 panels (panels 0–4)
# -------------------------------------------------------------------------

for idx, cfg in enumerate(configs):
    row, col = idx // 3, idx % 3
    wcs = cfg['wcs']

    ax = fig.add_subplot(gs[row, col], projection=wcs)

    im = ax.imshow(cfg['data'], norm=common_norm, cmap='inferno', origin='lower')

    # Crop to consistent FOV around the region centre
    cx, cy = wcs.world_to_pixel_values(
        reg_center.ra.deg, reg_center.dec.deg)
    fov_pix = fov_arcsec / cfg['cdelt_arcsec']
    ax.set_xlim(cx - fov_pix, cx + fov_pix)
    ax.set_ylim(cy - fov_pix, cy + fov_pix)

    # Region boundary
    reg_pix_r = reg_radius_arcsec / cfg['cdelt_arcsec']
    ax.add_patch(Circle((cx, cy), radius=reg_pix_r,
                         edgecolor='cyan', facecolor='none',
                         lw=1.2, ls='--', zorder=5))

    # Beam indicator (lower-right, 8 % from edge)
    margin = 0.08 * 2 * fov_pix
    bx = cx + fov_pix - margin
    by = cy - fov_pix + margin
    beam_pix_r = cfg['beam_avg'] / 2.0 / cfg['cdelt_arcsec']
    ax.add_patch(Circle((bx, by), radius=beam_pix_r,
                         color='white', zorder=6))

    # Scale bar (upper-right, 8 % from edge)
    scalebar_arcsec = np.degrees(5.0e3 / 8.6e6) * 3600  # 5 kpc at 8.6 Mpc
    scalebar_pix = scalebar_arcsec / cfg['cdelt_arcsec']
    sb_x1 = cx + fov_pix - margin
    sb_x0 = sb_x1 - scalebar_pix
    sb_y = cy + fov_pix - margin
    ax.plot([sb_x0, sb_x1], [sb_y, sb_y],
            color='white', lw=2.5, solid_capstyle='butt', zorder=7)
    ax.text(0.5 * (sb_x0 + sb_x1), sb_y - 0.25 * margin, '5 kpc',
            color='white', ha='center', va='top', fontsize=17, zorder=7)

    # Inset label (upper-left): colored background matching the spectrum line
    beam_pc = round(cfg['beam_avg'] / 206265 * 8.6e6, -1)
    label_text = (f"{cfg['taper']}'' taper\n"
                  f"{cfg['beam_avg']:.1f}'' beam\n"
                  f"({beam_pc:.0f} pc)")
    ax.text(cx - fov_pix + 0.02 * 2 * fov_pix, cy + fov_pix - 0.02 * 2 * fov_pix,
            label_text,
            color='white', fontsize=15, va='top', ha='left', zorder=7,
            linespacing=1.4,
            bbox=dict(facecolor=spec_colors[idx], edgecolor='none',
                      boxstyle='round,pad=0.3', alpha=1.0))

    # Single shared colourbar on the top-right panel only
    if idx == 2:
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cbar.set_label('K km s$^{-1}$')
        cbar.ax.tick_params(labelsize=17)

    # Coordinate axis labels: Dec on left column, RA on bottom row
    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm')
    lon.set_ticks(spacing=4 * u.arcmin)
    lat.set_ticks(spacing=2 * u.arcmin)
    lon.set_ticks_position('b')
    lat.set_ticks_position('l')

    if col == 0:
        lat.set_axislabel('Dec (J2000)', fontsize=19)
    else:
        lat.set_ticklabel_visible(False)
        lat.set_axislabel('')

    if row == 1 and col < 2:
        lon.set_axislabel('RA (J2000)', fontsize=19)
    else:
        lon.set_ticklabel_visible(False)
        lon.set_axislabel('')

    lon.set_ticklabel(size=17)
    lat.set_ticklabel(size=17)

# -------------------------------------------------------------------------
# Spectrum panel (bottom right)
# -------------------------------------------------------------------------

ax_spec = fig.add_subplot(gs[1, 2])

for cfg, color in zip(configs, spec_colors):
    if not cfg['spec_file'].exists():
        print(f"Missing spectrum: {cfg['spec_file'].name}")
        continue

    with fits.open(cfg['spec_file']) as hdul:
        spec = OneDSpectrum.from_hdu(hdul[0])

    # Exclude first and last channels
    vel = spec.spectral_axis.to(u.km / u.s).value[1:-1]
    tb = spec.value[1:-1]

    ax_spec.plot(vel, tb, color=color, lw=1.5,
                 label=f"{cfg["beam_avg"]:.1f}'' ({cfg["taper"]}'' taper)",
                 drawstyle='steps-mid')

ax_spec.axhline(0, color='0.5', lw=0.8, ls='--')
ax_spec.set_xlabel(r'Velocity (km s$^{-1}$)')
ax_spec.set_ylabel('Flux density (Jy)', labelpad=2)
ax_spec.yaxis.tick_right()
ax_spec.yaxis.set_label_position('right')

ax_spec.invert_xaxis()

# -------------------------------------------------------------------------
# Save
# -------------------------------------------------------------------------

for suffix in ('png', 'pdf'):
    out = data_path / f'hi_summary.{suffix}'
    fig.savefig(out, dpi=150 if suffix == 'png' else None,
                bbox_inches='tight')
    print(f"Saved {out}")

plt.close(fig)
