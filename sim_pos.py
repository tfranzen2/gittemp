#!/usr/bin/env python3

from optparse import OptionParser
import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table,Column
from astropy.io.votable import writeto as writetoVO

# Read input parameters
usage="Usage: %prog [options] <file>\n"
parser = OptionParser(usage=usage)
parser.add_option('--nsrc',type="int", dest="nsrc", default=100, help="Number of simulated sources [default=%default]")
parser.add_option('--flux',type="float", dest="flux", default=1.0, help="Flux density of simulated sources, in mJy [default=%default]")
parser.add_option('--rad',type="float", dest="rad", default=1000.0, help="Maximum distance between simulated sources and map centre, in deg [default=%default]")
parser.add_option('--min_sep',type="float", dest="min_sep", default=0.0, help="Minimum separation between simulated sources, in arcsec [default=%default]")
parser.add_option('--output_root',type="string", dest="output_root", default='sim_pos', help="Output root filename [default=%default]")
(options, args) = parser.parse_args()
input_map=args[0]
nsrc=options.nsrc
flux=options.flux
min_sep=options.min_sep
rad=options.rad
output_root=options.output_root

# Convert simulated source flux from mJy to Jy
flux=flux/1000.0

# Get image dimensions
xdim = fits.getheader(input_map)['NAXIS1']
ydim = fits.getheader(input_map)['NAXIS2']

# Calculate central pixel
px_centre=xdim/2.0
py_centre=ydim/2.0

# Get corresponding (RA, Dec) position in deg
hdulist = fits.open(input_map) # read image header
w = wcs.WCS(hdulist[0].header,naxis=2) # get WCS
ref_ra, ref_dec = w.wcs_pix2world(px_centre, py_centre, 1)
c2=SkyCoord(ref_ra, ref_dec, unit="deg") # SkyCoord representation, defaults to ICRS frame
print ('Map centre RA (deg) =', '%.2f' % ref_ra)
print ('Map centre Dec (deg) =', '%.2f' % ref_dec)

# Get beam major and minor axes, and position angle
bmaj = fits.getheader(input_map)['BMAJ'] # in deg
bmin = fits.getheader(input_map)['BMIN'] # in deg
bpa = fits.getheader(input_map)['BPA'] # in deg
# Convert bmaj and bmin from deg to arcsec
bmaj=bmaj*3600.0
bmin=bmin*3600.0
print ('Beam major axis (arcsec) =', '%.1f' % bmaj)
print ('Beam minor axis (arcsec) =', '%.1f' % bmin)
print ('Beam position angle (deg) =', '%.1f' % bpa)

# Get pixel size (assume square pixel)
pixel_size = fits.getheader(input_map)['CDELT2'] # in deg
pixel_size=pixel_size*3600.0 # convert pixel size from deg to arcsec

i=0
k=0
ra_list=[]
dec_list=[]
x_list=[]
y_list=[]
print ('Generating simulated source positions')
while i<nsrc:
    if k >= 20000:
        print('Aborting as field is overcrowded. Reduce number of simulated sources and/or minimum separation between simulated sources.')
        exit()
    x = np.random.uniform(3.0,xdim-3.0) # x pixel position drawn randomly, must lie at least 3 pixels from map edge
    y = np.random.uniform(3.0,ydim-3.0) # y pixel position drawn randomly, must lie at least 3 pixels from map edge
    ra, dec = w.wcs_pix2world(x, y, 1) # convert pixel position to (RA, Dec) in deg
    c1 = SkyCoord(ra, dec, unit="deg") # SkyCoord representation, defaults to ICRS frame
    if c1.separation(c2).deg <= rad: # source must lie within rad deg from map centre
        if min_sep==0.0:
            ra_list.append(ra)
            dec_list.append(dec)
            x_list.append(x)
            y_list.append(y)
            i=i+1
        else:
            s1=1000000.0 # s1 is minumum separation with all other simulated sources, in arcsec
            for j in range(0,len(x_list)):
                s2=np.sqrt((x-x_list[j])**2+(y-y_list[j])**2)*pixel_size # this is less accurate than great circle distance formula but faster
                if s2 < s1:
                    s1=s2
            if s1 >= min_sep: # source must lie at least min_sep arcsec from any other simulated source
                ra_list.append(ra)
                dec_list.append(dec)
                x_list.append(x)
                y_list.append(y)
                i=i+1
    k=k+1

# Convert source positions from deg to sexagesimal format
ra_str=[]
dec_str=[]
for i in range(0,nsrc):
    c=SkyCoord(ra_list[i], dec_list[i], unit="deg")
    pos=c.to_string('hmsdms', sep=':')
    ra_str.append(pos.split()[0])
    dec_str.append(pos.split()[1])

# Write kvis file
f=open(output_root+'.ann','w')
print("colour red", file=f)
for i in range(0,nsrc):
    print("cross W", '%.6f' % ra_list[i], '%.6f' % dec_list[i], "0.1 0.1", file=f)
    
# Generate votable in Aegean format
island=list(range(1,nsrc+1))
source=np.repeat('source', nsrc)
background=np.zeros(nsrc)
local_rms=np.zeros(nsrc)
err_ra=np.zeros(nsrc)
err_dec=np.zeros(nsrc)
peak_flux=np.ones(nsrc)*flux
err_peak_flux=np.zeros(nsrc)
int_flux=np.ones(nsrc)*flux
err_int_flux=np.zeros(nsrc)
a=np.ones(nsrc)*bmaj
err_a=np.zeros(nsrc)
b=np.ones(nsrc)*bmin
err_b=np.zeros(nsrc)
pa=np.ones(nsrc)*bpa
err_pa=np.zeros(nsrc)
flags=np.zeros(nsrc)
residual_mean=np.zeros(nsrc)
residual_std=np.zeros(nsrc)
uuid=np.repeat('uuid', nsrc)
psf_a=np.zeros(nsrc)
psf_b=np.zeros(nsrc)
psf_pa=np.zeros(nsrc)
newvot = Table( [island, source, background, local_rms, ra_str, dec_str, ra_list, err_ra, dec_list, err_dec, peak_flux, err_peak_flux, int_flux, err_int_flux, a, err_a, b, err_b, pa, err_pa, flags, residual_mean, residual_std, uuid, psf_a, psf_b, psf_pa], names=('island', 'source', 'background', 'local_rms', 'ra_str', 'dec_str', 'ra', 'err_ra', 'dec', 'err_dec', 'peak_flux', 'err_peak_flux', 'int_flux', 'err_int_flux', 'a', 'err_a', 'b', 'err_b', 'pa', 'err_pa', 'flags', 'residual_mean', 'residual_std', 'uuid', 'psf_a', 'psf_b', 'psf_pa'), meta={'name': 'first table'} )
writetoVO(newvot, output_root+'.vot')
