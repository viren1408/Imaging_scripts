import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import Angle
#This code was written to check the offsets between two radio images of same pointing but different epochs, I wanted to median stack them but there seemed to be an offset in the sources which I want to correct! 
def read_fits_image(file_path):
    """Open and read the header and data from a FITS file."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found.")
    with fits.open(file_path) as hdul:
        header = hdul[0].header
        data = hdul[0].data
    return header, data

def get_wcs_translation(file1, file2):
    """Calculate the pixel and world coordinate shift between two FITS files."""
    header1, data1 = read_fits_image(file1)
    header2, data2 = read_fits_image(file2)
   #currently for GMRT image cube 
    if data1.ndim == 4:
        data1 = data1[0,0,:,:]
        data2 = data2[0,0,:,:]

    # Load WCS
    wcs1 = WCS(header1, naxis=2)
    wcs2 = WCS(header2, naxis=2)

    # Get the reference pixel coordinates (CRPIX1, CRPIX2)
    crpix1_file1 = np.array([header1['CRPIX1'], header1['CRPIX2']])
    crpix1_file2 = np.array([header2['CRPIX1'], header2['CRPIX2']])

    pixel_shift = crpix1_file2 -crpix1_file1

    #world coordinate shift
    world_coords_file1 =wcs1.pixel_to_world(crpix1_file1[0],crpix1_file1[1])
    world_coords_file2 = wcs2.pixel_to_world(crpix1_file2[0],crpix1_file2[1])
    world_shift = world_coords_file2.spherical_offsets_to(world_coords_file1) #getting spherical offsets (look documentation!)

    world_shift_arcsec = (world_shift[0].to(Angle('arcsec')), world_shift[1].to(Angle('arcsec')))

    return pixel_shift, world_shift_arcsec

file1 = 'G030.0-6.8_17_AUGUST.fits'
file2 = 'G030.0-6.8_20_AUGUST.fits'
pixel_shift, world_shift_arcsec = get_wcs_translation(file1, file2)

print(f"Pixel shift: {pixel_shift}")
print(f"World shift (RA, Dec) in arcseconds: {world_shift_arcsec[0]}, {world_shift_arcsec[1]}")
