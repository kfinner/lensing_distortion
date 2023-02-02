# Distort an image according to the SIS mass profile.
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
import astropy.units as u
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
from PIL import Image
import argparse

parser = argparse.ArgumentParser(description='Download and drizzle HST observations')
parser.add_argument("--img", dest="img", action="store", default="none",
                        help="The image that you would like to distort.")
parser.add_argument("--mass", dest="mass", action="store", default=5e14,
                        help="The mass of the galaxy or cluster.")
parser.add_argument("--zsource", dest="zsource", action="store", default=0.6,
                        help="The redshift of the source.")
parser.add_argument("--zlens", dest="zlens", action="store", default=0.3,
                        help="The redshift of the lens.")
parser.add_argument("--pixscale", dest="pixscale", action="store", default=0.1,
                        help="The pixel scale of the image in arcsecs per pixel.")
args = parser.parse_args()

def gauss2d(A, cx, cy, sigx, sigy, theta, size):
    xx, yy = np.meshgrid(np.arange(size),np.arange(size))
    return A * np.exp(-0.5* ( ( ( ( (xx-cx) * np.cos(theta) + (yy-cy) * np.sin(theta) ) / sigx )**2 ) + ( ( ( (xx-cx) * np.sin(theta) - (yy-cy) * np.cos(theta) ) / sigy )**2 ) ) )

#input_img = Image.open('/Users/kyle/Work/NIRWL/CANDELS/mass_maps/GN/GN_color_cropped.jpeg')
#input_img = np.array(input_img).astype(float)
#input_img = input_img.T
#input_img = input_img[0,:,:]

#input_img = fits.open(args.img)[0].data
#import models
#pad_size = max(input_img.shape)//4
#if pad_size % 2 == 1: pad_size+=1
#input_img = np.pad(input_img, pad_size)

for i in range(0,500,5):
    input_img = gauss2d(1, 1000+i, 1000, 3, 3, 0, 2001)
#models.Gaussian(0, 100, 1000+i, 1000, 3, 3, 0).get_postage(2001)

    sigma_v = ( 100 * np.pi * const.G.to('Mpc^3 / (M_sun s^2)')**3 * (args.mass*u.solMass)**2 * cosmo.critical_density(args.zlens).to('M_sun/Mpc^3') / 3 )**(1./6.)
    sigma_v = sigma_v.to('km/s')

    beta = cosmo.angular_diameter_distance_z1z2(args.zlens, args.zsource) / cosmo.angular_diameter_distance(args.zsource) # Lensing efficiency.
    theta_e = 4 * np.pi * sigma_v**2 / const.c.to('km/s')**2 * beta # In radians.
    theta_e *= 206265 # in arcsecs
    theta_e /= args.pixscale # in pixels
    print(theta_e)


    lens_pos = input_img.shape[1]//2, input_img.shape[0]//2
    theta_img = np.zeros((input_img.shape[0],input_img.shape[1]))

    xx,yy = np.meshgrid(np.arange(input_img.shape[1]),np.arange(input_img.shape[0]))
    xx -= lens_pos[0]
    yy -= lens_pos[1]
    radius = np.sqrt(xx**2+yy**2)

    '''
    This is where all the magic happens.
    It takes the output positions and finds the input pixels
    that will contribute to its signal.
    e.g., for a new position (dxs,dys) find the pixel that contributes to its signal.
    Following the SIS profile, that is the distance from the Einstein radius.
    (theta_e / distance from the center)
    '''

    dxs = xx-(theta_e * xx / (radius + (radius == 0)))
    dys = yy-(theta_e * yy / (radius + (radius == 0)))

    xx += lens_pos[0]
    yy += lens_pos[1]
    dxs += lens_pos[0]
    dys += lens_pos[1]
    """
    # Trim the dxs/dys that are off the edge of the image.
    mask_dxs = (dxs.astype(int) >= 0) & (dxs.astype(int) < input_img.shape[1])
    mask_dys = (dys.astype(int) >= 0) & (dys.astype(int) < input_img.shape[0])

    xx = xx[mask_dxs]
    dxs = dxs[mask_dxs]
    yy = yy[mask_dys]
    dys = dys[mask_dys]
    """
    theta_img[yy,xx] = input_img[dys.astype(int),dxs.astype(int)]

    theta_img = theta_img#.T
    plt.imshow(theta_img, origin='lower')
    plt.axvline(theta_img.shape[1]//2, color='black')
    plt.axhline(theta_img.shape[0]//2, color='black')
    plt.xlim(800,1200)
    plt.ylim(800,1200)
    print(i)
    plt.pause(0.1)
    #plt.close('all')
#distorted_img = Image.fromarray(theta_img.astype(np.uint8))
#plt.imshow(distorted_img)
#distorted_img.save(f'/Users/kyle/Work/sis_model_demo/images/distorted_image_test.jpg')
