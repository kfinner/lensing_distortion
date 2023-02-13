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
from scipy.special import gamma, gammaincinv, gammainc

parser = argparse.ArgumentParser(description='Download and drizzle HST observations')
parser.add_argument("--img", dest="img", action="store", default="none",
                        help="The image that you would like to distort.")
parser.add_argument("--mass", dest="mass", action="store", default=1e12, type=float,
                        help="The mass of the galaxy or cluster.")
parser.add_argument("--zsource", dest="zsource", action="store", default=0.6, type=float,
                        help="The redshift of the source.")
parser.add_argument("--zlens", dest="zlens", action="store", default=0.3, type=float,
                        help="The redshift of the lens.")
parser.add_argument("--pixscale", dest="pixscale", action="store", default=0.001, type=float,
                        help="The pixel scale of the image in arcsecs per pixel.")
args = parser.parse_args()

def gauss2d(A, cx, cy, sigx, sigy, theta, size):
    xx, yy = np.meshgrid(np.arange(size),np.arange(size))
    return A * np.exp(-0.5* ( ( ( ( (xx-cx) * np.cos(theta) + (yy-cy) * np.sin(theta) ) / sigx )**2 ) + ( ( ( (xx-cx) * np.sin(theta) - (yy-cy) * np.cos(theta) ) / sigy )**2 ) ) )

def sersic2d(flux,cenx,ceny,hlr,e,theta,n,size=51):
    xx, yy = np.meshgrid(np.arange(size),np.arange(size))
    bn = gammaincinv(2*n, 0.5)#2*n - 1/3 + 4/(405*n) + 46 / (25515*n**2) # Capaccioli 1989 (good approximation for 0.5 < n < 10)

    a = hlr
    b = hlr * (1-e)

    A = (xx - cenx) * np.cos(theta) + (yy - ceny) * np.sin(theta)
    B = -(xx - cenx) * np.sin(theta) + (yy - ceny) * np.cos(theta)
    R = np.sqrt( (A / a)**2 + (B / b)**2)

    return flux * np.exp(-bn * ( R**(1./n) - 1) )  #/ (bn**2 * hlr**2)

#input_img = Image.open('/Users/kyle/Work/NIRWL/CANDELS/mass_maps/GN/GN_color_cropped.jpeg')
#input_img = np.array(input_img).astype(float)
#input_img = input_img.T
#input_img = input_img[0,:,:]

#input_img = fits.open(args.img)[0].data
#import models
#pad_size = max(input_img.shape)//4
#if pad_size % 2 == 1: pad_size+=1
#input_img = np.pad(input_img, pad_size)


input_img = gauss2d(10, 1000+10, 1000, 3, 3, 0, 2001)
#models.Gaussian(0, 100, 1000+i, 1000, 3, 3, 0).get_postage(2001)

sigma_v = ( 100 * np.pi * const.G.to('Mpc^3 / (M_sun s^2)')**3 * (args.mass*u.solMass)**2 * cosmo.critical_density(args.zlens).to('M_sun/Mpc^3') / 3 )**(1./6.)
sigma_v = sigma_v.to('km/s')

beta = cosmo.angular_diameter_distance_z1z2(args.zlens, args.zsource) / cosmo.angular_diameter_distance(args.zsource) # Lensing efficiency.
theta_e = 4 * np.pi * sigma_v**2 / const.c.to('km/s')**2 * beta # In radians.
theta_e *= 206265 # in arcsecs
print(theta_e)
theta_e /= float(args.pixscale) # in pixels
print('The Einstein radius in pixels is', theta_e, args.pixscale)


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

theta_img[yy,xx] = input_img[dys.astype(int),dxs.astype(int)]

# Change this into arcsecs.
sersic = sersic2d(10, 1000, 1000, 0.3/args.pixscale, 0, 0, 1, size=2001)

# We could easily implement an SED/filter scaling here.
sr = sersic * 1.4
sg = sersic * 1
sb = sersic * 0.4

tr = theta_img * 0.4
tg = theta_img * 1
tb = theta_img * 2.0

r = sr+tr
g = sg+tg
b = sb+tb

s_max = np.array([r.max(), g.max(), b.max()]).max()
r /= s_max
g /= s_max
b /= s_max


img = np.zeros((r.shape[0], r.shape[1], 3), dtype=float)
img[:,:,0] = r * 255
img[:,:,1] = g * 255
img[:,:,2] = b * 255
img = img.astype(np.uint8)
img = np.flip(img, 0)
img = Image.fromarray(img)


#theta_img = theta_img#.T
plt.imshow(img, origin='lower')
#plt.axvline(theta_img.shape[1]//2, color='black')
#plt.axhline(theta_img.shape[0]//2, color='black')
plt.xlim(1000-theta_e*3-50,1000+theta_e*3+51)
plt.ylim(1000-theta_e*3-50,1000+theta_e*3+51)

plt.show()
#distorted_img = Image.fromarray(theta_img.astype(np.uint8))
#plt.imshow(distorted_img)
#distorted_img.save(f'/Users/kyle/Work/sis_model_demo/images/distorted_image_test.jpg')
