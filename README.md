# lensing_distortion
Distorts an image according to an SIS profile.

to run: python sis_iamge_sersic_lens.py from terminal

dependencies are:
astropy
pillow
numpy
matplotlib
argparse
scipy

Current state:
1. Creates a Gaussian background galaxy and distorts it according to an SIS mass profile.
2. Adds a Sersic galaxy into the image.
3. Artificially sets the colors by non-scientifically (but quickly) scaling RGB channels.

Things that can be incluced in future versions:
- Galaxy SED / filter to determine observed flux.
- Create a second image that has been dispersed by a grism.
