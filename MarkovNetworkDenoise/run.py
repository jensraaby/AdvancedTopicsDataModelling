# Jens Raaby, May 2013

from scipy import misc
import icm

# This script performs noise removal with the given parameters, and saves the resulting image
noisy = misc.imread('noisyImage.png')
tcv = icm.ICMCleaner(noisy, 1, 1, 0)
tcv.optimise(20)
tcv.save('restored20.png')
