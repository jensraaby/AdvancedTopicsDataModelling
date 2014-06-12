# Jens Raaby, May 2013

from scipy import misc
import icm

stats = {}
def parameter_selection(max_iterations=20):
    """
    This is a method for evaluating different parameter values
    """

    noisy = misc.imread('noisyImage.png')
 
    
    # change these values if desired:
    h = [-3, -2, -1, 0, 1, 2, 3]
    eta = [0, 1, 2, 3, 4, 5]
    beta = 1
    for hv in h:
        for ev in eta:
            c = icm.ICMCleaner(noisy, beta, ev, hv)
            n, E_1, times = stats[(ev, hv)] = c.optimise(max_iterations)
            fname = "parametersearch/eta%d_h%d_%diter.png" % (ev, hv, n)
            c.save(fname)


# Simply run it with 20 iterations - will take a while
parameter_selection(20)

