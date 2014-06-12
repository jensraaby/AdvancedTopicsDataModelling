# Jens Raaby, May 2013

import numpy as np
import scipy as sp
import time
from scipy import misc


class ICMCleaner:

    def __init__(self, noisy, beta, eta, h):
        """
        Sets up an instance of the ICM noise-removal class.
        Reports initial system energy
        """
        # Create binary version of the image, as a byte array
        # n.b. A Python binary is also stored as a byte
        self.Y = np.int8(noisy)

        # Substitute the pixel intensities for our 2 values
        self.Y[noisy == 1] = -1
        self.Y[noisy == 0] = 1

        # Create the 'hidden' layer
        self.X = self.Y.copy()

        # store the energy parameters
        self.beta = beta
        self.eta = eta
        self.h = h

        # to save computation, neighbour indices are cached
        self.neighbour_cache = {}

        # Get the initial energy
        start = time.time()
        self.initial_E = self.energy()
        elapsed = (time.time() - start)
        print "Beta: %f\nEta: %f\nh: %f" % (beta, eta, h)
        print "Initial energy is ", self.initial_E

    def optimise(self, max_iterations=1):
        """
        Iteratively remove noise using ICM algorithm. Stops when energy converges.

        This method also times each iteration and prints statistics as it progresses
        """
        n = 1
        times = []
        energies = []
        E_0 = self.energy()
        energies.append(E_0)

        while n < (max_iterations + 1):
            print "Iteration %d:" % n

            start = time.time()
            self.iterate()
            elapsed = (time.time() - start)
            times.append(elapsed)

            E_1 = self.energy()
            energies.append(E_1)
            print "\tEnergy = %d" % E_1
            print "\tTime = %fs" % elapsed
            if E_1 == E_0:
                print "Converged after %d iterations " % n
                break
            E_0 = E_1
            n += 1

        return n, E_1, times

    def iterate(self):
        """
        Loop over all pixels and chooses the values which minimises the local energy
        """
        for (x, y), value in np.ndenumerate(self.X):
            E_now, E_flipped = self.eval(x, y)

            # now we need to see if flipping would reduce energy
            # E_flipped = self.eval(x,y,True)
            if E_now > E_flipped:
                self.X[(x, y)] *= -1

    def eval(self, x, y):
        """
        Compute the (local) energy for the given coordinates.

        This function returns the energy with the value in its current state
        and if it was flipped to the other value.
        """
        nb = self.neighbours(x, y)

        this_point = self.X[(x, y)]
        etaterm = self.eta * this_point * self.Y[(x, y)]
        hterm = self.h * this_point
        betaterm = 0
        betaterm_flipped = 0
        # loop over the neighbours - the slowest part of the program
        for n in nb:
            betaterm += this_point * self.X[n]
            betaterm_flipped += -this_point * self.X[n]
        # deferred multiplication:
        betaterm *= self.beta
        betaterm_flipped *= self.beta
        # evaluate the energies:
        E = hterm - betaterm - etaterm
        Eflipped = (-hterm) + etaterm - betaterm_flipped
        return E, Eflipped

    def energy(self):
        """
        Compute the energy of the entire system
        """
        beta_term = 0
        eta_term = 0
        h_term = 0

        # first loop over the "inner" nodes
        #  these all have the same edges (right and below)
        for x in xrange(0, self.X.shape[0] - 1):
            for y in xrange(0, self.X.shape[1] - 1):
                xy = self.X[(x, y)]
                # product with right neighbour
                beta_term += xy * self.X[(x + 1, y)]
                # product with bottom neighbour
                beta_term += xy * self.X[(x, y + 1)]

        # loop over the right hand edge
        right_x = self.X.shape[0] - 1
        for y in range(self.X.shape[1] - 1):
            this_value = self.X[(right_x, y)]
            beta_term += this_value * self.X[(right_x, y + 1)]

        # loop over the bottom edge
        bottom_y = self.X.shape[1] - 1
        for x in range(self.X.shape[0] - 1):
            this_value = self.X[(x, bottom_y)]
            beta_term += this_value * self.X[(x + 1, bottom_y)]

        # multiply through with the coefficient
        beta_term *= self.beta
        eta_term = self.eta * np.sum(self.X * self.Y)
        h_term = self.h * np.sum(self.X)
        return h_term - beta_term - eta_term

    def save(self, filename):
        """
        Saves the restored image to disk.
        This method handles converting the binary format back to 8-bit
        """

        restored = np.uint8(self.X)
        # this sets negatives to have 255
        restored[self.X == -1] = 0
        restored[self.X == 1] = 255

        sp.misc.imsave(filename, restored)

    def neighbours(self, x, y):
        """
        Returns array of neighbour coordinates. Results are cached in a dictionary
        since this function is called for every pixel on each iteration, with no change
        to the output.

        The actual logic is in a global method outside the class definition
        """
        if (x, y) in self.neighbour_cache:
            return self.neighbour_cache[(x, y)]
        else:
            ret = self.neighbour_cache[(x, y)] = neighbours(x, y, self.X.shape)
            return ret


def neighbours(x, y, size):
    """
    Basic method to find the indices for the neighbours of a given location in a matrix
    """
    max_x = size[0] - 1
    max_y = size[1] - 1

    try:
        if (x > max_x or y > max_y or x < 0 or x < 0):
            raise Exception("Invalid coordinates")

        # build a list of coordinates for the neighbours
        b = []
        if x > 0:
            b.extend([(x - 1, y)])
        if x < max_x:
            b.extend([(x + 1, y)])

        if y > 0:
            b.extend([(x, y - 1)])
        if y < max_y:
            b.extend([(x, y + 1)])

        return b

    except:
        print "Invalid coordinates: %s, %s" % (x, y)
