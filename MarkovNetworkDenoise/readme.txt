# Markov random field for image de-noising
# Jens Raaby, May 2013

Purpose
=======
This code can remove noise from a binary image using the
`Iterated conditional modes’ (ICM) technique.
Underlying this is an undirected graphical model (Markov random field).
Pixels in the noisy image are connected to corresponding pixels in the (unknown) original
image. The pixels in the original image are connected to their direct neighbours.

Using the fact that adjacent pixels are most likely to be the same, an
energy function is optimised over the entire Markov random field.

More details in ``Pattern Recognition and Machine Learning’’ section 8.3.3 (Bishop 2007)

Dependencies
============
Python 2.7 (tested with Python 2.7.4 on Mac OS X 10.8.3)
 
Python Libraries Used:
- numpy
- scipy
- matplotlib

Files
=====
icm.py:
The graphical model and ICM algorithm

run.py:
A simple run (with at most 20 iterations).
Saves the cleaned image to disk.

parameter_search.py:
Tries multiple values for eta and h

Running
=======

To run the ICM algorithm with the default parameters (h=0, eta=1, beta=1)
execute the run.py script, for example (in a terminal):
$ python run.py

To run grid search with various values of h and eta, execute:
$ python parameter_search.py

This will write images to the parametersearch sub-directory.