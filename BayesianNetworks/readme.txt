README 
Jens Raaby, May 2013

Dependencies:
Python 2.7 (tested on Mac OS 10.8)

Purpose:
A small library for (binary variable) Bayesian Networks.
Joint probabilities, ancestral sampling , markov blanket sampling
and rejection sampling are implemented.

Disclaimer: 
This is not intended to be used in production code.
There are many improvements which could be made before reusing.


-----------------
To run the demo code:

For the example network from figure 8.2 of ``Pattern Recognition and Machine Learning’’ (Bishop 2007), the code to run ancestral sampling and markov blanket sampling is in the file run.py

From a terminal, execute:
python run.py


The file wetgrass.py implements the simpler Wet Grass network. The script runs through several demonstrations of the sampling.py package: ancestral sampling, markov blanket sampling and rejection sampling.

Execute:
python wetgrass.py
