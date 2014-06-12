# ATDM Bayesian networks
# Jens Raaby, May 2013

# This script runs some examples using the popular wet grass network
import sampling

def main():

    # We create the network used in the Sampling and MCMC lecture slides
    bn = sampling.BN()
    bn.add_node("Cloudy",[],{(True,): 0.5})
    bn.add_node("Sprinkler",["Cloudy"],{(True,): 0.1, (False,): 0.5})
    bn.add_node("Rain",["Cloudy"],{(False,):(0.2),(True,): (0.8)})
    bn.add_node("WetGrass",["Sprinkler","Rain"],{(False,False):0.01, # otherwise we can get divide by zero
                                                 (False,True): 0.9,
                                                 (True,False): 0.9, 
                                                 (True,True):  0.99})
    bn.validate()
    
    # Try out the joint sample and take an ancestral sample
    bn.print_joint_probability()
    bn.print_ancestral_sample()
    bn.reset_sampling()
    
    print "Ancestral sampling for only Cloudy and Sprinkler: \n%s" % bn.ancestral_sample(["Cloudy","Sprinkler"])
    
    print "\n/////////////////\n"
    # Perform lots of samples, printing the results as a frequency table
    bn.print_multi_sample(10000)
    
    print "\n/////////////////\n"
    
    # Get the markov blanket for Rain
    mb = bn.markov_blanket("Rain")
    print "The markov blanket for Rain is \n\t%s" % mb
    
    # Try markov blanket sampling:
    blanket_nodes = {"Cloudy":False, "WetGrass": True, "Sprinkler": False}
    bn.markov_blanket_sampling("Rain",blanket_nodes)
    print "\n/////////////////\n"
    
    
    obs = {"Cloudy":True, "Rain": False}
    # try rejection sampling:
    bn.rejection_sampling(blanket_nodes,500)
        
if __name__ == '__main__':
    main()
