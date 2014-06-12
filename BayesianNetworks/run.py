# ATDM Assignment 1
# Jens Raaby, May 2013

import sampling
    
def main():
    
    
    bn_given = sampling.BN()
    bn_given.add_node("x1",[],{(True,): 0.3})
    bn_given.add_node("x2",[],{(True,): 0.2})
    bn_given.add_node("x3",[],{(True,): 0.5})
    bn_given.add_node("x4",["x1","x2","x3"],{(False,False,False): 0.05,
                                             (False,True, False): 0.7,
                                             (True, False,False): 0.3,
                                             (True, True, False): 0.9,
                                             (False,False,True):  0.5,
                                             (False,True, True):  0.75,
                                             (True, False,True):  0.7,
                                             (True, True, True):  0.95})
                                             
   
                                        
    bn_given.add_node("x5",["x1","x3"],{(False,False):0.05,
                                        (False,True): 0.8,
                                        (True,False): 0.07,
                                        (True,True):  0.8})

    bn_given.add_node("x6",["x4"],{(False,):0.2,(True,): 0.7})
    
    bn_given.add_node("x7",["x4","x5"],{(False,False):0.1,
                                        (False,True): 0.3,
                                        (True,False): 0.3,
                                        (True,True):  0.7})
                                        

    
    
    # print_graph.BNtoPNG(bn_given,'given.png')
    print bn_given.joint_probability()
    
    # Try out the joint sample and take an ancestral sample
    bn_given.print_joint_probability()
    bn_given.print_ancestral_sample()
    bn_given.reset_sampling()
    
    print "\n/////////////////\n"
    
    # Try markov blanket sampling:
    print "Markov Blanket Sampling X5"
    blanket_nodes = {"x1":False, "x3": True, "x4": False, "x7": True}
    v = bn_given.markov_blanket_sampling("x5",blanket_nodes)
    print "Sampled %s" % v
    print "\n/////////////////\n"
    
    
if __name__ == '__main__':
    main()
