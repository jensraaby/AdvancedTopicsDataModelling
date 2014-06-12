# ATDM Bayesian networks
# Jens Raaby, May 2013

""" 
This file contains the code for Bayesian Networks with Binary Random Variables
Ancestral sampling and markov blanket sampling are implemented
"""

# We use the built in random generator and ordered dictionary libraries
import random
from collections import OrderedDict
        
class BN:
    """ A simple Bayesian network class.
    Stores nodes in a dictionary (str name -> BinaryRandomVariable)
    Edges are encoded in the nodes themselves
    """
    
    def __init__(self):
       # This dictionary keeps nodes, ordered by insertion
       self._nodes = OrderedDict() 
       
       # This dictionary keeps the ancestral ordering for the nodes
       self.ordering = None

    
    def add_node(self, name, parents, pt):
        """Creates a random variable and inserts as a node"""
        if name in self._nodes:
            raise Exception("A random variable already exists with the name (%s)" % name)
        
        # we validate the parents, while adding references to list 'ps'
        # ps must be a list so the order is maintained
        ps = list()
        for p in parents:
            if p not in self._nodes:
                raise Exception("Error adding %s: Parent %s does not exist" % (name,p))
            else:
                ps.append(self._nodes[p])
        
        # create the node
        n = BinaryRandomVariable(name,ps,pt)
        
        # we add child references to the parents
        for p in ps:
            p.add_child(name)

        # insert the node
        self._nodes[name] = n
        
    
    def get_binary_random_variable(self, name):
        """Returns the node for a given random variable"""
        return self._nodes[name]
        
    
    def observe_node(self,node,value):
        """Sets the value of the given node"""
        self._nodes[node].observe(value)
        
    def observe_nodes(self,node_values={}):
        """Sets several values to the given observed values.
        Input is a dict from node name to observation (True or False)"""
        for n,v in node_values.items():
            self.observe_node(n,v)
            
    def sample_node(self,node):
        """Samples the given node using a random number"""
        return self._nodes[node].sample()
    
    def joint_probability(self):
        """
        Compute joint probability of all nodes.
          This is done using the equation: 
            p(\mathbf{X}) = \prod_{i=1}^{N} p(\mathbf{X}_i | pa_i)
        """
        if self.ordering is None:
            self.validate()
        px = 1.0
        # we iterate over the nodes in ancestral order
        for k in self.ordering:
            node = self._nodes[k]
            if node.is_root():
                # no ancestors involved
                px *= node.p()
                
            else:
                # generate the probability conditions for all parents being true:
                conditions = tuple([True for i in xrange(len(node.parents))])
                # get the probability of this node given all parents sampled True
                px *= node.p(conditions)
                
        return px
        
    def print_joint_probability(self):
        """Computes and prints the joint probability for the network"""
        jp = self.joint_probability()
        print "Joint probability: \n\tp(%s) = %s" % (','.join([n for n in self._nodes.keys()]), jp)
        return jp


    def ancestral_sample(self, selected=None):
        """Assigns values to all variables using ancestral sampling
        If selected contains a list of nodes, then only their assignments will be returned."""
        if self.ordering is None:
            self.validate()
        
        if selected is None:
            for k in self.ordering:
                node = self._nodes[k]
                node.sample()
            return [(name,n.sample()) for (name,n) in self._nodes.items() ]
        else:
            # if only interested in a subset of nodes,
            #  stop sampling process when they have been sampled
            remaining = list(selected)
            for k in self.ordering:
                node = self._nodes[k]
                node.sample()

                if k in remaining:
                    remaining.remove(k)
                    
                # check if further sampling needed
                if len(remaining) == 0:
                    return [(name, self.sample_node(name)) for name in selected ]
                
                
            
        
    def print_ancestral_sample(self):
        sampling = self.ancestral_sample()
        print  "Ancestral sample: \n \t %s" % '\n\t '.join(["%s: %s"%(n,v) for (n,v) in sampling])
        return sampling
    
    def print_multi_sample(self,N):
        """
        Performs N ancestral samples and computes the frequencies of each state.
        The results are printed (along with proportions of the total)
        """
        stats = {}
        for n in xrange(N):
            sample = self.ancestral_sample()
            count = stats.setdefault(tuple(sample),0)
            stats[tuple(sample)] += 1
            self.reset_sampling()

        stats = OrderedDict(sorted(stats.items(), key=lambda x: x[1],reverse=True))
        print "Frequencies after %s samples: \n\t" %N, "\t".join(["%80s: %4s (%4s)\n" % (sample,stats[sample],stats[sample]/float(N)) for sample in stats.keys() ])
        return stats
   
    def markov_blanket(self,node):
        """
        Identifies the markov blanket of a given variable.
        This is the set of the parents, the children and the co-parents
        of the children.
        """
        if self.ordering is None:
            self.validate()

        n = self._nodes[node] # raises exception for missing node
        mb = set()
        
        for p in n.parents:
            mb.add(p.name)
  
        for c in n.children:
            mb.add(c)
            cps = ([p.name for p in self._nodes[c].parents])
            for cp in cps:
                if not cp == node:
                    mb.add(cp)
            
        return mb
        
    def markov_blanket_sampling(self,node,assignments):
        """Generates a sample for a given variable given the 
        markov blanket assignments.
        
        Takes the variable name, and the assignments to variables.
        The blanket variables should be assigned, but this isn't checked
        """
        n = self._nodes[node]
        mb = self.markov_blanket(node)
        
        # set the interest node to true for the sampling
        n.observe(True)
        print n,n._set
        # Set the assignments
        self.observe_nodes(assignments)
        

        numerator = (n.get_probability()) 
        for p in n.children:
            # print p,": %s" % self.get_node(p).get_probability()
            numerator *= self.get_node(p).get_probability()
           
        # print numerator
        # set the variable to false
        n.observe(False)
        p_false = 1 - n.get_probability()
        # print "Prob false %s" % p_false
        p_not =  n.get_probability() * self.product([self.get_node(p).get_probability() for p in n.children])
        
        # set the variable to true
        n.observe(True)
        # print "Prob true %s" % n.get_probability()
        p = n.get_probability() * self.product([self.get_node(p).get_probability() for p in n.children])
        denominator = p_not + p

        
        p_n_given_mb = (numerator/denominator)
        
       
        print "p(%s=True | %s) = %s" % (node,assignments, p_n_given_mb)
        rp = random.uniform(0.0,1.0)
        val = rp < p_n_given_mb
        n.observe(val)
        return val
 
    def rejection_sampling(self,evidence={},N=100):
        """
        If any variables are observed, then any samples 
        which do not agree are ignored until we find one that does
        Note that if no variables are observed this is equivalent to 
        ancestral_sample()
        
        evidence is the dictionary of assignments (can be empty)
        N is the number of samples to generate
        
        """
        print "Rejection sampling given evidence: \n\t"
        print "\n".join(["%20s: %10s" % (v,a) for v,a in evidence.items()])
        e_nodes = evidence.keys()
        stats = {}
        N_orig = N
        N_attempts = 0
        
        while N>0:
            failed = False
            self.reset_sampling()
            s = self.ancestral_sample()
            
            # verify against evidence:
            for n in e_nodes:
                samples = dict((x,y) for x, y in s)
                if not samples[n] == evidence[n]:
                    failed = True

            if not failed:
                count = stats.setdefault(tuple(s),0)
                stats[tuple(s)] += 1
                N -= 1
            N_attempts +=1
            
        stats = OrderedDict(sorted(stats.items(), key=lambda x: x[1],reverse=True))
        print "Rejected %s samples out of %s" %(N_attempts-N_orig,N_attempts)
        print "Sample frequencies after %s samples: \n\t" %N_orig, "\t".join(["%80s: %4s (%4s)\n" % (sample,stats[sample],stats[sample]/float(N_orig)) for sample in stats.keys() ])
        return stats
           
    def reset_sampling(self):
        """
        Removes all observations from the network
        """
        for (name,rv) in self._nodes.items():
            rv.reset()
    
    def reset_node(self,node):
        self._nodes[node].reset()
        
    def print_ancestral_order(self):
        print "Ancestral ordering: \n\t", self.ordering
    
    # Utility functions
    def number_nodes(self):
        return len(self._nodes)
    
    def validate(self):
        """
        Sets the ancestral order for all the nodes
        """
        if not self.ordering is None:
            return
        
        print "Setting ancestral order for network"
        self.ordering = {}
        nextID = 1
        roots = []
        for i in self._nodes:
            if self._nodes[i].is_root():
                roots.append(i)
                self.ordering[i] = nextID
                self._nodes[i].set_order(nextID)
                nextID += 1
        
        # order the next level of nodes
        self.order_nodes(roots,nextID)
   
    def order_nodes(self,parents,nextID):
        """ Recursive method for setting ancestral order for
        all the descendant nodes for the given parents"""
        nextlevel = []
        for n in parents:
            for c in self._nodes[n].children:
                # only assign once for each node
                if not c in nextlevel:
                    nextlevel.append(c)
                    
        # order the nextlevel:
        for p in nextlevel:
            self.ordering[p] = nextID
            self._nodes[p].set_order(nextID)
            
            nextID += 1
            
        # recursive call:
        if len(nextlevel) > 1:
            self.order_nodes(nextlevel,nextID)
    
    def node_sampled(self,node):
        """
        Return true is the given node is sampled or observed
        """
        return not self._nodes[node]._set is None
    
    def get_node(self,node):
        return self._nodes[node]
            
    def product(self,numbers):
        """
        """
        r = 1.0
        for n in numbers:
            r *= n
        return r
        
    def __str__(self):
        return "BN (%s nodes)" % self.number_nodes()
        


class BinaryRandomVariable:
    """ A BinaryRandomVariable is a random variable that can take 2 states:
     - True or False - with some probability p(True)
     
     The variable can have N parents, in which case the probability table (PT)
     must have size 2^N. That is, you must enumerate the probability of this variable
     being true, given all the combinations of the parent variables
    """
    def __init__(self, name, parents=list(), pt={}):
        self.name = name
        self.parents = parents

        # number of children can vary, so make it an array
        self.children = []
        
        # verify the PT dimensions
        # Since this is a binary variable there are 2^N (N = num. parents)
        if not (len(pt) == 2**len(parents) ):
            raise Exception("Wrong size probability table for %s parents, should be %s" % (len(parents),2**len(parents)))
       
        # store the probability table, set the initial state to indicate unsampled
        self.pt = pt
        self._set = None


    def is_root(self):
        return len(self.parents) == 0
        
    def is_leaf(self):
        return len(self.children) == 0
    
    def set_order(self,order_index):
        """Stores the topographical order of this node locally
        We assume that all probability tuples are ordered using this order field.
        It's mainly for bookkeeping that we store this
        """
        self._order = order_index
        
        
    def observe(self,value):
        """
        Observes the value of this node with the given value.
        
        """
        if not isinstance(value,bool):
            raise Exception("Binary variables can only be observed as True or False")
        self._set = value
        
    
    def p(self,conditions={}):
        """
        Get the probability of this event, given conditions 
        (if this is a root node, then just returns p(true))
        
        The conditions should be a tuple of truth values ordered by parents.
        Note the order of the conditions is assumed to be the same as the order  
        used when creating the random variable.
        """
        if self.is_root():
            return self.pt[(True,)]
        else:
            # we do still have a risk that the parents were supplied out of order
            return self.pt[conditions]


    def sample(self):
        """Take a sample for this node. Generated using PRNG and ancestral sampling"""
        
        # generate a random probability
        rp = random.uniform(0.0,1.0)
        
        # if the node was already sampled, we just return that value
        if not self._set is None:
            return self._set

        if self.is_root():
            self._set = rp < self.pt[(True,)]
            return self._set
        
        # when there are parents:   
        samples = [None for p in self.parents]
        for i in xrange(len(self.parents)):
            samples[i] = self.parents[i].sample()
            
        # look up the probability based on the parents samples
        conditions = tuple(samples)
        self._set = rp < self.pt[conditions]
        return self._set

    def get_probability(self):
        """
        Similar to sample(), but just returns probability based on parents and current set state
        """
        if self.is_root():
            return self.pt[(True,)]

        # when there are parents:   
        samples = [None for p in self.parents]
        for i in xrange(len(self.parents)):
            samples[i] = self.parents[i].sample()
            
        # look up the probability based on the parents samples
        conditions = tuple(samples)
        
        # if this variable is set, then return the prob:
        if not self._set is None:
            if self._set:
                return self.pt[conditions]
            else:
                return 1-self.pt[conditions]
        # otherwise just the prob that it is true given the ancestors:
        return self.pt[conditions]
        
        
    def reset(self):
        """Clear any sampled observation"""
        self._set = None
        

    def parents_orders(self):
        """Get the ancestral ordering for the parents.
        Not currently used, but might be useful in future"""
        return [p._order for p in self.parents]
        
   
        
    def add_child(self, child):
        """Adds a child to the node. Never adds duplicates"""
        if not child in set(self.children):
            self.children.append(child)
         
             
    def __str__(self):
        if self.is_root():
            return "Root: p(%s) = %s" % (self.name,self.pt[(True,)])
        return "Node: p(%s | %s)" % (self.name, [p.name for p in self.parents])

  
