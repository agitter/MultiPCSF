# Copyright 2013 Massachusetts Institute of Technology
# BSD-2-Clause license https://github.com/agitter/MultiPCSF/blob/master/LICENSE

from optparse import OptionParser
import networkx
import sys
import os

__author__ = 'Anthony Gitter'

def main(argList):
    # Parse the arguments, which either come from the command line or a list
    # provided by the Python code calling this function
    parser = CreateParser()
    (opts, args) = parser.parse_args(argList)

    if opts.name == "None":
        raise RuntimeError("Must specify an output filename")

    # Create the network
    if opts.model == "ba":
        network = CreateBA(opts.n, opts.m)
    else:
        # Shouldn't be able to get to this case
        raise RuntimeError("%s is not a recognized model" % opts.model)

    # Save the network
    filename = os.path.join(opts.outPath,opts.name) + "_%s_n%d_m%d.txt" % (opts.model, opts.n, opts.m)
    print "Writing %s" % filename
    with open(filename, "w") as f:
        # Generate the edge weights as we iterate through the edges
        # TODO support non-uniform edge weights
        for n1, n2 in network.edges_iter():
            # Transform from 0-indexed to 1-indexed nodes
            f.write("N%d N%d %f\n" % (n1+1, n2+1, 0.9))

# Return a Barabasi-Albert network
def CreateBA(n, m):
    print "Creating Barabasi-Albert network with n=%d and m=%d" % (n, m)

    graph = networkx.barabasi_albert_graph(n,m)

    if not networkx.is_connected(graph):
        raise RuntimeError("Network is not connected")

    print "Network has %d edges" % graph.size()

    return graph


# Setup the option parser
def CreateParser():
    parser = OptionParser(description="Use Python's NetworkX package to generate random connected networks")
    parser.add_option("--name", type="string", dest="name", default="None", help="The prefix of the output filename")
    parser.add_option("--outpath", type="string", dest="outPath", default=".", help="The output directory")
    parser.add_option("--model", type="choice", dest="model", default="ba", choices=["ba"], help="The random graph generation model to use.  Currently only the Barabasi-Albert (ba) preferential attachment algorithm for scale-free networks is supported.")
    parser.add_option("--n", type="int", dest="n", default=100, help="The number of nodes in the graph")
    parser.add_option("--m", type="int", dest="m", default=10, help="Number of edges to attach from a new node to existing nodes.  The network will have m*(n-m) edges.")
    return parser

if __name__ == "__main__":
    # Use the command line arguments to setup the options (the same as the default OptionParser behavior)
    main(sys.argv[1:])
