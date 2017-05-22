# Copyright 2013 Massachusetts Institute of Technology
# BSD-2-Clause license https://github.com/agitter/MultiPCSF/blob/master/LICENSE

import optparse
import networkx
import sys
import os
import random
import math
import NetworkUtil

__author__ = 'Anthony Gitter'

def main(argList):
    # Parse the arguments, which either come from the command line or a list
    # provided by the Python code calling this function
    parser = CreateParser()
    (opts, args) = parser.parse_args(argList)

    print "Parameters: %s" % opts

    if opts.networkFile == "None":
        raise RuntimeError("Must specify an network filename")

    if opts.pathwaySource == "load" and (opts.pathwayPath == "None" or opts.pathwayListFile == "None"):
        raise RuntimeError("Must specify pathwayPath and pathwayListFile when loading pathways")

    # Create the output path if needed
    if not os.path.exists(opts.outPath):
        print "Creating output directory %s" % opts.outPath
        os.makedirs(opts.outPath)


    # Load the interaction network
    network = NetworkUtil.LoadNetwork(opts.networkFile, weight=True)

    # Load or generate the pathways
    if opts.pathwaySource == "load":
        pathways = LoadPathways(opts.pathwayPath, opts.pathwayListFile)
    elif opts.pathwaySource == "generate":
        pathways = GeneratePathways(network, opts.numPathways, opts.branching, opts.depth, opts.outPath, opts.name)
    else:
        # Shouldn't be able to get to this case
        raise RuntimeError("%s is not a recognized pathway source" % opts.pathwaySource)
    

    # Sample from the pathways
    CreateSamples(pathways, opts.samples, opts.fraction, opts.outPath, opts.name, opts.noise, opts.sampleGroups, set(network.nodes()))


# Generate and write samples.  Each sample receives a prize of 1
# (true prizes as well as noisy prizes).
# Writes a list of all samples and a map from samples to pathways.
# Sample groups will be written to the sample list if desired.
def CreateSamples(pathways, samples, fraction, outPath, name, noise, sampleGroups, networkNodes):
    if fraction <= 0 or fraction > 1:
        raise RuntimeError("fraction must be in (0, 1]")
    if noise < 0:
        raise RuntimeError("noise must be non-negative")

    with open(os.path.join(outPath,"%s_sampleList.txt" % name), "w") as listFile:
        with open(os.path.join(outPath,"%s_sampleMap.txt" % name), "w") as mapFile:

            for pathway in pathways:
                # The number of true prizes
                numPrizes = int(math.ceil(fraction * pathway.order()))

                # Calculate the set of possible noisy nodes for each pathway and
                # how many noisy prizes to sample
                potentialNoisy = networkNodes.difference(set(pathway.nodes()))
                numNoisy = int(math.ceil(noise * numPrizes))
                if len(potentialNoisy) < numNoisy:
                    raise RuntimeError("Cannot sample %d noisy prizes" % numNoisy)
                
                for i in range(1, samples+1):
                    prizes = random.sample(pathway.nodes(), numPrizes)
                    if numNoisy > 0:
                        noisyPrizes = random.sample(potentialNoisy, numNoisy)
                    else:
                        noisyPrizes = set()

                    # Write the prizes
                    filename = "%s_sample%d.txt" % (CleanString(pathway.graph["name"]), i)
                    fullFilename = os.path.join(outPath, filename)
                    with open(fullFilename, "w") as f:
                        for prize in prizes:
                            f.write("%s 1\n" % prize)
                        for noisyPrize in noisyPrizes:
                            f.write("%s 1\n" % noisyPrize)

                    # Add the sample to the summary output files
                    group = ""
                    if sampleGroups == "pathway":
                        group = "\tgroup_%s" % pathway.graph["name"]
                    listFile.write("%s%s\n" % (filename, group))
                    mapFile.write("%s\t%s\n" % (filename, pathway.graph["filename"]))


# Remove special characters from a string
def CleanString(str):
    out = str.replace("(", "-")
    out = out.replace(")", "-")
    return out


# Generate and write synthetic pathways using depth first search.
# Return a dictionary mapping pathway objects to filenames
def GeneratePathways(network, numPathways, branching, depth, outPath, runName):
    if depth < 1:
        raise RuntimeError("depth must be positive")
    if branching < 1:
        raise RuntimeError("branching must be positive")

    pathways = []
    for i in range(1, numPathways+1):
        name = "%s_syntheticPathway%d" % (runName, i)
        root = random.choice(network.nodes())
        print "Root is %s" % root
        pathway = DFS(root, networkx.Graph(), network, branching, depth)
        pathway.graph["name"]=name
        print "Generated %s with %d nodes and %d edges" % (name, pathway.order(), pathway.size())
        pathways.append(pathway)

        # Sanity checks to verify a valid pathway was generated
        # Make sure the pathway is a tree
        if not pathway.order() == pathway.size() + 1:
            raise RuntimeError("Error generating pathway.  Pathway is not a tree.")
        # The tree can't have more nodes than a perfect tree with the given depth and branching factor,
        # which is given by a geometric series starting at 0 see
        # (http://en.wikipedia.org/wiki/Summation#Some_summations_involving_exponential_terms) (formula
        # simplified because branching is always positive)
        if branching > 1:
            nodeBound = (branching**(depth+1)-1)/(branching-1)
        else:
            nodeBound = depth+1
        if pathway.order() > nodeBound:
            raise RuntimeError("Error generating pathway.  Too many nodes.")
        # No node can be more than the maximum depth away from the root
        if networkx.eccentricity(pathway, v=root) > depth:
            raise RuntimeError("Error generating pathway.  Maximum depth exceeds %d." % depth)

        # Store the pathway
        filename = os.path.join(outPath,"%s.txt" % name)
        pathway.graph["filename"]=filename
        with open(filename, "w") as f:
            for n1, n2 in pathway.edges_iter():
                f.write("%s %s\n" % (n1, n2))

    return pathways


# Helper function for the depth first search used to generate synthetic pathways
def DFS(curNode, pathway, network, branching, depth):
    visited = 0
    # Iterate through the neighbors of curNode in a random order
    neighbors = network.neighbors(curNode)
    random.shuffle(neighbors)
    for neighbor in neighbors:
        if visited >= branching:
            return pathway

        # To constrain the pathway to have a tree structure, do not
        # visit neighbors that are already in the pathway
        if not neighbor in pathway:
            # Add the edge to the pathway
            pathway.add_edge(curNode, neighbor)
            #print "Adding %s %s" % (curNode, neighbor) # Debugging
            visited += 1
            # See if we should continue the DFS from the neighbor
            if depth > 1:
                pathway = DFS(neighbor, pathway, network, branching, depth-1)

    # Return the pathway if there are no more neighbors to visit and the maximum
    # branching factor was not reached
    return pathway


# Load pathways from file as undirected graphs.  The file should not have a header and
# must have the two node names in the first two whitespace-separated columns.
# All columns besides the node names are ignored.  Does not check whether
# the node names match the node names in the network or if the pathway nodes/edges
# are present in the network.
def LoadPathways(pathwayPath, listFile):
    pathways = []
    with open(listFile) as inFile:
        for pathwayLine in inFile:
            pathwayLine = pathwayLine.strip()
            # Each line is a relative path to a pathway file
            pathway = NetworkUtil.LoadGraphiteNetwork(os.path.join(pathwayPath, pathwayLine))

            pathway.graph["filename"] = os.path.join(pathwayPath, pathwayLine)
            if pathwayLine.endswith(".txt"):
                pathwayLine = pathwayLine[0:-4] # Remove ".txt"
            pathway.graph["name"] = pathwayLine

            # Debugging
            print "Loaded %s with %d nodes and %d edges" % (pathway.graph["name"], pathway.order(), pathway.size())
            # Add the pathway to the list
            pathways.append(pathway)
    return pathways


# Setup the option parser
def CreateParser():
    parser = optparse.OptionParser(description="Generate prizes from canonical or synthetic pathways")
    parser.add_option("--outpath", type="string", dest="outPath", default=".", help="The output directory, which will be created if it does not exist.")
    parser.add_option("--name", type="string", dest="name", default="None", help="A name used to generate two output files in the output directory.  One lists all of the samples that were generated.  The other maps the samples to the pathway they were generated from.  Also used if generating pathways.")
    parser.add_option("--networkfile", type="string", dest="networkFile", default="None", help="The filename (including path) of the undirected interaction network.")
    parser.add_option("--pathwaysource", type="choice", dest="pathwaySource", default="generate", choices=["generate","load"], help="'load' or 'generate' the pathways that will be sampled from.")
    parser.add_option("--samples", type="int", dest="samples", default=10, help="The number of samples that will be generated from each pathway.")
    parser.add_option("--fraction", type="float", dest="fraction", default=0.5, help="The fraction of nodes in a pathway that will be written as prizes.")
    parser.add_option("--noise", type="float", dest="noise", default=0, help="Controls the number of noisy prizes, which are prizes that are not true pathway members.  This parameter is relative to the number of true prizes, which is controlled by the 'fraction' parameter, and can be greater than 1.")
    parser.add_option("--samplegroups", type="choice", dest="sampleGroups", default="none", choices=["none", "pathway"], help="Group samples and constrain forest learning so that only samples within the same group have similar forests.  'none' (default) places all samples in the same group and 'pathway' groups samples by the pathway they were drawn from.")

    # Options that are only needed when loading pathways
    loadGroup = optparse.OptionGroup(parser, "Load pathway options", "Options that are only required when the 'load' option was selected")
    loadGroup.add_option("--pathwaypath", type="string", dest="pathwayPath", default="None", help="The directory that contains the pathway files.")
    loadGroup.add_option("--pathwaylistfile", type="string", dest="pathwayListFile", default="None", help="The filename (including path) of a list of which pathways in the pathwaypath directory should be used to generate samples.")
    parser.add_option_group(loadGroup)

    # Options that are only needed when generating pathways
    genGroup = optparse.OptionGroup(parser, "Generate pathway options", "Options that are only required when the 'generate' option was selected")
    genGroup.add_option("--numpathways", type="int", dest="numPathways", default=1, help="The number of synthetic pathways to generate")
    genGroup.add_option("--branching", type="int", dest="branching", default=1, help="The maximum branching factor of nodes in the synthetic pathway.  Set to 1 to generate a linear cascade.")
    genGroup.add_option("--depth", type="int", dest="depth", default=5, help="The maximum depth from the root to leaves.  This depth may not be achieved due to the toplogy of the interaction network.  The root has depth 0.")
    parser.add_option_group(genGroup)
    return parser

if __name__ == "__main__":
    # Use the command line arguments to setup the options (the same as the default OptionParser behavior)
    main(sys.argv[1:])
