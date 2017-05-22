# Copyright 2013 Massachusetts Institute of Technology
# BSD-2-Clause license https://github.com/agitter/MultiPCSF/blob/master/LICENSE

import optparse
import sys
import os
import math
import random
import GenerateNetwork
import GenerateSamples
import ConstrainedMultiSample
import Pooled
import NetworkUtil

__author__ = 'Anthony Gitter'

def main(argList):
    # Parse the arguments, which either come from the command line or a list
    # provided by the Python code calling this function
    parser = CreateParser()
    (opts, args) = parser.parse_args(argList)

    print "Parameters: %s" % opts

    # Ensure that pooled learning and constrained multi-sample learning are not both enabled
    if (not opts.pool == "none") and opts.iterations > 1:
        raise RuntimeError("Cannot run pooled learning and constrained multi-sample learning simultaneously.  Set pool=none or iterations=1")

    # Create the output path if needed and write the parameters to a file
    if not os.path.exists(opts.outPath):
        print "Creating output directory %s" % opts.outPath
        os.makedirs(opts.outPath)
    with open(os.path.join(opts.outPath,"params.log"), "w") as logFile:
        logFile.write("%s\n" % opts)


    # Set a seed for the random number generator, which may be None if one wasn't specified
    random.seed(opts.seed)


    # Generate the interaction network if needed
    if opts.networkSource == "generate":
        # Create a pseudo command line argument list
        genNetArgs = ["--name", opts.networkName, "--outpath", opts.outPath, "--model", opts.networkModel, "--n", opts.n, "--m", opts.m]
        GenerateNetwork.main(map(str, genNetArgs))
        networkFile = os.path.join(opts.outPath,opts.networkName) + "_%s_n%d_m%d.txt" % (opts.networkModel, opts.n, opts.m)
    elif opts.networkSource == "load":
        networkFile = opts.networkName
    else:
        # Shouldn't be able to reach this case
        raise RuntimeError("%s is not a recognized networksource" % opts.networkSource)


    # Generate the samples
    genSamplesArgs = ["--outpath", opts.outPath, "--name", opts.sampleName, "--networkfile", networkFile, "--pathwaysource", opts.pathwaySource, "--samples", opts.samples, "--fraction", opts.fraction, "--noise", opts.noise, "--samplegroups", opts.sampleGroups]
    if opts.pathwaySource == "load":
        genSamplesArgs.extend(["--pathwaypath", opts.pathwayPath, "--pathwaylistfile", opts.pathwayListFile])
    elif opts.pathwaySource =="generate":
        genSamplesArgs.extend(["--numpathways", opts.numPathways, "--branching", opts.branching, "--depth", opts.pathwayDepth])
    else:
        # Shouldn't be able to reach this case
        raise RuntimeError("%s is not a recognized pathwaysource" % opts.pathwaySource)
    GenerateSamples.main(map(str, genSamplesArgs))


    # Run prize collecting Steiner forest, commented out the options that are not yet supported
    pcsfArgs = ["--interactomepath", os.path.dirname(networkFile)]
    pcsfArgs.extend(["--terminalpath", opts.outPath])
    pcsfArgs.extend(["--resultpath", opts.outPath])
    pcsfArgs.extend(["--undirectedfile", os.path.basename(networkFile)])
    pcsfArgs.extend(["--beta", opts.beta])
    pcsfArgs.extend(["--terminalfile", "%s_sampleList.txt" % opts.sampleName])
    #pcsfArgs.extend(["--directedfile", "None"])
    #pcsfArgs.extend(["--tfdnafile", "None"])
    #pcsfArgs.extend(["--mrnabeta, 1.0])
    pcsfArgs.extend(["--msgpath", opts.msgPath])
    pcsfArgs.extend(["--depth", opts.msgDepth])
    pcsfArgs.extend(["--W", opts.W])
    pcsfArgs.extend(["--workers", opts.workers])
    pcsfArgs.extend(["--dummyneighbors", opts.dummyNeighbors])
    # TODO Pooled.py needs to support --mu.  Will cause an error now.
    pcsfArgs.extend(["--mu", opts.m])

    # Use pooled learning or constrained multi-sample learning for the Steiner forests
    if opts.pool == "none":
        pcsfArgs.extend(["--iterations", opts.iterations])
        if opts.iterations > 1:
            pcsfArgs.extend(["--lambda1", opts.lambda1])
            pcsfArgs.extend(["--lambda2", opts.lambda2])
            pcsfArgs.extend(["--artificialprizes", opts.artificialPrizes])
            pcsfArgs.extend(["--itermode", opts.iterMode])

        ConstrainedMultiSample.main(map(str, pcsfArgs))
    else:
        # Each way of aggregating prizes is a separate call to Pooled.py
        if "average" in opts.pool:
            poolArgs = list(pcsfArgs)
            poolArgs.extend(["--pool", "average"])
            Pooled.main(map(str, poolArgs))
        if "max" in opts.pool:
            poolArgs = list(pcsfArgs)
            poolArgs.extend(["--pool", "max"])
            Pooled.main(map(str, poolArgs))


    # Evaluate the Steiner forests
    # Use the map from samples to pathways to load the Steiner forest built from the sample and then
    # compare it to the pathway it was generated from
    print "Evaluating Steiner forests"
    if opts.pool == "none":
        EvaluateForests("itr1", opts)
        # Also iterate the final iteration forests (with and without pruning) if multiple iterations were run
        if opts.iterations > 1:
            EvaluateForests("itr%s" % opts.iterations, opts)
            if opts.artificialPrizes == "positive" or opts.artificialPrizes == "positiveWeighted":
                EvaluateForests("final", opts)
    else:
        if "average" in opts.pool:
            EvaluatePooled("average", opts)
        if "max" in opts.pool:
            EvaluatePooled("max", opts)

    # Evaluate the entire interactome to update upper bounds on pathway recall
    EvaluateUpperBound(networkFile, opts)


# Evaluate the Steiner forests learned when pooling the prizes with the specified type of aggregation
def EvaluatePooled(poolType, opts):
    # Can only evaluate pooled runs if all groups contain samples generated from the same pathway
    if opts.sampleGroups == "pathway":
        steinerPath = os.path.join(opts.outPath, poolType)
        # The list of pooled Steiner forest output files and their corresponding pathway files
        forest2Pathway = LoadPooledTuples(opts.outPath, "%s_sampleList.txt" % opts.sampleName, "%s_sampleMap.txt" % opts.sampleName, steinerPath, opts.W, opts.msgDepth)
        Evaluate(forest2Pathway, os.path.join(opts.outPath, "evaluation_%s.txt" % poolType), opts.fraction, opts.noise, False)
    else:
        print "Must group samples by pathway in order to evaluate pooled runs"

# Calcualte an upper bound on pathway node and edge recall
# Only considers an undirected network for now
def EvaluateUpperBound(networkFile, opts):
    network2Pathway = []
    pathways = set()
    with open(os.path.join(opts.outPath, "%s_sampleMap.txt" % opts.sampleName)) as smFile:
        for line in smFile:
            parts = line.split()
            # parts[0] is the sample name, parts[1] is the full pathway path
            pathways.add(parts[1])
    for pathway in pathways:
        network2Pathway.append((networkFile, pathway))
    Evaluate(network2Pathway, os.path.join(opts.outPath, "evaluation_interactome.txt"), 0, 0, True)


# Evaluate the Steiner forests in the subdirectory in the outPath specified by itrName
def EvaluateForests(itrName, opts):
    steinerPath = os.path.join(opts.outPath, itrName)
    # The list of Steiner forest output files and their corresponding pathway files
    forest2Pathway = LoadForestTuples(opts.outPath, "%s_sampleMap.txt" % opts.sampleName, steinerPath, opts.W, opts.msgDepth)
    Evaluate(forest2Pathway, os.path.join(opts.outPath, "evaluation_%s.txt" % itrName), opts.fraction, opts.noise, False)


# The core of the evaluation that computes and writes the overlaps.
# Takes a list of tuples of networks to evaluate and the gold standard for each network as input
# as well as the output file name.  Optionally can evaluate weighted networks.
def Evaluate(network2Pathway, outFileName, fraction, noise, weightedNetworks=False):
    with open(outFileName, "w") as outFile:
        npSum = 0
        nrSum = 0
        epSum = 0
        erSum = 0
        outFile.write("Steiner forest\tPathway\tTrue prizes\tNoisy prizes\tForest nodes\tPathway nodes\tIntersection nodes\tNode precision\tNode recall\tForest edges\tPathway edges\tIntersection edges\tEdge precision\tEdge recall\n")
        # The name forestFile assumes the networks to evaluate are Steiner forests, but they can
        # be any network        
        for forestFile, pathwayFile in network2Pathway:
            # For each Steiner forest, compute the precision and recall with respect to the original pathway
            forest = NetworkUtil.LoadNetwork(forestFile, weight=weightedNetworks)
            # Remove the artificial node if the forest is not empty
            if "DUMMY" in forest:
                forest.remove_node("DUMMY")
            # NetworkUtil.LoadNetwork only works for the simple format used when writing synthetic
            # pathways.  LoadGraphiteNetwork works for the simple format and the graphite edge list.
            pathway = NetworkUtil.LoadGraphiteNetwork(pathwayFile)
            intersection = NetworkUtil.Intersection(forest, pathway)
            if forest.order() == 0:
                nPrecision = 0
            else:
                nPrecision = float(intersection.order())/forest.order()
            npSum += nPrecision
            nRecall = float(intersection.order())/pathway.order()
            nrSum += nRecall
            if forest.size() == 0:
                ePrecision = 0
            else:
                ePrecision = float(intersection.size())/forest.size()
            epSum += ePrecision
            eRecall = float(intersection.size())/pathway.size()
            erSum += eRecall
            truePrizes = int(math.ceil(fraction*pathway.order()))
            noisyPrizes = int(math.ceil(noise*truePrizes))
            outFile.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%f\t%f\n" % (os.path.basename(forestFile), os.path.basename(pathwayFile), truePrizes, noisyPrizes, forest.order(), pathway.order(), intersection.order(), nPrecision, nRecall, forest.size(), pathway.size(), intersection.size(), ePrecision, eRecall))
        # Write the average node/edge precision/recall
        outFile.write("Average\t\t\t\t\t\t\t%f\t%f\t\t\t\t%f\t%f\n" % (npSum/len(network2Pathway), nrSum/len(network2Pathway), epSum/len(network2Pathway), erSum/len(network2Pathway)))


# Load a list of tuples from Steiner forest output files on the pooled prizes and the pathway files
# from which those prizes were sampled.
# Checks that all samples within a group were generated from the same pathway
def LoadPooledTuples(samplePath, sampleList, sampleMapFilename, steinerPath, W, depth):
    # Load the group-to-sample mapping.  The other two maps aren't used.
    terminalMap, sampleMap, countMap = ConstrainedMultiSample.LoadTerminalFiles(samplePath, sampleList)
    # Load a map from samples to the pathways they were generated from
    sample2Pathway = dict()
    with open(os.path.join(samplePath, sampleMapFilename)) as smFile:
        for line in smFile:
            parts = line.split()
            # parts[0] is the sample name, parts[1] is the full pathway path
            # Clean the name the same way ConstrainedMultiSample does
            name = parts[0]
            if name.endswith(".txt"):
                name = name[0:-4] # Remove ".txt"
            sample2Pathway[name] = parts[1]
    # Create tuples of Steiner forests from the pooled data to the reference pathway
    # There will be one entry per group
    forest2Pathway = []
    for group in sampleMap.iterkeys():
        forestFile = os.path.join(steinerPath, "symbol_pooled_%s_%s_1.0_%d.txt" % (group, str(W), depth))
        # All samples in this group
        sampleNames = sampleMap[group]
        # Set a default pathway, then verify it is the same for all pathways
        pathway = sample2Pathway[sampleNames[0]]
        for sampleName in sampleNames:
            if not pathway == sample2Pathway[sampleName]:
                raise RuntimeError("Error evaluating pooled Steiner forest.  %s has different reference pathways" % group)
        forest2Pathway.append((forestFile, pathway))
    return forest2Pathway


# Load a list of tuples of Steiner forest output files and the pathway files from which their prizes were sampled
def LoadForestTuples(outPath, sampleMapName, steinerPath, W, depth):
    forest2Pathway = []
    with open(os.path.join(outPath, sampleMapName)) as smFile:
        for line in smFile:
            parts = line.split()
            # parts[0] is the sample name, parts[1] is the full pathway path
            sampleName = parts[0]
            if sampleName.endswith(".txt"):
                sampleName = sampleName[0:-4] # Remove ".txt"
            forestFile = "symbol_%s_%s_1.0_%d.txt" % (sampleName, str(W), depth)
            forest2Pathway.append((os.path.join(steinerPath, forestFile), parts[1]))
    return forest2Pathway


# Setup the option parser
def CreateParser():
    parser = optparse.OptionParser(description="Load or generate a PPI network, generate or load pathways, sample nodes from the pathways, run prize collecting Steiner forest algorithm on the samples, and evaluate how well the pathways are recovered.")

    # General options
    parser.add_option("--outpath", type="string", dest="outPath", default=".", help="The output directory, which will be created if it does not exist.")
    parser.add_option("--networksource", type="choice", dest="networkSource", default="generate", choices=["generate","load"], help="'load' or 'generate' (default) the interaction network.")
    parser.add_option("--networkname", type="string", dest="networkName", default="None", help="If networksource is 'generate', this is the prefix of the network file that will be created in the output directory.  If networksource is 'load', this is the filename (including path) of the interaction network.")
    parser.add_option("--pool", type="choice", dest="pool", default="none", choices=["none", "average", "max", "average|max"], help="Pool prizes for all samples and learn a single Steiner forest on the pool data.  Set to 'none' (default) to disable pooling.  Set to 'average', 'max', or 'average|max' (runs both options) to enable pooling.  Cannot be used with constrained multi-sample learning.")

    # Network generation options
    genNetGroup = optparse.OptionGroup(parser, "Generate network options", "Options that are only required when the 'generate' option was selected for networksource")
    genNetGroup.add_option("--networkmodel", type="choice", dest="networkModel", default="ba", choices=["ba"], help="The random graph generation model to use.  Currently only the Barabasi-Albert (ba) preferential attachment algorithm for scale-free networks is supported.")
    genNetGroup.add_option("--n", type="int", dest="n", default=100, help="The number of nodes in the graph")
    genNetGroup.add_option("--m", type="int", dest="m", default=10, help="Number of edges to attach from a new node to existing nodes.  The network will have m*(n-m) edges.")
    parser.add_option_group(genNetGroup)

    # Sample generation options
    parser.add_option("--samplename", type="string", dest="sampleName", default="None", help="A name used to generate two output files in the output directory.  One lists all of the samples that were generated.  The other maps the samples to the pathway they were generated from.  Also used if generating pathways.")
    # Will be obtained from the other options
    #parser.add_option("--networkfile", type="string", dest="networkFile", default="None", help="The filename (including path) of the undirected interaction network.")
    parser.add_option("--pathwaysource", type="choice", dest="pathwaySource", default="generate", choices=["generate","load"], help="'load' or 'generate' (default) the pathways that will be sampled from.")
    parser.add_option("--samples", type="int", dest="samples", default=10, help="The number of samples that will be generated from each pathway.")
    parser.add_option("--fraction", type="float", dest="fraction", default=0.5, help="The fraction of nodes in a pathway that will be written as prizes.")
    parser.add_option("--noise", type="float", dest="noise", default=0, help="Controls the number of noisy prizes, which are prizes that are not true pathway members.  This parameter is relative to the number of true prizes, which is controlled by the 'fraction' parameter, and can be greater than 1.")
    parser.add_option("--samplegroups", type="choice", dest="sampleGroups", default="none", choices=["none", "pathway"], help="Group samples and constrain forest learning so that only samples within the same group have similar forests.  'none' (default) places all samples in the same group and 'pathway' groups samples by the pathway they were drawn from.")

    # Options that are only needed when loading pathways
    loadGroup = optparse.OptionGroup(parser, "Load pathway options", "Options that are only required when the 'load' option was selected for pathwaysource")
    loadGroup.add_option("--pathwaypath", type="string", dest="pathwayPath", default="None", help="The directory that contains the pathway files.")
    loadGroup.add_option("--pathwaylistfile", type="string", dest="pathwayListFile", default="None", help="The filename (including path) of a list of which pathways in the pathwaypath directory should be used to generate samples.")
    parser.add_option_group(loadGroup)

    # Options that are only needed when generating pathways
    genGroup = optparse.OptionGroup(parser, "Generate pathway options", "Options that are only required when the 'generate' option was selected for pathwaysource")
    genGroup.add_option("--numpathways", type="int", dest="numPathways", default=1, help="The number of synthetic pathways to generate")
    genGroup.add_option("--branching", type="int", dest="branching", default=1, help="The maximum branching factor of nodes in the synthetic pathway.  Set to 1 to generate a linear cascade.")
    genGroup.add_option("--pathwaydepth", type="int", dest="pathwayDepth", default=5, help="The maximum depth from the root to leaves in the synthetic pathway.  This depth may not be achieved due to the toplogy of the interaction network.  The root has depth 0.")
    parser.add_option_group(genGroup)

    # Constrained Steiner forest algorithm options
    # Does not support post-translational modifications or TF-gene interactions at this time
    # Commented out options are not supported or obtained from prior steps in the pipeline
    #parser.add_option("--interactomepath", type="string", dest="interactomePath", help="This path points to the directory where all interaction files are deposited (i.e., protein-protein interactions, kinase-substrate interactions, transcription factor-DNA interactions)",default='None')
    #parser.add_option("--terminalpath",type="string",dest="terminalPath",help="This path points to the directory where the terminal files are deposited.",default='None')
    #parser.add_option("--resultpath",type="string",dest="resultPath",help="This path points to the directory where the outputs will be located.",default='None')
    #parser.add_option("--undirectedfile",type="string",dest="undirectedFile",help="The name of the interaction file where protein-protein interaction data with probabilistic weights (e.g in [0,1]) are available. Columns should be ordered [prot1 prot2 weight].",default='None')
    parser.add_option("--beta",type="float",dest="beta",help="Beta parameter given here is to scale node penalties of protein terminals to the edge costs.  This scaling is only performed once when the initial stp files are created.",default=1.0)
    #parser.add_option("--terminalfile",type="string",dest="masterTerminalFile",help="A file in terminalpath that lists the files that give the node prizes for each sample.  All listed filenames should be relative to terminal path.  If gene penalties are given in the terminal files, gene names should end with '_MRNA'",default='None')
    #parser.add_option("--directedfile",type="string",dest="directedFile",help="Optional: The name of the interaction file where directed interactions (i.e. kinase-substrate) with probabilistic weights are available.  Columns should be ordered [substrate kinase weight].",default='None')
    #parser.add_option("--tfdnafile",type="string",dest="tfdnaFile",help="Optional: The name of the interaction file where TF-DNA interactions with probabilistic weights are available. Columns should be ordered [TF gene weight].  Gene names should not end with '_MRNA' because '_MRNA' is automatically appended to them.",default='None')
    #parser.add_option("--mrnabeta",type="string",dest="mrnaBeta",help="Optional: The beta parameter for gene terminal nodes. Its default value is equal to the --beta value.  It is applied to all terminals whose name ends with '_MRNA'.  The scaling is only performed once when the initial stp files are created.",default='None')
    parser.add_option("--msgpath",type="string",dest="msgPath",help="This path points to the directory where the message-passing code is available",default='None')
    parser.add_option("--msgdepth",type="int",dest="msgDepth",help="Depth parameter for the Steiner tree message passing algorithm",default=10)
    parser.add_option("--W", type="float", dest="W",help="The cost of the edges from the artificial root node to its neighbors.",default=1.0)
    parser.add_option("--iterations", type="int", dest="iterations", help="The number of iterations to run if using constrained multi-sample prize collecting Steiner tree.  Set to 1 (the default) to learn standard prize collecting Steiner trees for each sample (set of prizes).", default=1)
    parser.add_option("--workers", type="int", dest="workers", help="The number of worker processes to use in the multiprocessing pool or threads to use in multi-threaded belief propagation.  Should not exceed the number of cores available.  Defaults to the number of CPUs.", default=-1)
    parser.add_option("--dummyneighbors", type="choice", dest="dummyNeighbors", default="prizes", choices=["prizes", "nonprizes"], help="Connect the dummy node to all 'prizes' (default) or 'nonprizes'.")
    parser.add_option("--seed", type="int", dest="seed", default=None, help="An integer seed to use for network generation, pathway generation, and sample generation as applicable.  Default is no seed, which gives pseudo-random behavior.")
    parser.add_option("--mu", type="float", dest="mu", default=0, help="A parameter used to penalize high-degree nodes from being selected as Steiner nodes.  Does not affect prize nodes but does affect artificial prizes.  The penalty is -mu*degree.  Set mu <= 0 to disable the penalty (default).  Does not yet work with pooled prizes.")

    multiGroup = optparse.OptionGroup(parser, "Constrained multi-sample options", "Options that are only required when the iterations option is greater than 1.")
    multiGroup.add_option("--lambda1", type="float", dest="lambda1",help="The tradeoff coefficient for the penalty incurred by nodes in the Steiner forests that are not in the set of common nodes.",default=1.0)
    multiGroup.add_option("--lambda2", type="float", dest="lambda2",help="The tradeoff coefficient for the reward on the size of the set of common nodes when using unweighted artificial prizes or the power to which the node frequency is taken (called alpha elsewhere) for weighted prizes.",default=1.0)
    multiGroup.add_option("--artificialprizes", type="choice", dest="artificialPrizes", default="negativeWeighted", choices=["positive", "negative", "positiveWeighted", "negativeWeighted"], help="Use 'positive' or 'negative' prizes to encourage trees to include common set proteins.  Use 'positiveWeighted' or 'negativeWeighted' (default) prizes to construct weighted artificial prizes based on the node frequency in the most recent forests.")
    multiGroup.add_option("--itermode", type="choice", dest="iterMode", default="batch", choices=["batch", "random"], help="Learn forests simultaneously in 'batch' (default) or sequentially in 'random' order.  Batch mode computes artificial prizes with respect to all forests at the previous iteration.  Random mode computes prizes for a specific sample given the most recent forests for all other samples.")
    parser.add_option_group(multiGroup)

    return parser

if __name__ == "__main__":
    # Use the command line arguments to setup the options (the same as the default OptionParser behavior)
    main(sys.argv[1:])
