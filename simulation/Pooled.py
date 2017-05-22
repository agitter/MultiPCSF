# Copyright 2013 Massachusetts Institute of Technology
# BSD-2-Clause license https://github.com/agitter/MultiPCSF/blob/master/LICENSE

from optparse import OptionParser
import sys
import os
import time
import ConstrainedMultiSample

__author__ = 'Anthony Gitter'


# Pool prizes across samples and learn a single forest
def main(argList):
    # Parse the arguments, which either come from the command line or a list
    # provided by the Python code calling this function
    parser = CreateParser()
    (opts, args) = parser.parse_args(argList)

    print "Starting pooled Steiner forest learning %s" % time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print "Parameters: %s" % opts

    # TODO Add error checking of inputs

    # Load all of the proteins in the interactome, ignoring
    # genes.  These are needed if connecting the dummy node to
    # non-prizes
    allProts = ConstrainedMultiSample.LoadProteins(opts.interactomePath, opts.undirectedFile, opts.directedFile, opts.tfdnaFile)
        

    # Create a directory to hold the intermediate and final output with
    # the same name as the pooling type
    outPath = os.path.join(opts.resultPath,opts.pool)
    if not os.path.exists(outPath):
        os.makedirs(outPath)

    # The pooled prize file has to be relative to the directory that the original terminal files are in
    # so update the terminal file parameter in opts and store the original separately
    origTerminalPath = opts.terminalPath
    opts.terminalPath = outPath


    # Load the list of terminal files and the group-to-sample mapping
    # The sample names aren't needed
    terminalMap, sampleMap, countMap = ConstrainedMultiSample.LoadTerminalFiles(origTerminalPath, opts.masterTerminalFile)

    # Pool prizes, create stp files, and learn forests for all groups
    for group in terminalMap.iterkeys():
        print "%d samples in group %s" % (countMap[group], group)
        name = "pooled_%s" % group

        # Pool samples within each group, creating one prize file per group
        print "Pooling prizes for group %s using %s" % (group, opts.pool)
        terminalFiles = terminalMap[group]
        prizeFile = "pooledPrizes_%s.txt" % group
        PoolPrizes(terminalFiles, origTerminalPath, opts.pool, os.path.join(outPath,prizeFile))

        # Create the stp file for this group
        stpFile = ConstrainedMultiSample.CreateStp(opts, outPath, prizeFile, name)

        # Setup how the dummy node will be connected
        # to the network, either all prizes or all non-prizes (potential Steiner nodes)
        dnFile = "%s_dummyNeighbors.txt" % name
        ConstrainedMultiSample.DummyNeighbors(allProts, outPath, stpFile, dnFile, opts.dummyNeighbors)

        # Learn the Steiner forest for this group
        ConstrainedMultiSample.LearnSteiner(opts, outPath, outPath, name, dnFile, opts.workers)

    print "Finishing pooled Steiner forest learning %s" % time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())


# Pool prizes within a group using the specified aggregation type
# and write them to a file
def PoolPrizes(terminalFiles, terminalPath, poolType, outFile):
    if not (poolType == "average" or poolType == "max" or poolType == "sum"):
        raise RuntimeError("%s is not a recognized pooling strategy" % poolType)

    # Track which proteins have a prize in any terminal file
    prizeProts = set()
    # Store a list of prize dictionaries, one for each terminal file
    prizeDicts = []
    for terminalFile in terminalFiles:
        prizeDict = dict()
        with open(os.path.join(terminalPath,terminalFile)) as f:
            for line in f:
                parts = line.strip().upper().split()
                prizeProts.add(parts[0])
                prizeDict[parts[0]] = float(parts[1])
        prizeDicts.append(prizeDict)

    # For each protein, query the prize dictionaries to build a list of all prizes
    # for that protein (prizes default to 0) and perform the specified aggregation
    with open(outFile, "w") as out:
        for prot in prizeProts:
            prizeList = []
            for prizeDict in prizeDicts:
                if prot in prizeDict:
                    prizeList.append(prizeDict[prot])
                else:
                    prizeList.append(0.0)
            # Perform the aggregation
            if poolType == "average":
                prize = sum(prizeList)/len(prizeList)
            elif poolType == "max":
                prize = max(prizeList)
            elif poolType == "sum":
                prize = sum(prizeList)
            out.write("%s\t%f\n" % (prot, prize))


# Setup the option parser
def CreateParser():
    # Options that are fixed for the time being:
    # Connection type for the artificial node / exclude / targetfile
    # Output file names
    # Beta (only allowed to be set once initially, range of Beta is not supported)
    # Species
    # W cannot be given as a range of values
    parser = OptionParser()
    parser.add_option("--interactomepath", type="string", dest="interactomePath", help="This path points to the directory where all interaction files are deposited (i.e., protein-protein interactions, kinase-substrate interactions, transcription factor-DNA interactions)",default='None')
    parser.add_option("--terminalpath",type="string",dest="terminalPath",help="This path points to the directory where the terminal files are deposited.",default='None')
    parser.add_option("--resultpath",type="string",dest="resultPath",help="This path points to the directory where the outputs will be located.",default='None')
    parser.add_option("--undirectedfile",type="string",dest="undirectedFile",help="The name of the interaction file where protein-protein interaction data with probabilistic weights (e.g in [0,1]) are available. Columns should be ordered [prot1 prot2 weight].",default='None')
    parser.add_option("--beta",type="float",dest="beta",help="Beta parameter given here is to scale node penalties of protein terminals to the edge costs.  This scaling is only performed once when the initial stp files are created.",default=1.0)
    parser.add_option("--terminalfile",type="string",dest="masterTerminalFile",help="A file in terminalpath that lists the files that give the node prizes for each sample.  All listed filenames should be relative to terminal path.  If gene penalties are given in the terminal files, gene names should end with '_MRNA'.  Optinonally can include a tab-separated second column that assigns each sample to a group so the prizes are only pooled within the same group.",default='None')
    #parser.add_option("--resultfilename",type="string",dest="resultfilename",help="The name of the file where the combined information will be written to be used as input in the message passing tool.",default='None')
    parser.add_option("--directedfile",type="string",dest="directedFile",help="Optional: The name of the interaction file where directed interactions (i.e. kinase-substrate) with probabilistic weights are available.  Columns should be ordered [substrate kinase weight].",default='None')
    parser.add_option("--tfdnafile",type="string",dest="tfdnaFile",help="Optional: The name of the interaction file where TF-DNA interactions with probabilistic weights are available. Columns should be ordered [TF gene weight].  Gene names should not end with '_MRNA' because '_MRNA' is automatically appended to them.",default='None')
    parser.add_option("--mrnabeta",type="string",dest="mrnaBeta",help="Optional: The beta parameter for gene terminal nodes. Its default value is equal to the --beta value.  It is applied to all terminals whose name ends with '_MRNA'.  The scaling is only performed once when the initial stp files are created.",default='None')


    #parser.add_option("--outputpath", type="string", dest="outputpath", help="This path points to the directory where the output files will be written",default='None')
    parser.add_option("--msgpath",type="string",dest="msgPath",help="This path points to the directory where the message-passing code is available",default='None')
    parser.add_option("--depth",type="int",dest="depth",help="Depth parameter",default=10)
    #parser.add_option("--conn",type="string", dest="conn", help="How to connect the artificial node to the interactome: 1 to all nodes in the interactome, 2 to all node in the interactome, except the terminals, 3 to a given set of nodes in the interactome, 4 to a given set of nodes in the interactome, except the terminals.")
    #parser.add_option("--stppath", type="string", dest="stppath", help="This path points to the directory where the stp file is available.")
    #parser.add_option("--stpfile", type="string", dest="stpfile", help="The name of the stp file (without the file extension)")
    parser.add_option("--W", type="float", dest="W",help="The cost of the edges from the artificial root node to its neighbors.",default=1.0)
    #parser.add_option("--beta", type="string", dest="beta",help="the range and increment for the beta parameter. format: start_end_increment")
    #parser.add_option("--exclude", type="string", dest="exclude", help="If set to 1 and the connection type is 2, 3, or 4 then remove prizes from all nodes that are connected to the artificial node")
    #parser.add_option("--targetfile", type="string", dest="targetfile", help="a subset of the proteins in the interactome that will be connected to the artificial root in connection types 3 and 4",default='None')
    #parser.add_option("--species", type="string", dest="species", help="the organism that you are working on",default='None')


    parser.add_option("--workers", type="int", dest="workers", help="The number threads to use in multi-threaded belief propagation.  Should not exceed the number of cores available.  Defaults to the number of CPUs.", default=-1)
    # Determine how to connect the dummy node
    parser.add_option("--dummyneighbors", type="choice", dest="dummyNeighbors", default="prizes", choices=["prizes", "nonprizes"], help="Connect the dummy node to all 'prizes' (default) or 'nonprizes'.")
    parser.add_option("--pool", type="choice", dest="pool", default="average", choices=["average", "max", "sum"], help="How to pool prizes, either 'average' (default), 'max', or 'sum'.")
    return parser


if __name__=="__main__":
    # Use the command line arguments to setup the options (the same as the default OptionParser behavior)
    main(sys.argv[1:])
