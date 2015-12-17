from optparse import OptionParser
import sys
import os
import shutil
import itertools
import multiprocessing
import time
import random
import NetworkUtil

__author__ = 'Anthony Gitter'


# Iterate between finding a common set of nodes among the samples
# and learning Steiner forests constrained to use those
# common nodes whenever possible
def main(argList):
    # Parse the arguments, which either come from the command line or a list
    # provided by the Python code calling this function
    parser = CreateParser()
    (opts, args) = parser.parse_args(argList)

    print "Starting constrained multi-patient Steiner forest %s" % time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print "Parameters: %s" % opts


    # TODO Add error checking of inputs
    if opts.iterations < 1:
        raise RuntimeError("Must have at least 1 iteration")
    # TODO Should allow the option to run serially without the pool because a
    # pool with 1 worker is not efficient
    if opts.workers < 1:
        opts.workers = multiprocessing.cpu_count()

    # Assume negative prizes to implement the common set
    # and change if using positive common set prizes
    negativePrizes = True
    if "positive" in opts.artificialPrizes:
        negativePrizes = False

    # Assume unweighted prizes
    weightedPrizes = False
    if "Weighted" in opts.artificialPrizes:
        weightedPrizes = True

    # Assume batch mode
    batchMode = True
    if opts.iterMode == "random":
        batchMode = False

    # Load all of the proteins in the interactome, ignoring
    # genes.  The artificial prizes will be created for a subset of these nodes.
    allProts = LoadProteins(opts.interactomePath, opts.undirectedFile, opts.directedFile, opts.tfdnaFile)

    # Load the negative prizes for the degree penalties or an empty dictionary
    # if they aren't being used
    directedFile = "None"
    if opts.directedFile != "None":
        directedFile = os.path.join(opts.interactomePath, opts.directedFile)
    degPenalties = NetworkUtil.DegreePenalties(opts.mu, os.path.join(opts.interactomePath, opts.undirectedFile), directedFile)

    # Create the initial stp files
    # New directory to hold the original data before the iterations begin
    # These stp files will be read and updated at subsequent iterations
    initPath = os.path.join(opts.resultPath,"initial")
    if not os.path.exists(initPath):
        os.makedirs(initPath)


    # Load the list of terminal files and the sample-to-group mapping
    terminalMap, sampleMap, countMap = LoadTerminalFiles(opts.terminalPath, opts.masterTerminalFile)
    # Store the groups in a fixed order
    groups = sorted(terminalMap.iterkeys())
    for group in groups:
        print "%d samples in group %s" % (countMap[group], group)


    # Create a pool for creating .stp files and learning Steiner forests in parallel
    # using the specified number of workers.  Use it to create the initial
    # .stp files.  Even when running the subsequent iterations in random sequential
    # order, create a pool to learn the initial trees and final pruned trees (if applicable).
    print "Creating a pool with %d workers" % opts.workers
    pool = multiprocessing.Pool(opts.workers)
    initialStpMap = dict()
    for group in groups:
        terminalFiles = terminalMap[group]
        sampleNames = sampleMap[group]
        # opts and initPath are invariant arguments for each sample
        zippedArgs = itertools.izip(itertools.repeat(opts), itertools.repeat(initPath), terminalFiles, sampleNames)
        initialStpMap[group] = pool.map(CreateStpHelper, zippedArgs) # Blocks until all are finished


    # Store which proteins don't have prizes for each patient.
    # These are the nodes that could potentially be Steiner nodes for
    # each sample.  This can't be recovered from the stp files at later
    # iterations because both original prizes and artificial prizes will exist.
    # Also track how the dummy node will be connected
    # to the networks, either all prizes or all non-prizes (potential Steiner nodes)
    potentialSteinerMap = dict()
    dummyNeighborMap = dict()
    for group in groups:
        numSamples = countMap[group]
        sampleNames = sampleMap[group]
        initialStps = initialStpMap[group]
        potentialSteiner = [] # A list of sets
        dummyNeighborFiles = [] # A list of filenames
        for i in range(numSamples):
            dnFile = sampleNames[i] + "_dummyNeighbors.txt"
            dummyNeighborFiles.append(dnFile)
            potentialSteiner.append(DummyNeighbors(allProts, initPath, initialStps[i], dnFile, opts.dummyNeighbors))
        potentialSteinerMap[group] = potentialSteiner
        dummyNeighborMap[group] = dummyNeighborFiles


    itrPath = os.path.join(opts.resultPath,"itr1")
    if not os.path.exists(itrPath):
        os.makedirs(itrPath)

    # Initialize the artificial prizes to be an empty dictionary so that
    # we learn the initial trees independently
    artificialPrizes = dict()
    # Write the unused itr1 artificial prizes so that the files exist for post-processing
    for group in groups:
        NetworkUtil.WriteDict(os.path.join(itrPath,"artificialPrizes_%s.txt" % group), artificialPrizes)
    print "%d artificial prizes at iteration 1" % len(artificialPrizes)


    # Add the degree penalties to the initial stp files.  Pass in the empty artificial prize
    # dictionary, which won't have an effect.
    for group in groups:
        sampleNames = sampleMap[group]
        numSamples = countMap[group]
        potentialSteiner = potentialSteinerMap[group]
        dummyNeighborFiles = dummyNeighborMap[group]
        for i in range(numSamples):
            # Copy the dummy neighbors, which must be in the same directory as the stp file
            UpdateStp(artificialPrizes, degPenalties, potentialSteiner[i], initPath, itrPath, sampleNames[i])
            shutil.copyfile(os.path.join(initPath,dummyNeighborFiles[i]), os.path.join(itrPath,dummyNeighborFiles[i]))


    # Learn the first iteration Steiner forests in parallel
    # Run single-threaded belief propagation when using the worker pool
    lastForestMap = dict()
    for group in groups:
        numSamples = countMap[group]
        sampleNames = sampleMap[group]
        dummyNeighborFiles = dummyNeighborMap[group]
        zippedArgs = itertools.izip(itertools.repeat(opts), itertools.repeat(itrPath), itertools.repeat(itrPath),  sampleNames, dummyNeighborFiles, itertools.repeat(1))
        pool.map(LearnSteinerHelper, zippedArgs)
        lastForests = [] # A list of sets, where each set contains the Steiner forest nodes
        for i in range(numSamples):
            lastForests.append(LoadForestNodes("%s/symbol_%s_%s_1.0_%d.txt" % (itrPath, sampleNames[i], str(opts.W), opts.depth)))
        lastForestMap[group] = lastForests


    # Learn the forests at all remaining iterations and return the directory
    # that contains the forests from the last iteration.
    if opts.iterations > 1:
        if batchMode:
            itrPath = Batch(opts, pool, initPath, allProts, sampleMap, potentialSteinerMap, dummyNeighborMap, lastForestMap, countMap, weightedPrizes, negativePrizes, degPenalties)
        else:
            itrPath = RandSequential(opts, initPath, allProts, sampleMap, potentialSteinerMap, dummyNeighborMap, lastForestMap, countMap, weightedPrizes, negativePrizes, degPenalties)


    # Prune Steiner nodes from the forests that are not used to reach any prizes and
    # are only present because they were in the common set.
    # This is not necessary if only 1 iteration was run because in that case there
    # is no common set.
    # It is also not necessary if negative prizes were used.
    if opts.iterations > 1 and (not negativePrizes):
        print "Learning final forests"
        print "Pruning forests from %s" % itrPath
        finalPath = os.path.join(opts.resultPath,"final")
        if not os.path.exists(finalPath):
            os.makedirs(finalPath)

        # Nothing is returned by these operations so they can be performed
        # simultaneously independent of the groupings
        sampleNames = FlattenDict(sampleMap, groups)
        dummyNeighborFiles = FlattenDict(dummyNeighborMap, groups)
        potentialSteiner = FlattenDict(potentialSteinerMap, groups)

        for i in range(len(sampleNames)):
            forestFile = "%s/symbol_%s_%s_1.0_%d.txt" % (itrPath, sampleNames[i], str(opts.W), opts.depth)
            FilterStpEdges(forestFile, initPath, finalPath, sampleNames[i], degPenalties, potentialSteiner[i])
            shutil.copyfile(os.path.join(initPath,dummyNeighborFiles[i]), os.path.join(finalPath,dummyNeighborFiles[i]))

        zippedArgs = itertools.izip(itertools.repeat(opts), itertools.repeat(finalPath), itertools.repeat(finalPath),  sampleNames, dummyNeighborFiles, itertools.repeat(1))
        pool.map(LearnSteinerHelper, zippedArgs)

    print "Finishing constrained multi-patient Steiner forest %s" % time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())

    pool.close()


# Learn all forests in an iteration sequentially in a random order, asssigning
# prizes for a fixed sample given all of the forests for the other samples
# in the same group.
# Those other forests will be the most recently forest learned for the sample
# regardless of whether it was learned at iteration i or i-1.
def RandSequential(opts, initPath, allProts, sampleMap, potentialSteinerMap, dummyNeighborMap, lastForestMap, countMap, weightedPrizes, negativePrizes, degPenalties):
    print "Learning forests in random sequential mode"

    # Iterate (rounds 2+)
    itrPath = initPath
    for itr in range(2,opts.iterations+1):
        #lastPath = itrPath
        itrPath = os.path.join(opts.resultPath,"itr%d" % itr)
        if not os.path.exists(itrPath):
            os.makedirs(itrPath)

        # Only constrain the Steiner forests to be similar to other samples in the same group
        for group in sampleMap.iterkeys():
            sampleNames = sampleMap[group]
            numSamples = countMap[group]
            potentialSteiner = potentialSteinerMap[group]
            dummyNeighborFiles = dummyNeighborMap[group]
            lastForests = lastForestMap[group]

            if len(sampleNames) != numSamples or len(potentialSteiner) != numSamples or len(dummyNeighborFiles) != numSamples or len(lastForests) != numSamples:
                raise RuntimeError("Must have the same number of samples in group %s" % group)

            # Randomly choose the order in which to learn forests at this iteration
            order = range(numSamples)
            random.shuffle(order)

            # Write the order to a file
            with open(os.path.join(itrPath, "sampleOrder_%s.txt" % group), "w") as f:
                for index in order:
                    f.write("%d\t%s\n" % (index, sampleNames[index]))

            # Iterate over all samples in the random order
            for index in order:
                # Create artificial prizes for this sample using all N-1 lastForests
                otherLastForests = list(lastForests)
                otherLastForests.pop(index)
                if weightedPrizes:
                    # lambda2 is used as the alpha parameter
                    artificialPrizes = CreateWgtPrizes(allProts, otherLastForests, opts.lambda1, opts.lambda2, negativePrizes)
                else:
                    # Use all N-1 other sets of potential Steiner nodes
                    otherPotentialSteiner = list(potentialSteiner)
                    otherPotentialSteiner.pop(index)
                    artificialPrizes = CreateUnwgtPrizes(allProts, otherPotentialSteiner, otherLastForests, opts.lambda1, opts.lambda2, negativePrizes)
                NetworkUtil.WriteDict(os.path.join(itrPath,"%s_artificialPrizes.txt" % sampleNames[index]), artificialPrizes)

                # Update the stp file based on the artificial prizes and degree penalties and copy the dummy neighbors
                UpdateStp(artificialPrizes, degPenalties, potentialSteiner[index], initPath, itrPath, sampleNames[index])
                shutil.copyfile(os.path.join(initPath,dummyNeighborFiles[index]), os.path.join(itrPath,dummyNeighborFiles[index]))

                # Learn a new forest for this sample and update lastForests
                # All samples (besides the first and last in the random order) will use last forests
                # that are a mix of forests from this iteration and the previous iteration
                LearnSteiner(opts, itrPath, itrPath, sampleNames[index], dummyNeighborFiles[index], opts.workers)
                lastForests[index] = LoadForestNodes("%s/symbol_%s_%s_1.0_%d.txt" % (itrPath, sampleNames[index], str(opts.W), opts.depth))
            # Store all forests learned for this group at this iteration so they can be
            # retreived at the next iteration
            lastForestMap[group] = lastForests

    return itrPath


# Learn all forests in an iteration in batch in parallel using the multiprocessing pool.
# Update the artificial prizes once per iteration using the forests from the previous
# iteration in the same group
def Batch(opts, pool, initPath, allProts, sampleMap, potentialSteinerMap, dummyNeighborMap, lastForestMap, countMap, weightedPrizes, negativePrizes, degPenalties):
    print "Learning forests in parallel batch mode"

    # Iterate (rounds 2+)
    itrPath = initPath
    for itr in range(2,opts.iterations+1):
        #lastPath = itrPath
        itrPath = os.path.join(opts.resultPath,"itr%d" % itr)
        if not os.path.exists(itrPath):
            os.makedirs(itrPath)

        # Only constrain the Steiner forests to be similar to other samples in the same group
        for group in sampleMap.iterkeys():
            sampleNames = sampleMap[group]
            numSamples = countMap[group]
            potentialSteiner = potentialSteinerMap[group]
            dummyNeighborFiles = dummyNeighborMap[group]
            lastForests = lastForestMap[group]

            if len(sampleNames) != numSamples or len(potentialSteiner) != numSamples or len(dummyNeighborFiles) != numSamples or len(lastForests) != numSamples:
                raise RuntimeError("Must have the same number of samples in group %s" % group)

            # Update artificial prizes based on the forests from the previous iteration
            if weightedPrizes:
                # lambda2 is used as the alpha parameter
                artificialPrizes = CreateWgtPrizes(allProts, lastForests, opts.lambda1, opts.lambda2, negativePrizes)
            else:
                artificialPrizes = CreateUnwgtPrizes(allProts, potentialSteiner, lastForests, opts.lambda1, opts.lambda2, negativePrizes)
            NetworkUtil.WriteDict(os.path.join(itrPath,"artificialPrizes_%s.txt" % group), artificialPrizes)
            print "%d artificial prizes in group %s at iteration %d" % (len(artificialPrizes), group, itr)

            # Update the stp files based on the new artificial prizes and degree penalties
            # and copy the potential Steiner node files, which need to be in itrPath
            for i in range(numSamples):
                UpdateStp(artificialPrizes, degPenalties, potentialSteiner[i], initPath, itrPath, sampleNames[i])
                shutil.copyfile(os.path.join(initPath,dummyNeighborFiles[i]), os.path.join(itrPath,dummyNeighborFiles[i]))

            # Learn new Steiner forests in parallel
            zippedArgs = itertools.izip(itertools.repeat(opts), itertools.repeat(itrPath), itertools.repeat(itrPath),  sampleNames, dummyNeighborFiles, itertools.repeat(1))
            pool.map(LearnSteinerHelper, zippedArgs)
            lastForests = []
            for i in range(numSamples):
                lastForests.append(LoadForestNodes("%s/symbol_%s_%s_1.0_%d.txt" % (itrPath, sampleNames[i], str(opts.W), opts.depth)))
            lastForestMap[group] = lastForests

    return itrPath


# TODO Could reduce the edge costs so that subtrees that have some high weight
# artificial prizes and low weight real prizes are not pruned
# Filter the initial stp file so that it includes only edges
# that were selected in a previous run of the Steiner forest
# algorithm.  Add the node degree penalties.
def FilterStpEdges(forestFile, inPath, outPath, sampleName, degreePenalties, potentialSteiner):
    edges = LoadForestEdges(forestFile)
    with open(os.path.join(inPath,sampleName + ".stp")) as inStp:
        with open(os.path.join(outPath,sampleName + ".stp"), "w") as outStp:
            for line in inStp:
                parts = line.split()
                if parts[0] == "D":
                    # If the edge is directed it must have been used in
                    # the same direction in the forest.  Both the stp
                    # file and the forest output follow the <biologicalTarget> <biologicalParent>
                    # convention, i.e. the edge points toward the root in the
                    # belief propagation not the biological parent.
                    if "%s %s" % (parts[1], parts[2]) in edges:
                        outStp.write(line)
                elif parts[0] == "E":
                    # If the edge is undirected allow it to be used in
                    # either direction in the forest
                    if ("%s %s" % (parts[1], parts[2]) in edges) or ("%s %s" % (parts[2], parts[1]) in edges):
                        outStp.write(line)
                else:
                    outStp.write(line)
            # Write the degree penalties
            for node in degreePenalties.iterkeys():
                if node in potentialSteiner:
                    outStp.write("W %s %f\n" % (node, degreePenalties[node]))


# Update the initial stp file by adding artificial prizes present in the
# dictionary.  Add the negative prizes used to implement penalties on
# node degree as well.
def UpdateStp(artificialPrizes, degreePenalties, potentialSteiner, inPath, outPath, sampleName):
    newStpFile = os.path.join(outPath,sampleName + ".stp")
    shutil.copyfile(os.path.join(inPath,sampleName + ".stp"), newStpFile)
    with open(newStpFile, "a") as outStp:
        # Create a set of nodes with an artificial prize or a degree penalty
        keys = set(artificialPrizes.iterkeys()).union(degreePenalties.iterkeys())

        # Append the new prizes to the end of the file summing the two types of
        # artificial prizes if both are present
        for node in keys:
            if node in potentialSteiner:
                artificialPrize = 0
                degreePenalty = 0
                if node in artificialPrizes:
                    artificialPrize = artificialPrizes[node]
                if node in degreePenalties:
                    degreePenalty = degreePenalties[node]
                outStp.write("W %s %f\n" % (node, artificialPrize + degreePenalty))


# Created weighted artificial prizes that are derived from the frequency
# of a node in the collection of Steiner forests.
def CreateWgtPrizes(allProts, lastForests, lambda1, alpha, negativePrizes):
    forestFreq = NetworkUtil.SetFrequency(lastForests)
    artificialPrizes = {}
    if negativePrizes:
        # Need to iterate over all proteins when creating negative prizes
        for node in allProts:
            freq = 0
            if node in forestFreq:
                freq = forestFreq[node]
            # Only create non-zero prizes, i.e. for nodes that are not in all
            # forests
            if freq < 1:
                artificialPrizes[node] = -lambda1 * ((1-freq) ** alpha)
    else:
        # For positive prizes only need to iterate over the nodes that appear
        # in some forest
        for node in forestFreq.iterkeys():
            freq = forestFreq[node]
            # Frequently is guaranteed to be > 0 because the keys are only
            # the union of all forest nodes
            artificialPrizes[node] = lambda1 * (freq ** alpha)
                
    return artificialPrizes

# Create unweighted artificial prizes by explicitly constructing a common
# set and then generating prizes for all proteins accordingly depending
# on whether positive or negative prizes are being used.
def CreateUnwgtPrizes(allProts, potentialSteiner, lastForests, lambda1, lambda2, negativePrizes):
    # Construct the common set
    commonSet = set()    
    for prot in allProts:
        # The same scoring function works for both positive and negative prizes
        inScore = Penalty(prot, potentialSteiner, lastForests)
        if lambda1 * inScore < lambda2:
            commonSet.add(prot)
    
    # Build the prizes based on the common set.  With positive artificial prizes
    # add prizes when a node is in the common set.  With negative artificial
    # prizes add a prize when a node is not in the common set.
    if negativePrizes:
        artificialPrizes = dict.fromkeys(allProts.difference(commonSet), -lambda1)
    else:
        artificialPrizes = dict.fromkeys(commonSet, lambda1)
    return artificialPrizes


# This penalty function charges a cost of 1 for each sample in which the
# protein is a potential Steiner node (i.e. not a prize node) but
# not in the Steiner forest.  The same function is used to calculate
# the reward when using negative prizes.
def Penalty(prot, potentialSteiner, lastForests):
    if len(potentialSteiner) != len(lastForests):
        raise RuntimeError("Must have the same number of samples")

    count = 0
    for i in range(len(lastForests)):
        if (prot not in lastForests[i]) and (prot in potentialSteiner[i]):
            count += 1

    return count


# Load all edges of the Steiner forest as set whose members of of the form
# "<protein1> <protien2>".  Expects the "symbol_*" output from the
# message passing algorithm
def LoadForestEdges(forestFile):
    forest = set()
    with open(forestFile) as f:
        for line in f:
            forest.add(line.strip())
    return forest


# Load all members of the Steiner forest, both Steiner nodes and prizes,
# and return them as a set.  Expects the "symbol_*" output from the
# message passing algorithm
def LoadForestNodes(forestFile):
    forest = set()
    with open(forestFile) as f:
        for line in f:
            forest.update(line.split())
    return forest


# A helper function that takes all of the LearnSteiner functions as a tuple
def LearnSteinerHelper(args):
    return LearnSteiner(*args)


# Learn a Steiner forest
# TODO call the Python function directly instead of using a system call
def LearnSteiner(opts, stpPath, outPath, sampleName, targetFile, threads):
    print "Learning Steiner forest for %s" % sampleName
    pyFile = os.path.join(os.path.dirname(__file__),"PCSF.py")
    os.system("python %s --outputpath=%s --msgpath=%s --depth=%d --conn=3 --targetfile=%s --stppath=%s/ --stpfile=%s --W=%s_%s_%s --beta=1.0_1.0_1.0 --exclude=2 --species=human --threads=%s" % (pyFile, outPath, opts.msgPath, opts.depth, targetFile, stpPath, sampleName, str(opts.W), str(opts.W), str(opts.W), str(threads)))


# Find and return a set of the potential Steiner nodes, the non-prize
# proteins (not mRNAs) that are present in the network.  Also
# write a list of dummy node neighbors, which can be either prizes or
# non-prizes.
def DummyNeighbors(allProts, path, stpFile, dnFile, neighborType):
    prizes = set()
    with open(os.path.join(path,stpFile)) as f:
        for line in f:
            parts = line.split()
            # mRNAs can't be Steiner nodes
            if parts[0] == "W" and not parts[1].endswith("_MRNA"):
                prizes.add(parts[1])

    psNodes = allProts.difference(prizes)    

    if neighborType == "prizes":
        NetworkUtil.WriteCollection(os.path.join(path,dnFile), prizes)
    elif neighborType == "nonprizes":
        NetworkUtil.WriteCollection(os.path.join(path,dnFile), psNodes)
    else:
        raise RuntimeError("%s is not a valid type of dummy node neighbor connection" % neighborType)

    return psNodes


# A helper function that takes all of the CreateStp args as a tuple
def CreateStpHelper(args):
    return CreateStp(*args)


# Create the initial .stp file
# TODO call the Python function directly instead of using a system call
def CreateStp(opts, initPath, terminalFile, sampleName):
    stpFile = sampleName + ".stp"
    print "Creating %s" % stpFile
    pyFile = os.path.join(os.path.dirname(__file__),"MessagePassingInput.py")
    os.system("python %s --interactomepath=%s --terminalpath=%s --resultpath=%s --indirectedfilename=%s --beta=%f --terminalfile=%s --resultfilename=%s --directedfilename=%s --tfdnafilename=%s --mrnabeta=%s" % (pyFile, opts.interactomePath, opts.terminalPath, initPath, opts.undirectedFile, opts.beta, terminalFile, stpFile, opts.directedFile, opts.tfdnaFile, opts.mrnaBeta))
    return stpFile


# Create a single list of all items in a dictionary of iterables.  Takes a
# list of keys so that the order of items in the returned list
# is deterministic.
def FlattenDict(dictionary, keys):
    if not len(keys) == len(dictionary.keys()):
        raise RuntimeError("Dictionary does not have %d keys" % len(keys))

    flat = []
    for key in keys:
        flat.extend(dictionary[key])

    return flat


# Load all terminal files from the master file and create a map from
# groups to member samples.  Default to placing all samples in a single
# group.  Returns a map from groups to terminal files, a map from
# groups to sample names, and a map from groups to sample counts
def LoadTerminalFiles(terminalPath, masterTerminalFile):
    terminalMap = dict()
    sampleMap = dict()
    defaultCount = 0
    with open(os.path.join(terminalPath,masterTerminalFile)) as f:
        for line in f:
            curSample = line.strip().split("\t")
            if len(curSample) > 2:
                raise RuntimeError("%s has a line with %d columns" % (masterTerminalFile, len(curSample)))
            elif len(curSample) == 1:
                group = "defaultGroup"
                defaultCount += 1
            else:
                group = curSample[1]

            # Create lists for the group if it does not yet exist
            if group not in terminalMap:
                terminalMap[group] = []
                sampleMap[group] = []
            terminalMap[group].append(curSample[0])

            name = curSample[0]
            if name.endswith(".txt"):
                name = name[0:-4] # Remove ".txt"
            sampleMap[group].append(name)

    # Count the number of samples in each group
    countMap = dict()
    for group in terminalMap.iterkeys():
        countMap[group] = len(terminalMap[group])

    print "Assigned %d samples to the default group" % defaultCount
    return terminalMap, sampleMap, countMap


# Load all proteins in the interaction network but not the genes
# (the targets of protein-DNA interactions)
def LoadProteins(interactomePath, undirectedFile, directedFile, tfdnaFile):
    prots = set()
    # Load a proteins in the PPI
    with open(os.path.join(interactomePath,undirectedFile)) as f:
        for line in f:
            parts = line.strip().upper().split()
            if len(parts) == 3:
                prots.add(parts[0])
                prots.add(parts[1])

    # Load the kinases and substrates if these interactions exist
    if directedFile != "None":
        with open(os.path.join(interactomePath,directedFile)) as f:
            for line in f:
                parts = line.strip().upper().split()
                if len(parts) == 3:
                    prots.add(parts[0])
                    prots.add(parts[1])

    # Load the TFs
    if tfdnaFile != "None":
        with open(os.path.join(interactomePath,tfdnaFile)) as f:
            for line in f:
                parts = line.strip().upper().split()
                if len(parts) == 3:
                    # Load only the TF, not the gene
                    prots.add(parts[0])

    return prots


# Setup the option parser
def CreateParser():
    # Options that are fixed for the time being:
    # Connection type for the artificial node / exclude / targetfile
    # Output file names
    # Beta (only allowed to be set once initially, range of Beta is not supported)
    # Species
    # W cannot be given as a range of values
    parser = OptionParser()
    parser.add_option("--interactomepath", type="string", dest="interactomePath", help="This path points to the directory that contains the interaction network files",default='None')
    parser.add_option("--terminalpath",type="string",dest="terminalPath",help="This path points to the directory that contains the terminal (node prize) files",default='None')
    parser.add_option("--resultpath",type="string",dest="resultPath",help="This path points to the directory where the output files will be written.",default='None')
    parser.add_option("--undirectedfile",type="string",dest="undirectedFile",help="The name of the protein-protein interaction file in the interactomepath directory.  The file is expected to contain undirected interactions with probabilistic weights (e.g in [0,1]). Columns should be ordered [prot1 prot2 weight].",default='None')
    parser.add_option("--terminalfile",type="string",dest="masterTerminalFile",help="A file in the terminalpath directory that lists the files that give the node prizes for each sample.  All listed filenames should be relative to terminal path.  If gene penalties are given in the terminal files, gene names should end with '_MRNA'.  Optionally can include a tab-separated second column that assigns each sample to a group so the forests are only constrained to be similar to other samples in the same group.",default='None')
    #parser.add_option("--resultfilename",type="string",dest="resultfilename",help="The name of the file where the combined information will be written to be used as input in the message passing tool.",default='None')
    parser.add_option("--directedfile",type="string",dest="directedFile",help="Optional: The name of the interaction file where directed interactions (e.g. kinase-substrate) with probabilistic weights are available.  Columns should be ordered [substrate kinase weight].",default='None')
    parser.add_option("--tfdnafile",type="string",dest="tfdnaFile",help="Optional: The name of the interaction file where TF-DNA interactions with probabilistic weights are available. Columns should be ordered [TF gene weight].  Gene names should not end with '_MRNA' because '_MRNA' is automatically appended to them.",default='None')
    parser.add_option("--mrnabeta",type="string",dest="mrnaBeta",help="Optional: The beta parameter for gene terminal nodes. Its default value is equal to the --beta value.  It is applied to all terminals whose name ends with '_MRNA'.  The scaling is only performed once when the initial stp files are created.",default='None')


    #parser.add_option("--outputpath", type="string", dest="outputpath", help="This path points to the directory where the output files will be written",default='None')
    parser.add_option("--msgpath",type="string",dest="msgPath",help="The path and file name of the msgsteiner executable",default='./msgsteiner')
    parser.add_option("--depth",type="int",dest="depth",help="Depth parameter that limits the maximum depth from the Steiner tree root to the leaves",default=10)
    #parser.add_option("--conn",type="string", dest="conn", help="How to connect the artificial node to the interactome: 1 to all nodes in the interactome, 2 to all node in the interactome, except the terminals, 3 to a given set of nodes in the interactome, 4 to a given set of nodes in the interactome, except the terminals.")
    #parser.add_option("--stppath", type="string", dest="stppath", help="This path points to the directory where the stp file is available.")
    #parser.add_option("--stpfile", type="string", dest="stpfile", help="The name of the stp file (without the file extension)")
    parser.add_option("--W", type="float", dest="W",help="The cost of the edges from the artificial root node to its neighbors.",default=1.0)
    #parser.add_option("--beta", type="string", dest="beta",help="the range and increment for the beta parameter. format: start_end_increment")
    #parser.add_option("--exclude", type="string", dest="exclude", help="If set to 1 and the connection type is 2, 3, or 4 then remove prizes from all nodes that are connected to the artificial node")
    #parser.add_option("--targetfile", type="string", dest="targetfile", help="a subset of the proteins in the interactome that will be connected to the artificial root in connection types 3 and 4",default='None')
    #parser.add_option("--species", type="string", dest="species", help="the organism that you are working on",default='None')

    parser.add_option("--beta",type="float",dest="beta",help="The scaling factor applied to the node prizes, which is used to control the relative strength of node prizes and edge costs.  This scaling is only performed once when the initial stp files are created.",default=1.0)
    parser.add_option("--lambda", type="float", dest="lambda1",help="The tradeoff coefficient for the penalty incurred by nodes in the Steiner forests that are not in the set of common nodes.",default=1.0)
    # Redefine lambda2 in terms of how many trees should contain a node before you include it in the common set?
    parser.add_option("--alpha", type="float", dest="lambda2",help="The tradeoff coefficient for the reward on the size of the set of common nodes when using unweighted artificial prizes or the power to which the node frequency is taken for weighted prizes.",default=1.0)
    parser.add_option("--mu", type="float", dest="mu", default=0, help="A parameter used to penalize high-degree nodes from being selected as Steiner nodes.  Does not affect prize nodes but does affect artificial prizes.  The penalty is -mu*degree.  Set mu <= 0 to disable the penalty (default).")
    # Fix the number of iterations until other convergence criteria is implemented
    parser.add_option("--iterations", type="int", dest="iterations", help="The number of iterations to run", default=10)
    # Non-positive values will map to cpu_count()
    parser.add_option("--workers", type="int", dest="workers", help="The number of worker processes to use in the multiprocessing pool or threads to use in multi-threaded belief propagation.  Should not exceed the number of cores available.  Defaults to the number of CPUs.", default=-1)
    # Use positive or negative artificial prizes to encourage use of the common set or weighted prizes to avoid
    # constructing a common set
    parser.add_option("--artificialprizes", type="choice", dest="artificialPrizes", default="negativeWeighted", choices=["positive", "negative", "positiveWeighted", "negativeWeighted"], help="Use 'positive' or 'negative' prizes to encourage trees to include common set proteins.  Use 'positiveWeighted' or 'negativeWeighted' (default) prizes to construct weighted artificial prizes based on the node frequency in the most recent forests.")
    # Determine how to connect the dummy node
    parser.add_option("--dummyneighbors", type="choice", dest="dummyNeighbors", default="prizes", choices=["prizes", "nonprizes"], help="Connect the dummy node to all 'prizes' (default) or 'nonprizes'.")
    parser.add_option("--itermode", type="choice", dest="iterMode", default="batch", choices=["batch", "random"], help="Learn forests simultaneously in 'batch' (default) or sequentially in 'random' order.  Batch mode computes artificial prizes with respect to all forests at the previous iteration.  Random mode computes prizes for a specific sample given the most recent forests for all other samples.")
    return parser


if __name__=="__main__":
    # Use the command line arguments to setup the options (the same as the default OptionParser behavior)
    main(sys.argv[1:])
