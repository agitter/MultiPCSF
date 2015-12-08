import networkx

__author__ = 'Anthony Gitter'

# Load an interaction network as an undirected graph.
# Each line should have the format "node1 node2 weight" where the presence of
# the weight column is specified by the weight argument.
# There should be no header line in the file.  Can also be used to load
# the "symbol_*" output from the Steiner tree message passing algorithm,
# but all edges will be stored as undirected edges and the dummy node
# will be included.
# TODO Is this function needed anymore?
def LoadNetwork(networkFile, weight=False):
    if weight:
        graph = networkx.read_edgelist(networkFile, delimiter=None, data=(("weight",float),))
    else:
        graph = networkx.read_edgelist(networkFile, delimiter=None)
    #print "Loaded network %s with %d nodes and %d edges" % (networkFile, graph.order(), graph.size())
    return graph


# Load an interaction network as an undirected graph from a graphite
# edge list.  Only the first two columns, the node names, are used.  Edge
# attributes are ignored.  Multiple instances of the same edge are collapsed
# and self edges are ignored.
# TODO This is partly redundant with LoadNetwork
def LoadGraphiteNetwork(networkFile):
    graph = networkx.Graph()
    with open(networkFile) as f:
        # Read all edges in this pathway and build a graph
        # If an edge is listed twice it will only be stored once because this is not a MultiGraph
        for edgeLine in f:
            edgeParts = edgeLine.split()
            # Ignore self edges
            if not edgeParts[0] == edgeParts[1]:
                graph.add_edge(edgeParts[0], edgeParts[1])
    return graph

# Load an interaciton network as an undirected graph and calculate
# a penalty (negative prize) for each node.  The penalty is
# -mu * degree.  If mu <= 0 then return an empty dictionary.
def DegreePenalties(mu, undirNetworkFile, dirNetworkFile):
    penalties = dict()

    # mu <= 0 means no penalties are desired
    if mu > 0:
        # TODO need to support directed edges as well
        if dirNetworkFile != None and dirNetworkFile != "None":
            raise RuntimeError("Degree penalties for directed networks are not yet supported")
        # Interaction networks are assumed to be weighted
        network = LoadNetwork(undirNetworkFile, True)
        for node in network:
            penalty = -mu * network.degree(node)
            penalties[node] = penalty

    return penalties


# Return the node and edge intersection of two undirected graphs
# Unlike the networkx intersection, the graphs may have different node sets.
def Intersection(graph1, graph2):
    copy = graph1.copy()
    # Iterate through all graph1 edges and remove those that are not in graph2
    for edge in graph1.edges():
        # *edge unpacks the edge tuple
        if not graph2.has_edge(*edge):
            copy.remove_edge(*edge)
    # Iterate through all graph1 nodes and remove those that are not in graph2
    for node in graph1.nodes():
        if not node in graph2:
            copy.remove_node(node)
    return copy


# Take a list of sets and return a dict that maps elements in the sets
# to the number of sets they appear in as ints.  This is not written for
# maximum efficiency but rather readability.
def SetCounts(setList):
    keys = set()
    for curSet in setList:
        keys.update(curSet)

    # Initialize the dictionary that will be used to store the counts of each element
    countDict = dict.fromkeys(keys, 0)
    # Populate the counts
    for curSet in setList:
        for element in curSet:
            countDict[element] += 1

    return countDict


# Take a list of sets and return a dict that maps elements in the sets
# to the fraction of sets they appear in.  This is not written for
# maximum efficiency but rather readability.
def SetFrequency(setList):
    n = float(len(setList)) # Want floating point division
    countDict = SetCounts(setList)

    # Transform the counts into frequencies
    freqDict = {}
    for key in countDict.keys():
        freqDict[key] = countDict[key] / n

    return freqDict


# Write a collection to a file with one item per line
def WriteCollection(filename, collection):
    with open(filename, "w") as f:
            for item in collection:
                f.write(str(item) + "\n")


# Write a dictionary to a tab-delimited file with one key-value pair
# per line using str() to format the keys and values
def WriteDict(filename, dictionary):
    with open(filename, "w") as f:
            for item in dictionary.iterkeys():
                f.write("%s\t%s\n" % (str(item), str(dictionary[item])))
