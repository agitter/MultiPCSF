__author__ = 'Nurcan Tuncbag' # Modified by Anthony Gitter
__email__ = 'ntuncbag@mit.edu'

import os, math
import networkx
import pickle, operator
import string
import subprocess
from optparse import OptionParser

def sort_table(table, cols):
    """ sort a table by multiple columns
        table: a list of lists (or tuple of tuples) where each inner list 
               represents a row
        cols:  a list (or tuple) specifying the column numbers to sort by
               e.g. (1,0) would sort by column 1, then by column 0
    """
    for col in reversed(cols):
        table = sorted(table, key=operator.itemgetter(col))
    return table

def idmatchHuman(species):
#    species = species.lower()
#    if species == 'human': 
#         idfile=pickle.load(open('/nfs/data/ntuncbag/GeneSymbolMapping/humGeneToEnsPep.pkl','r'))
#    if species == 'mouse':
#         idfile=pickle.load(open('/nfs/data/ntuncbag/GeneSymbolMapping/mm9SymbolToEnsPep_pep2gene.pkl','r'))
#    return idfile
    return dict()

# Return a list of all unique receptors in the file
def givenset(stppath,targetfile):
    artificialTargets = set()
    file = open(os.path.join(stppath, targetfile),"r")
    #file = open(stppath+''+targetfile,"r")
    while 1:
        line = file.readline()
        if line == "": break
        receptor = line.strip().split()[0]#.replace("9606.","")
        artificialTargets.add(receptor)
    file.close()
    artificialTargets = list(artificialTargets)
    return artificialTargets


# Return a list of all unique receptors in the file, minus terminal nodes
# in the STP file
def givenset_terminalexcluded(stppath, targetfile, stpfile):
    artificialTargets = set()
    file = open(os.path.join(stppath, targetfile),"r")
    #file = open(stppath+''+targetfile,"r")
    while 1:
        line = file.readline()
        if line == "": break
        receptor = line.strip().split()[0]#.replace("9606.","")
        artificialTargets.add(receptor)
    file.close()
    artificialTargets = list(artificialTargets)
    file = open(os.path.join(stppath,"%s.stp" % stpfile), "r")
    #file = open(stppath+'/'+stpfile+'.stp',"r")
    nodelist = set()
    terminalnodelist = set()
    while 1:
        line = file.readline()
        if line == "":
            break
#        if line.startswith("E") or line.startswith("D"):
#            temp = line.split()
#            node1 = temp[1]
#            node2 = temp[2]
#            nodelist.add(node1)
#            nodelist.add(node2)
        if line.startswith("W"):
            temp = line.strip().split()
            node = temp[1]
            terminalnodelist.add(node)
    
    for item in terminalnodelist:
        try:
            artificialTargets.remove(item)
        except ValueError:
            continue
    return artificialTargets

# Rank receptors (named in targetfile) by their degree in the network
# in stpfile
def ReceptomeRanking(stppath,stpfile,targetfile,species):
    G = networkx.Graph()
    receptorDegree = []
    idfile = idmatchHuman(species)
    file = open(os.path.join(stppath,"%s.stp" % stpfile), "r")
    #file = open(stppath+'/'+stpfile+'.stp',"r")
    nodelist = givenset(stppath, targetfile) #'diseases.txt') #'artificialTargetFile')
    artificial = "DUMMY"
    edgelist = []
    #cost = []
    terminals = []
    terminalset = set()
    while 1:
        line = file.readline()
        if line == "": break
        if line.startswith("E") or line.startswith("D"):
            temp = line.split()
            node1 = temp[1]
            node2 = temp[2]
#            nodelist.add(node1)
#            nodelist.add(node2)
            #cost.append(float(temp[3]))
            G.add_edge(node1,node2)

    for node in nodelist:
        deg = networkx.degree(G,node)
        try:
            nodename = idfile[node]
        except KeyError:
            nodename = [node]
        if deg != {} and deg != 0 and deg != []:
            receptorDegree.append([nodename, deg])
    sortedlist = sort_table(receptorDegree, [1,1])
    index1 = -1
    #print sortedlist[-10:-1]
    dreceptomeRank = {}
    for i in range(len(sortedlist)):
        index = index1-i
        proteins = sortedlist[index][0]
        for protein in proteins:
#            print protein, i
            dreceptomeRank[protein] = i
    print len(dreceptomeRank)
    return dreceptomeRank

# Return a set of all nodes in the network
def allnodes(stppath, stpfile):
    file = open(os.path.join(stppath,"%s.stp" % stpfile), "r")
    #file = open(stppath+'/'+stpfile+'.stp',"r")
    nodelist = set()
    while 1:
        line = file.readline()
        if line == "":
            break
        if line.startswith("E") or line.startswith("D"):
            temp = line.split()
            node1 = temp[1]
            node2 = temp[2]
#            if node1.endswith('mrna') == False:
            nodelist.add(node1)
#            if node2.endswith('mrna') == False:
            nodelist.add(node2)
    print len(nodelist)
    return nodelist

# Return a list of all nodes in the network that are not terminals
def terminalexcluded(stppath, stpfile):
    file = open(os.path.join(stppath,"%s.stp" % stpfile), "r")
    #file = open(stppath+'/'+stpfile+'.stp',"r")
    nodelist = set()
    terminalnodelist = set()
    while 1:
        line = file.readline()
        if line == "":
            break
        if line.startswith("E") or line.startswith("D"):
            temp = line.split()
            node1 = temp[1]
            node2 = temp[2]
            nodelist.add(node1)
            nodelist.add(node2)
        if line.startswith("W"):
            temp = line.strip().split()
            node = temp[1]
            terminalnodelist.add(node)
    nodelist = list(nodelist)
    for item in terminalnodelist:
        try:
            nodelist.remove(item)
        except ValueError:
            continue
    return nodelist

# Return all edges for which neither endpoint is a member of knockoutlist
def knockout(stppath,stpfile,knockoutlist):
    file = open(os.path.join(stppath,"%s.stp" % stpfile), "r")
    #file = open(stppath+'/'+stpfile+'.stp',"r")
    edgelist = []
    while 1:
        line = file.readline()
        if line == "": break
        if line.startswith("E") or line.startswith("D"):
            temp = line.split()
            node1 = temp[1]
            node2 = temp[2]
#            nodelist.add(node1)
#            nodelist.add(node2)
            #cost.append(float(temp[3]))
            if (node1 not in knockoutlist) and (node2 not in knockoutlist):
                edgelist.append([temp[0], node1,node2,float(temp[3])])
    return edgelist


# Modified to return a lists of all the lines in the input to the message passing
# code instead of writing them to disk
def PrepareInputFile(stppath,stpfile,connectiontype, outputpath, w, b, D, knockoutlist, exclude, targetfile,species):
    # nodelist is the list of nodes that will be connected to the dummy node
    if connectiontype == '1':
        nodelist = allnodes(stppath, stpfile)
    if connectiontype == '2':
        nodelist = terminalexcluded(stppath, stpfile)
    if connectiontype == '3':
        nodelist = givenset(stppath, targetfile)#'diseases.txt')#'artificialTargetFile')
    if connectiontype == '4':
        nodelist = givenset_terminalexcluded(stppath, targetfile, stpfile)

    ind = 0
    file = open(os.path.join(stppath,"%s.stp" % stpfile), "r")
    #file = open(stppath+'/'+stpfile+'.stp',"r")
    artificial = "DUMMY"
    edgelist = []
    #cost = []
    terminals = []
    terminalset = set()
    while 1:
        line = file.readline()
        if line == "": break
        if line.startswith("E") or line.startswith("D"):
            temp = line.split()
            node1 = temp[1]
            node2 = temp[2]
            ind += 1
#            nodelist.add(node1)
#            nodelist.add(node2)
            #cost.append(float(temp[3]))
            #print ind
#            if connectiontype == '3':
#                if (node1 in nodelist or node2 in nodelist) and float(temp[3])<=0.225:
#                    edgelist.append([temp[0], node1,node2, 0.01])
#                elif (node1 not in nodelist and node2 not in nodelist):                
#                    edgelist.append([temp[0], node1,node2,float(temp[3])])
#            else:
            edgelist.append([temp[0], node1,node2,float(temp[3])])

        if line.startswith("W"):
            temp = line.strip().split()
            node = temp[1]
            prize = float(temp[2])#1.0 #float(temp[2])/4.0 #1.0 #float(temp[2])/4.0
            if connectiontype == '1':
                if node.endswith('_MRNA'):
                #if node.endswith('mrna'):
                    terminals.append("W %s %.4f\n" % (node, b*prize))
                    terminalset.add(node)
                else:
                    terminals.append("W %s %.4f\n" % (node, b*prize))
                    terminalset.add(node)
            else:
                if exclude == True:
                    if node not in nodelist:
                        if node.endswith('_MRNA'):
                        #if node.endswith('mrna'):
                            terminals.append("W %s %.4f\n" % (node, b*prize))
                            terminalset.add(node)
                        else:
                            terminals.append("W %s %.4f\n" % (node, b*prize))
                            terminalset.add(node)
                if exclude == False:
                    if node.endswith('_MRNA'):
                    #if node.endswith('mrna'):
                        terminals.append("W %s %.4f\n" % (node, b*prize))
                        terminalset.add(node)
                    else:
                        terminals.append("W %s %f\n" % (node, b*prize))
                        terminalset.add(node)


    file.close()
    # If there are nodes to be knocked-out, update the edge list to remove
    # edges involved a KO'd node
    if knockoutlist != []:
        edgelist = knockout(stppath,stpfile,knockoutlist)
    # No longer write the input to disk, instead build the fild in memory
    #file = open("%s/%s_%s_%s_%s.dat" % (outputpath, stpfile, str(w), str(b), str(D)),"w")
    inputData = []
    for item in edgelist:
        #file.writelines("%s %s %s %f\n" %(item[0], item[1], item[2], item[3]))
        inputData.append("%s %s %s %f\n" %(item[0], item[1], item[2], item[3]))
    for item in nodelist:
        #file.writelines("D %s %s %.4f\n" %(item, artificial, w))
        inputData.append("D %s %s %.4f\n" %(item, artificial, w))
    #for item in terminals:
        #file.writelines(item)#"W %s %f\n" % (item, b*1.0/len(terminalnodelist)))
    inputData.extend(terminals)
    #file.writelines("W %s 100.0\n" % artificial)
    inputData.append("W %s 100.0\n" % artificial)
    #file.writelines("R %s\n\n" % artificial)
    inputData.append("R %s\n\n" % artificial)
    #file.close()
    idfile = idmatchHuman(species)
    file = open("%s/gbm_cell_line.attr" % outputpath,"w")
    for t in terminalset:
        try: 
            sym_tlist = idfile[t]
            for sym_t  in sym_tlist:
                file.writelines(sym_t+' = 1\n')
        except KeyError:
            file.writelines(t+' = 1\n')
    file.close()
    print "Input Data Prepared: w = %s, b = %s, connectiontype = %s" % (w, b, connectiontype)
    return inputData

# Start the message passing algorithm
# Takes the input data as a list of lines instead of reading them from disk
def RunMSGAlgorithm(msgpath, D, w, b, outputpath, stpfile, connectiontype, inputData, threads):
    # No longer have an input file.  Open the output files.
    #inputfile = "%s/%s_%s_%s_%s.dat" % (outputpath, stpfile, str(w), str(b),str(D))
    resultfilename = "%s/%s_%s_%s_%s.txt" % (outputpath, stpfile, str(w), str(b),str(D))
    objectivefilename = "%s/%s_%s_%s_%s.objective" % (outputpath, stpfile, str(w), str(b),str(D))
    with open(resultfilename, "w") as resultFile:
        with open(objectivefilename, "w") as objFile:
            #os.system( "%s/msgsteiner -d %s -t 1000000 -o -r 1e-5 -g 1e-3 < %s > %s 2> %s" % (msgpath, D, inputfile, resultfile, objectivefile))
            # Start a subprocess with a 1 line buffer size
            #subprocArgs = ["%s/msgsteiner" % msgpath, "-d", D, "-t", "1000000", "-o", "-r", "1e-5", "-g", "1e-3"]
            subprocArgs = ["%s/msgsteiner9" % msgpath, "-d", D, "-t", "1000000", "-o", "-r", "1e-5", "-g", "1e-3", "-j", str(threads)]
            subproc = subprocess.Popen(subprocArgs, bufsize=1, stdin=subprocess.PIPE, stdout=resultFile, stderr=objFile)
            for line in inputData:
                subproc.stdin.write(line) # Lines are already newline terminated
            subproc.stdin.close()
            subproc.wait()
    print "MSG Run Finished with the parameters: w = %s, b = %s, D = %s, connectiontype = %s" % (w, b,D, connectiontype)
    return 1

# Write statistics about the forest including the sizes of the individual trees
# Returns the connected components
def OutputPCSFCheck(w, b, D, outputpath, stpfile,species):
    resultfile = "%s/%s_%s_%s_%s.txt" % (outputpath, stpfile, str(w), str(b),str(D))
    outputfile = "%s/%s_%s_%s_%s.output" % (outputpath, stpfile, str(w), str(b),str(D))
    file = open(resultfile,"r")
    outfile = open(outputfile,"w")
    # G is the graph in which the dummy node has been removed
    G = networkx.Graph()
    # H is the graph returned from the message passing algorithm
    H = networkx.Graph()
    dummyPartners = []
    while 1:
        line = file.readline()
        if line == "": break
        temp = line.strip().split()
        node1 = temp[0]
        node2 = temp[1]
        H.add_edge(node1,node2)
        if node1 != "DUMMY" and node2 != "DUMMY":
            G.add_edge(node1,node2)
            #edgeweights.append(float(temp[3]))
        if node1 == 'DUMMY':
            dummyPartners.append(node2)
        if node2 == 'DUMMY':
            dummyPartners.append(node1)
#    print len(dummyPartners)
    file.close()
    # A list, where each element in the list is itself a list of nodes in the component
    nodeList = networkx.connected_components(G) #component.connected_components(G)
    treesize = []
    nonredundanttreesize = []
    tempnodelist = []
    nodeClusterDict = {}
    for nodes in nodeList:
        for node in nodes:
            nodeClusterDict[node] = len(nodes)
        if len(nodes) >= 10:
             #for n in nodes:
            nonredundanttreesize.append(len(nodes))
            for n in nodes:
#                    if n not in terminalnodelist:
#                        print n
                tempnodelist.append(n)
#                print 'batch'
        treesize.append(len(nodes))
    outfile.writelines("Parameters that are used: w = %s, b = %s, D = %s\n\n" % (str(w), str(b), str(D)))
    # Special case if the Steiner forest that was output was empty
    if len(treesize) == 0:
        outfile.writelines("PCSF characteristics:\n-----\nmin_tree size = %d\nmax_tree size = %d\nmean_tree size = %f\nnumber of trees = %d\ntotal size = %d\nnumber of singletons = %d\n" % (0, 0, 0, 0, 0, 0))
        print "Steiner forest %s is empty" % resultfile
    else:
        outfile.writelines("PCSF characteristics:\n-----\nmin_tree size = %d\nmax_tree size = %d\nmean_tree size = %f\nnumber of trees = %d\ntotal size = %d\nnumber of singletons = %d\n" % (min(treesize), max(treesize), sum(treesize)/float(len(treesize)), len(treesize), len(H.nodes())-1, len(H.nodes())-sum(treesize)-1))
    outfile.writelines('----------------\n')
    outfile.writelines('Size of the trees in the forest:\n')
    ind = 0
    for ts in treesize:
        ind += 1
        outfile.writelines('T'+str(ind)+'\t'+str(ts)+'\n')
    outfile.writelines('-------------\n')
    #idfile = idmatchHuman(species)
##    receptomeDegDict = ReceptomeRanking(stppath,stpfile)
    for item in dummyPartners:
        try:
            subTsize = nodeClusterDict[item]
            if subTsize >= 2:
##                #print item, subTsize
                #try:
                #    for s in idfile[item]:
                #outfile.writelines( s + "\t" + str(subTsize) + '\n') #'ranking = ', receptomeDegDict[s], 'TreeSize',subTsize
                #except KeyError:
                outfile.writelines(item + "\t" + str(subTsize) + '\n') #'ranking = ', receptomeDegDict[s], 'TreeSize', subTsize
        except KeyError:
            continue
    outfile.close()
    return nodeList

# Return the set of nodes and an unweighted graph object from the message passing results
def SteinerTree(resultpath, resultfilename):
    G = networkx.Graph()
    nodeSet = set()
    file = open(resultpath+'/'+resultfilename,'r')
    while 1:
        line = file.readline()
        if line == "": break
        temp = line.strip().split()
        nodeSet.add(temp[0])
        nodeSet.add(temp[1])
        #if temp[0].endswith('mrna') == False and temp[2].endswith('mrna') == False:
        G.add_edge(temp[0], temp[1])#, weight=float(temp[2]))
    file.close()
    return nodeSet, G

# Returns all, undirected, and directed edges between members of the connected components
# of the forest
def NetworkMapping(outputpath, stppath, stpfile, w, b, D,species):
    idfile = idmatchHuman(species)
    # The directed and undirected edges in the orginal stp file, with weights
    G1 = networkx.Graph()
    G2 = networkx.DiGraph()
    # All, undirected, and directed edges between members of the connected components
    # of the forest
    H = networkx.Graph()
    H1 = networkx.Graph()
    H2 = networkx.DiGraph()
    # An undirected, weighted network containing all original edges in the stp file
    # between nodes used in the forest
    I = networkx.Graph()
    nodeList = OutputPCSFCheck(w, b, D, outputpath, stpfile, species)
    file = open(stppath+"/"+stpfile+'.stp', "r")
    while 1:
        line = file.readline()
        if line == "": break
        if line.startswith("E"): # or line.startswith("D"):
            temp = line.strip().split()
            G1.add_edge(temp[1], temp[2], weight=float(temp[3]))
        if line.startswith("D"): # or line.startswith("D"):
            temp = line.strip().split()
            G2.add_edge(temp[1], temp[2], weight=float(temp[3]))
    resultfilename = "%s_%s_%s_%s.txt" % (stpfile, str(w), str(b), str(D))
    # The undirected, unweighted version of the graph returned by the message passing algorithm
    S = SteinerTree(outputpath, resultfilename)[1]
    # Iterate over connected components
    for nodeSet in nodeList:
        for node1 in nodeSet:
            for node2 in nodeSet:
                if G1.has_edge(node1,node2):
                #if G1.has_edge(node1,node2) == True:
                    weight1 = G1.get_edge_data(node1, node2)
                    H1.add_edge(node1, node2, weight=weight1['weight'])
                    H.add_edge(node1, node2, weight=weight1['weight'])
                if G2.has_edge(node1,node2):
                #if G2.has_edge(node1,node2) == True:
                    weight1 = G2.get_edge_data(node1, node2)
                    H2.add_edge(node1, node2, weight=weight1['weight'])
                    H.add_edge(node1, node2, weight=weight1['weight'])
    file.close()
    # The nodes used in the forest
    nodeSet = SteinerTree(outputpath, resultfilename)[0]
    for node1 in nodeSet:
        for node2 in nodeSet:
            if G1.has_edge(node1,node2):
            #if G1.has_edge(node1,node2) == True:
                weight1 = G1.get_edge_data(node1, node2)
                #H1.add_edge(node1, node2, weight=weight1['weight'])
                I.add_edge(node1, node2, weight=weight1['weight'])
            if G2.has_edge(node1,node2):
            #if G2.has_edge(node1,node2) == True:
                weight1 = G2.get_edge_data(node1, node2)
                #H2.add_edge(node1, node2, weight=weight1['weight'])
                I.add_edge(node1, node2, weight=weight1['weight'])
    symbolfilename = "%s/symbol_fullnetwork_%s_%s_%s_%s.txt" % (outputpath, stpfile, str(w), str(b), str(D))
    symbolfile = open(symbolfilename,'w')
    #print len(H)
    for edge in I.edges():
        node1, node2 = edge
        #try:
        #node1 = idfile[node1]
        #except KeyError:
        node1 = [node1]
        #try:
        #node2 = idfile[node2]
        #except KeyError:
        node2 = [node2]
        if S.has_edge(edge[0], edge[1]) == True:
            for n1 in node1:
                for n2 in node2:
                    symbolfile.writelines('%s\tsteiner\t%s\n'%(n1,n2))
        if S.has_edge(edge[0], edge[1]) == False and I.has_edge(edge[0], edge[1]) == True:
            for n1 in node1:
                for n2 in node2:
                    symbolfile.writelines('%s\tintra\t%s\n'%(n1,n2))
        if I.has_edge(edge[0], edge[1]) == False:
            for n1 in node1:
                for n2 in node2:
                    symbolfile.writelines('%s\tinter\t%s\n'%(n1,n2))
    symbolfile.close()

    return H, H1, H2

# Unused
def NetworkMapping2(outputpath, stppath, stpfile, w, b, D):
    idfile = idmatchHuman(species)
    G1 = networkx.Graph()
    G2 = networkx.DiGraph()
    H = networkx.Graph()
    H1 = networkx.Graph()
    H2 = networkx.DiGraph()
    file = open(stppath+"/"+stpfile+'.stp', "r")
    while 1:
        line = file.readline()
        if line == "": break
        if line.startswith("E"): # or line.startswith("D"):
            temp = line.strip().split()
            G1.add_edge(temp[1], temp[2], weight=float(temp[3]))
        if line.startswith("D"): # or line.startswith("D"):
            temp = line.strip().split()
            G2.add_edge(temp[1], temp[2], weight=float(temp[3]))
    resultfilename = "%s_%s_%s_%s.txt" % (stpfile, str(w), str(b), str(D))
    nodeSet = SteinerTree(outputpath, resultfilename)[0]
    for node1 in nodeSet:
        for node2 in nodeSet:
            #if G1.has_edge(node1,node2) == True:
            if G1.has_edge(node1,node2):
                weight1 = G1.get_edge_data(node1, node2)
                H1.add_edge(node1, node2, weight=weight1['weight'])
                H.add_edge(node1, node2, weight=weight1['weight'])
            #if G2.has_edge(node1,node2) == True:
            if G2.has_edge(node1,node2):
                weight1 = G2.get_edge_data(node1, node2)
                H2.add_edge(node1, node2, weight=weight1['weight'])
                H.add_edge(node1, node2, weight=weight1['weight'])
    file.close
    symbolfilename = "%s/symbol_fullnetwork_%s_%s_%s_%s.txt" % (outputpath, stpfile, str(w), str(b), str(D))
    symbolfile = open(symbolfilename,'w')
    print len(H)
    for edge in H.edges():
        node1, node2 = edge
        try:
            node1 = idfile[node1]
        except KeyError:
            node1 = [node1]
        try:
            node2 = idfile[node2]
        except KeyError:
            node2 = [node2]
        for n1 in node1:
            for n2 in node2:
                symbolfile.writelines('%s %s\n'%(n1,n2))
    symbolfile.close()
    
    return H, H1, H2

# Writes the pairs of nodes in the message passing output
# Does not do the mapping presently
def sifidConverter(outputpath, stpfile, w, b, D, species):
    idfile = idmatchHuman(species)
    inputfile = "%s/%s_%s_%s_%s.txt" % (outputpath, stpfile, str(w), str(b),str(D))
    symbolfilename = "%s/symbol_%s_%s_%s_%s.txt" % (outputpath, stpfile, str(w), str(b), str(D))
    file = open(inputfile,'r')
    symbolfile = open(symbolfilename,'w')
#    print idfile.keys()[0]
    while 1:
        line = file.readline()
        if line == "": break
        node1, node2 = line.split()[0:2]
        #try:
        #    node1 = idfile[node1]
        #except KeyError:
        node1 = [node1]
        #try:
        #    node2 = idfile[node2]
        #except KeyError:
        node2 = [node2]
        for n1 in node1:
            for n2 in node2:
                symbolfile.writelines('%s %s\n'%(n1,n2))
    symbolfile.close()
    file.close()

# Unused
def interclusterconnectivity(T1,T2,H):
    inter = 0
    tempG = networkx.Graph()
    for node1 in T1:
        for node2 in T2:
            tempG.add_edge(node1, node2)
    for edges in tempG.edges():
        if H.has_edge(edges[0],edges[1]):
   #             edges = list(edges)
   #             if (node1 in edges) and (node2 in edges):
            inter += 1
             #       print edges
    return inter

# Unused
def intraclusterconnectivity(T, H):
    intra = 0
    for i in range(len(T)):
        for j in range(len(T)):
            if i < j:
                node1 = T[i]
                node2 = T[j]
                for edges in H.edges():
                    edges = list(edges)
                    if (node1 in edges) and (node2 in edges):
                        intra += 1
    return intra, len(T)

# Returns a list of numbers from start to end, non-inclusive of the end value
def frange(start, end=None, inc=None):
    "A range function, that does accept float increments..."

    if end == None:
        end = start + 0.0
        start = 0.0

    if inc == None:
        inc = 1.0
    if start == end:
        return [start]
    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)
        
    return L


def main():
    #PARSING PARAMETERS
    parser = OptionParser()

    parser.add_option("--outputpath", type="string", dest="outputpath", help="This path points to the directory where the output files will be written",default='None')
    parser.add_option("--msgpath",type="string",dest="msgpath",help="This path points to the directory where the message-passing code is available",default='None')
    parser.add_option("--depth",type="string",dest="depth",help="Depth parameter",default='None')
    parser.add_option("--conn",type="string", dest="conn", help="How to connect the artificial node to the interactome: 1 to all nodes in the interactome, 2 to all nodes in the interactome, except the terminals, 3 to a given set of nodes in the interactome, 4 to a given set of nodes in the interactome, except the terminals.")
    parser.add_option("--stppath", type="string", dest="stppath", help="This path points to the directory where the stp file is available.")
    parser.add_option("--stpfile", type="string", dest="stpfile", help="The name of the stp file (without the file extension)")
    parser.add_option("--W", type="string", dest="W",help="the range and increment for the artificial edge costs. format: start_end_increment")
    parser.add_option("--beta", type="string", dest="beta",help="the range and increment for the beta parameter. format: start_end_increment")
    parser.add_option("--exclude", type="string", dest="exclude", help="If set to 1 and the connection type is 2, 3, or 4 then remove prizes from all nodes that are connected to the artificial node")
    parser.add_option("--targetfile", type="string", dest="targetfile", help="a subset of the proteins in the interactome that will be connected to the artificial root in connection types 3 and 4",default='None')
    parser.add_option("--species", type="string", dest="species", help="the organism that you are working on",default='None')
    # New option to support multi-threaded message passing in msgsteiner9
    parser.add_option("--threads", type="int", dest="threads", help="The number of threads used to run the message passing algorithm.  Defaults to 1.", default=1)
    (options, args) = parser.parse_args()
    outputpath = options.outputpath
    msgpath = options.msgpath
    D = options.depth
    connectiontype = options.conn
    stppath = options.stppath
    stpfile = options.stpfile
    W = options.W
    beta = options.beta
    exclude = (options.exclude == '1')
    #exclude = options.exclude
    #if exclude == '1':
    #    exclude = True
    #if exclude == '2':
    #    exclude = False
    targetfile = options.targetfile
    species = options.species
    betainit, betaend, betaincrement = beta.split('_')
    Winit, Wend, Wincrement = W.split('_')
    betalist = frange(float(betainit), float(betaend), float(betaincrement))
    Wlist = frange(float(Winit), float(Wend), float(Wincrement))
    for w in Wlist:
        for b in betalist:
           
            # PrepareInputFile now returns a list of the lines in the input file instead of writing them to disk
            inputData = PrepareInputFile(stppath,stpfile,connectiontype, outputpath, w, b, D, [], exclude, targetfile, species)
            RunMSGAlgorithm(msgpath, D, w, b, outputpath, stpfile, connectiontype, inputData, options.threads)
            nodeList = OutputPCSFCheck(w, b, D, outputpath, stpfile,species)
            sifidConverter(outputpath, stpfile, w, b, str(D), species)
            H = NetworkMapping(outputpath, stppath, stpfile, w, b, D, species)[0]
            # Using subprocess now to send input to the message passing code, there is no file to remove
            #inputfile = "%s/%s_%s_%s_%s.dat" % (outputpath, stpfile, str(w), str(b),str(D))
            #os.system('rm %s' % inputfile)

main()

