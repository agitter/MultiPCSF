__author__ = 'Nurcan Tuncbag' # Modified by Anthony Gitter

import sys
from optparse import OptionParser

def Kinase_Substrate(interactomepath, filename):
    kipdict = {}
    file = open(interactomepath+"/"+filename,"r")
    while 1:
        line = file.readline()
        if line == "": break
        # Make the protein names upper case
        temp = line.strip().upper().split()
        kinase = temp[1]
        substrate = temp[0]
        weight = temp[2]
        kippair = [kinase, substrate]
        kippair.sort()
        # Use the same number of underscores as the undirected dictionary
        kipdict[kippair[0]+"____"+kippair[1]] = substrate + "_" + kinase + "_" + weight
    file.close()
    return kipdict


# Change the default directedfilename from None to "None"
def combine_directed_indirected_PPI(interactomepath, indirectedfilename, directedfilename="None"):
    file = open(interactomepath+'/'+indirectedfilename,"r")
    interactiondict = {}
    combined = []
    interactome = set()
    while 1:
        line = file.readline()
        if line == "": break
        # Make upper case
        temp = line.strip().upper().split()
        node1 = temp[0]
        node2 = temp[1]
        # Skip lines that don't have exactly 3 parts, ignoring
        # lines without weights and nodes with spaces in their names
        if len(temp) != 3:
            print "Skipping: " + line
        else:
            weight = float(temp[2])

            weight = weight2Cost(weight)
            # Don't filter high cost edges, or make the threshold a parameter
            #if node1 != node2 and float(weight) <= 0.5:
            if node1 != node2:
                pair = [node1,node2]
                pair.sort()
                interactome.add(node1)
                interactome.add(node2)
                interactiondict[pair[0]+"____"+pair[1]] = str(weight)
    file.close()
    if directedfilename == "None":
        for interaction in interactiondict.keys():
            n1 = interaction.split("____")[0]
            n2 = interaction.split("____")[1]
            combined.append(["E", n1, n2, interactiondict[interaction]])
    else:
        kipdict = Kinase_Substrate(interactomepath, directedfilename)
        for interaction in interactiondict.keys():
            try:
                kipint = kipdict[interaction]
                kinaseweight = float(kipint.split("_")[2])

                kinaseweight = weight2Cost(kinaseweight)
                combined.append(["D", kipint.split("_")[0], kipint.split("_")[1], kinaseweight])
            except KeyError:
                n1 = interaction.split("____")[0]
                n2 = interaction.split("____")[1]
                combined.append(["E", n1, n2, interactiondict[interaction]])
        for item in kipdict.keys():
            node1 = kipdict[item].split("_")[0]
            node2 = kipdict[item].split("_")[1]
            kinaseweight = float(kipdict[item].split("_")[2])
            kinaseweight = weight2Cost(kinaseweight)
            if item not in interactiondict.keys():
                combined.append(["D", node1, node2, kinaseweight])
                interactome.add(node1)
                interactome.add(node2)
    interactome = list(interactome)
    return combined, interactome

def connectTF_DNA(interactomepath, tfdnafilename, mrnaterminals, interactome):
    tfdna = set()
    if tfdnafilename == 'None':
       return tfdna
    else:
        file = open(interactomepath+'/'+tfdnafilename,"r")
        while 1:
            line = file.readline()
            if line == "": break
            # Make names upper case
            temp = line.upper().split()
            node1 = temp[0]
            node2 = temp[1]
            weight = float(temp[2])
            # Check for _MRNA
            if (node1 in interactome) and (node2+"_MRNA" in mrnaterminals):
                # Use the weight in the file
                tfdna.add("D %s %s %f\n" % (node2+"_MRNA", node1, weight2Cost(weight)))
        print "%d TF-DNA edges" % len(tfdna)
        file.close()
        return tfdna

def terminalNodes(terminalpath, terminalfile, beta, mrnabeta):
    file = open(terminalpath+'/'+terminalfile,"r")
    terminals = []
    mrnaterminals = []
    while 1:
        line = file.readline()
        if line == "": break
        # Make names upper case
        temp = line.strip().upper().split()
        nodename =  temp[0]
        weight = float(temp[1])
        # Check for _MRNA
        if nodename.endswith("_MRNA"):
            terminals.append(["W", nodename, mrnabeta*weight])
            mrnaterminals.append(nodename)
        else:
            terminals.append(["W", nodename, beta*weight])
    file.close()
    return terminals, mrnaterminals

# Force the weight into the interval [0, 0.99] and then
# transform to a cost by taking 1-weight
def weight2Cost(weight):
    if weight > 0.99:
        weight = 0.99
    elif weight < 0:
        weight = 0
    return 1 - weight

def main():
    #PARSING PARAMETERS
    parser = OptionParser()

    # Updated the parameter descriptions to match the code change
    parser.add_option("--interactomepath", type="string", dest="interactomepath", help="This path points to the directory where all interaction files are deposited (i.e., protein-protein interactions, kinase-substrate interactions, transcription factor-DNA interactions)",default='None')
    parser.add_option("--terminalpath",type="string",dest="terminalpath",help="This path points to the directory where the terminal file is deposited.",default='None')
    parser.add_option("--resultpath",type="string",dest="resultpath",help="This path points to the directory where the outputs will be located.",default='None')
    parser.add_option("--indirectedfilename",type="string",dest="indirectedfilename",help="The name of the interaction file where protein-protein interaction data with probabilistic weights (e.g in [0,1]) are available.",default='None')
    parser.add_option("--beta",type="string",dest="beta",help="Beta parameter given here is to scale node penalties of protein terminals to the edge costs",default='None')
    parser.add_option("--terminalfile",type="string",dest="terminalfile",help="The name of the terminal file where the set of the terminal nodes with the node penalties are available.  If gene penalties are given, gene names should end with '_MRNA'",default='None')
    parser.add_option("--resultfilename",type="string",dest="resultfilename",help="The name of the file where the combined information will be written to be used as input in the message passing tool.",default='None')

    # Change the default from None to 'None'
    parser.add_option("--directedfilename",type="string",dest="directedfilename",help="Optional: The name of the interaction file where directed interactions (i.e. kinase-substrate) with probabilistic weights are available.  Columns should be ordered [substrate kinase weight].",default='None')
    parser.add_option("--tfdnafilename",type="string",dest="tfdnafilename",help="Optional: The name of the interaction file where TF-DNA interactions with probabilistic weights are available. Columns should be ordered [TF gene weight] and gene names should not end with '_MRNA' because '_MRNA' is automatically appended to them.",default='None')
    parser.add_option("--mrnabeta",type="string",dest="mrnabeta",help="Optional: The beta parameter for gene terminal nodes. Its default value is equal to the --beta value.  It is applied to all terminals whose name ends with '_MRNA'.",default='None')
    parser.add_option("--root",type="string",dest="root",help="Optional: The name of the root node.",default='None')

    (options, args) = parser.parse_args()
    interactomepath = options.interactomepath
    terminalpath = options.terminalpath
    indirectedfilename = options.indirectedfilename
    directedfilename = options.directedfilename
    tfdnafilename = options.tfdnafilename
    mrnabeta = options.mrnabeta
    beta = options.beta
    terminalfile = options.terminalfile
    resultfilename = options.resultfilename
    resultpath = options.resultpath
    root = options.root
    if interactomepath == 'None' or terminalpath == 'None' or indirectedfilename == 'None' or beta == 'None' or resultpath == 'None' or terminalfile == 'None' or resultfilename == 'None':
        print "***"
        print "***ERROR: There are not enough arguments given by the user. Please enter at least the required arguments!!!"
        print "***"
        return -1
    else:
        beta = float(beta)
        if mrnabeta == 'None':
            mrnabeta = beta
        else:
            mrnabeta = float(mrnabeta)
        combined, interactome = combine_directed_indirected_PPI(interactomepath, indirectedfilename, directedfilename)
        terminal, mrnaterminals = terminalNodes(terminalpath, terminalfile, beta, mrnabeta)
        tfdna = connectTF_DNA(interactomepath, tfdnafilename, mrnaterminals,interactome)
        file = open(resultpath+"/"+ resultfilename, "w")
        for item in combined:
            file.writelines('%s %s %s %s\n' % (item[0], item[1], item[2], item[3]))
        if tfdna != set():
            for item in tfdna:
                file.writelines(item)
        for item in terminal:
            # Don't print any gene prizes if there are no TF-gene edges
            if tfdna == set():
                # Check for "_MRNA"
                if not item[1].endswith("_MRNA"):
                    file.writelines('%s %s %s\n' % (item[0], item[1], item[2]))
            if tfdna != set():
                file.writelines('%s %s %s\n' % (item[0], item[1], item[2]))
        if root != 'None':
            # Capitalize the root
            file.writelines('R %s\n' % root.upper())
        file.close()

        file = open(resultpath + "/" + resultfilename.split(".")[0] + ".attr","w")
        file.writelines("terminalAttr\n")
        for item in terminal:
            file.writelines(item[1] + " = 1\n")
        file.close()

        # Log the arguments provided
        file = open(resultpath + "/" + resultfilename.split(".")[0] + ".log","w")
        for arg in sys.argv:
            file.write(arg + " ")
        file.write("\n")
        file.close()

        return 1

if __name__=="__main__":
    main()

