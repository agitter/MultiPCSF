[Gitter et al 2014]: http://www.worldscientific.com/doi/abs/10.1142/9789814583220_0005
[Omics Integrator]: http://fraenkel.mit.edu/omicsintegrator
[msgsteiner]: http://areeweb.polito.it/ricerca/cmp/code/bpsteiner
[TCGA 2012]: http://www.nature.com/nature/journal/v490/n7418/full/nature11412.html
[Szklarczyk et al 2011]: http://nar.oxfordjournals.org/content/39/suppl_1/D561.long

# Multi-PCSF
An implementation of the Multi-PCSF algorithm described in [Gitter et al 2014].
This is an early, rough version of the PCSF code and is no longer under active
development.  Please see [Omics Integrator] for the current version of PCSF
from Ernest Fraenkel's lab.  However, Omics Integrator does not yet support the
multi-sample feature introduced by Multi-PCSF.  The Omics Integrator website
describes how to install the [msgsteiner] dependency, which is also required
by Multi-PCSF.

## Example
BreastCancer.sh in the scripts subdirectory provides an example of how to run
Multi-PCSF.  Before running the script, the `msgpath` variable must be set to
the location of the msgsteiner executable, including the file name.

## Data
The breast cancer tumor sample data and protein-protein interaction network
data described in the Multi-PCSF manuscript are provided as an example
dataset.  If you use these data in a manuscript, cite [TCGA 2012] for
the breast cancer data and [Szklarczyk et al 2011] for the STRING
protein-protein interaction network and see their respective websites
for the terms of use.

## Usage
```
Usage: ConstrainedMultiSample.py [options]

Options:
  -h, --help            show this help message and exit
  --interactomepath=INTERACTOMEPATH
                        This path points to the directory where all
                        interaction files are deposited (i.e., protein-protein
                        interactions, kinase-substrate interactions,
                        transcription factor-DNA interactions)
  --terminalpath=TERMINALPATH
                        This path points to the directory where the terminal
                        files are deposited.
  --resultpath=RESULTPATH
                        This path points to the directory where the outputs
                        will be located.
  --undirectedfile=UNDIRECTEDFILE
                        The name of the interaction file where protein-protein
                        interaction data with probabilistic weights (e.g in
                        [0,1]) are available. Columns should be ordered [prot1
                        prot2 weight].
  --beta=BETA           Beta parameter given here is to scale node penalties
                        of protein terminals to the edge costs.  This scaling
                        is only performed once when the initial stp files are
                        created.
  --terminalfile=MASTERTERMINALFILE
                        A file in terminalpath that lists the files that give
                        the node prizes for each sample.  All listed filenames
                        should be relative to terminal path.  If gene
                        penalties are given in the terminal files, gene names
                        should end with '_MRNA'.  Optinonally can include a
                        tab-separated second column that assigns each sample
                        to a group so the forests are only constrained to be
                        similar to other samples in the same group.
  --directedfile=DIRECTEDFILE
                        Optional: The name of the interaction file where
                        directed interactions (i.e. kinase-substrate) with
                        probabilistic weights are available.  Columns should
                        be ordered [substrate kinase weight].
  --tfdnafile=TFDNAFILE
                        Optional: The name of the interaction file where TF-
                        DNA interactions with probabilistic weights are
                        available. Columns should be ordered [TF gene weight].
                        Gene names should not end with '_MRNA' because '_MRNA'
                        is automatically appended to them.
  --mrnabeta=MRNABETA   Optional: The beta parameter for gene terminal nodes.
                        Its default value is equal to the --beta value.  It is
                        applied to all terminals whose name ends with '_MRNA'.
                        The scaling is only performed once when the initial
                        stp files are created.
  --msgpath=MSGPATH     The path and file name of the msgsteiner executable
  --depth=DEPTH         Depth parameter
  --W=W                 The cost of the edges from the artificial root node to
                        its neighbors.
  --iterations=ITERATIONS
                        The number of iterations to run
  --lambda=LAMBDA1      The tradeoff coefficient for the penalty incurred by
                        nodes in the Steiner forests that are not in the set
                        of common nodes.
  --alpha=LAMBDA2       The tradeoff coefficient for the reward on the size of
                        the set of common nodes when using unweighted
                        artificial prizes or the power to which the node
                        frequency is taken for weighted prizes.
  --workers=WORKERS     The number of worker processes to use in the
                        multiprocessing pool or threads to use in multi-
                        threaded belief propagation.  Should not exceed the
                        number of cores available.  Defaults to the number of
                        CPUs.
  --artificialprizes=ARTIFICIALPRIZES
                        Use 'positive' or 'negative' prizes to encourage trees
                        to include common set proteins.  Use
                        'positiveWeighted' or 'negativeWeighted' (default)
                        prizes to construct weighted artificial prizes based
                        on the node frequency in the most recent forests.
  --dummyneighbors=DUMMYNEIGHBORS
                        Connect the dummy node to all 'prizes' (default) or
                        'nonprizes'.
  --itermode=ITERMODE   Learn forests simultaneously in 'batch' (default) or
                        sequentially in 'random' order.  Batch mode computes
                        artificial prizes with respect to all forests at the
                        previous iteration.  Random mode computes prizes for a
                        specific sample given the most recent forests for all
                        other samples.
  --mu=MU               A parameter used to penalize high-degree nodes from
                        being selected as Steiner nodes.  Does not affect
                        prize nodes but does affect artificial prizes.  The
                        penalty is -mu*degree.  Set mu <= 0 to disable the
                        penalty (default).
```

## Developers
* Nurcan Tuncbag
* Anthony Gitter