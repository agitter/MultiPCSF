[Gitter et al 2014]: http://www.worldscientific.com/doi/abs/10.1142/9789814583220_0005
[Omics Integrator]: http://fraenkel.mit.edu/omicsintegrator
[msgsteiner]: http://areeweb.polito.it/ricerca/cmp/code/bpsteiner
[TCGA 2012]: http://www.nature.com/nature/journal/v490/n7418/full/nature11412.html
[Szklarczyk et al 2011]: http://nar.oxfordjournals.org/content/39/suppl_1/D561.long
[TCGA]: http://cancergenome.nih.gov/publications/publicationguidelines
[STRING]: http://string-db.org/newstring_cgi/show_download_page.pl

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
([TCGA], [STRING]) for the terms of use.

## Usage
Only the most commonly used options are described below.  Use
`python ConstrainedMultiSample.py -h` to view the complete usage message.
See the provided example data for file formatting guidelines.
```
Usage: ConstrainedMultiSample.py [options]

Options:
  -h, --help            show this help message and exit
  --interactomepath=INTERACTOMEPATH
                        This path points to the directory that contains the
                        interaction network files
  --terminalpath=TERMINALPATH
                        This path points to the directory that contains the
                        terminal (node prize) files
  --resultpath=RESULTPATH
                        This path points to the directory where the output
                        files will be written.
  --undirectedfile=UNDIRECTEDFILE
                        The name of the protein-protein interaction file in
                        the interactomepath directory.  The file is expected
                        to contain undirected interactions with probabilistic
                        weights (e.g in [0,1]). Columns should be ordered
                        [prot1 prot2 weight].
  --terminalfile=MASTERTERMINALFILE
                        A file in the terminalpath directory that lists the
                        files that give the node prizes for each sample.  All
                        listed filenames should be relative to terminal path.
                        If gene penalties are given in the terminal files,
                        gene names should end with '_MRNA'.  Optionally can
                        include a tab-separated second column that assigns
                        each sample to a group so the forests are only
                        constrained to be similar to other samples in the same
                        group.
  --msgpath=MSGPATH     The path and file name of the msgsteiner executable
  --depth=DEPTH         Depth parameter that limits the maximum depth from the
                        Steiner tree root to the leaves
  --W=W                 The cost of the edges from the artificial root node to
                        its neighbors.
  --beta=BETA           The scaling factor applied to the node prizes, which
                        is used to control the relative strength of node
                        prizes and edge costs.  This scaling is only performed
                        once when the initial stp files are created.
  --lambda=LAMBDA1      The tradeoff coefficient for the penalty incurred by
                        nodes in the Steiner forests that are not in the set
                        of common nodes.
  --alpha=LAMBDA2       The tradeoff coefficient for the reward on the size of
                        the set of common nodes when using unweighted
                        artificial prizes or the power to which the node
                        frequency is taken for weighted prizes.
  --mu=MU               A parameter used to penalize high-degree nodes from
                        being selected as Steiner nodes.  Does not affect
                        prize nodes but does affect artificial prizes.  The
                        penalty is -mu*degree.  Set mu <= 0 to disable the
                        penalty (default).
  --iterations=ITERATIONS
                        The number of iterations to run
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
```

## Developers
* Nurcan Tuncbag
* Anthony Gitter