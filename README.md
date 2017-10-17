[Gitter et al 2014]: http://www.worldscientific.com/doi/abs/10.1142/9789814583220_0005
[Omics Integrator]: https://github.com/fraenkel-lab/OmicsIntegrator
[msgsteiner]: http://areeweb.polito.it/ricerca/cmp/code/bpsteiner
[TCGA 2012]: http://www.nature.com/nature/journal/v490/n7418/full/nature11412.html
[Szklarczyk et al 2011]: http://nar.oxfordjournals.org/content/39/suppl_1/D561.long
[TCGA]: http://cancergenome.nih.gov/publications/publicationguidelines
[STRING]: http://string-db.org/cgi/access.pl?footer_active_subpage=licensing
[Database of Cell Signaling]: http://stke.sciencemag.org/about/help/cm
[Gough 2002]: https://doi.org/10.1111/j.1749-6632.2002.tb04532.x
[Microsoft Research]: https://www.microsoft.com/en-us/research/lab/microsoft-research-new-england/

# Multi-PCSF
[![DOI](https://zenodo.org/badge/47654267.svg)](https://zenodo.org/badge/latestdoi/47654267)

This repository contains an implementation of the multi-sample prize-collecting
Steiner forest (Multi-PCSF) algorithm described in [Gitter et al 2014]. This
version of the PCSF code is provided for reproducibility of the results in the
manuscript but is no longer under active development.  Please see [Omics
Integrator] for the current version of PCSF from Ernest Fraenkel's lab.
However, Omics Integrator does not yet support the multi-sample feature
introduced by Multi-PCSF.  The Omics Integrator website describes how to install
the [msgsteiner] dependency, which is also required by Multi-PCSF.

## Example
`BreastCancer.sh` in the scripts subdirectory provides an example of how to run
Multi-PCSF.  Before running the script, the `msgpath` variable must be set to
the location of the msgsteiner executable, including the file name.

## Data
The breast cancer tumor sample data and protein-protein interaction network data
described in the Multi-PCSF manuscript are provided as an example dataset.  If
you use these data in a manuscript, cite [TCGA 2012] for the breast cancer data
and [Szklarczyk et al 2011] for the STRING protein-protein interaction network
and see their respective websites ([TCGA], [STRING]) for the terms of use.

The *Science Signaling* Database of Cell Signaling EGFR pathway that was used to
simulate samples is also provided in the `data` subdirectory.  If you use this
pathway in a manuscript, cite [Gough 2002] and see the [Database of Cell
Signaling] website for the terms of use.

## Usage
Only the most commonly used options are described below.  Use `python
ConstrainedMultiSample.py -h` to view the complete usage message. See the
provided example data for file formatting guidelines.  Please open an issue with
any usage questions.
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

## Output
Several subdirectories are created in the directory specified by the
`--resultpath` argument.  The `initial` and `itr*` directories (one for each of
the iterations specified by the `--iterations` argument) provide detailed
information about intermediate results.  Except for the last `itr*` directory,
these can typically be deleted after Multi-PCSF terminates.

The location of the final Multi-PCSF networks depends on the settings. If
`--artificialprizes` was set to one of the negative prize options or only one
iteration was run, the output networks are in the last `itr*` directory.  If
positive artificial prizes were used, a post-processing pruning step is
executed.  This runs the Steiner forest algorithm once more for each sample to
prune nodes in the network that do not connect real prize nodes to the forest
but rather were included only due to their positive artificial prizes.  In this
case, the output networks are in the `final` directory.

The output directory contains intermediate files and the following files that
are most useful for interpreting and visualizing the networks.  For each input
file `<sample>` listed in the `--terminalfile` input file, several output files
will be created:
* `symbol_<sample>_<options>.txt`: `<sample>` is the input sample name and
`<options>` are the values of the `W`, `beta`, and `depth` arguments. This
space-separated file contains a line for each edge in the output network, where
each line provides the names of the interacting proteins. The artificial root
node `DUMMY` is still present. This is typically the most relevant
representation of the output network. The edges are the same as the edges in the
msgsteiner output file `<sample>_<options>.txt`.
* `symbol_fullnetwork_<sample>_<options>.txt`: `<sample>` is the input sample name
and `<options>` are the values of the `W`, `beta`, and `depth` arguments. This
tab-separated file contains a line for each edge in the output network. The
artificial root node has been removed.  The `steiner` edges are the edges from
the optimal Steiner forest.  The `intra` edges are additional edges that have
been added back to the Steiner forest, which are sometimes useful for
identifying alternative pathway connections.
* `<sample>_<options>.output`: Summary statistics of the Steiner forest produced.
* `<sample>_<options>.objective`: Output messages from the msgsteiner program,
including optimization progress.

The other files are intermediate files used to create the input for msgsteiner
or prepare the output network file from the msgsteiner output.

## Simulation
The `simulation` subdirectory contains the code that was used to simulate
input samples from synthetic or real pathways.  This code currently serves as
extended documentation and is not runnable.  It uses an old version of
`ConstrainedMultiSample.py` and needs to be updated to use the refactored
version, which accepts different command line arguments.

## Roadmap
* Implement multi-sample functionality in Omics Integrator
* Refactor simulation code to generate prizes from known or synthetic pathways
* Document support for distinct groups of samples

## Developers
* Nurcan Tuncbag
* Anthony Gitter

## Acknowledgements
Portions of the Multi-PCSF software were developed with support from [Microsoft
Research] while Anthony Gitter was a postdoctoral researcher there.  We thank
Microsoft for granting permission to release the code as open source under the
Sample Code Exception and Paul Oka in particular for coordinating the release.
We acknowledge all authors of [Gitter et al 2014] for their role in the
algorithm development.
