#!/bin/bash

# Set to the path to the msgsteiner executable
msgpath=./msgsteiner

interactomepath=../data
undirectedfile=STRING_PPI_GeneSymbol.txt
terminalpath=../data/TCGA_BreastCancer
terminalfile=SampleList.txt
resultpath=../results/BreastCancer
depth=10
W=1.0
iterations=5
alpha=2
beta=0.5
lambda=1.0
workers=1
artificialprizes=positiveWeighted
dummyneighbors=prizes
itermode=random

# Create the output directory if it does not exist
mkdir -p $resultpath

CMD="python ../ConstrainedMultiSample.py \
            --interactomepath=$interactomepath \
            --terminalpath=$terminalpath
            --resultpath=$resultpath \
            --undirectedfile=$undirectedfile \
            --beta=$beta \
            --terminalfile=$terminalfile \
            --msgpath=$msgpath \
            --depth=$depth \
            --W=$W \
            --iterations=$iterations \
            --lambda1=$lambda \
            --lambda2=$alpha \
            --workers=$workers \
            --artificialprizes=$artificialprizes \
            --dummyneighbors=$dummyneighbors \
            --itermode=$itermode \
            --workers=$workers"

echo $CMD
$CMD
