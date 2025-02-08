# GOTHiC: Gene Ontology Topology from Hi-C

GOTHiC is a python tool that allows for the creation, annotation, and analysis of Hi-C networks. In essence, it leverages graph theory techniques to allow for the assessment of structure-function relationships in nuclear chromatin.

![A summary of GOTHiC's primary workflow](https://github.com/LavalleeAdamLab/GOTHiC/assets/61287366/7d13adaa-5a70-4c1b-b0a7-7116a9736b1b)

## Using GOTHiC:
To use GOTHiC, copy `gothic.py` into a directory within your python PATH, and simply 'from gothic import *'. Individual functions are available for specific tasks involving the conversion of .sam or .hic files to Hi-C graphs, graph functional annotation with Gene Ontology terms or REACTOME metabolic pathways, clustering analysis, and visualization. Function documentation can be found in the Docs/ directory, and scripts are available in the Scripts/ directory for running the most typical workflow. These scripts can be used for Hi-C network creation and normalization from SAM file, graph annotation with Gene Ontology Terms, clustering analysis, A/B compartment analysis, and visualization.  

Dependencies:
graph-tool (install [here](https://graph-tool.skewed.de/installation.html))

required packages: pandas numpy matplotlib scipy tqdm plotly goatools seaborn fanc iced gprofiler sklearn markov_clustering 

![Heatmap visualization of A/B compartments and Network visualization of GO terms with GOTHiC](https://github.com/user-attachments/assets/45a5e6b3-aba6-4c50-97c8-3112fcdb055e)





