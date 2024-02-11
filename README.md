# Visualizing Next Generation Sequencing Data

This tool is an interactive visualisation of phenotypes (human and mouse) that have been associates with SNPs or other locations along the genome. The visualisation includes the annotation of each SNP and also the density distribution of all SNPs on the chromosomes. Using the density distribution clusters can be automatically defined by giving extra prerequisites, such as minimum cluster size and minimum number of included SNPs and a threshold value of the kde-score. The density between different datasets and experiments can also be visually compared to each other.


**Configurations**

Two files of datasets (tab-delimited csv-files) with each up to six different experiments can be uploaded. Configurations such as the header names for the datasets, genome type or genome reference can be set. To store the configurations for the next session, they can be exported in a text file and uploaded the next time.


**Simple Data Visualisation**

The SNPs and other locations on the genome are visualised as coloured lines on the chromosomes. Different experiments and datasets can be selected via checkboxes. Interactions such as zooming helps to get better insights in the dataset.


**Density Visualisation**

Through the use of kernel density estimation with the Epanechnikov kernel the density distribution from the given SNPs is defined. The calculation can be adapted to personal requirements by changing different parameters. The density is then used to define clusters. For details of the cluster a link to the UCSC genome browser is included.


**Density Comparison**

The density between two different dataset or experiments can be visual compared by setting one feature to green and the other one to red. The more the area is marked as yellow, the features are aligned.
