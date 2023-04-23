
# MUPI - a mutational-network analysis pipeline

MUPI is an open source network analysis pipeline for understanding relations between mutations in cancer.


## Installation
To run this pipeline use python 3+. You also need to install required libraries. Install them using command below:

```
pip install -m requiermments.txt
```
## Running
### Quick Run

To run all analysis run run_all.py folder. Generate a result and data folder and set their path in the main function of this file. Put all of the data listed below in the data folder:

#### Somatic mutation file (simple_somatic_mutation.open.tsv)
This will be used to generate the mutational network. It should be a tsv file with rows showing mutations in patients' genes. This file should be in icgc format with at least two columns. gene_affected for the gene_id of mutated gene and icgc_specimen_id for the id of the patient with that mutation. The default name for this file is simple_somatic_mutation.open.tsv and it should be located in the data folder.

#### Gene list (gene_list_hg19.xlsx)
This file contains a mapping between gene-id and gene-name. It also contains information about coding or non-coding genes and the location and strand of each gene.

#### Ppi files (9606.protein.links.v11.5 and 9606.protein.info.v11.5)
Download 9606.protein.links.v11.5 and 9606.protein.info.v11.5 from string website [(this link)](https://string-db.org/cgi/download?sessionId=b0rskhtvOySj&species_text=Homo+sapiens).
Ppi interactions are in the 9606.protein.links.v11.5 file. Genes related to each protein is in the 9606.protein.info.v11.5 file.

### Running different modules
You can also run each module separately. The description is as follow:

#### Network generation

To generate the mutational network run the network_generator.py file. Generate an instance of a MutationalNetwork object. Input data folder and data file name when initiating this object. Call the generate_network function with desired threshold as input. The output of this function would be an instance of networkx graph object. It can be saved using save_network function. The input data should be in the format of icgc-somatic mutation data. More specifically it should be a table formatted file (like csv, tsv or xlsx files) with a gene column named gene_affected and patient column named icgc_specimen_id. These values can be saved in gene_key and patient_key and it can be changed after initiating the MutationalNetwork object.


#### Threshold selection


#### Generated network evaluation
To evaluate the generated network run the network/evaluator.py file. In the main function it will generate an instance of NetworkEvaluator and run the evaluate function with the desired network (input with networkx type). You can change the sample network path to load your network. It will give you a report of the number of coding and non-coding nodes and links in the network. It uses the gene_list_hg19.xlsx file in the data folder to distinguish between coding and non-coding genes.


#### Community detection
All of the community detection algorithms are in the community_detection/community_detection.py file. You can load your networkx network and run this file with the desired community detection algorithm. By default in the main function results/network0.15.pkl network will be loaded and communities will be found using hierarchical ensemble algorithm. The results will be saved in the results folder.

#### Community evaluation
This module will make a report of patient coverage of founded communities. It will use the gene-patient.pkl file which is made by the data_preprocessor.py module. It will generate two excel files coms-patient-coverage.csv and patient-coms-covered.csv. The first one includes the number and ratio of patients covered by each community. The second one includes the number of communities that cover each patient. You can check the plots of these two reports by calling the plot method of this module.

#### Biological analysis
**Ppi interactions**

First run the biological_analysis/ppi/generate_network.py file to generate the ppi-interaction network. You can change the confidence min score for putting an edge between two proteins in the network. The network will be saved in the result folder as ppi-net-{confidence-level}.pkl. 
Then run the biological_analysis/ppi/ppi_analysis.py file to analyze the network and the communities. This module will use the ppi-net-{confidence-level}.pkl and gene-name-dict.pkl files generated in the result folder. evaluate_coms_ppi_links function will generate many random subnetworks for each community (with the same number of nodes and edges) and compute the ratio of the ppi-links for each generated random network. Then it will check if the average of ppi links in a random network is greater in our communities or not. It will calculate the enrichment value for each community which is the ratio of ppi-links in our community to the average number of ppi-links in all random communities. The results will be saved in ppi-analysis-rand{random-cnt}.csv in the results folder.

**Gene set enrichment (gene ontology and kegg pathway)**

This module will use stringdb api to find enriched (with fdr <= 0.05) gene ontology and kegg pathways related to each community. It will use the gene-name-dict.pkl file to convert gene_id to gene name for genes in each community, the results will be saved in the result folder as gene-set-enrichment.csv file. 

**Cancer pathways analysis**


**Gene location analysis** 

Using this module you can get a report of number cis and trans links and number of links with genes on the same strand in the network and each community. It will also give you a report of avg distance between genes of cis links. This module will use the gene-location-dic.pkl file which is generated by the preprocess_data.py module. 

**Literature study**
