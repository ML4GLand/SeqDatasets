# Starting with a new dataset

# Preparing for your data

## ML4GLand datasets

- If you are working with an IGVF dataset, I recommend you just generate a new repository from the https://github.com/IGVF-UCSD/dataset_template repo and clone that to wherever you are working

## Non-ML4GLand datasets

- Create a new folder on the cluster that will house your data
- Create a GitHub repo and link this to the new folder
- Follow the directory structure outlined here: https://github.com/ML4GLand/dataset_template

## Tips

- Name your dataset something descriptive and machine readable (no spaces!)

# Data acquisition

The data you are about to start working with inevitably exists somewhere other than our cluster, so we will need to move it over from where it currently is. This is sometimes non-trivial, so I’ve put together some code to help with different scenarios depending on where the data starts out (https://github.com/IGVF-UCSD/single_cell_utilities/tree/main/data_acquisition). You will likely need to write some custom scripts based on these to transfer the some combination of the follwowing

1. **Raw fastqs** —> these should be deposited within a `fastq/`subdirectory in your main dataset directory. Include the date of acquisition
2. **Processed** —> the exact form of this can vary substantially, but the outputs should be placed in a `processed/` subdirectory within the main directory. Include the date of acquisition
3. ******************Annotated****************** —> the exact form of this can vary substantially, but the outputs should be placed in an `annotation/` subdirectory within the main directory. Include the date of acquisition
4. **Analysis** products —> the exact form of this can vary substantially, but the outputs should be placed in an `analysis/` subdirectory within the main directory. Include the date of acquisition
5. ****************Metadata**************** —> the exact form of this can vary substantially, but the outputs should be placed in a `metadata/` subdirectory within the main directory. Include the date of acquisition

Before you move on, add a well-documented, short notebook describing the process of acquiring it to the `bin/` directory.

# Data processing

Data will usually come in either as raw sequencing **fastq** files, or in already processed and analyzed formats. If it is the former, you will need to process the data yourself. There are many ways to do this, and I have provided some code and pipelines you can run to get to the processed state. A couple things to note:

- Make sure that you version your processing pipeline, similar to data acquisition. Stick the results in a subdirectory within the`processed/` subdirectory with a date (e.g. `08Sep23`).
- Make sure the processing is reproducible in some way…TODO

# Data annotation

Most of the current downstream analyses performed on single cell data rely on being able to group cells together by some sort of biologically relevant property. This is most commonly cell-type, but can be the results of a more nuanced clustering of data as well (or a known covariate). The goal of this section is to get to those annotations through a sort of preliminary analysis of the data. This usually includes:

1. Computing quality control metrics on each sample within a dataset
2. Filtering low quality cells based on the computed metrics
3. Normalizing the data to deal with unwanted properties of the dataset (e.g. unstable variance across features)
4. Integrating samples to get all cells in the same feature space and to correct for unwanted technical variation (horizontal integration)
5. Integrating modalities (vertical or diagonal integration) to project a cell into a single feature space that captures multiple experimental observations
6. Clustering based on these feature spaces
7. Automatic annotation of clusters as a first pass
8. Manual refinement of annotation by experts

It is also possible that the data received has already undergone one or more of these steps already. **If this is the case, I recommend that you put together a few notebooks where you skeptically investigate the data in as many ways as possible.**

Whatever route you end up taking, create the following in a new `annotation/` subdirectory:

1. `{DatasetID}_CellClusterAssignment.tsv` — a metadata table with **one row per barcode** that includes:
    - bc
    - CellClusterID
    - RNAUMI
    - RNAGenes
    - ATACFragments (if multiome)
    - Any other additional columns you want (e.g. %mito, %ribo etc.)
2. `{DatasetID}_ClusterMetadata.tsv` — a separate table with cluster level metadata. This one should be **one row per cluster** and include:
    - CellClusterID — should match above for barcodes (e.g. the unique IDs here should exactly match the unique IDs in the CellCluster column above)
    - ManualAnnotationLabel —
    - nCells
    - MeanRNAUMIsPerCell
    - MeanATACFragmentsPerCell
3. `{DatasetID}_thresholds.tsv` — a separate table of per sample thresholds used for filtering cells (if different ones were used per sample)

Other useful files to generate during this stage include:

1. `reductions/` —> any dimensionality reductions
2. `loadings/` —> any feature loadings
3. `qc/` —> if per sample qc was done first

# Data analysis

Once you have a confident set of annotations, there are set of relatively general analysis steps you can perform. These are distinguished from more specific downstream analysis tasks that will depend on experimental design (e.g. these will likely be applicable to almost any single cell dataset). These include:

1. Identifying cluster markers for features at single cell resolution
2. Create pseudobulk representations of the dataset for each metadata field of interest. Put these in a subdirectory called `pseudobulk`. Note that this will depend on the assay, but the concept should be the same. Some examples:
    - A `tagAlign` for ATAC
    - A `tsv` for RNA
3. Summarizing features at pseudobulk resolution:
    - Identifying cluster markers
    - Calling peaks for ATAC
4. Creating feature matrices from called summarized features
    - fragment x cell matrices for ATAC

A couple reminders:

- **Version this data as well!** You may end up modifiying cluster assignments, or using a different tool for generating pseudobulks

# Data submission

Finally, if the dataset gets published, or you are a part of a consortium, you will need to submit the data. This usually involves creating a metadata table that points to all relevant files and running some kind of command line submission script. I will be adding resources to this as the pipeline for it becomes more robust, but a couple tips:

- Make sure you fully understand the experimental set-up prior to generating a metadata table. You should understand the nuances, including but not limited to the number of samples, number of biological replicates per sample, the number of technical replicates per sample, the sequencer used, the library protocol used, etc.
- Come up with a unique identifier **per file** that captures the
    - biosample type (e.g. condition, timepoint etc.) — 0hr_control (if relevant)
    - the sample ID (e.g. patient ID or internal sample ID) — DM0b
    - the modality used — scRNA, scATAC
    - the technical replicate number — _1 (shallow), _3 (deep)
    - the lane/plate of sequencing - 1, 2, 3, 4, 13G etc.
    - the read type - R1, R2, I1, R3
    
    An example could be `scATAC_24hr_3-cyt_dm21a_fastq_2_1_R1` which is the modality, biosample type, sample ID, file type (optional), lane, tehcnical replicate number and read type.
    

# Some final TODOs

- Create a README that documents the major steps followed for this dataset, and highlights any nuances worth noting
- You likely won’t be able to host most of the actual data files on GitHub due to their size. So make sure that you provide a link to the publicly available data files prior to publishing this repo. Also provide some code for someone to be able to download and hit the ground running

# A note on data wrangling

Inevitably, analysis tools require data to be written to disk in different formats to be usable. If you started from raw data, this usually isn’t a problem. If you were given some processed version of the data, it will likely be written to disk in a specific format (e.g. an special hdf5 file called `h5ad` for [Scanpy](https://www.notion.so/Scanpy-cc1d7c2ba4f447f79cbc43aebf991742?pvs=21) ). If you want to do an analysis in R, you’ll have to convert things. Unfortunately, this is very non-trivial and improper conversions can lead to a lot confusion and improper downstream analysis. To mitigate this, I’ve put together some code that covers many data wrangling cases in a reproducible fashion: https://github.com/IGVF-UCSD/single_cell_utilities/tree/main/data_wrangling. Note, that because many formats are tool specific, this code is mainly organized by tool.
