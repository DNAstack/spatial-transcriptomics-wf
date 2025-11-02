# spatial-transcriptomics-wf
Repo for testing and developing a common postmortem-derived brain sequencing (PMDBS) and non-human workflow harmonized across ASAP with human and mouse spatial transcriptomics data for both Nanostring GeoMx and 10x Visium platforms. The main goal is to uncover spatially distinct gene expression profiles across tissue samples, enabling insights into tissue architecture, cell type composition, and disease-related molecular patterns.

Common workflows, tasks, utility scripts, and docker images reused across harmonized ASAP workflows are defined in [the wf-common repository](https://github.com/ASAP-CRN/wf-common).


# Table of contents

- [Workflows](#workflows)
- [Inputs](#inputs)
- [Outputs](#outputs)
	- [Output structure](#output-structure)
- [Docker images](#docker-images)


# Workflows

Worfklows are defined in [the `workflows` directory](workflows). There is [the `spatial_geomx` workflow directory](workflows/spatial_geomx) and [the `spatial_visium` workflow directory](workflows/spatial_visium).

These workflows are set up to analyze spatial transcriptomics data: Nanostring GeoMx in WDL using command line, R, and Python scripts and 10x Visium in WDL using command line and Python scripts.

_Note: Unlike our other workflows (e.g., [pmdbs-sc-rnaseq-wf](https://github.com/ASAP-CRN/pmdbs-sc-rnaseq-wf/tree/main)), the spatial workflows do not perform cross-team cohort analysis, as integrating spatial coordinates and mapping dimensions across datasets adds significant complexity._

## Nanostring GeoMx workflow overview

For the Nanostring GeoMx workflow, we start with raw output files from the instrument and convert them into counts. Then, we clean the data by removing unreliable segments and genes, adjust for technical noise, and combine data from different slides. Finally, we cluster based on segment (which may contain many cells) and visualize their transcriptional profiles in a UMAP space.

**Nanostring GeoMx workflow diagram:**
![Nanostring GeoMx workflow diagram](workflows/spatial_geomx/workflow_diagram.svg "Workflow diagram")

**Nanostring GeoMx entrypoint**: [workflows/spatial_geomx/main.wdl](workflows/spatial_geomx/main.wdl)

**Nanostring GeoMx input template**: [workflows/spatial_geomx/inputs.json](workflows/spatial_geomx/inputs.json)

The Nanostring GeoMx workflow is broken up into three main chunks:

1. [Preprocessing](#preprocessing)
1. [Process to adata](#process-to-adata)
1. [Cohort analysis](#cohort-analysis)

### Preprocessing

Run once per slide; only rerun when the preprocessing workflow version is updated. Preprocessing outputs are stored in the originating team's raw and staging data buckets.

### Process to adata
Run once per slide. Intermediate files from previous runs are not reused and are stored in timestamped directories.

### Cohort analysis

Run once per team (all slide from a single team) if `project.run_project_cohort_analysis` is set to `true`. Additional slides requires this entire analysis to be rerun. Intermediate files from previous runs are not reused and are stored in timestamped directories.

## 10x Visium workflow overview

For the 10x Visium workflow, we start with raw output files from the instrument and convert them into counts. Then, we clean the data by removing unreliable spots and genes, adjust for technical noise, and combine data from different samples. We cluster based on spots (which may contain many cells) and visualize their transcriptional profiles in a UMAP space. We also identify spatially variable genes to uncover expression patterns that are location-specific within the tissue.

**10x Visium workflow diagram:**
![10x Visium workflow diagram](workflows/spatial_visium/workflow_diagram.svg "Workflow diagram")

**10x Visium entrypoint**: [workflows/spatial_visium/main.wdl](workflows/spatial_visium/main.wdl)

**10x Visium input template**: [workflows/spatial_visium/inputs.json](workflows/spatial_visium/inputs.json)

The 10x Visium workflow is broken up into two main chunks:

1. [Preprocessing](#preprocessing)
2. [Cohort analysis](#cohort-analysis)

### Preprocessing

Run once per sample; only rerun when the preprocessing workflow version is updated. Preprocessing outputs are stored in the originating team's raw and staging data buckets.

### Cohort analysis

Run once per team (all samples from a single team) if `project.run_project_cohort_analysis` is set to `true`. This can be rerun using different sample subsets; including additional samples requires this entire analysis to be rerun. Intermediate files from previous runs are not reused and are stored in timestamped directories.


# Inputs

## Nanostring GeoMx inputs

An input template file can be found at [workflows/spatial_geomx/inputs.json](workflows/spatial_geomx/inputs.json).

| Type | Name | Description |
| :- | :- | :- |
| Array[[Project](#nanostring-geomx-project)] | projects | The project ID, set of slides and their associated samples, reads and metadata, output bucket locations, and whether or not to run project-level cohort analysis. |
| File | geomxngs_config_pkc | The GeoMx DSP configuration file to associate assay targets with GeoMx HybCode barcodes and Seq Code primers; see https://nanostring.com/products/geomx-digital-spatial-profiler/geomx-dsp-configuration-files/. |
| Int? | min_segment_reads | Minimum number of segment reads. [1000] |
| Int? | min_percent_reads_trimmed | Minimum % of reads trimmed. [80] |
| Int? | min_percent_reads_stitched | Minimum % of reads stitched. [80] |
| Int? | min_percent_reads_aligned | Minimum % of reads aligned. [80] |
| Int? | min_saturation | Minimum sequencing saturation. [50] |
| Int? | min_neg_ctrl_count | Minimum negative control counts. [1] |
| Int? | max_ntc_count | Maximum counts observed in NTC well. [1000] |
| Int? | min_nuclei | Minimum # of nuclei estimated. [100] |
| Int? | min_segment_area | Minimum segment area. [5000] |
| File? | cell_type_markers_list | CSV file containing a list of major cell type markers; used for detecting genes of interest. |
| Float? | min_genes_detected_in_percent_segment | Minimum % of segments that detect the genes. [0.01] |
| Int? | n_comps | Number of principal components to compute. [30] |
| String? | batch_key | Key in AnnData object for batch information. ['batch_id'] |
| Float? | leiden_resolution | Value controlling the coarseness of the Leiden clustering. [0.4] |
| String | container_registry | Container registry where workflow Docker images are hosted. |
| String? | zones | Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f'] |

## 10x Visium inputs

An input template file can be found at [workflows/spatial_visium/inputs.json](workflows/spatial_visium/inputs.json).

| Type | Name | Description |
| :- | :- | :- |
| Array[[Project](#10x-visium-project)] | projects | The project ID, set of samples and their associated reads and metadata, output bucket locations, and whether or not to run project-level cohort analysis. |
| File | spaceranger_reference_data | Space Ranger transcriptome reference data; see https://www.10xgenomics.com/support/software/space-ranger/downloads and [10x Visium notes](#10x-visium-notes). |
| File | visium_probe_set_csv | Visium probe-based assays target genes in Space Ranger transcriptome; see https://www.10xgenomics.com/support/software/space-ranger/downloads and [10x Visium notes](#10x-visium-notes). |
| Int? | filter_cells_min_counts | Minimum number of counts required for a cell to pass filtering. [5000] |
| Int? | filter_cells_min_genes | Minimum number of genes required for a cell to pass filtering. [3000] |
| Int? | filter_genes_min_cells | Minimum number of cells expressed required for a gene to pass filtering. [10] |
| Float? | filter_mt_max_percent | Maximum percentage of mitochondrial read counts for a cell to pass filtering. [0.2] |
| Float? | normalize_target_sum | The total count to which each cell's gene expression values will be normalized. [10000] |
| Int? | n_top_genes | Number of highly-variable genes to keep. [3000] |
| Int? | n_comps | Number of principal components to compute. [30] |
| String? | batch_key | Key in AnnData object for batch information. ['batch_id'] |
| Float? | leiden_resolution | Value controlling the coarseness of the Leiden clustering. [0.4] |
| String | container_registry | Container registry where workflow Docker images are hosted. |
| String? | zones | Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f'] |

## Structs

### Nanostring GeoMx structs

#### Nanostring GeoMx Project

| Type | Name | Description |
| :- | :- | :- |
| String | team_id | Unique identifier for team; used for naming output files. |
| String | dataset_id | Unique identifier for dataset; used for naming output files. |
| String | dataset_doi_url | Generated Zenodo DOI URL referencing the dataset. |
| Array[[Slide](#nanostring-geomx-slide)] | slides | The set of slides associated with this project. |
| File | project_sample_metadata_csv | CSV containing all sample information including batch, condition, etc. |
| File | geomx_config_ini | The configuration (.ini) file, containing pipeline processing parameters that is used by the GeoMx NGS pipeline to assist in converting the FASTQ files to DCC files. It is from the GeoMx DSP readout package. Sections can include `[Sequencing]`, `[Processing_v2]`, `[AOI_List]`, and `[Targets]`; see [GeoMx configuration (.ini) files notes](#geomx-configuration-(.ini)-files). |
| Boolean | run_project_cohort_analysis | Whether or not to run cohort analysis within the project. |
| String | raw_data_bucket | Raw data bucket; intermediate output files that are not final workflow outputs are stored here. |
| String | staging_data_bucket | Staging data bucket; final project-level outputs are stored here. |

#### Nanostring GeoMx Slide

| Type | Name | Description |
| :- | :- | :- |
| String | slide_id | Unique identifier for the slide within the project; used for naming output files. |
| File | geomx_lab_annotation_xlsx | The annotation (.xlsx) file/lab worksheet, containing phenotypic data from the GeoMx DSP readout package; see [GeoMx Lab Worksheet notes](#geomx-lab-worksheet). |
| Array[[Sample](#nanostring-geomx-sample)] | samples | The set of samples associated with this project. |

#### Nanostring GeoMx Sample

| Type | Name | Description |
| :- | :- | :- |
| String | sample_id | Unique identifier for the sample within the project. |
| String? | batch | The sample's batch. |
| File | fastq_R1 | Path to the sample's read 1 FASTQ file. |
| File | fastq_R2 | Path to the sample's read 2 FASTQ file. |
| File? | fastq_I1 | Optional FASTQ index 1. |
| File? | fastq_I2 | Optional FASTQ index 2. |

### 10x Visium structs

#### 10x Visium Project

| Type | Name | Description |
| :- | :- | :- |
| String | team_id | Unique identifier for team; used for naming output files. |
| String | dataset_id | Unique identifier for dataset; used for naming output files. |
| String | dataset_doi_url | Generated Zenodo DOI URL referencing the dataset. |
| Array[[Sample](#10x-visium-sample)] | samples | The set of samples associated with this project. |
| Boolean | run_project_cohort_analysis | Whether or not to run cohort analysis within the project. |
| String | raw_data_bucket | Raw data bucket; intermediate output files that are not final workflow outputs are stored here. |
| String | staging_data_bucket | Staging data bucket; final project-level outputs are stored here. |

#### 10x Visium Sample

| Type | Name | Description |
| :- | :- | :- |
| String | sample_id | Unique identifier for the sample within the project. |
| String? | batch | The sample's batch. |
| File | fastq_R1 | Path to the sample's read 1 FASTQ file. |
| File | fastq_R2 | Path to the sample's read 2 FASTQ file. |
| File? | fastq_I1 | Optional FASTQ index 1. |
| File? | fastq_I2 | Optional FASTQ index 2. |
| File | visium_brightfield_image | The 10x Visium brightfield image, which is a high-resolution image of a tissue section and used for plotting spatial coordinates. |
| String | visium_slide_serial_number | The 10x Visium slide serial number obtained from the ASAP sample metadata. The unique identifier printed on the label of each Visium slide; see https://www.10xgenomics.com/support/software/space-ranger/3.0/getting-started/space-ranger-glossary. |
| String | visium_capture_area | The 10x Visium slide capture area obtained from the ASAP sample metadata. Active regions for capturing expression data on a Visium slide; see https://www.10xgenomics.com/support/software/space-ranger/3.0/getting-started/space-ranger-glossary. |

## Generating the inputs JSON

The inputs JSON may be generated manually, however when running a large number of samples, this can become unwieldly. The [`generate_inputs` utility script](https://github.com/ASAP-CRN/wf-common/blob/main/util/generate_inputs) may be used to automatically generate the inputs JSON (`inputs.{staging_env}.{source}-{cohort_dataset}.{date}.json`) and a sample list TSV (`{team_id}.{source}-{cohort_dataset}.sample_list.{date}.tsv`); same as the one generated in [the write_cohort_sample_list task](https://github.com/ASAP-CRN/wf-common/wdl/tasks/write_cohort_sample_list.wdl)). The script requires the libraries outlined in [the requirements.txt file](https://github.com/ASAP-CRN/wf-common/util/requirements.txt) and the following inputs:

- `project-tsv`: One or more project TSVs with one row per sample and columns team_id, ASAP_dataset_id, ASAP_sample_id, batch, fastq_R1s, fastq_R2s, fastq_I1s, fastq_I2s, embargoed, source, dataset, dataset_DOI_url, and SPATIAL columns if applicable: geomx_config, geomx_dsp_config, geomx_annotation_file, visium_cytassist, visium_probe_set, visium_slide_ref, and visium_capture_area. All samples from all projects may be included in the same project TSV, or multiple project TSVs may be provided.
	- `team_id`: A unique identifier for the team from which the sample(s) arose.
	- `ASAP_dataset_id`: A generated unique identifier for the dataset from which the sample(s) arose.
	- `ASAP_sample_id`: A generated unique identifier for the sample within the project.
	- `batch`: The sample's batch.
	- `fastq_R1s`: The gs uri to read 1 of sample FASTQ.
		- This is appended to the `project-tsv` from the `fastq-locs-txt`: FASTQ locations for all samples provided in the `project-tsv`. Each sample is expected to have one set of paired fastqs located at `${fastq_path}/${sample_id}*`. The read 1 file should include 'R1' somewhere in the filename. Generate this file e.g. by running `gcloud storage ls gs://fastq_bucket/some/path/**.fastq.gz >> fastq_locs.txt`.
	- `fastq_R2s`: The gs uri to read 2 of sample FASTQ.
		- This is appended to the `project-tsv` from the `fastq-locs-txt`: FASTQ locations for all samples provided in the `project-tsv`. Each sample is expected to have one set of paired fastqs located at `${fastq_path}/${sample_id}*`. The read 2 file should include 'R2' somewhere in the filename. Generate this file e.g. by running `gcloud storage ls gs://fastq_bucket/some/path/**.fastq.gz >> fastq_locs.txt`.
	- `fastq_I1s`: The gs uri to sample FASTQ index 1.
	- `fastq_I2s`: The gs uri to sample FASTQ index 2.
	- `embargoed`: The internal QC/embargo status of dataset.
	- `source`: The source of dataset (e.g. 'pmdbs').
	- `dataset`: The assigned dataset name without the source (e.g. 'sn-rnaseq')
	- `dataset_DOI_url`: Generated Zenodo DOI URL referencing the dataset.
	- `geomx_config`, `geomx_dsp_config`, `geomx_annotation_file`: go to [Nanostring GeoMx inputs](#nanostring-geomx-inputs)
	- `geomx_slide`: The GeoMx DSP identifier for the slide from which the sample(s) arose.
	- `ASAP_geomx_slide_id`: A unique identifier for the slide per batch/run from which the sample(s) arose.
	- `visium_cytassist`, `visium_probe_set`, `visium_slide_ref`, `visium_capture_area`: go to [10x Visium inputs](#10x-visium-inputs)
- `inputs-template`: The inputs template JSON file into which the `projects` information derived from the `project-tsv` will be inserted. Must have a key ending in `*.projects`. Other default values filled out in the inputs template will be written to the output inputs.json file.
- `run-project-cohort-analysis`: Optionally run project-level cohort analysis for provided projects. This value will apply to all projects. [false]
- `workflow_name`: WDL workflow name.
- `cohort-dataset`: Dataset name in cohort bucket name (e.g. 'sc-rnaseq').

Example usage:

```bash
./wf-common/util/generate_inputs \
	--project-tsv metadata.tsv \
	--inputs-template workflows/spatial_geomx/inputs.json \
	--run-project-cohort-analysis \
	--workflow-name spatial_geomx_analysis \
	--cohort-dataset spatial-geomx

./wf-common/util/generate_inputs \
	--project-tsv metadata.tsv \
	--inputs-template workflows/spatial_visium/inputs.json \
	--run-project-cohort-analysis \
	--workflow-name spatial_visium_analysis \
	--cohort-dataset spatial-visium
```

# Outputs

## Output structure

- `cohort_id`: either the `team_id` for project-level downstream analysis, or the `cohort_id` for the full cohort
- `workflow_run_timestamp`: format: `%Y-%m-%dT%H-%M-%SZ`
- The list of samples used to generate the cohort analysis will be output alongside other cohort analysis outputs in the staging data bucket (`${cohort_id}.sample_list.tsv`)
- The `MANIFEST.tsv` file in the staging data bucket describes the file name, md5 hash, timestamp, workflow version, workflow name, and workflow release for the run used to generate each file in that directory

### Raw data (intermediate files and final outputs for all runs of the workflow)

The raw data bucket will contain *some* artifacts generated as part of workflow execution. Following successful workflow execution, the artifacts will also be copied into the staging bucket as final outputs.

In the workflow, task outputs are either specified as `String` (final outputs, which will be copied in order to live in raw data buckets and staging buckets) or `File` (intermediate outputs that are periodically cleaned up, which will live in the cromwell-output bucket). This was implemented to reduce storage costs.

```bash
asap-raw-{cohort,team-xxyy}-{source}-{dataset}
└── workflow_execution
    └── spatial_geomx
        ├── cohort_analysis
        │   └──${cohort_analysis_workflow_version}
        │      └── ${workflow_run_timestamp}
        │          └── <cohort_analysis outputs>
        ├── process_to_adata
        │   └──${process_to_adata_workflow_version}
        │      └── ${workflow_run_timestamp}
        │          └── <process_to_adata outputs>
        └── preprocess
            ├── fastq_to_dcc
            │   └── ${fastq_to_dcc_task_version}
            │       └── <fastq_to_dcc output>
            ├── dcc_to_rds
            │   └── ${dcc_to_rds_task_version}
            │       └── <dcc_to_rds output>
            └── qc
                └── ${qc_task_version}
                    └── <qc output>

asap-raw-{cohort,team-xxyy}-{source}-{dataset}
└── workflow_execution
    └── spatial_visium
        ├── cohort_analysis
        │   └──${cohort_analysis_workflow_version}
        │      └── ${workflow_run_timestamp}
        │          └── <cohort_analysis outputs>
        └── preprocess
            ├── spaceranger_count
            │   └── ${spaceranger_count_task_version}
            │       └── <spaceranger_count output>
            ├── counts_to_adata
            │   └── ${counts_to_adata_task_version}
            │       └── <counts_to_adata output>
            └── qc
                └── ${qc_task_version}
                    └── <qc output>
```

### Staging data (intermediate workflow objects and final workflow outputs for the latest run of the workflow)

Following QC by researchers, the objects in the dev or uat bucket are synced into the curated data buckets, maintaining the same file structure. Curated data buckets are named `asap-curated-{team-xxyy}-{source}-{dataset}`.

Data may be synced using [the `promote_staging_data` script](#promoting-staging-data).

```bash
asap-dev-{team-xxyy}-{source}-{dataset}
└── spatial_geomx
    ├── cohort_analysis
    │   ├── ${team_id}.sample_list.tsv
    │   ├── ${team_id}.merged_metadata.csv
    │   ├── ${team_id}.merged_processed.h5ad
    │   ├── ${team_id}.all_genes.csv
    │   ├── ${team_id}.hvg_genes.csv
    │   ├── ${team_id}.hvg_dispersion.png
    │   ├── ${team_id}.umap_cluster.png
    │   ├── ${team_id}.final.h5ad
    │   ├── ${team_id}.final_metadata.csv
    │   └── MANIFEST.tsv
    ├── process_to_adata
    │   ├── ${slideN_id}.segment_gene_detection_plot.png
    │   ├── ${slideN_id}.gene_detection_rate.csv
    │   ├── ${slideN_id}.q3_negprobe_plot.png
    │   ├── ${slideN_id}.normalization_plot.png
    │   └── MANIFEST.tsv
    └── preprocess
        ├── ${slideA_id}.DCC.zip
        ├── ${slideA_id}.geomxngs_out_dir.tar.gz
        ├── ${slideA_id}.NanoStringGeoMxSet.rds
        ├── ${slideA_id}.qc.rds
        ├── ${slideA_id}.segment_qc_summary.csv
        ├── ${slideA_id}.probe_qc_summary.csv
        ├── ${slideA_id}.gene_count.csv
        ├── MANIFEST.tsv
        ├── ...
        ├── ${slideN_id}.DCC.zip
        ├── ${slideN_id}.geomxngs_out_dir.tar.gz
        ├── ${slideN_id}.NanoStringGeoMxSet.rds
        ├── ${slideN_id}.qc.rds
        ├── ${slideN_id}.segment_qc_summary.csv
        ├── ${slideN_id}.probe_qc_summary.csv
        ├── ${slideN_id}.gene_count.csv
        └── MANIFEST.tsv

asap-dev-{team-xxyy}-{source}-{dataset}
└── spatial_visium
    ├── cohort_analysis
    │   ├── ${team_id}.sample_list.tsv
    │   ├── ${team_id}.merged_cleaned_unfiltered.h5ad
    │   ├── ${team_id}.merged_metadata.csv
    │   ├── ${team_id}.all_genes.csv
    │   ├── ${team_id}.hvg_genes.csv
    │   ├── ${team_id}.qc_violin.png
    │   ├── ${team_id}.qc_dist.png
    │   ├── ${team_id}.hvg_dispersion.png
    │   ├── ${team_id}.umap_cluster.png
    │   ├── ${team_id}.spatial_scatter.png
    │   ├── ${team_id}.final.h5ad
    │   ├── ${team_id}.final_metadata.csv
    │   ├── ${team_id}.moran_top_10_variable_genes.csv
    │   ├── ${team_id}.moran_top_4_variable_genes_spatial_scatter.png
    │   └── MANIFEST.tsv
    └── preprocess
        ├── ${sampleA_id}.raw_feature_bc_matrix.h5
        ├── ${sampleA_id}.filtered_feature_bc_matrix.h5
        ├── ${sampleA_id}.molecule_info.h5
        ├── ${sampleA_id}.metrics_summary.csv
        ├── ${sampleA_id}.spaceranger_spatial_outputs.tar.gz
        ├── ${sampleA_id}.aligned_fiducials.jpg
        ├── ${sampleA_id}.detected_tissue_image.jpg
        ├── ${sampleA_id}.tissue_hires_image.png
        ├── ${sampleA_id}.tissue_lowres_image.png
        ├── ${sampleA_id}.scalefactors_json.json
        ├── ${sampleA_id}.tissue_positions.csv
        ├── ${sampleA_id}.spatial_enrichment.csv
        ├── ${sampleA_id}.cleaned_unfiltered.h5ad
        ├── ${sampleA_id}.qc.h5ad
        ├── MANIFEST.tsv
        ├── ...
        ├── ${sampleN_id}.raw_feature_bc_matrix.h5
        ├── ${sampleN_id}.filtered_feature_bc_matrix.h5
        ├── ${sampleN_id}.molecule_info.h5
        ├── ${sampleN_id}.metrics_summary.csv
        ├── ${sampleN_id}.spaceranger_spatial_outputs.tar.gz
        ├── ${sampleN_id}.aligned_fiducials.jpg
        ├── ${sampleN_id}.detected_tissue_image.jpg
        ├── ${sampleN_id}.tissue_hires_image.png
        ├── ${sampleN_id}.tissue_lowres_image.png
        ├── ${sampleN_id}.scalefactors_json.json
        ├── ${sampleN_id}.tissue_positions.csv
        ├── ${sampleN_id}.spatial_enrichment.csv
        ├── ${sampleN_id}.cleaned_unfiltered.h5ad
        ├── ${sampleN_id}.qc.h5ad
        └── MANIFEST.tsv
```

## Promoting staging data

The [`promote_staging_data` script](https://github.com/ASAP-CRN/wf-common/util/promote_staging_data) can be used to promote staging data that has been approved to the curated data bucket for a team or set of teams.

This script compiles bucket and file information for both the initial (staging) and target (prod) environment. It also runs data integrity tests to ensure staging data can be promoted and generates a Markdown report. It (1) checks that files are not empty and are not less than or equal to 10 bytes (factoring in white space) and (2) checks that files have associated metadata and is present in MANIFEST.tsv.

If data integrity tests pass, this script will upload a combined MANIFEST.tsv and the data promotion Markdown report under a metadata/{timestamp} directory in the staging bucket. Previous manifest files and reports will be kept. Next, it will rsync all files in the staging bucket to the curated bucket's upstream, downstream, cohort_analysis, and metadata directories. **Exercise caution when using this script**; files that are not present in the source (staging) bucket will be deleted at the destination (curated) bucket.

If data integrity tests fail, staging data cannot be promoted. The combined `MANIFEST.tsv`, Markdown report, and `promote_staging_data_script.log` will be locally available.

The script defaults to a dry run, printing out the files that would be copied or deleted for each selected team.

### Options

```
-h  Display this message and exit
-t  Space-delimited team(s) to promote data for
-l  List available teams
-s  Source name in bucket name
-d  Space-delimited dataset name(s) in team bucket name, must follow the same order as {team}
-w  Workflow name used as a directory in bucket
-p  Promote data. If this option is not selected, data that would be copied or deleted is printed out, but files are not actually changed (dry run)
```

### Usage

```bash
# List available teams
./wf-common/util/promote_staging_data -t cohort -l -s pmdbs -d spatial-geomx -w spatial_geomx
./wf-common/util/promote_staging_data -t cohort -l -s mouse -d spatial-visium -w spatial_visium

# Print out the files that would be copied or deleted from the staging bucket to the curated bucket for teams team-edwards and team-vila
./wf-common/util/promote_staging_data -t team-edwards team-vila -s pmdbs -d spatial-geomx-th spatial-geomx-thlc -w spatial_geomx

# Promote data for team-edwards and team-vila
./wf-common/util/promote_staging_data -t team-edwards team-vila -s pmdbs -d spatial-geomx-th spatial-geomx-thlc -w spatial_geomx -p
```

# Docker images

Docker images are defined in [the `docker` directory](docker). Each image must minimally define a `build.env` file and a `Dockerfile`.

Example directory structure:
```bash
docker
├── geomxngs
│   ├── build.env
│   └── Dockerfile
└── spatial_py
    ├── build.env
    ├── Dockerfile
    ├── requirements.txt
    └── scripts
        ├── visium_counts_to_adata.py
        ├── visium_qc.py
        ├── visium_merge_and_plot_qc.py
        ├── visium_process.py
        ├── integrate_harmony.py
        ├── cluster.py
        ├── visium_plot_spatial.py
        ├── visium_spatially_variable_genes.py
        └── ...
```

## The `build.env` file

Each target image is defined using the `build.env` file, which is used to specify the name and version tag for the corresponding Docker image. It must contain at minimum the following variables:

- `IMAGE_NAME`
- `IMAGE_TAG`

All variables defined in the `build.env` file will be made available as build arguments during Docker image build.

The `DOCKERFILE` variable may be used to specify the path to a Dockerfile if that file is not found alongside the `build.env` file, for example when multiple images use the same base Dockerfile definition.

## Building Docker images

Docker images can be build using the [`build_docker_images`](https://github.com/DNAstack/bioinformatics-scripts/blob/main/scripts/build_docker_images) script.

```bash
# Build a single image
./build_docker_images -d docker/geomxngs

# Build all images in the `docker` directory
./build_docker_images -d docker

# Build and push all images in the docker directory, using the `dnastack` container registry
./build_docker_images -d docker -c dnastack -p
```

## Tool and library versions

| Image | Major tool versions | Links | Workflow |
| :- | :- | :- | :- |
| geomxngs | <ul><li>[geomxngs v3.1.1.6](https://nanostring.app.box.com/v/GeoMxSW3-1-0/folder/233772026049)</li></ul> | [Dockerfile](https://github.com/ASAP-CRN/spatial-transcriptomics-wf/tree/main/docker/geomxngs) | spatial_geomx |
| spaceranger | <ul><li>[spaceranger v3.1.2](https://www.10xgenomics.com/support/software/space-ranger/latest/release-notes/release-notes-for-SR#v-3-1-2)</li></ul> | [Dockerfile](https://github.com/ASAP-CRN/spatial-transcriptomics-wf/tree/main/docker/spaceranger) | spatial_visium |
| util | <ul><li>[google-cloud-cli 524.0.0-slim](https://cloud.google.com/sdk/docs/release-notes#52400_2025-05-28)</li></ul> | [Dockerfile](https://github.com/ASAP-CRN/wf-common/tree/main/docker/util) | both |
| spatial_r | R (v4.4.2) packages: <ul><li>[NanoStringNCTools v1.14.0](https://github.com/Nanostring-Biostats/NanoStringNCTools)</li><li>[Seurat v5.2.1](https://github.com/satijalab/seurat/releases/tag/v5.2.1)</li><li>[SeuratData v0.2.2.9001](https://github.com/satijalab/seurat-data/tags)</li><li>[SeuratDisk v0.0.0.9021](https://github.com/mojaveazure/seurat-disk)</li></ul> | [Dockerfile](https://github.com/ASAP-CRN/spatial-transcriptomics-wf/tree/main/docker/spatial_r) | spatial_geomx |
| spatial_py | Python (v3.12.5) libraries: <ul><li>[squidpy v1.6.2](https://github.com/scverse/squidpy/releases/tag/v1.6.2)</li><li>[matplotlib v3.10.0](https://github.com/matplotlib/matplotlib/releases/tag/v3.10.0)</li><li>[seaborn v0.13.2](https://github.com/mwaskom/seaborn/releases/tag/v0.13.2)</li><li>[harmonypy v0.0.10](https://github.com/slowkow/harmonypy/releases/tag/v0.0.10)</li><li>[scanpy v1.10.4](https://github.com/scverse/scanpy/releases/tag/1.10.4)</li></ul> | [Dockerfile](https://github.com/ASAP-CRN/spatial-transcriptomics-wf/tree/main/docker/spatial_py) | both |


# wdl-ci

[`wdl-ci`](https://github.com/DNAstack/wdl-ci) provides tools to validate and test workflows and tasks written in [Workflow Description Language (WDL)](https://github.com/openwdl/wdl). In addition to the tests packaged in `wdl-ci`, the [spatial-transcriptomics-wdl-ci-custom-test-dir](./spatial-transcriptomics-wdl-ci-custom-test-dir) is a directory containing custom WDL-based tests that are used to test workflow tasks. `wdl-ci` in this repository is set up to run on pull request.

In general, `wdl-ci` will use inputs provided in the [wdl-ci.config.json](./wdl-ci.config.json) and compare current outputs and validated outputs based on changed tasks/workflows to ensure outputs are still valid by meeting the critera in the specified tests. For example, if the rds to adata task in our workflow was changed, then this task would be submitted and that output would be considered the "current output". When inspecting the converted adata object, there is a test specified in the [wdl-ci.config.json](./wdl-ci.config.json) called, "check_hdf5". The test will compare the "current output" and "validated output" (provided in the [wdl-ci.config.json](./wdl-ci.config.json)) to make sure that the .h5ad file is still a valid HDF5 file.


# Notes

## Nanostring GeoMx notes

The [GeoMx DSP Instrument User Manual](https://nanostring.com/wp-content/uploads/2022/06/MAN-10152-01-GeoMx-DSP-Instrument-User-Manual.pdf) and [GeoMx DSP Data Analysis User Manual](https://nanostring.com/wp-content/uploads/2022/06/MAN-10154-01-GeoMx-DSP-Data-Analysis-User-Manual.pdf) provide a comprehensive overview of the GeoMx DSP workflow- from sample preparation and slide scanning to sequencing and data analysis. The required input files listed below for our `spatial_geomx` pipeline are described in detail in these manuals. These files are included as part of the GeoMx Readout Package, which are used during the experimental runs.

### GeoMx Lab Worksheet

The GeoMx Lab Worksheet contains information on the contents of each well of each plate and can be useful when doing library preparation and pooling. Some columns include:
- `Sample_ID`
- `slide name`
- `scan name`
- `panel`
- `roi`
- `segment`
- `aoi`
- `area`
- ...

### GeoMx configuration (.ini) files

An example of the configuration (.ini) file can be found in the [GeoMx-NGS-Pipeline-Dataset](https://nanostring.app.box.com/v/GeoMxSW3-1-0/file/1385968928681).

| Section         | Description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| `[Sequencing]`   | Specifies sequencing platform, read pattern, library prep, etc.            |
| `[Processing_v2]` | Core processing parameters including adapters and filters.                |
| `[AOI_List]`      | Area of Interest, which refers to the specific regions selected on a GeoMx slide for spatial transcriptomic profiling. |
| `[Targets]`       | Maps readout tag sequence identifiers to their corresponding nucleotide sequences. |

### GeoMx DSP configuration (.pkc) files

The Nanostring GeoMx configuration (.pkc) files were obtained from https://nanostring.com/products/geomx-digital-spatial-profiler/geomx-dsp-configuration-files/.
- [Human_WTA_v1.0](https://nanostring.com/wp-content/uploads/Hs_R_NGS_WTA_v1.0.pkc_.zip) for Human Whole Transcriptome Atlas
- [Mouse_WTA_v1.0](https://nanostring.com/wp-content/uploads/Mm_R_NGS_WTA_v1.0.zip) or [Mouse_WTA_v2.0](https://nanostring.com/wp-content/uploads/2024/06/Mm_R_NGS_WTA_v2.0.zip) for Mouse Whole Transcriptome Atlas

## 10x Visium notes

The Space Ranger reference data were obtained from https://www.10xgenomics.com/support/software/space-ranger/downloads.
- [Human reference (GRCh38)](https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz)
- [Mouse reference (mm10)](https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-mm10-2020-A.tar.gz)

The Space Ranger probe set data were obtained from https://www.10xgenomics.com/support/software/space-ranger/downloads.
- [Human Transcriptome v2](https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv)
- [Mouse Transcriptome v2](https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv)
