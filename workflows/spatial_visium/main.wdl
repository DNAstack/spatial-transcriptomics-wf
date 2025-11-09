version 1.0

# Harmonized human PMDBS and non-human spatial transcriptomics workflow entrypoint for 10x Visium data

import "structs.wdl"
import "../../wf-common/wdl/tasks/get_workflow_metadata.wdl" as GetWorkflowMetadata
import "preprocess/preprocess.wdl" as Preprocess
import "cohort_analysis/cohort_analysis.wdl" as CohortAnalysis

workflow spatial_visium_analysis {
	input {
		Array[Project] projects

		File spaceranger_reference_data
		File? visium_probe_set_csv

		# Processing parameters
		Int filter_cells_min_counts = 5000
		Int filter_cells_min_genes = 3000
		Int filter_genes_min_cells = 10
		Float filter_mt_max_percent = 0.2
		Float normalize_target_sum = 10000
		Int n_top_genes = 3000
		Int n_comps = 30
		String batch_key = "batch_id"
		Float leiden_resolution = 0.4

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

	String workflow_execution_path = "workflow_execution"
	String workflow_name = "spatial_visium"
	String workflow_version = "v1.0.1"
	String workflow_release = "https://github.com/ASAP-CRN/spatial-transcriptomics-wf/releases/tag/spatial_visium_analysis-~{workflow_version}"

	call GetWorkflowMetadata.get_workflow_metadata {
		input:
			zones = zones
	}

	scatter (project in projects) {
		String project_raw_data_path_prefix = "~{project.raw_data_bucket}/~{workflow_execution_path}/~{workflow_name}"

		call Preprocess.preprocess {
			input:
				team_id = project.team_id,
				dataset_id = project.dataset_id,
				dataset_doi_url = project.dataset_doi_url,
				samples = project.samples,
				spaceranger_reference_data = spaceranger_reference_data,
				visium_probe_set_csv = visium_probe_set_csv,
				workflow_name = workflow_name,
				workflow_version = workflow_version,
				workflow_release = workflow_release,
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = project_raw_data_path_prefix,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}

		Array[String] preprocessing_output_file_paths = flatten([
			preprocess.raw_counts,
			preprocess.filtered_counts,
			preprocess.molecule_info,
			preprocess.metrics_summary_csv,
			preprocess.spatial_outputs_tar_gz,
			flatten(preprocess.spatial_images),
			preprocess.scalefactors_json,
			preprocess.tissue_positions_csv,
			preprocess.spatial_enrichment_csv,
			preprocess.initial_adata_object,
			preprocess.qc_adata_object
		]) #!StringCoercion

		if (project.run_project_cohort_analysis) {
			call CohortAnalysis.cohort_analysis as project_cohort_analysis {
				input:
					cohort_id = project.team_id,
					project_sample_ids = preprocess.project_sample_ids,
					preprocessed_adata_objects = preprocess.qc_adata_object,
					preprocessing_output_file_paths = preprocessing_output_file_paths,
					filter_cells_min_counts = filter_cells_min_counts,
					filter_cells_min_genes = filter_cells_min_genes,
					filter_genes_min_cells = filter_genes_min_cells,
					filter_mt_max_percent = filter_mt_max_percent,
					normalize_target_sum = normalize_target_sum,
					n_top_genes = n_top_genes,
					n_comps = n_comps,
					batch_key = batch_key,
					leiden_resolution = leiden_resolution,
					workflow_name = workflow_name,
					workflow_version = workflow_version,
					workflow_release = workflow_release,
					run_timestamp = get_workflow_metadata.timestamp,
					raw_data_path_prefix = project_raw_data_path_prefix,
					staging_data_buckets = project.staging_data_buckets,
					billing_project = get_workflow_metadata.billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}
	}

	output {
		# Sample-level outputs
		## Sample list
		Array[Array[Array[String]]] project_sample_ids = preprocess.project_sample_ids

		## Preprocess
		Array[Array[File]] raw_counts = preprocess.raw_counts
		Array[Array[File]] filtered_counts = preprocess.filtered_counts
		Array[Array[File]] molecule_info = preprocess.molecule_info
		Array[Array[File]] metrics_summary_csv = preprocess.metrics_summary_csv
		Array[Array[File]] spatial_outputs_tar_gz = preprocess.spatial_outputs_tar_gz
		Array[Array[Array[File]]] spatial_images = preprocess.spatial_images
		Array[Array[File]] scalefactors_json = preprocess.scalefactors_json
		Array[Array[File]] tissue_positions_csv = preprocess.tissue_positions_csv
		Array[Array[File]] spatial_enrichment_csv = preprocess.spatial_enrichment_csv
		Array[Array[File]] initial_adata_object = preprocess.initial_adata_object
		Array[Array[File]] qc_adata_object = preprocess.qc_adata_object

		# Project cohort analysis outputs
		## List of samples included in the cohort
		Array[File?] project_cohort_sample_list = project_cohort_analysis.cohort_sample_list

		# Merged, processed (filtered, normalized, dimensionality reduced), integrated, and clustered adata objects, and plots
		Array[File?] project_merged_adata_object = project_cohort_analysis.merged_adata_object
		Array[File?] project_merged_metadata_csv = project_cohort_analysis.merged_metadata_csv
		Array[File?] project_all_genes_csv = project_cohort_analysis.all_genes_csv
		Array[File?] project_hvg_genes_csv = project_cohort_analysis.hvg_genes_csv
		Array[Array[File]?] project_qc_plots_png = project_cohort_analysis.qc_plots_png
		Array[File?] project_processed_adata_object = project_cohort_analysis.processed_adata_object
		Array[File?] project_hvg_plot_png = project_cohort_analysis.hvg_plot_png
		Array[File?] project_integrated_adata_object = project_cohort_analysis.integrated_adata_object
		Array[File?] project_clustered_adata_object = project_cohort_analysis.clustered_adata_object
		Array[File?] project_umap_cluster_plots_png = project_cohort_analysis.umap_cluster_plots_png

		# Image features outputs
		Array[File?] project_spatial_scatter_plot_png = project_cohort_analysis.spatial_scatter_plot_png

		# Spatial statistics outputs
		Array[File?] project_final_adata_object = project_cohort_analysis.final_adata_object
		Array[File?] project_final_metadata_csv = project_cohort_analysis.final_metadata_csv
		Array[File?] project_moran_top_10_variable_genes_csv = project_cohort_analysis.moran_top_10_variable_genes_csv
		Array[File?] project_moran_top_4_variable_genes_spatial_scatter_plot_png = project_cohort_analysis.moran_top_4_variable_genes_spatial_scatter_plot_png

		Array[Array[File]?] preprocess_manifests = project_cohort_analysis.preprocess_manifest_tsvs
		Array[Array[File]?] project_manifests = project_cohort_analysis.cohort_analysis_manifest_tsvs
	}

	meta {
		description: "Harmonized human postmortem-derived brain sequencing (PMDBS) and non-human spatial transcriptomics workflow for 10x Visium data"
	}

	parameter_meta {
		projects: {help: "The project ID, set of samples and their associated reads and metadata, output bucket locations, and whether or not to run project-level downstream analysis."}
		spaceranger_reference_data: {help: "Space Ranger transcriptome reference data; see https://www.10xgenomics.com/support/software/space-ranger/downloads."}
		visium_probe_set_csv: {help: "Visium probe-based assays target genes in Space Ranger transcriptome; see https://www.10xgenomics.com/support/software/space-ranger/downloads."}
		filter_cells_min_counts: {help: "Minimum number of counts required for a cell to pass filtering. [5000]"}
		filter_cells_min_genes: {help: "Minimum number of genes required for a cell to pass filtering. [3000]"}
		filter_genes_min_cells: {help: "Minimum number of cells expressed required for a gene to pass filtering. [10]"}
		filter_mt_max_percent: {help: "Maximum percentage of mitochondrial read counts for a cell to pass filtering. [0.2]"}
		normalize_target_sum: {help: "The total count to which each cell's gene expression values will be normalized. [10000]"}
		n_top_genes: {help: "Number of highly-variable genes to keep. [3000]"}
		n_comps: {help: "Number of principal components to compute. [30]"}
		batch_key: {help: "Key in AnnData object for batch information. ['batch_id']"}
		leiden_resolution: {help: "Value controlling the coarseness of the Leiden clustering. [0.4]"}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}
