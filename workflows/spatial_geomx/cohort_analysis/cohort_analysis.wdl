version 1.0

# Merge processed adata object, integrate and cluster

import "../../../wf-common/wdl/tasks/write_cohort_sample_list.wdl" as WriteCohortSampleList
import "../../integrate_data/integrate_data.wdl" as IntegrateData
import "../../../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs

workflow cohort_analysis {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids
		Array[File] processed_adata_objects

		# If provided, these files will be uploaded to the staging bucket alongside other intermediate files made by this workflow
		Array[String] preprocessing_output_file_paths = []
		Array[String] processing_output_file_paths = []

		# Integrate and cluster parameters
		Int n_top_genes
		Int n_comps
		String batch_key
		Float leiden_resolution

		String workflow_name
		String workflow_version
		String workflow_release
		String crn_release_version
		String run_timestamp
		String raw_data_path_prefix
		Array[String] staging_data_buckets
		String billing_project
		String container_registry
		String zones
	}

	String sub_workflow_name = "cohort_analysis"
	String sub_workflow_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String raw_data_path = "~{raw_data_path_prefix}/~{sub_workflow_name}/~{sub_workflow_version}/~{run_timestamp}"

	call WriteCohortSampleList.write_cohort_sample_list {
		input:
			cohort_id = cohort_id,
			project_sample_ids = project_sample_ids,
			billing_project = billing_project,
			workflow_info = workflow_info,
			raw_data_path = raw_data_path,
			container_registry = container_registry,
			zones = zones
	}

	call merge_and_prep {
		input:
			cohort_id = cohort_id,
			processed_adata_objects = processed_adata_objects,
			n_top_genes = n_top_genes,
			n_comps = n_comps,
			batch_key = batch_key,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call IntegrateData.integrate_data {
		input:
			cohort_id = cohort_id,
			merged_and_processed_adata_object = merge_and_prep.merged_and_processed_adata_object, #!FileCoercion
			n_comps = n_comps,
			batch_key = batch_key,
			leiden_resolution = leiden_resolution,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call export_final_artifacts {
		input:
			cohort_id = cohort_id,
			clustered_adata_object = integrate_data.clustered_adata_object,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call UploadFinalOutputs.upload_final_outputs as upload_preprocess_files {
		input:
			output_file_paths = preprocessing_output_file_paths,
			staging_data_buckets = staging_data_buckets,
			staging_data_path = "~{workflow_name}/release/~{crn_release_version}/preprocess",
			billing_project = billing_project,
			zones = zones
	}

	call UploadFinalOutputs.upload_final_outputs as upload_process_to_adata_files {
		input:
			output_file_paths = processing_output_file_paths,
			staging_data_buckets = staging_data_buckets,
			staging_data_path = "~{workflow_name}/release/~{crn_release_version}/process_to_adata",
			billing_project = billing_project,
			zones = zones
	}

	Array[String] cohort_analysis_final_output_paths = flatten([
		[
			write_cohort_sample_list.cohort_sample_list
		],
		[
			merge_and_prep.merged_metadata_csv,
			merge_and_prep.merged_and_processed_adata_object,
			merge_and_prep.all_genes_csv,
			merge_and_prep.hvg_genes_csv,
			merge_and_prep.hvg_plot_png
		],
		[
			integrate_data.umap_cluster_plots_png
		],
		[
			export_final_artifacts.final_adata_object,
			export_final_artifacts.final_metadata_csv
		]
	]) #!StringCoercion

	call UploadFinalOutputs.upload_final_outputs as upload_cohort_analysis_files {
		input:
			output_file_paths = cohort_analysis_final_output_paths,
			staging_data_buckets = staging_data_buckets,
			staging_data_path = "~{workflow_name}/release/~{crn_release_version}/~{sub_workflow_name}",
			billing_project = billing_project,
			zones = zones
	}

	output {
		File cohort_sample_list = write_cohort_sample_list.cohort_sample_list #!FileCoercion

		# Merged and prepped AnnData object
		File merged_metadata_csv = merge_and_prep.merged_metadata_csv #!FileCoercion
		File merged_and_processed_adata_object = merge_and_prep.merged_and_processed_adata_object #!FileCoercion
		File all_genes_csv = merge_and_prep.all_genes_csv #!FileCoercion
		File hvg_genes_csv = merge_and_prep.hvg_genes_csv #!FileCoercion
		File hvg_plot_png = merge_and_prep.hvg_plot_png #!FileCoercion

		# Integrate data outputs
		File integrated_adata_object = integrate_data.integrated_adata_object
		File clustered_adata_object = integrate_data.clustered_adata_object
		File umap_cluster_plots_png = integrate_data.umap_cluster_plots_png

		# Export final outputs
		File final_adata_object = export_final_artifacts.final_adata_object #!FileCoercion
		File final_metadata_csv = export_final_artifacts.final_metadata_csv #!FileCoercion

		Array[File] preprocess_manifest_tsvs = upload_preprocess_files.manifests #!FileCoercion
		Array[File] process_to_adata_manifest_tsvs = upload_process_to_adata_files.manifests #!FileCoercion
		Array[File] cohort_analysis_manifest_tsvs = upload_cohort_analysis_files.manifests #!FileCoercion
	}

	meta {
		description: "Run team-level and/or cross-team cohort analysis on the Nanostring GeoMx data by dimensionality reduction, sample integration, and clustering."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		project_sample_ids: {help: "Associated team ID and sample ID; used to generate a sample list."}
		processed_adata_objects: {help: "An array of processed adata objects to run cohort analysis on."}
		preprocessing_output_file_paths: {help: "Selected preprocessed output files to upload to the staging bucket alongside selected cohort analysis output files."}
		processing_output_file_paths: {help: "Selected processed output files to upload to the staging bucket alongside selected cohort analysis output files."}
		n_top_genes: {help: "Number of highly-variable genes to keep. [3000]"}
		n_comps: {help: "Number of principal components to compute. [30]"}
		batch_key: {help: "Key in AnnData object for batch information. ['batch_id']"}
		leiden_resolution: {help: "Value controlling the coarseness of the Leiden clustering. [0.4]"}
		workflow_name: {help: "Workflow name; stored in the file-level manifest and final manifest with all saved files."}
		workflow_version: {help: "Workflow version; stored in the file-level manifest and final manifest with all saved files."}
		workflow_release: {help: "GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		crn_release_version: {help: "CRN Cloud release version; used to organize outputs and for the CRN Cloud release."}
		run_timestamp: {help: "UTC timestamp; stored in the file-level manifest and final manifest with all saved files."}
		raw_data_path_prefix: {help: "Raw data bucket path prefix; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis`)."}
		staging_data_buckets: {help: "Array of staging data buckets to upload intermediate files to (i.e., DEV or UAT buckets depending on internal QC status)."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}

task merge_and_prep {
	input {
		String cohort_id
		Array[File] processed_adata_objects

		Int n_top_genes
		Int n_comps
		String batch_key

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(processed_adata_objects, "GB") * 2 + 20)
	Int disk_size = ceil(size(processed_adata_objects, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		geomx_merge_and_prep \
			--adata-paths-input ~{sep=' ' processed_adata_objects} \
			--n-top-genes ~{n_top_genes} \
			--n-comps ~{n_comps} \
			--batch-key ~{batch_key} \
			--output-prefix ~{cohort_id} \
			--adata-output ~{cohort_id}.merged_processed.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.merged_metadata.csv" \
			-o "~{cohort_id}.merged_processed.h5ad" \
			-o "~{cohort_id}.all_genes.csv" \
			-o "~{cohort_id}.hvg_genes.csv" \
			-o "~{cohort_id}.hvg_dispersion.png"

	>>>

	output {
		String merged_metadata_csv = "~{raw_data_path}/~{cohort_id}.merged_metadata.csv"
		String merged_and_processed_adata_object = "~{raw_data_path}/~{cohort_id}.merged_processed.h5ad"
		String all_genes_csv = "~{raw_data_path}/~{cohort_id}.all_genes.csv"
		String hvg_genes_csv = "~{raw_data_path}/~{cohort_id}.hvg_genes.csv"
		String hvg_plot_png = "~{raw_data_path}/~{cohort_id}.hvg_dispersion.png"
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.0"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 15
		zones: zones
	}

	meta {
		description: "Merge sample-level AnnData objects to a single cohort-level AnnData object, and prepare for downstream analysis by annotating highly-variable genes (HVG) and performing PCA."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		processed_adata_objects: {help: "An array of processed AnnData object to merge."}
		n_top_genes: {help: "Number of highly-variable genes to keep. [3000]"}
		n_comps: {help: "Number of principal components to compute. [30]"}
		batch_key: {help: "Key in AnnData object for batch information. ['batch_id']"}
		raw_data_path: {help: "Raw data bucket path for merged adata and HVG plot outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}

task export_final_artifacts {
	input {
		String cohort_id
		File clustered_adata_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(clustered_adata_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(clustered_adata_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		geomx_export_final_artifacts \
			--cohort-id ~{cohort_id} \
			--adata-input ~{clustered_adata_object} \
			--adata-output ~{cohort_id}.final.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.final.h5ad" \
			-o "~{cohort_id}.final_metadata.csv"

	>>>

	output {
		String final_adata_object = "~{raw_data_path}/~{cohort_id}.final.h5ad"
		String final_metadata_csv = "~{raw_data_path}/~{cohort_id}.final_metadata.csv"
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.0"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 15
		zones: zones
	}

	meta {
		description: "Export clustered AnnData object as final and grab the metadata."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		clustered_adata_object: {help: "Clustered AnnData object."}
		raw_data_path: {help: "Raw data bucket path for merged adata and HVG plot outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}
