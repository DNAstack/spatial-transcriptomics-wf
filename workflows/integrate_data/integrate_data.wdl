version 1.0

# Perform dataset integration, UMAP, and clustering steps

workflow integrate_data {
	input {
		String cohort_id
		File merged_and_processed_adata_object

		Int n_comps
		String batch_key
		Float leiden_resolution

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	call integrate_sample_data {
		input:
			cohort_id = cohort_id,
			merged_and_processed_adata_object = merged_and_processed_adata_object,
			batch_key = batch_key,
			container_registry = container_registry,
			zones = zones
	}

	call cluster {
		input:
			cohort_id = cohort_id,
			integrated_adata_object = integrate_sample_data.integrated_adata_object,
			n_comps = n_comps,
			leiden_resolution = leiden_resolution,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		File integrated_adata_object = integrate_sample_data.integrated_adata_object
		File clustered_adata_object = cluster.clustered_adata_object
		File umap_cluster_plots_png = cluster.umap_cluster_plots_png #!FileCoercion
	}

	meta {
		description: "Perform sample integration with Harmony and cluster the 10x Visium data based on transcriptonal similarity."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		merged_and_processed_adata_object: {help: "Merged and processed AnnData object to run integrate data workflow on."}
		n_comps: {help: "Number of principal components to compute. [30]"}
		batch_key: {help: "Key in AnnData object for batch information. ['batch_id']"}
		leiden_resolution: {help: "Value controlling the coarseness of the Leiden clustering. [0.4]"}
		raw_data_path: {help: "Raw data bucket path for integrate data workflow outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}

task integrate_sample_data {
	input {
		String cohort_id
		File merged_and_processed_adata_object

		String batch_key

		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(merged_and_processed_adata_object, "GB") * 5 + 20)
	Int disk_size = ceil(size(merged_and_processed_adata_object, "GB") * 3 + 50)

	command <<<
		set -euo pipefail

		integrate_harmony \
			--adata-input ~{merged_and_processed_adata_object} \
			--batch-key ~{batch_key} \
			--adata-output ~{cohort_id}.harmony_integrated.h5ad
	>>>

	output {
		File integrated_adata_object = "~{cohort_id}.harmony_integrated.h5ad"
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.1"
		cpu: 4
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 10
		zones: zones
	}

	meta {
		description: "Perform batch correction on the processed AnnData object using the Harmony algorithm."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		merged_and_processed_adata_object: {help: "Merged and processed AnnData object to run integrate data workflow on."}
		batch_key: {help: "Key in AnnData object for batch information. ['batch_id']"}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}

task cluster {
	input {
		String cohort_id
		File integrated_adata_object

		Int n_comps
		Float leiden_resolution

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(integrated_adata_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(integrated_adata_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		cluster \
			--adata-input ~{integrated_adata_object} \
			--n-comps ~{n_comps} \
			--resolution ~{leiden_resolution} \
			--plots-prefix ~{cohort_id} \
			--adata-output ~{cohort_id}.clustered.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.umap_cluster.png"
	>>>

	output {
		File clustered_adata_object = "~{cohort_id}.clustered.h5ad"
		String umap_cluster_plots_png = "~{raw_data_path}/~{cohort_id}.umap_cluster.png"
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.1"
		cpu: 4
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 10
		zones: zones
	}

	meta {
		description: "Compute the nearest neighbors distance matrix and UMAP embedding based on the PCA-reduced spatial data, and perform Leiden clustering."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		integrated_adata_object: {help: "Integrated AnnData object."}
		n_comps: {help: "Number of principal components to compute. [30]"}
		leiden_resolution: {help: "Value controlling the coarseness of the Leiden clustering. [0.4]"}
		raw_data_path: {help: "Raw data bucket path for clustered outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}
