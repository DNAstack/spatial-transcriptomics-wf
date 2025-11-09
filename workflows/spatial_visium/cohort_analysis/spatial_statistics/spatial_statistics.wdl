version 1.0

# Run spatial statistics

workflow spatial_statistics {
	input {
		String cohort_id
		File clustered_adata_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	call spatially_variable_gene_analysis {
		input:
			cohort_id = cohort_id,
			clustered_adata_object = clustered_adata_object,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		# Moran’s I global spatial auto-correlation statistics
		File final_adata_object = spatially_variable_gene_analysis.final_adata_object #!FileCoercion
		File final_metadata_csv = spatially_variable_gene_analysis.final_metadata_csv #!FileCoercion
		File moran_top_10_variable_genes_csv = spatially_variable_gene_analysis.moran_top_10_variable_genes_csv #!FileCoercion
		File moran_top_4_variable_genes_spatial_scatter_plot_png = spatially_variable_gene_analysis.moran_top_4_variable_genes_spatial_scatter_plot_png #!FileCoercion
	}

	meta {
		description: "Perform spatial statistics for downstream analysis."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		clustered_adata_object: {help: "Leiden clustered AnnData object."}
		raw_data_path: {help: "Raw data bucket path for spatial statistics workflow outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}

task spatially_variable_gene_analysis {
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

		visium_spatially_variable_genes \
			--cohort-id ~{cohort_id} \
			--adata-input ~{clustered_adata_object} \
			--adata-output ~{cohort_id}.final.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.final.h5ad" \
			-o "~{cohort_id}.final_metadata.csv" \
			-o "~{cohort_id}.moran_top_10_variable_genes.csv" \
			-o "~{cohort_id}.moran_top_4_variable_genes_spatial_scatter.png"
	>>>

	output {
		String final_adata_object = "~{raw_data_path}/~{cohort_id}.final.h5ad"
		String final_metadata_csv = "~{raw_data_path}/~{cohort_id}.final_metadata.csv"
		String moran_top_10_variable_genes_csv = "~{raw_data_path}/~{cohort_id}.moran_top_10_variable_genes.csv"
		String moran_top_4_variable_genes_spatial_scatter_plot_png = "~{raw_data_path}/~{cohort_id}.moran_top_4_variable_genes_spatial_scatter.png"
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.1"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 15
		zones: zones
	}

	meta {
		description: "Identify spatially variable genes by computing Moran's I Score."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		clustered_adata_object: {help: "Leiden clustered AnnData object."}
		raw_data_path: {help: "Raw data bucket path for spatial statistics workflow outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}
