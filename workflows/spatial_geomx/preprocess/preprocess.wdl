version 1.0

# Generate a QC'ed RDS object by converting FASTQ files to DCC (digital count conversion) files to a NanoStringGeoMxSet object

import "../structs.wdl"

workflow preprocess {
	input {
		String team_id
		String dataset_id
		String dataset_doi_url
		Array[Slide] slides

		File project_sample_metadata_csv
		File geomx_config_ini
		File geomxngs_config_pkc

		Int min_segment_reads
		Int min_percent_reads_trimmed
		Int min_percent_reads_stitched
		Int min_percent_reads_aligned
		Int min_saturation
		Int min_neg_ctrl_count
		Int max_ntc_count
		Int min_nuclei
		Int min_segment_area

		String workflow_name
		String workflow_version
		String workflow_release
		String run_timestamp
		String raw_data_path_prefix
		String billing_project
		String container_registry
		String zones
	}

	# Task and subworkflow versions
	String sub_workflow_name = "preprocess"
	String fastq_to_dcc_task_version = "1.0.0"
	String dcc_to_rds_task_version = "1.0.0"
	String qc_task_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String dcc_raw_data_path = "~{workflow_raw_data_path_prefix}/fastq_to_dcc/~{fastq_to_dcc_task_version}"
	String rds_raw_data_path = "~{workflow_raw_data_path_prefix}/dcc_to_rds/~{dcc_to_rds_task_version}"
	String qc_raw_data_path = "~{workflow_raw_data_path_prefix}/qc/~{qc_task_version}"

	scatter (slide_object in slides) {
		String fastq_to_dcc_output = "~{dcc_raw_data_path}/~{slide_object.slide_id}.geomxngs_out_dir.tar.gz"
		String dcc_to_rds_output = "~{rds_raw_data_path}/~{slide_object.slide_id}.NanoStringGeoMxSet.rds"
		String qc_output = "~{qc_raw_data_path}/~{slide_object.slide_id}.qc.rds"
	}

	# For each sample, outputs an array of true/false: [fastq_to_dcc_complete, dcc_to_rds_complete, qc_complete]
	call check_output_files_exist {
		input:
			fastq_to_dcc_output_files = fastq_to_dcc_output,
			dcc_to_rds_output_files = dcc_to_rds_output,
			qc_output_files = qc_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (slide_index in range(length(slides))) {
		Slide slide = slides[slide_index]

		scatter (sample_index in range(length(slide.samples))) {
			Sample sample = slide.samples[sample_index]
			Array[String] project_sample_id = [team_id, sample.sample_id, dataset_doi_url]
			Array[File] fastq_R1s = sample.fastq_R1s
			Array[File] fastq_R2s = sample.fastq_R2s
		}

		String fastq_to_dcc_complete = check_output_files_exist.sample_preprocessing_complete[slide_index][0]
		String dcc_to_rds_complete = check_output_files_exist.sample_preprocessing_complete[slide_index][1]
		String qc_complete = check_output_files_exist.sample_preprocessing_complete[slide_index][2]

		String fastq_to_dcc_geomxngs_dcc_zip = "~{dcc_raw_data_path}/~{slide.slide_id}.DCC.zip"
		String fastq_to_dcc_geomxngs_output_tar_gz = "~{dcc_raw_data_path}/~{slide.slide_id}.geomxngs_out_dir.tar.gz"

		if (fastq_to_dcc_complete == "false") {
			call fastq_to_dcc {
				input:
					slide_id = slide.slide_id,
					fastq_R1s = flatten(fastq_R1s),
					fastq_R2s = flatten(fastq_R2s),
					geomx_config_ini = geomx_config_ini,
					raw_data_path = dcc_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File geomxngs_dcc_zip_output = select_first([fastq_to_dcc.geomxngs_dcc_zip, fastq_to_dcc_geomxngs_dcc_zip]) #!FileCoercion
		File geomxngs_output_tar_gz_output = select_first([fastq_to_dcc.geomxngs_output_tar_gz, fastq_to_dcc_geomxngs_output_tar_gz]) #!FileCoercion

		String dcc_to_rds_object = "~{rds_raw_data_path}/~{slide.slide_id}.NanoStringGeoMxSet.rds"

		if (dcc_to_rds_complete == "false") {
			call dcc_to_rds {
				input:
					team_id = team_id,
					dataset_id = dataset_id,
					slide_id = slide.slide_id,
					project_sample_metadata_csv = project_sample_metadata_csv,
					geomxngs_dcc_zip = geomxngs_dcc_zip_output,
					geomx_lab_annotation_xlsx = slide.geomx_lab_annotation_xlsx,
					geomxngs_config_pkc = geomxngs_config_pkc,
					raw_data_path = rds_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File initial_rds_object_output = select_first([dcc_to_rds.initial_rds_object, dcc_to_rds_object]) #!FileCoercion

		String qc_metrics_rds_object = "~{qc_raw_data_path}/~{slide.slide_id}.qc.rds"
		String qc_segment_summary_csv = "~{qc_raw_data_path}/~{slide.slide_id}.segment_qc_summary.csv"
		String qc_probe_summary_csv = "~{qc_raw_data_path}/~{slide.slide_id}.probe_qc_summary.csv"
		String qc_gene_count_csv = "~{qc_raw_data_path}/~{slide.slide_id}.gene_count.csv"

		if (qc_complete == "false") {
			call qc {
				input:
					slide_id = slide.slide_id,
					initial_rds_object = initial_rds_object_output,
					min_segment_reads = min_segment_reads,
					min_percent_reads_trimmed = min_percent_reads_trimmed,
					min_percent_reads_stitched = min_percent_reads_stitched,
					min_percent_reads_aligned = min_percent_reads_aligned,
					min_saturation = min_saturation,
					min_neg_ctrl_count = min_neg_ctrl_count,
					max_ntc_count = max_ntc_count,
					min_nuclei = min_nuclei,
					min_segment_area = min_segment_area,
					raw_data_path = qc_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File qc_rds_object_output = select_first([qc.qc_rds_object, qc_metrics_rds_object]) #!FileCoercion
		File segment_qc_summary_csv_output = select_first([qc.segment_qc_summary_csv, qc_segment_summary_csv]) #!FileCoercion
		File probe_qc_summary_csv_output = select_first([qc.probe_qc_summary_csv, qc_probe_summary_csv]) #!FileCoercion
		File gene_count_csv_output = select_first([qc.gene_count_csv, qc_gene_count_csv]) #!FileCoercion
	}

	output {
		# Sample list
		Array[Array[String]] project_sample_ids = flatten(project_sample_id)

		# GeoMxNGSPipeline outputs including converted DCC files
		Array[File] geomxngs_dcc_zip = geomxngs_dcc_zip_output #!FileCoercion
		Array[File] geomxngs_output_tar_gz = geomxngs_output_tar_gz_output #!FileCoercion

		# Initial RDS object
		Array[File] initial_rds_object = initial_rds_object_output #!FileCoercion

		# QC RDS object and tables
		Array[File] qc_rds_object = qc_rds_object_output #!FileCoercion
		Array[File] segment_qc_summary_csv = segment_qc_summary_csv_output #!FileCoercion
		Array[File] probe_qc_summary_csv = probe_qc_summary_csv_output #!FileCoercion
		Array[File] gene_count_csv = gene_count_csv_output #!FileCoercion
	}

	meta {
		description: "Preprocess the Nanostring GeoMx data by running the GeoMxNGSPipeline, convert DCC files to NanoStringGeoMxSet objects, and QC with GeoMx R libraries."
	}

	parameter_meta {
		team_id: {help: "Name of the CRN Team; stored in the NanoStringGeoMxSet objects."}
		dataset_id: {help: "Generated ASAP dataset ID; stored in the NanoStringGeoMxSet objects."}
		dataset_doi_url: {help: "Generated Zenodo DOI URL referencing the dataset."}
		slides: {help: "An array of Slide struct, set of slides and their Samples and GeoMx experimental information including the annotation (.xlsx) file/lab worksheet, containing phenotypic data from the GeoMx DSP readout package."}
		samples: {help: "An array of Sample struct, set of samples within a slide and their associated reads."}
		geomx_config_ini: {help: "The configuration (.ini) file, containing pipeline processing parameters that is used by the GeoMx NGS pipeline to assist in converting the FASTQ files to DCC files. It is from the GeoMx DSP readout package."}
		geomxngs_config_pkc: {help: "The GeoMx DSP configuration file to associate assay targets with GeoMx HybCode barcodes and Seq Code primers."}
		min_segment_reads: {help: "Minimum number of segment reads. [1000]"}
		min_percent_reads_trimmed: {help: "Minimum % of reads trimmed. [80]"}
		min_percent_reads_stitched: {help: "Minimum % of reads stitched. [80]"}
		min_percent_reads_aligned: {help: "Minimum % of reads aligned. [80]"}
		min_saturation: {help: "Minimum sequencing saturation. [50]"}
		min_neg_ctrl_count: {help: "Minimum negative control counts. [1]"}
		max_ntc_count: {help: "Maximum counts observed in NTC well. [1000]"}
		min_nuclei: {help: "Minimum # of nuclei estimated. [100]"}
		min_segment_area: {help: "Minimum segment area. [5000]"}
		workflow_name: {help: "Workflow name; stored in the file-level manifest and final manifest with all saved files."}
		workflow_version: {help: "Workflow version; stored in the file-level manifest and final manifest with all saved files."}
		workflow_release: {help: "GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		run_timestamp: {help: "UTC timestamp; stored in the file-level manifest and final manifest with all saved files."}
		raw_data_path_prefix: {help: "Raw data bucket path prefix; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess`)."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}

task check_output_files_exist {
	input {
		Array[String] fastq_to_dcc_output_files
		Array[String] dcc_to_rds_output_files
		Array[String] qc_output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			dcc_file=$(echo "${output_files}" | cut -f 1)
			rds_file=$(echo "${output_files}" | cut -f 2)
			qc_file=$(echo "${output_files}" | cut -f 3)

			if gcloud storage ls --billing-project=~{billing_project} "${dcc_file}"; then
				if gcloud storage ls --billing-project=~{billing_project} "${rds_file}"; then
					if gcloud storage ls --billing-project=~{billing_project} "${qc_file}"; then
						# If we find all outputs, don't rerun anything
						echo -e "true\ttrue\ttrue" >> sample_preprocessing_complete.tsv
					else
						# If we find fastq_to_dcc and dcc_to_rds outputs, then run (or rerun) qc
						echo -e "true\ttrue\tfalse" >> sample_preprocessing_complete.tsv
					fi
				else
					# If we only find fastq_to_dcc, then run (or rerun) dcc_to_rds and qc
					echo -e "true\tfalse\tfalse" >> sample_preprocessing_complete.tsv
				fi
			else
				# If we don't find fast_to_dcc output, we must need to run (or rerun) preprocessing
				echo -e "false\tfalse\tfalse" >> sample_preprocessing_complete.tsv
			fi
		done < <(paste ~{write_lines(fastq_to_dcc_output_files)} ~{write_lines(dcc_to_rds_output_files)} ~{write_lines(qc_output_files)})
	>>>

	output {
		Array[Array[String]] sample_preprocessing_complete = read_tsv("sample_preprocessing_complete.tsv")
	}

	runtime {
		docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:524.0.0-slim"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk 20 HDD"
		preemptible: 3
		zones: zones
	}

	meta {
		description: "Checks for existing preprocessing files per sample and skips certain preprocessing steps if they exist."
	}

	parameter_meta {
		fastq_to_dcc_output_files: {help: "Converted DCC output file to detect (`<sample>.geomxngs_out_dir.tar.gz`)."}
		dcc_to_rds_output_files: {help: "Converted NanoStringGeoMxSet (.RDS) object output file to detect (`<sample>.NanoStringGeoMxSet.rds`)."}
		qc_output_files: {help: "QC'ed output file to detect (`<sample>.qc.rds`)."}
		billing_project: {help: "Billing project to charge GCP costs."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}

task fastq_to_dcc {
	input {
		String slide_id

		Array[File] fastq_R1s
		Array[File] fastq_R2s

		File geomx_config_ini

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 8
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size(flatten([fastq_R1s, fastq_R2s]), "GB") * 2 + 50)

	command <<<
		set -euo pipefail
	
		# Ensure fastqs are in the same directory
		mkdir fastqs
		while read -r fastq || [[ -n "${fastq}" ]]; do
			if [[ -n "${fastq}" ]]; then
				fastq_basename=$(basename "${fastq}")
				fastq_sample_id=$(cut -d'_' -f1-4 <<< "${fastq_basename}")
				validated_fastq_name=$(fix_fastq_names --fastq "${fastq}" --sample-id "${fastq_sample_id}")
				if [[ -e "fastqs/${validated_fastq_name}" ]]; then
					echo "[WARNING] Skipping fastq renaming; already in proper format [${validated_fastq_name}]."
				else
					ln -s "${fastq}" "fastqs/${validated_fastq_name}"
				fi
			fi
		done < <(cat \
			~{write_lines(fastq_R1s)} \
			~{write_lines(fastq_R2s)})

		expect <<EOF
		set timeout -1
		spawn geomxngspipeline \
			--ini=~{geomx_config_ini} \
			--in="$(pwd)/fastqs" \
			--out="~{slide_id}_geomxngs_out_dir" \
			--check-illumina-naming=false
		send -- "2"
		expect eof
		EOF

		# Only keep the sample that was processed or else all samples in the config (.ini) file is processed with zero RTS_ID count causing errors downstream
		mkdir ~{slide_id}.DCC
		touch sample_names.txt
		while read -r file || [[ -n "${file}" ]]; do
			basename "$file" | cut -d '_' -f 1-4 | grep "^DSP"
		done < <(ls fastqs) > sample_names.txt
		sed 's/$/.dcc/' sample_names.txt > sample_names_with_dcc.txt
		while read -r file || [[ -n "${file}" ]]; do
			cp ~{slide_id}_geomxngs_out_dir/"$file" ./~{slide_id}.DCC/
		done < sample_names_with_dcc.txt
		zip -r ~{slide_id}.DCC.zip ~{slide_id}.DCC

		random_dcc=$(head -1 sample_names_with_dcc.txt)
		detect_empty_counts=$(grep "^Raw,0$" ./~{slide_id}.DCC/"${random_dcc}" || [[ $? == 1 ]])
		if [[ -n "${detect_empty_counts}" ]]; then
			echo "[ERROR] Checked a DCC file generated from present fastqs and there are no counts [${random_dcc}]. Exiting."
			exit 1
		fi

		tar -czvf "~{slide_id}.geomxngs_out_dir.tar.gz" "~{slide_id}_geomxngs_out_dir"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{slide_id}.DCC.zip" \
			-o "~{slide_id}.geomxngs_out_dir.tar.gz"
	>>>

	output {
		String geomxngs_dcc_zip = "~{raw_data_path}/~{slide_id}.DCC.zip"
		String geomxngs_output_tar_gz = "~{raw_data_path}/~{slide_id}.geomxngs_out_dir.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/geomx_ngs:3.1.1.6"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 15
		zones: zones
	}

	meta {
		description: "Processes raw sequencing data from Nanostring GeoMx experiments using Nanostring's GeoMxNGSPipeline to demultiplex raw reads using probe and sample information, align and count unique molecular identifiers (UMIs), and generate gene expression count matrices along with QC metrics."
	}

	parameter_meta {
		slide_id: {help: "Generated slide ID; used to name output files."}
		fastq_R1s: {help: "Sample's read 1 FASTQ file."}
		fastq_R2s: {help: "Sample's read 2 FASTQ file."}
		geomx_config_ini: {help: "The configuration (.ini) file, containing pipeline processing parameters that is used by the GeoMx NGS pipeline to assist in converting the FASTQ files to DCC files. It is from the GeoMx DSP readout package."}
		raw_data_path: {help: "Raw data bucket path for GeoMxNGSPipeline outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess/fastq_to_dcc/<fastq_to_dcc_task_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}

task dcc_to_rds {
	input {
		String team_id
		String dataset_id
		String slide_id

		File project_sample_metadata_csv

		File geomxngs_dcc_zip
		File geomx_lab_annotation_xlsx
		File geomxngs_config_pkc

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 4
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([project_sample_metadata_csv, geomxngs_dcc_zip, geomx_lab_annotation_xlsx, geomxngs_config_pkc], "GB") * 4 + 30)

	command <<<
		set -euo pipefail

		unzip -d ./dcc_files_dir -j ~{geomxngs_dcc_zip}

		geomx_counts_to_rds \
			--team-id ~{team_id} \
			--dataset-id ~{dataset_id} \
			--slide-id ~{slide_id} \
			--sample-metadata ~{project_sample_metadata_csv} \
			--dcc-dir ./dcc_files_dir \
			--pkc-file ~{geomxngs_config_pkc} \
			--annotation-file ~{geomx_lab_annotation_xlsx} \
			--output ~{slide_id}.NanoStringGeoMxSet.rds

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{slide_id}.NanoStringGeoMxSet.rds"
	>>>

	output {
		String initial_rds_object = "~{raw_data_path}/~{slide_id}.NanoStringGeoMxSet.rds"
	}

	runtime {
		docker: "~{container_registry}/spatial_r:1.0.0"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 15
		zones: zones
	}

	meta {
		description: "Convert digital count conversion (DCC) files to a NanoStringGeoMxSet (.RDS) object."
	}

	parameter_meta {
		team_id: {help: "Name of the CRN Team; stored in the NanoStringGeoMxSet objects."}
		dataset_id: {help: "Generated ASAP dataset ID; stored in the NanoStringGeoMxSet objects."}
		slide_id: {help: "Generated slide ID; used to name output files."}
		project_sample_metadata_csv: {help: "QC'ed SAMPLE.csv metadata used to match and store sample_id, ASAP_sample_id, and batch in the NanoStringGeoMxSet object."}
		geomxngs_dcc_zip: {help: "DCC files for each sample compressed in a ZIP file."}
		geomx_lab_annotation_xlsx: {help: "The annotation (.xlsx) file/lab worksheet, containing phenotypic data from the GeoMx DSP readout package."}
		geomxngs_config_pkc: {help: "The GeoMx DSP configuration file to associate assay targets with GeoMx HybCode barcodes and Seq Code primers."}
		raw_data_path: {help: "Raw data bucket path for converted RDS file outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess/dcc_to_rds/<dcc_to_rds_task_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}

task qc {
	input {
		String slide_id

		File initial_rds_object

		Int min_segment_reads
		Int min_percent_reads_trimmed
		Int min_percent_reads_stitched
		Int min_percent_reads_aligned
		Int min_saturation
		Int min_neg_ctrl_count
		Int max_ntc_count
		Int min_nuclei
		Int min_segment_area

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 4
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size(initial_rds_object, "GB") * 2 + 30)

	command <<<
		set -euo pipefail

		geomx_qc \
			--slide-id ~{slide_id} \
			--input ~{initial_rds_object} \
			--min-reads ~{min_segment_reads} \
			--percent-trimmed ~{min_percent_reads_trimmed} \
			--percent-stitched ~{min_percent_reads_stitched} \
			--percent-aligned ~{min_percent_reads_aligned} \
			--percent-saturation ~{min_saturation} \
			--min-neg-count ~{min_neg_ctrl_count} \
			--max-ntc-count ~{max_ntc_count} \
			--min-nuclei ~{min_nuclei} \
			--min-area ~{min_segment_area} \
			--output ~{slide_id}.qc.rds

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{slide_id}.qc.rds" \
			-o "~{slide_id}.segment_qc_summary.csv" \
			-o "~{slide_id}.probe_qc_summary.csv" \
			-o "~{slide_id}.gene_count.csv"
	>>>

	output {
		String qc_rds_object = "~{raw_data_path}/~{slide_id}.qc.rds"
		String segment_qc_summary_csv = "~{raw_data_path}/~{slide_id}.segment_qc_summary.csv"
		String probe_qc_summary_csv = "~{raw_data_path}/~{slide_id}.probe_qc_summary.csv"
		String gene_count_csv = "~{raw_data_path}/~{slide_id}.gene_count.csv"
	}

	runtime {
		docker: "~{container_registry}/spatial_r:1.0.0"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 15
		zones: zones
	}

	meta {
		description: "Calculate QC metrics on Nanostring GeoMx data with segments and probes, and generate gene-level count data."
	}

	parameter_meta {
		slide_id: {help: "Generated slide ID; used to name output files."}
		initial_rds_object: {help: "The initial RDS object converted from counts."}
		min_segment_reads: {help: "Minimum number of segment reads. [1000]"}
		min_percent_reads_trimmed: {help: "Minimum % of reads trimmed. [80]"}
		min_percent_reads_stitched: {help: "Minimum % of reads stitched. [80]"}
		min_percent_reads_aligned: {help: "Minimum % of reads aligned. [80]"}
		min_saturation: {help: "Minimum sequencing saturation. [50]"}
		min_neg_ctrl_count: {help: "Minimum negative control counts. [1]"}
		max_ntc_count: {help: "Maximum counts observed in NTC well. [1000]"}
		min_nuclei: {help: "Minimum # of nuclei estimated. [100]"}
		min_segment_area: {help: "Minimum segment area. [5000]"}
		raw_data_path: {help: "Raw data bucket path for converted RDS file outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess/dcc_to_rds/<dcc_to_rds_task_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}
