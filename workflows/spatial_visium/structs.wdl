version 1.0

struct Sample {
	String sample_id
	String? batch

	Array[File]+ fastq_R1s
	Array[File]+ fastq_R2s
	Array[File] fastq_I1s
	Array[File] fastq_I2s

	File visium_brightfield_image
	String visium_slide_serial_number
	String visium_capture_area
}

struct Project {
	String asap_team_id
	String asap_dataset_id
	String asap_dataset_doi_url
	Array[Sample] samples

	Boolean run_project_cohort_analysis

	String raw_data_bucket
	Array[String] staging_data_buckets
}
