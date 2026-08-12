version 1.0

# Generate a preprocessed AnnData object

import "../../wf-common/wdl/structs.wdl"
import "../../sc-rnaseq-wf/workflows/preprocess/preprocess.wdl" as ScRnaPreprocess
	alias Project as RnaProject
	alias Sample as RnaSample

workflow preprocess {
	input {
		String team_id
		String dataset_id
		String dataset_doi_url
		Array[Sample] samples

		File cellranger_arc_reference_data
		Float cellbender_fpr

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
	String cellranger_arc_task_version = "1.0.0"
	# TODO This must match https://github.com/ASAP-CRN/sc-rnaseq-wf/blob/main/workflows/preprocess/preprocess.wdl
	String cellbender_task_version = "1.0.1"
	String rna_adata_task_version = "1.2.0"
	String atac_adata_task_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String cellranger_arc_raw_data_path = "~{workflow_raw_data_path_prefix}/cellranger_arc/~{cellranger_arc_task_version}"
	String rna_cellbender_raw_data_path = "~{workflow_raw_data_path_prefix}/rna/remove_technical_artifacts/~{cellbender_task_version}"
	String rna_adata_raw_data_path = "~{workflow_raw_data_path_prefix}/rna/counts_to_adata/~{rna_adata_task_version}"
	String atac_adata_raw_data_path = "~{workflow_raw_data_path_prefix}/atac/counts_to_adata/~{atac_adata_task_version}"

	scatter (sample_object in samples) {
		String cellranger_arc_count_output = "~{cellranger_arc_raw_data_path}/~{dataset_id}.~{sample_object.sample_id}.raw_feature_bc_matrix.h5"
		String rna_cellbender_count_output = "~{rna_cellbender_raw_data_path}/~{dataset_id}.~{sample_object.sample_id}.cellbender.h5"
		String rna_initial_adata_object_output = "~{rna_adata_raw_data_path}/~{dataset_id}.~{sample_object.sample_id}.cleaned_unfiltered.h5ad"
		String atac_initial_adata_object_output = "~{atac_adata_raw_data_path}/~{dataset_id}.~{sample_object.sample_id}.cleaned_unfiltered.h5ad"
	}

	# For each sample, outputs an array of true/false: [cellranger_arc_counts_complete, rna_remove_technical_artifacts_complete, rna_initial_adata_object_complete, atac_initial_adata_object_complete]
	call check_output_files_exist {
		input:
			cellranger_arc_count_output_files = cellranger_arc_count_output,
			rna_remove_technical_artifacts_output_files = rna_cellbender_count_output,
			rna_initial_adata_object_output_files = rna_initial_adata_object_output,
			atac_initial_adata_object_output_files = atac_initial_adata_object_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (index in range(length(samples))) {
		Sample sample = samples[index]

		String dataset_sample_id = "~{dataset_id}.~{sample.sample_id}"
		Array[String] project_sample_id = [team_id, sample.sample_id, dataset_doi_url]

		String cellranger_arc_count_complete = check_output_files_exist.sample_preprocessing_complete[index][0]
		String rna_cellbender_remove_background_complete = check_output_files_exist.sample_preprocessing_complete[index][1]
		String rna_initial_adata_object_complete = check_output_files_exist.sample_preprocessing_complete[index][2]
		String atac_initial_adata_object_complete = check_output_files_exist.sample_preprocessing_complete[index][3]

		String cellranger_arc_multiome_outputs_tar_gz = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.cellranger_arc_multiome_outputs.tar.gz"
		String cellranger_arc_raw_counts = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.raw_feature_bc_matrix.h5"
		String cellranger_arc_filtered_counts = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.filtered_feature_bc_matrix.h5"
		String cellranger_arc_molecule_info = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.molecule_info.h5"
		String cellranger_arc_metrics_summary_csv = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.metrics_summary.csv"
		String cellranger_arc_per_barcode_metrics_summary_csv = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.per_barcode_metrics.csv"
		String cellranger_arc_gex_possorted_bam = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.gex_possorted_bam.bam"
		String cellranger_arc_gex_possorted_bam_index = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.gex_possorted_bam.bam.bai"
		String cellranger_arc_atac_possorted_bam = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.atac_possorted_bam.bam"
		String cellranger_arc_atac_possorted_bam_index = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.atac_possorted_bam.bam.bai"
		String cellranger_arc_atac_fragments_tsv_gz = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.atac_fragments.tsv.gz"
		String cellranger_arc_atac_fragments_tsv_gz_tbi = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.atac_fragments.tsv.gz.tbi"
		String cellranger_arc_atac_peaks_bed = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.atac_peaks.bed"
		String cellranger_arc_atac_cut_sites_bigwig = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.atac_cut_sites.bigwig"
		String cellranger_arc_atac_peak_annotation_tsv = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.atac_peak_annotation.tsv"
		String cellranger_arc_feature_linkage_bedpe = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.feature_linkage.bedpe"

		if (cellranger_arc_count_complete == "false") {
			call cellranger_arc_count {
				input:
					dataset_sample_id = dataset_sample_id,
					sample_id = sample.sample_id,
					fastq_R1s = sample.fastq_R1s,
					fastq_R2s = sample.fastq_R2s,
					fastq_I1s = sample.fastq_I1s,
					fastq_I2s = sample.fastq_I2s,
					cellranger_arc_reference_data = cellranger_arc_reference_data,
					raw_data_path = cellranger_arc_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File multiome_outputs_tar_gz_output = select_first([cellranger_arc_count.multiome_outputs_tar_gz, cellranger_arc_multiome_outputs_tar_gz]) #!FileCoercion
		File raw_counts_output = select_first([cellranger_arc_count.raw_counts, cellranger_arc_raw_counts]) #!FileCoercion
		File filtered_counts_output = select_first([cellranger_arc_count.filtered_counts, cellranger_arc_filtered_counts]) #!FileCoercion
		File molecule_info_output = select_first([cellranger_arc_count.molecule_info, cellranger_arc_molecule_info]) #!FileCoercion
		File metrics_summary_csv_output = select_first([cellranger_arc_count.metrics_summary_csv, cellranger_arc_metrics_summary_csv]) #!FileCoercion
		File per_barcode_metrics_summary_csv_output = select_first([cellranger_arc_count.per_barcode_metrics_summary_csv, cellranger_arc_per_barcode_metrics_summary_csv]) #!FileCoercion
		File gex_possorted_bam_output = select_first([cellranger_arc_count.gex_possorted_bam, cellranger_arc_gex_possorted_bam]) #!FileCoercion
		File gex_possorted_bam_index_output = select_first([cellranger_arc_count.gex_possorted_bam_index, cellranger_arc_gex_possorted_bam_index]) #!FileCoercion
		File atac_possorted_bam_output = select_first([cellranger_arc_count.atac_possorted_bam, cellranger_arc_atac_possorted_bam]) #!FileCoercion
		File atac_possorted_bam_index_output = select_first([cellranger_arc_count.atac_possorted_bam_index, cellranger_arc_atac_possorted_bam_index]) #!FileCoercion
		File atac_fragments_tsv_gz_output = select_first([cellranger_arc_count.atac_fragments_tsv_gz, cellranger_arc_atac_fragments_tsv_gz]) #!FileCoercion
		File atac_fragments_tsv_gz_tbi_output = select_first([cellranger_arc_count.atac_fragments_tsv_gz_tbi, cellranger_arc_atac_fragments_tsv_gz_tbi]) #!FileCoercion
		File atac_peaks_bed_output = select_first([cellranger_arc_count.atac_peaks_bed, cellranger_arc_atac_peaks_bed]) #!FileCoercion
		File atac_cut_sites_bigwig_output = select_first([cellranger_arc_count.atac_cut_sites_bigwig, cellranger_arc_atac_cut_sites_bigwig]) #!FileCoercion
		File atac_peak_annotation_tsv_output = select_first([cellranger_arc_count.atac_peak_annotation_tsv, cellranger_arc_atac_peak_annotation_tsv]) #!FileCoercion
		File feature_linkage_bedpe_output = select_first([cellranger_arc_count.feature_linkage_bedpe, cellranger_arc_feature_linkage_bedpe]) #!FileCoercion
	
		String cellbender_report_html = "~{rna_cellbender_raw_data_path}/~{dataset_sample_id}.cellbender_report.html"
		String cellbender_removed_background_counts = "~{rna_cellbender_raw_data_path}/~{dataset_sample_id}.cellbender.h5"
		String cellbender_filtered_removed_background_counts = "~{rna_cellbender_raw_data_path}/~{dataset_sample_id}.cellbender_filtered.h5"
		String cellbender_cell_barcodes_csv = "~{rna_cellbender_raw_data_path}/~{dataset_sample_id}.cellbender_cell_barcodes.csv"
		String cellbender_graph_pdf = "~{rna_cellbender_raw_data_path}/~{dataset_sample_id}.cellbender.pdf"
		String cellbender_log = "~{rna_cellbender_raw_data_path}/~{dataset_sample_id}.cellbender.log"
		String cellbender_metrics_csv = "~{rna_cellbender_raw_data_path}/~{dataset_sample_id}.cellbender_metrics.csv"
		String cellbender_posterior_probability = "~{rna_cellbender_raw_data_path}/~{dataset_sample_id}.cellbend_posterior.h5"

		if (rna_cellbender_remove_background_complete == "false") {
			call ScRnaPreprocess.remove_technical_artifacts as rna_remove_technical_artifacts {
				input:
					dataset_sample_id = dataset_sample_id,
					raw_counts = raw_counts_output,
					cellbender_fpr = cellbender_fpr,
					raw_data_path = rna_cellbender_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File report_html_output = select_first([rna_remove_technical_artifacts.report_html, cellbender_report_html]) #!FileCoercion
		File removed_background_counts_output = select_first([rna_remove_technical_artifacts.removed_background_counts, cellbender_removed_background_counts]) #!FileCoercion
		File filtered_removed_background_counts_output = select_first([rna_remove_technical_artifacts.filtered_removed_background_counts, cellbender_filtered_removed_background_counts]) #!FileCoercion
		File cell_barcodes_csv_output = select_first([rna_remove_technical_artifacts.cell_barcodes_csv, cellbender_cell_barcodes_csv]) #!FileCoercion
		File graph_pdf_output = select_first([rna_remove_technical_artifacts.graph_pdf, cellbender_graph_pdf]) #!FileCoercion
		File log_output = select_first([rna_remove_technical_artifacts.log, cellbender_log]) #!FileCoercion
		File metrics_csv_output = select_first([rna_remove_technical_artifacts.metrics_csv, cellbender_metrics_csv]) #!FileCoercion
		File posterior_probability_output = select_first([rna_remove_technical_artifacts.posterior_probability, cellbender_posterior_probability]) #!FileCoercion

		String preprocessed_rna_adata_object = "~{rna_adata_raw_data_path}/~{dataset_sample_id}.cleaned_unfiltered.h5ad"

		if (rna_initial_adata_object_complete == "false") {
			call ScRnaPreprocess.counts_to_adata as rna_counts_to_adata {
				input:
					dataset_sample_id = dataset_sample_id,
					sample_id = sample.sample_id,
					batch = select_first([sample.batch]),
					sex = select_first([sample.sex]),
					team_id = team_id,
					dataset_id = dataset_id,
					cellbender_counts = removed_background_counts_output,
					raw_data_path = rna_adata_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File preprocessed_rna_adata_object_output = select_first([rna_counts_to_adata.initial_adata_object, preprocessed_rna_adata_object]) #!FileCoercion
		
		String preprocessed_atac_adata_object = "~{atac_adata_raw_data_path}/~{dataset_sample_id}.cleaned_unfiltered.h5ad"

		if (atac_initial_adata_object_complete == "false") {
			call atac_counts_to_adata {
				input:
					dataset_sample_id = dataset_sample_id,
					sample_id = sample.sample_id,
					batch = select_first([sample.batch]),
					sex = select_first([sample.sex]),
					team_id = team_id,
					dataset_id = dataset_id,
					fragments_tsv_gz = atac_fragments_tsv_gz_output,
					raw_data_path = atac_adata_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File preprocessed_atac_adata_object_output = select_first([atac_counts_to_adata.initial_adata_object, preprocessed_atac_adata_object]) #!FileCoercion
	}

	output {
		# Sample list
		Array[Array[String]] project_sample_ids = project_sample_id

		# Cellranger ARC
		Array[File] multiome_outputs_tar_gz = multiome_outputs_tar_gz_output #!FileCoercion
		Array[File] raw_counts = raw_counts_output #!FileCoercion
		Array[File] filtered_counts = filtered_counts_output #!FileCoercion
		Array[File] molecule_info = molecule_info_output #!FileCoercion
		Array[File] metrics_summary_csv = metrics_summary_csv_output #!FileCoercion
		Array[File] per_barcode_metrics_summary_csv = per_barcode_metrics_summary_csv_output #!FileCoercion
		Array[File] gex_possorted_bam = gex_possorted_bam_output #!FileCoercion
		Array[File] gex_possorted_bam_index = gex_possorted_bam_index_output #!FileCoercion
		Array[File] atac_possorted_bam = atac_possorted_bam_output #!FileCoercion
		Array[File] atac_possorted_bam_index = atac_possorted_bam_index_output #!FileCoercion
		Array[File] atac_fragments_tsv_gz = atac_fragments_tsv_gz_output #!FileCoercion
		Array[File] atac_fragments_tsv_gz_tbi = atac_fragments_tsv_gz_tbi_output #!FileCoercion
		Array[File] atac_peaks_bed = atac_peaks_bed_output #!FileCoercion
		Array[File] atac_cut_sites_bigwig = atac_cut_sites_bigwig_output #!FileCoercion
		Array[File] atac_peak_annotation_tsv = atac_peak_annotation_tsv_output #!FileCoercion
		Array[File] feature_linkage_bedpe = feature_linkage_bedpe_output #!FileCoercion

		# RNA - Remove technical artifacts with Cellbender
		Array[File] report_html = report_html_output
		Array[File] removed_background_counts = removed_background_counts_output #!FileCoercion
		Array[File] filtered_removed_background_counts = filtered_removed_background_counts_output #!FileCoercion
		Array[File] cell_barcodes_csv = cell_barcodes_csv_output #!FileCoercion
		Array[File] graph_pdf = graph_pdf_output #!FileCoercion
		Array[File] log = log_output #!FileCoercion
		Array[File] metrics_csv = metrics_csv_output #!FileCoercion
		Array[File] posterior_probability = posterior_probability_output #!FileCoercion

		# RNA - AnnData counts
		Array[File] rna_initial_adata_object = preprocessed_rna_adata_object_output #!FileCoercion

		# ATAC - AnnData counts
		Array[File] atac_initial_adata_object = preprocessed_atac_adata_object_output #!FileCoercion
	}

	meta {
		description: "Preprocess the 10x Chromium Single Cell Epi Multiome data by running Cell Ranger ARC count and converting counts to AnnData object."
	}

	parameter_meta {
		team_id: {help: "Name of the CRN Team; stored in the AnnData objects."}
		dataset_id: {help: "Generated ASAP dataset ID; stored in the AnnData objects."}
		dataset_doi_url: {help: "Generated Zenodo DOI URL referencing the dataset."}
		samples: {help: "An array of Sample struct, set of samples and their associated reads and metadata information."}
		cellranger_arc_reference_data: {help: "CellRanger ARC reference data; see https://www.10xgenomics.com/support/software/cell-ranger-arc/downloads."}
		cellbender_fpr: {help: "Cellbender false positive rate. [0.0]"}
		workflow_name: {help: "Workflow name; stored in the file-level manifest and final manifest with all saved files."}
		workflow_version: {help: "Workflow version; stored in the file-level manifest and final manifest with all saved files."}
		workflow_release: {help: "GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		run_timestamp: {help: "UTC timestamp; stored in the file-level manifest and final manifest with all saved files."}
		raw_data_path_prefix: {help: "Raw data bucket path prefix; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess`)."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}

task check_output_files_exist {
	input {
		Array[String] cellranger_arc_count_output_files
		Array[String] rna_remove_technical_artifacts_output_files
		Array[String] rna_initial_adata_object_output_files
		Array[String] atac_initial_adata_object_output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			cellranger_arc_counts_file=$(echo "${output_files}" | cut -f 1)
			rna_cellbender_counts_file=$(echo "${output_files}" | cut -f 2)
			rna_initial_adata_object_file=$(echo "${output_files}" | cut -f 3)
			atac_initial_adata_object_file=$(echo "${output_files}" | cut -f 4)

			if gcloud storage ls --billing-project=~{billing_project} "${cellranger_arc_counts_file}"; then
				# Check RNA outputs
				if gcloud storage ls --billing-project=~{billing_project} "${rna_cellbender_counts_file}"; then
					if gcloud storage ls --billing-project=~{billing_project} "${rna_initial_adata_object_file}"; then
						# If we find all RNA outputs, don't rerun anything on the RNA side
						rna_result="true\ttrue\ttrue"
					else
						# If we find cellranger and cellbender outputs, but don't find RNA adata outputs, just rerun counts_to_adata
						rna_result="true\ttrue\tfalse"
					fi
				else
					# If we find cellranger, but not cellbender outputs, it does not matter if adata objects exist, so run (or rerun) both cellbender and counts_to_adata
					rna_result="true\tfalse\tfalse"
				fi

				# Check ATAC outputs
				if gcloud storage ls --billing-project=~{billing_project} "${atac_initial_adata_object_file}"; then
					atac_result="true"
				else
					atac_result="false"
				fi

				echo -e "${rna_result}\t${atac_result}" >> sample_preprocessing_complete.tsv
			else
				# If we don't find cellranger output, we must also need to run (or rerun) preprocessing for both RNA and ATAC
				echo -e "false\tfalse\tfalse\tfalse" >> sample_preprocessing_complete.tsv
			fi
		done < <(paste ~{write_lines(cellranger_arc_count_output_files)} ~{write_lines(rna_remove_technical_artifacts_output_files)} ~{write_lines(rna_initial_adata_object_output_files)} ~{write_lines(atac_initial_adata_object_output_files)})
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
}

task cellranger_arc_count {
	input {
		String dataset_sample_id
		String sample_id

		Array[File] fastq_R1s
		Array[File] fastq_R2s
		Array[File] fastq_I1s
		Array[File] fastq_I2s

		File cellranger_arc_reference_data

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 16
	Int mem_gb = 48
	Int disk_size = ceil((size(cellranger_arc_reference_data, "GB") + size(flatten([fastq_R1s, fastq_R2s, fastq_I1s, fastq_I2s]), "GB")) * 4 + 50)

	command <<<
		set -euo pipefail

		# Unpack refdata
		mkdir cellranger_arc_refdata
		tar \
			-zxvf ~{cellranger_arc_reference_data} \
			-C cellranger_arc_refdata \
			--strip-components 1

		# Ensure fastqs are in the same directory
		mkdir fastqs
		while read -r fastq || [[ -n "${fastq}" ]]; do
			if [[ -n "${fastq}" ]]; then
				validated_fastq_name=$(fix_fastq_names --fastq "${fastq}" --sample-id "~{sample_id}")
				if [[ -e "fastqs/${validated_fastq_name}" ]]; then
					echo "[ERROR] Something's gone wrong with fastq renaming; trying to create fastq [${validated_fastq_name}] but it already exists. Exiting."
					exit 1
				else
					ln -s "${fastq}" "fastqs/${validated_fastq_name}"
				fi
			fi
		done < <(cat \
			~{write_lines(fastq_R1s)} \
			~{write_lines(fastq_R2s)} \
			~{write_lines(fastq_I1s)} \
			~{write_lines(fastq_I2s)})

		create_cellranger_libraries_csv \
			--sample-id ~{sample_id} \
			--output ~{sample_id}.libraries.csv

		cellranger-arc --version

		/usr/bin/time -v \
		cellranger-arc count \
			--id=~{sample_id} \
			--reference="$(pwd)/cellranger_arc_refdata" \
			--libraries="~{sample_id}.libraries.csv" \
			--create-bam=true \
			--localcores ~{threads} \
			--localmem ~{mem_gb - 4}

		# TMP
		ls -lR

		# Save Cell Ranger ARC outs
		cp -r ~{sample_id}/outs multiome_outputs
		tar -czvf "~{dataset_sample_id}.cellranger_arc_multiome_outputs.tar.gz" multiome_outputs

		# Rename outputs to include dataset and sample ID
		mv ~{sample_id}/outs/raw_feature_bc_matrix.h5 ~{dataset_sample_id}.raw_feature_bc_matrix.h5
		mv ~{sample_id}/outs/filtered_feature_bc_matrix.h5 ~{dataset_sample_id}.filtered_feature_bc_matrix.h5
		mv ~{sample_id}/outs/molecule_info.h5 ~{dataset_sample_id}.molecule_info.h5
		mv ~{sample_id}/outs/summary.csv ~{dataset_sample_id}.summary.csv
		mv ~{sample_id}/outs/per_barcode_metrics.csv ~{dataset_sample_id}.per_barcode_metrics.csv
		mv ~{sample_id}/outs/gex_possorted_bam.bam ~{dataset_sample_id}.gex_possorted_bam.bam
		mv ~{sample_id}/outs/gex_possorted_bam.bam.bai ~{dataset_sample_id}.gex_possorted_bam.bam.bai
		mv ~{sample_id}/outs/atac_possorted_bam.bam ~{dataset_sample_id}.atac_possorted_bam.bam
		mv ~{sample_id}/outs/atac_possorted_bam.bam.bai ~{dataset_sample_id}.atac_possorted_bam.bam.bai
		mv ~{sample_id}/outs/atac_fragments.tsv.gz ~{dataset_sample_id}.atac_fragments.tsv.gz
		mv ~{sample_id}/outs/atac_fragments.tsv.gz.tbi ~{dataset_sample_id}.atac_fragments.tsv.gz.tbi
		mv ~{sample_id}/outs/atac_peaks.bed ~{dataset_sample_id}.atac_peaks.bed
		mv ~{sample_id}/outs/atac_cut_sites.bigwig ~{dataset_sample_id}.atac_cut_sites.bigwig
		mv ~{sample_id}/outs/atac_peak_annotation.tsv ~{dataset_sample_id}.atac_peak_annotation.tsv
		mv ~{sample_id}/outs/analysis/feature_linkage/feature_linkage.bedpe ~{dataset_sample_id}.feature_linkage.bedpe

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{dataset_sample_id}.cellranger_arc_multiome_outputs.tar.gz" \
			-o "~{dataset_sample_id}.raw_feature_bc_matrix.h5" \
			-o "~{dataset_sample_id}.filtered_feature_bc_matrix.h5" \
			-o "~{dataset_sample_id}.molecule_info.h5" \
			-o "~{dataset_sample_id}.summary.csv" \
			-o "~{dataset_sample_id}.per_barcode_metrics.csv" \
			-o "~{dataset_sample_id}.gex_possorted_bam.bam" \
			-o "~{dataset_sample_id}.gex_possorted_bam.bam.bai" \
			-o "~{dataset_sample_id}.atac_possorted_bam.bam" \
			-o "~{dataset_sample_id}.atac_possorted_bam.bam.bai" \
			-o "~{dataset_sample_id}.atac_fragments.tsv.gz" \
			-o "~{dataset_sample_id}.atac_fragments.tsv.gz.tbi" \
			-o "~{dataset_sample_id}.atac_peaks.bed" \
			-o "~{dataset_sample_id}.atac_cut_sites.bigwig" \
			-o "~{dataset_sample_id}.atac_peak_annotation.tsv" \
			-o "~{dataset_sample_id}.feature_linkage.bedpe"
	>>>

	output {
		String multiome_outputs_tar_gz = "~{raw_data_path}/~{dataset_sample_id}.cellranger_arc_multiome_outputs.tar.gz"
		String raw_counts = "~{raw_data_path}/~{dataset_sample_id}.raw_feature_bc_matrix.h5"
		String filtered_counts = "~{raw_data_path}/~{dataset_sample_id}.filtered_feature_bc_matrix.h5"
		String molecule_info = "~{raw_data_path}/~{dataset_sample_id}.molecule_info.h5"
		String metrics_summary_csv = "~{raw_data_path}/~{dataset_sample_id}.summary.csv"
		String per_barcode_metrics_summary_csv = "~{raw_data_path}/~{dataset_sample_id}.per_barcode_metrics.csv"
		String gex_possorted_bam = "~{raw_data_path}/~{dataset_sample_id}.gex_possorted_bam.bam"
		String gex_possorted_bam_index = "~{raw_data_path}/~{dataset_sample_id}.gex_possorted_bam.bam.bai"
		String atac_possorted_bam = "~{raw_data_path}/~{dataset_sample_id}.atac_possorted_bam.bam"
		String atac_possorted_bam_index = "~{raw_data_path}/~{dataset_sample_id}.atac_possorted_bam.bam.bai"
		String atac_fragments_tsv_gz = "~{raw_data_path}/~{dataset_sample_id}.atac_fragments.tsv.gz"
		String atac_fragments_tsv_gz_tbi = "~{raw_data_path}/~{dataset_sample_id}.atac_fragments.tsv.gz.tbi"
		String atac_peaks_bed = "~{raw_data_path}/~{dataset_sample_id}.atac_peaks.bed"
		String atac_cut_sites_bigwig = "~{raw_data_path}/~{dataset_sample_id}.atac_cut_sites.bigwig"
		String atac_peak_annotation_tsv = "~{raw_data_path}/~{dataset_sample_id}.atac_peak_annotation.tsv"
		String feature_linkage_bedpe = "~{raw_data_path}/~{dataset_sample_id}.feature_linkage.bedpe"
	}

	runtime {
		docker: "~{container_registry}/cellranger_arc:2.2.0"
		cpu: threads
		cpuPlatform: "Intel Cascade Lake"
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}

	meta {
		description: "Processes raw sequencing data including same-cell open chromatin and gene expression profiling from 10x Chromium Single Cell Epi Multiome experiments to generate raw and filtered feature-barcode count matrices, and ATAC peaks and fragments."
	}

	parameter_meta {
		dataset_sample_id: {help: "Generated ASAP dataset ID and sample ID; stored in the AnnData objects."}
		sample_id: {help: "Generated ASAP sample ID; used to name output files."}
		fastq_R1s: {help: "Sample's read 1 FASTQ file."}
		fastq_R2s: {help: "Sample's read 2 FASTQ file."}
		fastq_I1s: {help: "Optional FASTQ index 1."}
		fastq_I2s: {help: "Optional FASTQ index 2."}
		cellranger_arc_reference_data: {help: "CellRanger ARC reference data; see https://www.10xgenomics.com/support/software/cell-ranger-arc/downloads."}
		raw_data_path: {help: "Raw data bucket path for cellranger-arc count outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess/cellranger_arc/<cellranger_arc_task_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}

task atac_counts_to_adata {
	input {
		String dataset_sample_id
		String sample_id
		String batch
		String sex
		String team_id
		String dataset_id

		File fragments_tsv_gz

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int disk_size = ceil(size(fragments_tsv_gz, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		/usr/bin/time -v \
		counts_to_adata \
			--cellranger-arc-fragments ~{fragments_tsv_gz} \
			--sample-id ~{sample_id} \
			--batch ~{batch} \
			--sex ~{sex} \
			--team ~{team_id} \
			--dataset ~{dataset_id} \
			--adata-output ~{dataset_sample_id}.cleaned_unfiltered.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{dataset_sample_id}.cleaned_unfiltered.h5ad"
	>>>

	output {
		String initial_adata_object = "~{raw_data_path}/~{dataset_sample_id}.cleaned_unfiltered.h5ad"
	}

	runtime {
		docker: "~{container_registry}/sc_multiome_tools:1.0.0"
		cpu: 4
		cpuPlatform: "Intel Cascade Lake"
		memory: "32 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}

	meta {
		description: "Converts Cell Ranger ARC counts into AnnData objects using SnapATAC2."
	}

	parameter_meta {
		dataset_sample_id: {help: "Generated ASAP dataset ID and sample ID; stored in the AnnData objects."}
		sample_id: {help: "Generated ASAP sample ID; stored in the AnnData objects and used to name output files."}
		batch: {help: "The sample's batch; stored in the AnnData objects."}
		sex: {help: "The sample's sex; stored in the AnnData objects."}
		team_id: {help: "Name of the CRN Team; stored in the AnnData objects."}
		dataset_id: {help: "Generated ASAP dataset ID; stored in the AnnData objects."}
		fragments_tsv_gz: {help: "Count and barcode information for every ATAC fragment observed in the experiment in TSV format output by Cell Ranger ARC."}
		raw_data_path: {help: "Raw data bucket path for counts to adata outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess/atac/counts_to_adata/<atac_adata_task_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}
