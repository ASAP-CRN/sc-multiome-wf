version 1.0

# Generate a preprocessed AnnData object

import "../../wf-common/wdl/structs.wdl"

workflow preprocess {
	input {
		String team_id
		String dataset_id
		String dataset_doi_url
		Array[Sample] samples

		File cellranger_reference_data

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

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String cellranger_arc_raw_data_path = "~{workflow_raw_data_path_prefix}/cellranger_arc/~{cellranger_arc_task_version}"

	scatter (index in range(length(samples))) {
		Sample sample = samples[index]

		Array[String] project_sample_id = [team_id, sample.sample_id, dataset_doi_url]

		String cellranger_arc_multiome_outputs_tar_gz = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.cellranger_arc_multiome_outputs.tar.gz"
		String cellranger_arc_raw_counts = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.raw_feature_bc_matrix.h5"
		String cellranger_arc_filtered_counts = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.filtered_feature_bc_matrix.h5"
		String cellranger_arc_molecule_info = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.molecule_info.h5"
		String cellranger_arc_metrics_summary_csv = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.metrics_summary.csv"
		String cellranger_arc_per_barcode_metrics_summary_csv = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.per_barcode_metrics.csv"
		String cellranger_arc_gex_possorted_bam = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.gex_possorted_bam.bam"
		String cellranger_arc_gex_possorted_bam_index = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.gex_possorted_bam.bam.bai"
		String cellranger_arc_atac_possorted_bam = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.atac_possorted_bam.bam"
		String cellranger_arc_atac_possorted_bam_index = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.atac_possorted_bam.bam.bai"
		String cellranger_arc_atac_fragments_tsv_gz = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.atac_fragments.tsv.gz"
		String cellranger_arc_atac_fragments_tsv_gz_tbi = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.atac_fragments.tsv.gz.tbi"
		String cellranger_arc_atac_peaks_bed = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.atac_peaks.bed"
		String cellranger_arc_atac_cut_sites_bigwig = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.atac_cut_sites.bigwig"
		String cellranger_arc_atac_peak_annotation_tsv = "~{cellranger_arc_raw_data_path}/~{sample.sample_id}.atac_peak_annotation.tsv"

		call cellranger_arc_count {
			input:
				sample_id = sample.sample_id,
				fastq_R1s = sample.fastq_R1s,
				fastq_R2s = sample.fastq_R2s,
				fastq_I1s = sample.fastq_I1s,
				fastq_I2s = sample.fastq_I2s,
				cellranger_arc_reference_data = cellranger_reference_data,
				raw_data_path = cellranger_arc_raw_data_path,
				workflow_info = workflow_info,
				billing_project = billing_project,
				container_registry = container_registry,
				zones = zones
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
	}

	output {
		# Sample list
		Array[Array[String]] project_sample_ids = project_sample_id

		# Cellranger ARC
		Array[File?] multiome_outputs_tar_gz = multiome_outputs_tar_gz_output #!FileCoercion
		Array[File?] raw_counts = raw_counts_output #!FileCoercion
		Array[File?] filtered_counts = filtered_counts_output #!FileCoercion
		Array[File?] molecule_info = molecule_info_output #!FileCoercion
		Array[File?] metrics_summary_csv = metrics_summary_csv_output #!FileCoercion
		Array[File?] per_barcode_metrics_summary_csv = per_barcode_metrics_summary_csv_output #!FileCoercion
		Array[File?] gex_possorted_bam = gex_possorted_bam_output #!FileCoercion
		Array[File?] gex_possorted_bam_index = gex_possorted_bam_index_output #!FileCoercion
		Array[File?] atac_possorted_bam = atac_possorted_bam_output #!FileCoercion
		Array[File?] atac_possorted_bam_index = atac_possorted_bam_index_output #!FileCoercion
		Array[File?] atac_fragments_tsv_gz = atac_fragments_tsv_gz_output #!FileCoercion
		Array[File?] atac_fragments_tsv_gz_tbi = atac_fragments_tsv_gz_tbi_output #!FileCoercion
		Array[File?] atac_peaks_bed = atac_peaks_bed_output #!FileCoercion
		Array[File?] atac_cut_sites_bigwig = atac_cut_sites_bigwig_output #!FileCoercion
		Array[File?] atac_peak_annotation_tsv = atac_peak_annotation_tsv_output #!FileCoercion
	}
}

task check_output_files_exist {
	input {
		Array[String] cellranger_count_output_files
		Array[String] remove_technical_artifacts_output_files
		Array[String] initial_adata_object_output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			cellranger_counts_file=$(echo "${output_files}" | cut -f 1)
			cellbender_counts_file=$(echo "${output_files}" | cut -f 2)
			initial_adata_object_file=$(echo "${output_files}" | cut -f 3)

			if gsutil -u ~{billing_project} ls "${cellranger_counts_file}"; then
				if gsutil -u ~{billing_project} ls "${cellbender_counts_file}"; then
					if gsutil -u ~{billing_project} ls "${initial_adata_object_file}"; then
						# If we find all outputs, don't rerun anything
						echo -e "true\ttrue\ttrue" >> sample_preprocessing_complete.tsv
					else
						# If we find cellranger and cellbender outputs, but don't find adata outputs, just rerun counts_to_adata
						echo -e "true\ttrue\tfalse" >> sample_preprocessing_complete.tsv
					fi
				else
					# If we find cellranger, but not cellbender outputs, it does not matter if adata objects exist, so run (or rerun) both cellbender and counts_to_adata
					echo -e "true\tfalse\tfalse" >> sample_preprocessing_complete.tsv
				fi
			else
				# If we don't find cellranger output, we must also need to run (or rerun) preprocessing
				echo -e "false\tfalse\tfalse" >> sample_preprocessing_complete.tsv
			fi
		done < <(paste ~{write_lines(cellranger_count_output_files)} ~{write_lines(remove_technical_artifacts_output_files)} ~{write_lines(initial_adata_object_output_files)})
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
	Int mem_gb = 24
	Int disk_size = ceil((size(fastq_R1s, "GB") + size(fastq_R2s, "GB") + size(fastq_I1s, "GB") + size(fastq_I2s, "GB") + size(cellranger_arc_reference_data, "GB")) * 4 + 50)

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

		/usr/bin/time \
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
		tar -czvf "~{sample_id}.cellranger_arc_multiome_outputs.tar.gz" multiome_outputs

		# Rename outputs to include sample ID
		mv ~{sample_id}/outs/raw_feature_bc_matrix.h5 ~{sample_id}.raw_feature_bc_matrix.h5
		mv ~{sample_id}/outs/filtered_feature_bc_matrix.h5 ~{sample_id}.filtered_feature_bc_matrix.h5
		mv ~{sample_id}/outs/molecule_info.h5 ~{sample_id}.molecule_info.h5
		mv ~{sample_id}/outs/summary.csv ~{sample_id}.summary.csv
		mv ~{sample_id}/outs/per_barcode_metrics.csv ~{sample_id}.per_barcode_metrics.csv
		mv ~{sample_id}/outs/gex_possorted_bam.bam ~{sample_id}.gex_possorted_bam.bam
		mv ~{sample_id}/outs/gex_possorted_bam.bam.bai ~{sample_id}.gex_possorted_bam.bam.bai
		mv ~{sample_id}/outs/atac_possorted_bam.bam ~{sample_id}.atac_possorted_bam.bam
		mv ~{sample_id}/outs/atac_possorted_bam.bam.bai ~{sample_id}.atac_possorted_bam.bam.bai
		mv ~{sample_id}/outs/atac_fragments.tsv.gz ~{sample_id}.atac_fragments.tsv.gz
		mv ~{sample_id}/outs/atac_fragments.tsv.gz.tbi ~{sample_id}.atac_fragments.tsv.gz.tbi
		mv ~{sample_id}/outs/atac_peaks.bed ~{sample_id}.atac_peaks.bed
		mv ~{sample_id}/outs/atac_cut_sites.bigwig ~{sample_id}.atac_cut_sites.bigwig
		mv ~{sample_id}/outs/atac_peak_annotation.tsv ~{sample_id}.atac_peak_annotation.tsv

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.cellranger_arc_multiome_outputs.tar.gz" \
			-o "~{sample_id}.raw_feature_bc_matrix.h5" \
			-o "~{sample_id}.filtered_feature_bc_matrix.h5" \
			-o "~{sample_id}.molecule_info.h5" \
			-o "~{sample_id}.summary.csv" \
			-o "~{sample_id}.per_barcode_metrics.csv" \
			-o "~{sample_id}.gex_possorted_bam.bam" \
			-o "~{sample_id}.gex_possorted_bam.bam.bai" \
			-o "~{sample_id}.atac_possorted_bam.bam" \
			-o "~{sample_id}.atac_possorted_bam.bam.bai" \
			-o "~{sample_id}.atac_fragments.tsv.gz" \
			-o "~{sample_id}.atac_fragments.tsv.gz.tbi" \
			-o "~{sample_id}.atac_peaks.bed" \
			-o "~{sample_id}.atac_cut_sites.bigwig" \
			-o "~{sample_id}.atac_peak_annotation.tsv"
	>>>

	output {
		String multiome_outputs_tar_gz = "~{raw_data_path}/~{sample_id}.cellranger_arc_multiome_outputs.tar.gz"
		String raw_counts = "~{raw_data_path}/~{sample_id}.raw_feature_bc_matrix.h5"
		String filtered_counts = "~{raw_data_path}/~{sample_id}.filtered_feature_bc_matrix.h5"
		String molecule_info = "~{raw_data_path}/~{sample_id}.molecule_info.h5"
		String metrics_summary_csv = "~{raw_data_path}/~{sample_id}.summary.csv"
		String per_barcode_metrics_summary_csv = "~{raw_data_path}/~{sample_id}.per_barcode_metrics.csv"
		String gex_possorted_bam = "~{raw_data_path}/~{sample_id}.gex_possorted_bam.bam"
		String gex_possorted_bam_index = "~{raw_data_path}/~{sample_id}.gex_possorted_bam.bam.bai"
		String atac_possorted_bam = "~{raw_data_path}/~{sample_id}.atac_possorted_bam.bam"
		String atac_possorted_bam_index = "~{raw_data_path}/~{sample_id}.atac_possorted_bam.bam.bai"
		String atac_fragments_tsv_gz = "~{raw_data_path}/~{sample_id}.atac_fragments.tsv.gz"
		String atac_fragments_tsv_gz_tbi = "~{raw_data_path}/~{sample_id}.atac_fragments.tsv.gz.tbi"
		String atac_peaks_bed = "~{raw_data_path}/~{sample_id}.atac_peaks.bed"
		String atac_cut_sites_bigwig = "~{raw_data_path}/~{sample_id}.atac_cut_sites.bigwig"
		String atac_peak_annotation_tsv = "~{raw_data_path}/~{sample_id}.atac_peak_annotation.tsv"
	}

	runtime {
		docker: "~{container_registry}/cellranger_arc:2.2.0"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}
}
