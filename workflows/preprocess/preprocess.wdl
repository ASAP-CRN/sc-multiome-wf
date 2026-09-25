version 1.0

# Generate a preprocessed AnnData object

import "../structs.wdl"
import "../../sc-atacseq-wf/workflows/preprocess/preprocess.wdl" as ScAtacPreprocess
	alias Project as AtacProject
	alias Sample as AtacSample

workflow preprocess {
	input {
		String team_id
		String dataset_id
		String dataset_doi_url
		Array[Sample] samples

		File cellranger_arc_reference_data
		File cellranger_atac_reference_chrom_sizes
		Float cellbender_fpr
		Array[File]? vireo_assignment_files

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
	# TODO This must match https://github.com/ASAP-CRN/sc-rnaseq-wf/blob/main/workflows/preprocess/preprocess.wdl and https://github.com/ASAP-CRN/sc-atacseq-wf/blob/main/workflows/preprocess/preprocess.wdl
	String cellbender_task_version = "1.0.1"
	String rna_demux_counts_task_version = "1.0.0"
	String rna_demux_adata_task_version = "1.0.0"
	String atac_demux_fragments_task_version = "1.0.0"
	String atac_demux_adata_task_version = "1.0.0"
	String rna_adata_task_version = "1.2.0"
	String atac_adata_task_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String cellranger_arc_raw_data_path = "~{workflow_raw_data_path_prefix}/multiome/cellranger_arc/~{cellranger_arc_task_version}"
	String rna_cellbender_raw_data_path = "~{workflow_raw_data_path_prefix}/rna/remove_technical_artifacts/~{cellbender_task_version}"
	String rna_demux_counts_raw_data_path = "~{workflow_raw_data_path_prefix}/rna/split_demux_counts/~{rna_demux_counts_task_version}"
	String rna_demux_adata_raw_data_path = "~{workflow_raw_data_path_prefix}/rna/counts_to_adata/~{rna_demux_adata_task_version}"
	String atac_demux_fragments_raw_data_path = "~{workflow_raw_data_path_prefix}/atac/split_demux_fragments/~{atac_demux_fragments_task_version}"
	String atac_demux_adata_raw_data_path = "~{workflow_raw_data_path_prefix}/atac/fragments_to_adata/~{atac_demux_adata_task_version}"
	String rna_adata_raw_data_path = "~{workflow_raw_data_path_prefix}/rna/counts_to_adata/~{rna_adata_task_version}"
	String atac_adata_raw_data_path = "~{workflow_raw_data_path_prefix}/atac/fragments_to_adata/~{atac_adata_task_version}"

	scatter (sample_object in samples) {
		Boolean is_multiplexed = length(sample_object.multiplexed_samples) > 0

		String cellranger_arc_counts_output = "~{cellranger_arc_raw_data_path}/~{dataset_id}.~{sample_object.sample_id}.raw_feature_bc_matrix.h5"
		String rna_cellbender_counts_output = "~{rna_cellbender_raw_data_path}/~{dataset_id}.~{sample_object.sample_id}.cellbender.h5"
		String rna_initial_adata_object_output = "~{rna_adata_raw_data_path}/~{dataset_id}.~{sample_object.sample_id}.cleaned_unfiltered.h5ad"
		String atac_initial_adata_object_output = "~{atac_adata_raw_data_path}/~{dataset_id}.~{sample_object.sample_id}.cleaned_unfiltered.h5ad"
	
		scatter (mux_sample_object in sample_object.multiplexed_samples) {
			String mux_dataset_sample_pool_id = "~{dataset_id}.~{mux_sample_object.sample_id}.~{sample_object.sample_id}"
			String source_subject_id = mux_sample_object.source_subject_id
			String asap_subject_id = mux_sample_object.asap_subject_id
			String mux_sample_id = mux_sample_object.sample_id
			String mux_batch = select_first([mux_sample_object.batch])
			String mux_sex = select_first([mux_sample_object.sex])
			String mux_brain_region_level_1 = select_first([mux_sample_object.brain_region_level_1, "NA"])
			String mux_brain_region_level_2 = select_first([mux_sample_object.brain_region_level_2, "NA"])
			String mux_brain_region_level_3 = select_first([mux_sample_object.brain_region_level_3, "NA"])

			String rna_demux_counts_output = if is_multiplexed then "~{rna_demux_counts_raw_data_path}/~{mux_dataset_sample_pool_id}.demux_counts.h5ad" else "NA"
			String rna_demux_initial_adata_object_output = if is_multiplexed then "~{rna_demux_adata_raw_data_path}/~{mux_dataset_sample_pool_id}.cleaned_unfiltered.h5ad" else "NA"
			String atac_demux_fragments_output = if is_multiplexed then "~{atac_demux_fragments_raw_data_path}/~{mux_dataset_sample_pool_id}.fragments.tsv.gz" else "NA"
			String atac_demux_initial_adata_object_output = if is_multiplexed then "~{atac_demux_adata_raw_data_path}/~{mux_dataset_sample_pool_id}.cleaned_unfiltered.h5ad" else "NA"
		}
	}

	call check_output_files_exist as check_cellranger_arc_output_files_exist {
		input:
			output_files = cellranger_arc_counts_output,
			billing_project = billing_project,
			zones = zones
	}

	call check_output_files_exist as check_cellbender_output_files_exist {
		input:
			output_files = rna_cellbender_counts_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (index in range(length(samples))) {
		Sample sample = samples[index]

		String dataset_sample_id = "~{dataset_id}.~{sample.sample_id}"
		Array[String] project_sample_id = [team_id, sample.sample_id, dataset_doi_url]

		String cellranger_arc_counts_complete = check_cellranger_arc_output_files_exist.sample_preprocessing_complete[index][0]
		Boolean run_rna_cellbender = if cellranger_arc_counts_complete == "false" then true else !check_cellbender_output_files_exist.all_exist

		String cellranger_arc_multiome_outputs_tar_gz = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.cellranger_arc_multiome_outputs.tar.gz"
		String cellranger_arc_raw_counts = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.raw_feature_bc_matrix.h5"
		String cellranger_arc_filtered_counts = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.filtered_feature_bc_matrix.h5"
		String cellranger_arc_gex_molecule_info = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.gex_molecule_info.h5"
		String cellranger_arc_metrics_summary_csv = "~{cellranger_arc_raw_data_path}/~{dataset_sample_id}.summary.csv"
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

		if (cellranger_arc_counts_complete == "false") {
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
		File gex_molecule_info_output = select_first([cellranger_arc_count.gex_molecule_info, cellranger_arc_gex_molecule_info]) #!FileCoercion
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

		if (run_rna_cellbender) {
			call rna_remove_technical_artifacts {
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

		if (is_multiplexed[index]) {
			call check_output_files_exist as check_rna_demux_counts_output_files_exist {
				input:
					output_files = rna_demux_counts_output[index],
					billing_project = billing_project,
					zones = zones
			}

			Boolean run_rna_demux_counts = if run_rna_cellbender then true else !check_rna_demux_counts_output_files_exist.all_exist

			call check_output_files_exist as check_atac_split_demux_fragments_output_files_exist {
				input:
					output_files = atac_demux_fragments_output[index],
					billing_project = billing_project,
					zones = zones
			}

			Boolean run_atac_split_demux_fragments = if cellranger_arc_counts_complete == "false" then true else !check_atac_split_demux_fragments_output_files_exist.all_exist

			if (run_rna_demux_counts) {
				call rna_demux_counts {
					input:
						dataset_id = dataset_id,
						pool_id = sample.sample_id,
						source_subject_ids = source_subject_id[index],
						sample_ids = mux_sample_id[index],
						cellbender_counts = removed_background_counts_output,
						vireo_assignment_files = flatten(select_all([vireo_assignment_files])),
						raw_data_path = rna_demux_counts_raw_data_path,
						workflow_info = workflow_info,
						billing_project = billing_project,
						container_registry = container_registry,
						zones = zones
				}
			}

			Array[File] sample_demux_cellbender_counts_output = select_first([rna_demux_counts.sample_demux_cellbender_counts, rna_demux_counts_output[index]]) #!FileCoercion

			if (run_atac_split_demux_fragments) {
				call ScAtacPreprocess.split_demux_fragments as atac_split_demux_fragments {
					input:
						dataset_id = dataset_id,
						pool_id = sample.sample_id,
						source_subject_ids = source_subject_id[index],
						sample_ids = mux_sample_id[index],
						cellranger_atac_fragments = cellranger_arc_atac_fragments_tsv_gz, #!FileCoercion
						cellranger_atac_reference_chrom_sizes = cellranger_atac_reference_chrom_sizes,
						vireo_assignment_files = flatten(select_all([vireo_assignment_files])),
						raw_data_path = atac_demux_fragments_raw_data_path,
						workflow_info = workflow_info,
						billing_project = billing_project,
						container_registry = container_registry,
						zones = zones
				}
			}

			Array[File] sample_demux_fragments_tsv_gz_output = select_first([atac_split_demux_fragments.sample_split_fragments_tsv_gz, atac_demux_fragments_output[index]]) #!FileCoercion

			call check_output_files_exist as check_rna_demux_adata_output_files_exist {
				input:
					output_files = rna_demux_initial_adata_object_output[index],
					billing_project = billing_project,
					zones = zones
			}

			Boolean run_rna_demux_adata = if run_rna_demux_counts then true else !check_rna_demux_adata_output_files_exist.all_exist

			call check_output_files_exist as check_atac_demux_adata_output_files_exist {
				input:
					output_files = atac_demux_initial_adata_object_output[index],
					billing_project = billing_project,
					zones = zones
			}

			Boolean run_atac_demux_adata = if run_atac_split_demux_fragments then true else !check_atac_demux_adata_output_files_exist.all_exist

			scatter (mux_sample_index in range(length(sample.multiplexed_samples))) {
				if (run_rna_demux_adata) {
					call counts_to_adata as rna_demux_counts_to_adata {
						input:
							team_id = team_id,
							dataset_id = dataset_id,
							sample_id = mux_sample_id[index][mux_sample_index],
							batch = mux_batch[index][mux_sample_index],
							sex = mux_sex[index][mux_sample_index],
							brain_region_level_1 = mux_brain_region_level_1[index][mux_sample_index],
							brain_region_level_2 = mux_brain_region_level_2[index][mux_sample_index],
							brain_region_level_3 = mux_brain_region_level_3[index][mux_sample_index],
							pool_id = sample.sample_id,
							source_subject_id = source_subject_id[index][mux_sample_index],
							subject_id = asap_subject_id[index][mux_sample_index],
							cellbender_counts = sample_demux_cellbender_counts_output[mux_sample_index],
							raw_data_path = rna_demux_adata_raw_data_path,
							workflow_info = workflow_info,
							billing_project = billing_project,
							container_registry = container_registry,
							zones = zones
					}
				}

				File preprocessed_rna_demux_adata_object_output = select_first([rna_demux_counts_to_adata.initial_adata_object, rna_demux_initial_adata_object_output[index][mux_sample_index]]) #!FileCoercion

				if (run_atac_demux_adata) {
					call fragments_to_adata as atac_demux_fragments_to_adata {
						input:
							team_id = team_id,
							dataset_id = dataset_id,
							sample_id = mux_sample_id[index][mux_sample_index],
							batch = mux_batch[index][mux_sample_index],
							sex = mux_sex[index][mux_sample_index],
							brain_region_level_1 = mux_brain_region_level_1[index][mux_sample_index],
							brain_region_level_2 = mux_brain_region_level_2[index][mux_sample_index],
							brain_region_level_3 = mux_brain_region_level_3[index][mux_sample_index],
							pool_id = sample.sample_id,
							source_subject_id = source_subject_id[index][mux_sample_index],
							subject_id = asap_subject_id[index][mux_sample_index],
							fragments_tsv_gz = sample_demux_fragments_tsv_gz_output[mux_sample_index],
							raw_data_path = atac_demux_adata_raw_data_path,
							workflow_info = workflow_info,
							billing_project = billing_project,
							container_registry = container_registry,
							zones = zones
					}
				}

				File preprocessed_atac_demux_adata_object_output = select_first([atac_demux_fragments_to_adata.initial_adata_object, atac_demux_initial_adata_object_output[index][mux_sample_index]]) #!FileCoercion
			}
		}

		if (! is_multiplexed[index]) {
			call check_output_files_exist as check_rna_adata_output_files_exist {
				input:
					output_files = rna_initial_adata_object_output,
					billing_project = billing_project,
					zones = zones
			}

			String rna_adata_complete = check_rna_adata_output_files_exist.sample_preprocessing_complete[index][0]
			Boolean run_rna_adata = if run_rna_cellbender then true else (rna_adata_complete == "false")

			call check_output_files_exist as check_atac_adata_output_files_exist {
				input:
					output_files = atac_initial_adata_object_output,
					billing_project = billing_project,
					zones = zones
			}

			Boolean run_atac_adata = if cellranger_arc_counts_complete == "false" then true else !check_atac_adata_output_files_exist.all_exist

			if (run_rna_adata) {
				call counts_to_adata as rna_counts_to_adata {
					input:
						team_id = team_id,
						dataset_id = dataset_id,
						sample_id = sample.sample_id,
						batch = select_first([sample.batch]),
						sex = select_first([sample.sex]),
						brain_region_level_1 = select_first([sample.brain_region_level_1, "NA"]),
						brain_region_level_2 = select_first([sample.brain_region_level_2, "NA"]),
						brain_region_level_3 = select_first([sample.brain_region_level_3, "NA"]),
						cellbender_counts = removed_background_counts_output,
						raw_data_path = rna_adata_raw_data_path,
						workflow_info = workflow_info,
						billing_project = billing_project,
						container_registry = container_registry,
						zones = zones
				}
			}

			File preprocessed_rna_adata_object_output = select_first([rna_counts_to_adata.initial_adata_object, rna_initial_adata_object_output[index]]) #!FileCoercion
		
			if (run_atac_adata) {
				call fragments_to_adata as atac_fragments_to_adata {
					input:
						team_id = team_id,
						dataset_id = dataset_id,
						sample_id = sample.sample_id,
						batch = select_first([sample.batch]),
						sex = select_first([sample.sex]),
						brain_region_level_1 = select_first([sample.brain_region_level_1, "NA"]),
						brain_region_level_2 = select_first([sample.brain_region_level_2, "NA"]),
						brain_region_level_3 = select_first([sample.brain_region_level_3, "NA"]),
						fragments_tsv_gz = atac_fragments_tsv_gz_output,
						raw_data_path = atac_adata_raw_data_path,
						workflow_info = workflow_info,
						billing_project = billing_project,
						container_registry = container_registry,
						zones = zones
				}
			}

			File preprocessed_atac_adata_object_output = select_first([atac_fragments_to_adata.initial_adata_object, atac_initial_adata_object_output[index]]) #!FileCoercion
			
		}
	}

	output {
		# Sample list
		Array[Array[String]] project_sample_ids = project_sample_id

		# Cellranger ARC
		Array[File] multiome_outputs_tar_gz = multiome_outputs_tar_gz_output #!FileCoercion
		Array[File] raw_counts = raw_counts_output #!FileCoercion
		Array[File] filtered_counts = filtered_counts_output #!FileCoercion
		Array[File] gex_molecule_info = gex_molecule_info_output #!FileCoercion
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

		# RNA - Demultiplexed counts and AnnData object per-sample data
		Array[File] rna_subject_demux_cellbender_counts = flatten(select_all(rna_demux_counts.subject_demux_cellbender_counts))
		Array[File] rna_sample_demux_cellbender_counts = flatten(select_all(sample_demux_cellbender_counts_output)) #!FileCoercion
		Array[File] rna_initial_demux_adata_object = flatten(select_all(preprocessed_rna_demux_adata_object_output)) #!FileCoercion

		# ATAC - Demultiplexed fragments and AnnData object per-sample data
		Array[File] atac_subject_split_fragments_tsv_gz = flatten(select_all(atac_split_demux_fragments.subject_split_fragments_tsv_gz))
		Array[File] atac_sample_split_fragments_tsv_gz = flatten(select_all(sample_demux_fragments_tsv_gz_output)) #!FileCoercion
		Array[File] atac_initial_demux_adata_object = flatten(select_all(preprocessed_atac_demux_adata_object_output)) #!FileCoercion

		# RNA - AnnData object per-sample data
		Array[File] rna_initial_adata_object = select_all(preprocessed_rna_adata_object_output) #!FileCoercion

		# ATAC - AnnData object per-sample data
		Array[File] atac_initial_adata_object = select_all(preprocessed_atac_adata_object_output) #!FileCoercion
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
		cellranger_atac_reference_chrom_sizes: {help: "Chromosome sizes file (.chrom.sizes or .fa.fai) from the Cell Ranger ATAC reference (same as ARC), used to validate fragment coordinates during splitting."}
		cellbender_fpr: {help: "Cellbender false positive rate. [0.0]"}
		vireo_assignment_files: {help: "Vireo donor assignment CSV with columns: donor_id, sample, raw_bc. Covers all donors in the pool."}
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
		Array[String] output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		all_exist="true"

		while read -r file || [[ -n "${file}" ]]; do
			if [[ "${file}" == "NA" ]] || gcloud storage ls --billing-project=~{billing_project} "${file}"; then
				echo -e "true" >> sample_preprocessing_complete.tsv
			else
				echo -e "false" >> sample_preprocessing_complete.tsv
				all_exist="false"
			fi
		done < ~{write_lines(output_files)}

		echo "${all_exist}" > all_exist.txt
	>>>

	output {
		Array[Array[String]] sample_preprocessing_complete = read_tsv("sample_preprocessing_complete.tsv")
		Boolean all_exist = read_boolean("all_exist.txt")
	}

	runtime {
		docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:584.0.0-slim"
		cpu: 2
		cpuPlatform: "Intel Cascade Lake"
		memory: "4 GB"
		disks: "local-disk 20 HDD"
		preemptible: 3
		zones: zones
	}

	meta {
		description: "Checks for existing preprocessing files per pool or sample and skips certain preprocessing steps if they exist."
	}

	parameter_meta {
		output_files: {help: "Output file to detect."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
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
		mv ~{sample_id}/outs/gex_molecule_info.h5 ~{dataset_sample_id}.gex_molecule_info.h5
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
			-o "~{dataset_sample_id}.gex_molecule_info.h5" \
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
		String gex_molecule_info = "~{raw_data_path}/~{dataset_sample_id}.gex_molecule_info.h5"
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

task rna_remove_technical_artifacts {
	input {
		String dataset_sample_id

		File raw_counts

		Float cellbender_fpr

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int disk_size = ceil(size(raw_counts, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		nvidia-smi

		export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:512

		# The cellranger-arc raw matrix holds both Gene Expression and Peaks features
		/usr/bin/time \
		cellbender remove-background \
			--cuda \
			--input ~{raw_counts} \
			--output ~{dataset_sample_id}.cellbender. \
			--exclude-feature-types Peaks \
			--fpr ~{cellbender_fpr}

		mv ckpt.tar.gz "~{dataset_sample_id}.cellbender_ckpt.tar.gz"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{dataset_sample_id}.cellbender_report.html" \
			-o "~{dataset_sample_id}.cellbender.h5" \
			-o "~{dataset_sample_id}.cellbender_filtered.h5" \
			-o "~{dataset_sample_id}.cellbender_cell_barcodes.csv" \
			-o "~{dataset_sample_id}.cellbender.pdf" \
			-o "~{dataset_sample_id}.cellbender.log" \
			-o "~{dataset_sample_id}.cellbender_metrics.csv" \
			-o "~{dataset_sample_id}.cellbend_posterior.h5"
	>>>

	output {
		String report_html = "~{raw_data_path}/~{dataset_sample_id}.cellbender_report.html"
		String removed_background_counts = "~{raw_data_path}/~{dataset_sample_id}.cellbender.h5"
		String filtered_removed_background_counts = "~{raw_data_path}/~{dataset_sample_id}.cellbender_filtered.h5"
		String cell_barcodes_csv = "~{raw_data_path}/~{dataset_sample_id}.cellbender_cell_barcodes.csv"
		String graph_pdf = "~{raw_data_path}/~{dataset_sample_id}.cellbender.pdf"
		String log = "~{raw_data_path}/~{dataset_sample_id}.cellbender.log"
		String metrics_csv = "~{raw_data_path}/~{dataset_sample_id}.cellbender_metrics.csv"
		String posterior_probability = "~{raw_data_path}/~{dataset_sample_id}.cellbend_posterior.h5"
	}

	runtime {
		docker: "~{container_registry}/cellbender:0.3.0"
		cpu: 4
		memory: "64 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
		gpuType: "nvidia-tesla-t4"
		gpuCount: 1
	}

	meta {
		description: "Removes ambient RNA and technical noise from the Gene Expression features of raw cellranger-arc count matrices using CellBender."
	}

	parameter_meta {
		dataset_sample_id: {help: "Generated ASAP dataset ID and sample ID; stored in the AnnData objects."}
		raw_counts: {help: "Raw cellranger-arc count matrix (H5 format) to process; Peaks features are excluded and only Gene Expression features are processed."}
		cellbender_fpr: {help: "Cellbender false positive rate. [0.0]"}
		raw_data_path: {help: "Raw data bucket path for cellbender outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess/cellbender/<cellbender_task_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}

task rna_demux_counts {
	input {
		String dataset_id
		String pool_id
		Array[String] source_subject_ids
		Array[String] sample_ids

		File cellbender_counts
		Array[File] vireo_assignment_files

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int disk_size = ceil(size([cellbender_counts], "GB") + size(vireo_assignment_files, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		# Generate mapping of samples (pool), cell types (subject), and cell barcodes TSV
		## Filter out 'doublet' and 'unassigned' cell types (subject)
		echo -e "sample\tcell_type\tcell_barcode" > cell_barcodes.tsv
		for f in ~{sep=' ' vireo_assignment_files}; do
			awk -F ',' -v pool="~{pool_id}" '$2 ~ /^ASA/ && $3 ~ pool { OFS="\t"; print $3, $2, $4 }' "${f}" >> cell_barcodes.tsv
		done

		paste \
			~{write_lines(source_subject_ids)} \
			~{write_lines(sample_ids)} \
		> metadata.tsv

		duplicates=$(cut -f2 metadata.tsv | sort | uniq -d)
		if [[ -n "${duplicates}" ]]; then
			echo "[ERROR] Duplicate sample_id found for pool ~{pool_id}: ${duplicates}" >&2
			exit 1
		fi

		/usr/bin/time -v \
		demux_cellbender_counts \
			--cellbender-counts ~{cellbender_counts} \
			--cell-barcodes cell_barcodes.tsv \
			--metadata metadata.tsv \
			--pool-id ~{pool_id} \
			--dataset-id ~{dataset_id}

		files_to_upload=(*~{pool_id}.demux_counts.h5ad)
		upload_args=()
		for f in "${files_to_upload[@]}"; do
			upload_args+=(-o "${f}")
		done

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			"${upload_args[@]}"

		for f in "${files_to_upload[@]}"; do
			echo "~{raw_data_path}/${f}"
		done > sample_cellbender_counts_filenames.txt
	>>>

	output {
		Array[File] subject_demux_cellbender_counts = glob("*~{pool_id}.demux_counts.h5ad")
		Array[String] sample_demux_cellbender_counts = read_lines("sample_cellbender_counts_filenames.txt")
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
		description: "Demultiplex Cellbender counts and export as AnnData objects"
	}

	parameter_meta {
		dataset_id: {help: "Generated ASAP dataset ID; stored in the AnnData objects."}
		pool_id: {help: "Generated ASAP pool ID; used to name output files."}
		source_subject_ids: {help: "An array of provided source subject IDs; used for mapping."}
		sample_ids: {help: "An array of generated ASAP sample ID; used for mapping."}
		cellbender_counts: {help: "CellBender-cleaned count matrix (H5 format)."}
		vireo_assignment_files: {help: "Vireo donor assignment CSV with columns: donor_id, sample, raw_bc. Covers all donors in the pool."}
		raw_data_path: {help: "Raw data bucket path for counts to adata outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess/rna/demux/<rna_demux_task_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}

task counts_to_adata {
	input {
		String team_id
		String dataset_id
		String sample_id
		String batch
		String sex
		String brain_region_level_1
		String brain_region_level_2
		String brain_region_level_3
		String? pool_id
		String? source_subject_id
		String? subject_id

		File cellbender_counts

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	String adata_output_filename = if defined(pool_id) then "~{dataset_id}.~{sample_id}.~{pool_id}.cleaned_unfiltered.h5ad" else "~{dataset_id}.~{sample_id}.cleaned_unfiltered.h5ad"

	Int disk_size = ceil(size(cellbender_counts, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		/usr/bin/time -v \
		counts_to_adata \
			--cellbender-counts ~{cellbender_counts} \
			--team ~{team_id} \
			--dataset ~{dataset_id} \
			--sample-id ~{sample_id} \
			--batch ~{batch} \
			--sex ~{sex} \
			--brain-region-level-1 "~{brain_region_level_1}" \
			--brain-region-level-2 "~{brain_region_level_2}" \
			--brain-region-level-3 "~{brain_region_level_3}" \
			~{"--pool-id " + pool_id} \
			~{"--source-subject-id " + source_subject_id} \
			~{"--subject-id " + subject_id} \
			--adata-output ~{adata_output_filename}

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o ~{adata_output_filename}
	>>>

	output {
		String initial_adata_object = "~{raw_data_path}/~{adata_output_filename}"
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
		description: "Converts Cell Ranger ARC counts into AnnData objects using scanpy for demultiplexed samples and one per sample otherwise."
	}

	parameter_meta {
		team_id: {help: "Name of the CRN Team; stored in the AnnData objects."}
		dataset_id: {help: "Generated ASAP dataset ID; stored in the AnnData objects."}
		sample_id: {help: "Generated ASAP sample ID; stored in the AnnData objects and used to name output files."}
		batch: {help: "The sample's batch; stored in the AnnData objects."}
		sex: {help: "The sample's sex; stored in the AnnData objects."}
		brain_region_level_1: {help: "Abbreviation of most granular anatomical region (Level 1)."}
		brain_region_level_2: {help: "Abbreviation of intermediate level anatomical region (Level 2)."}
		brain_region_level_3: {help: "Abbreviation of coarse level anatomical region (Level 3)."}
		pool_id: {help: "Generated ASAP pool ID; stored in the AnnData objects."}
		source_subject_id: {help: "Source subject ID; stored in the AnnData objects."}
		subject_id: {help: "Generated ASAP subject ID; stored in the AnnData objects."}
		cellbender_counts: {help: "Cellbender counts (.h5) or demultiplexed Cellbender counts (.h5ad)."}
		raw_data_path: {help: "Raw data bucket path for counts to adata outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess/rna/counts_to_adata/<rna_adata_task_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}

task fragments_to_adata {
	input {
		String team_id
		String dataset_id
		String sample_id
		String batch
		String sex
		String brain_region_level_1
		String brain_region_level_2
		String brain_region_level_3
		String? pool_id
		String? source_subject_id
		String? subject_id

		File fragments_tsv_gz

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	String adata_output_filename = if defined(pool_id) then "~{dataset_id}.~{sample_id}.~{pool_id}.cleaned_unfiltered.h5ad" else "~{dataset_id}.~{sample_id}.cleaned_unfiltered.h5ad"

	Int disk_size = ceil(size(fragments_tsv_gz, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		/usr/bin/time -v \
		fragments_to_adata \
			--cellranger-arc-fragments ~{fragments_tsv_gz} \
			--team ~{team_id} \
			--dataset ~{dataset_id} \
			--sample-id ~{sample_id} \
			--batch ~{batch} \
			--sex ~{sex} \
			--brain-region-level-1 "~{brain_region_level_1}" \
			--brain-region-level-2 "~{brain_region_level_2}" \
			--brain-region-level-3 "~{brain_region_level_3}" \
			~{"--pool-id " + pool_id} \
			~{"--source-subject-id " + source_subject_id} \
			~{"--subject-id " + subject_id} \
			--adata-output ~{adata_output_filename}

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o ~{adata_output_filename}
	>>>

	output {
		String initial_adata_object = "~{raw_data_path}/~{adata_output_filename}"
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
		description: "Converts Cell Ranger ARC fragments into AnnData objects using SnapATAC2 for demultiplexed samples and one per sample otherwise."
	}

	parameter_meta {
		team_id: {help: "Name of the CRN Team; stored in the AnnData objects."}
		dataset_id: {help: "Generated ASAP dataset ID; stored in the AnnData objects."}
		sample_id: {help: "Generated ASAP sample ID; stored in the AnnData objects and used to name output files."}
		batch: {help: "The sample's batch; stored in the AnnData objects."}
		sex: {help: "The sample's sex; stored in the AnnData objects."}
		brain_region_level_1: {help: "Abbreviation of most granular anatomical region (Level 1)."}
		brain_region_level_2: {help: "Abbreviation of intermediate level anatomical region (Level 2)."}
		brain_region_level_3: {help: "Abbreviation of coarse level anatomical region (Level 3)."}
		pool_id: {help: "Generated ASAP pool ID; stored in the AnnData objects."}
		source_subject_id: {help: "Source subject ID; stored in the AnnData objects."}
		subject_id: {help: "Generated ASAP subject ID; stored in the AnnData objects."}
		fragments_tsv_gz: {help: "Count and barcode information for every ATAC fragment observed in the experiment in TSV format output by Cell Ranger ARC."}
		raw_data_path: {help: "Raw data bucket path for counts to adata outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/preprocess/atac/fragments_to_adata/<atac_adata_task_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}
