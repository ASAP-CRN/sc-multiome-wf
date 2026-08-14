version 1.0

# Harmonized human PMDBS brain sc/sn Multiome (RNA + ATAC) workflow entrypoint

import "../wf-common/wdl/structs.wdl"
import "../wf-common/wdl/tasks/get_workflow_metadata.wdl" as GetWorkflowMetadata
import "preprocess/preprocess.wdl" as Preprocess
import "../sc-rnaseq-wf/workflows/cohort_analysis/cohort_analysis.wdl" as ScRnaCohortAnalysis
import "../sc-atacseq-wf/workflows/cohort_analysis/cohort_analysis.wdl" as ScAtacCohortAnalysis
import "cohort_analysis/cohort_analysis.wdl" as MultiomeCohortAnalysis

workflow sc_multiome_analysis {
	input {
		String cohort_id
		Array[Project] projects

		# Preprocess
		File cellranger_arc_reference_data
		Float cellbender_fpr = 0.0

		# Filtering parameters
		Int pct_counts_mt_max = 10
		Float doublet_score_max = 0.2
		Array[Int] total_counts_limits = [100, 100000]
		Array[Int] n_genes_by_counts_limits = [100, 10000]

		# Allen Institute's Map My Cells
		File allen_brain_mmc_precomputed_stats_h5

		# Normalization parameters
		Int norm_target_sum = 10000 
		Int n_top_genes = 3000
		Int n_comps = 30

		# RNA - Sample integration
		String batch_key = "batch_id"
		String scvi_latent_key = "X_scVI"
		String scanvi_latent_key = "X_scANVI"
		String scanvi_predictions_key = "C_scANVI"

		# RNA and Multiome - Clustering parameters
		Int n_neighbors = 15
		Array[Float] leiden_res = [0.2, 0.5, 1.0]

		# RNA - Groups and features
		Array[String] rna_groups = ["sample", "batch", "cell_type", "leiden_res_0.05", "leiden_res_0.10", "leiden_res_0.20", "leiden_res_0.40"]
		Array[String] rna_features = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_rb", "doublet_score", "S_score", "G2M_score"]

		# ATAC - Sample integration
		String peakvi_latent_key = "X_peakVI"
		Int peakvi_max_epochs = 300
		Int peakvi_n_hidden = 128
		Int peakvi_batch_size = 128

		# ATAC - Groups and features
		Array[String] atac_groups = ["sample", "batch", "team", "dataset", "batch_id", "leiden"]
		Array[String] atac_features = ["n_fragment", "tsse", "frac_dup", "frac_mito", "doublet_score", "doublet_probability"]

		# Multiome - QC
		Float min_overlap_frac

		# Multiome - Sample integration
		String multivi_latent_key = "X_multiVI"
		Int multivi_max_epochs = 300
		Int multivi_n_hidden = 128
		Int multivi_batch_size = 128

		# Multiome - Cohort analysis
		Boolean peak_gene_only = true
		Float min_linkage_corr = 0.1
		Float min_linkage_sig = 0.0

		# Cohort analysis
		Boolean run_cross_team_cohort_analysis = false
		String cohort_raw_data_bucket
		Array[String] cohort_staging_data_buckets

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

	String workflow_execution_path = "workflow_execution"
	String workflow_name = "pmdbs_sc_multiome"
	String workflow_version = "v1.0.0"
	String workflow_release = "https://github.com/ASAP-CRN/sc-multiome-wf/releases/tag/sc_multiome_analysis-~{workflow_version}"
	String crn_release_version = "v6.0.0"

	call GetWorkflowMetadata.get_workflow_metadata {
		input:
			zones = zones
	}

	scatter (project in projects) {
		String project_raw_data_path_prefix = "~{project.raw_data_bucket}/~{workflow_execution_path}/~{workflow_name}"

		call Preprocess.preprocess {
			input:
				team_id = project.asap_team_id,
				dataset_id = project.asap_dataset_id,
				dataset_doi_url = project.asap_dataset_doi_url,
				samples = project.samples,
				cellranger_arc_reference_data = cellranger_arc_reference_data,
				cellbender_fpr = cellbender_fpr,
				workflow_name = workflow_name,
				workflow_version = workflow_version,
				workflow_release = workflow_release,
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = project_raw_data_path_prefix,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}

		Array[String] multiome_preprocessing_output_file_paths = flatten([
			preprocess.multiome_outputs_tar_gz,
			preprocess.raw_counts,
			preprocess.filtered_counts,
			preprocess.molecule_info,
			preprocess.metrics_summary_csv,
			preprocess.per_barcode_metrics_summary_csv,
			preprocess.gex_possorted_bam,
			preprocess.gex_possorted_bam_index,
			preprocess.atac_possorted_bam,
			preprocess.atac_possorted_bam_index,
			preprocess.atac_fragments_tsv_gz,
			preprocess.atac_fragments_tsv_gz_tbi,
			preprocess.atac_peaks_bed,
			preprocess.atac_cut_sites_bigwig,
			preprocess.atac_peak_annotation_tsv,
			preprocess.feature_linkage_bedpe
		]) #!StringCoercion

		Array[String] rna_preprocessing_output_file_paths = flatten([
			preprocess.report_html,
			preprocess.removed_background_counts,
			preprocess.filtered_removed_background_counts,
			preprocess.cell_barcodes_csv,
			preprocess.graph_pdf,
			preprocess.log,
			preprocess.metrics_csv,
			preprocess.posterior_probability,
			preprocess.rna_initial_adata_object
		]) #!StringCoercion

		Array[String] atac_preprocessing_output_file_paths = flatten([
			preprocess.atac_initial_adata_object
		]) #!StringCoercion

		if (project.run_project_cohort_analysis) {
			call ScRnaCohortAnalysis.cohort_analysis as rna_project_cohort_analysis {
				input:
					cohort_id = project.asap_team_id,
					project_sample_ids = preprocess.project_sample_ids,
					preprocessed_adata_objects = preprocess.rna_initial_adata_object,
					preprocessing_output_file_paths = rna_preprocessing_output_file_paths,
					pct_counts_mt_max = pct_counts_mt_max,
					doublet_score_max = doublet_score_max,
					total_counts_limits = total_counts_limits,
					n_genes_by_counts_limits = n_genes_by_counts_limits,
					allen_brain_mmc_precomputed_stats_h5 = allen_brain_mmc_precomputed_stats_h5,
					norm_target_sum = norm_target_sum,
					n_top_genes = n_top_genes,
					n_comps = n_comps,
					scvi_latent_key = scvi_latent_key,
					scanvi_latent_key = scanvi_latent_key,
					scanvi_predictions_key = scanvi_predictions_key,
					batch_key = batch_key,
					n_neighbors = n_neighbors,
					leiden_res = leiden_res,
					groups = rna_groups,
					features = rna_features,
					workflow_name = workflow_name,
					workflow_version = workflow_version,
					workflow_release = workflow_release,
					crn_release_version = "~{crn_release_version}/rna",
					run_timestamp = get_workflow_metadata.timestamp,
					raw_data_path_prefix = project_raw_data_path_prefix,
					staging_data_buckets = project.staging_data_buckets,
					billing_project = get_workflow_metadata.billing_project,
					container_registry = container_registry,
					zones = zones
			}

			call ScAtacCohortAnalysis.cohort_analysis as atac_project_cohort_analysis {
				input:
					cohort_id = project.asap_team_id,
					project_sample_ids = preprocess.project_sample_ids,
					preprocessed_adata_objects = preprocess.atac_initial_adata_object,
					preprocessing_output_file_paths = atac_preprocessing_output_file_paths,
					allen_brain_mmc_precomputed_stats_h5 = allen_brain_mmc_precomputed_stats_h5,
					n_top_genes = n_top_genes,
					n_comps = n_comps,
					batch_key = batch_key,
					peakvi_latent_key = peakvi_latent_key,
					peakvi_max_epochs = peakvi_max_epochs,
					peakvi_n_hidden = peakvi_n_hidden,
					peakvi_batch_size = peakvi_batch_size,
					groups = atac_groups,
					features = atac_features,
					workflow_name = workflow_name,
					workflow_version = workflow_version,
					workflow_release = workflow_release,
					crn_release_version = "~{crn_release_version}/atac",
					run_timestamp = get_workflow_metadata.timestamp,
					raw_data_path_prefix = project_raw_data_path_prefix,
					staging_data_buckets = project.staging_data_buckets,
					billing_project = get_workflow_metadata.billing_project,
					container_registry = container_registry,
					zones = zones
			}

			call MultiomeCohortAnalysis.cohort_analysis as multiome_project_cohort_analysis {
				input:
					cohort_id = project.asap_team_id,
					project_sample_ids = preprocess.project_sample_ids,
					preprocessing_output_file_paths = multiome_preprocessing_output_file_paths,
					cellranger_arc_feature_linkage_bedpe = preprocess.feature_linkage_bedpe,
					filtered_rna_adata_object = rna_project_cohort_analysis.filtered_adata_object,
					filtered_atac_adata_object = atac_project_cohort_analysis.merged_adata_object,
					labeled_cells_rna_adata_object = rna_project_cohort_analysis.labeled_cells_adata_object,
					min_overlap_frac = min_overlap_frac,
					scanvi_predictions_key = scanvi_predictions_key,
					batch_key = batch_key,
					multivi_latent_key = multivi_latent_key,
					multivi_max_epochs = multivi_max_epochs,
					multivi_n_hidden = multivi_n_hidden,
					multivi_batch_size = multivi_batch_size,
					n_neighbors = n_neighbors,
					leiden_res = leiden_res,
					peak_gene_only = peak_gene_only,
					min_linkage_corr = min_linkage_corr,
					min_linkage_sig = min_linkage_sig,
					workflow_name = workflow_name,
					workflow_version = workflow_version,
					workflow_release = workflow_release,
					crn_release_version = "~{crn_release_version}/multiome",
					run_timestamp = get_workflow_metadata.timestamp,
					raw_data_path_prefix = project_raw_data_path_prefix,
					staging_data_buckets = project.staging_data_buckets,
					billing_project = get_workflow_metadata.billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}
	}

	if (run_cross_team_cohort_analysis) {
		String cohort_raw_data_path_prefix = "~{cohort_raw_data_bucket}/~{workflow_execution_path}/~{workflow_name}"

		call ScRnaCohortAnalysis.cohort_analysis as rna_cross_team_cohort_analysis {
			input:
				cohort_id = cohort_id,
				project_sample_ids = flatten(preprocess.project_sample_ids),
				preprocessed_adata_objects = flatten(preprocess.rna_initial_adata_object),
				preprocessing_output_file_paths = flatten(rna_preprocessing_output_file_paths),
				pct_counts_mt_max = pct_counts_mt_max,
				doublet_score_max = doublet_score_max,
				total_counts_limits = total_counts_limits,
				n_genes_by_counts_limits = n_genes_by_counts_limits,
				allen_brain_mmc_precomputed_stats_h5 = allen_brain_mmc_precomputed_stats_h5,
				norm_target_sum = norm_target_sum,
				n_top_genes = n_top_genes,
				n_comps = n_comps,
				scvi_latent_key =scvi_latent_key,
				scanvi_latent_key = scanvi_latent_key,
				scanvi_predictions_key = scanvi_predictions_key,
				batch_key = batch_key,
				n_neighbors = n_neighbors,
				leiden_res = leiden_res,
				groups = rna_groups,
				features = rna_features,
				workflow_name = workflow_name,
				workflow_version = workflow_version,
				workflow_release = workflow_release,
				crn_release_version = "~{crn_release_version}/rna",
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = cohort_raw_data_path_prefix,
				staging_data_buckets = cohort_staging_data_buckets,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}

		call ScAtacCohortAnalysis.cohort_analysis as atac_cross_team_cohort_analysis {
			input:
				cohort_id = cohort_id,
				project_sample_ids = flatten(preprocess.project_sample_ids),
				preprocessed_adata_objects = flatten(preprocess.atac_initial_adata_object),
				preprocessing_output_file_paths = flatten(atac_preprocessing_output_file_paths),
				allen_brain_mmc_precomputed_stats_h5 = allen_brain_mmc_precomputed_stats_h5,
				n_top_genes = n_top_genes,
				n_comps = n_comps,
				batch_key = batch_key,
				peakvi_latent_key = peakvi_latent_key,
				peakvi_max_epochs = peakvi_max_epochs,
				peakvi_n_hidden = peakvi_n_hidden,
				peakvi_batch_size = peakvi_batch_size,
				groups = atac_groups,
				features = atac_features,
				workflow_name = workflow_name,
				workflow_version = workflow_version,
				workflow_release = workflow_release,
				crn_release_version = "~{crn_release_version}/atac",
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = cohort_raw_data_path_prefix,
				staging_data_buckets = cohort_staging_data_buckets,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}

		call MultiomeCohortAnalysis.cohort_analysis as multiome_cross_team_cohort_analysis {
			input:
				cohort_id = cohort_id,
				project_sample_ids = flatten(preprocess.project_sample_ids),
				preprocessing_output_file_paths = flatten(multiome_preprocessing_output_file_paths),
				cellranger_arc_feature_linkage_bedpe = flatten(preprocess.feature_linkage_bedpe),
				filtered_rna_adata_object = rna_cross_team_cohort_analysis.filtered_adata_object,
				filtered_atac_adata_object = atac_cross_team_cohort_analysis.merged_adata_object,
				labeled_cells_rna_adata_object = rna_cross_team_cohort_analysis.labeled_cells_adata_object,
				min_overlap_frac = min_overlap_frac,
				scanvi_predictions_key = scanvi_predictions_key,
				batch_key = batch_key,
				multivi_latent_key = multivi_latent_key,
				multivi_max_epochs = multivi_max_epochs,
				multivi_n_hidden = multivi_n_hidden,
				multivi_batch_size = multivi_batch_size,
				n_neighbors = n_neighbors,
				leiden_res = leiden_res,
				peak_gene_only = peak_gene_only,
				min_linkage_corr = min_linkage_corr,
				min_linkage_sig = min_linkage_sig,
				workflow_name = workflow_name,
				workflow_version = workflow_version,
				workflow_release = workflow_release,
				crn_release_version = "~{crn_release_version}/multiome",
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = cohort_raw_data_path_prefix,
				staging_data_buckets = cohort_staging_data_buckets,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}
	}

	output {
		# Sample-level outputs
		# Sample list
		Array[Array[Array[String]]] project_sample_ids = preprocess.project_sample_ids

		# Cell Ranger ARC
		Array[Array[File]] cellranger_arc_multiome_outputs_tar_gz = preprocess.multiome_outputs_tar_gz
		Array[Array[File]] cellranger_arc_raw_counts = preprocess.raw_counts
		Array[Array[File]] cellranger_arc_filtered_counts = preprocess.filtered_counts
		Array[Array[File]] cellranger_arc_molecule_info = preprocess.molecule_info
		Array[Array[File]] cellranger_arc_metrics_summary_csv = preprocess.metrics_summary_csv
		Array[Array[File]] cellranger_arc_per_barcode_metrics_summary_csv = preprocess.per_barcode_metrics_summary_csv
		Array[Array[File]] cellranger_arc_gex_possorted_bam = preprocess.gex_possorted_bam
		Array[Array[File]] cellranger_arc_gex_possorted_bam_index = preprocess.gex_possorted_bam_index
		Array[Array[File]] cellranger_arc_atac_possorted_bam = preprocess.atac_possorted_bam
		Array[Array[File]] cellranger_arc_atac_possorted_bam_index = preprocess.atac_possorted_bam_index
		Array[Array[File]] cellranger_arc_atac_fragments_tsv_gz = preprocess.atac_fragments_tsv_gz
		Array[Array[File]] cellranger_arc_atac_fragments_tsv_gz_tbi = preprocess.atac_fragments_tsv_gz_tbi
		Array[Array[File]] cellranger_arc_atac_peaks_bed = preprocess.atac_peaks_bed
		Array[Array[File]] cellranger_arc_atac_cut_sites_bigwig = preprocess.atac_cut_sites_bigwig
		Array[Array[File]] cellranger_arc_atac_peak_annotation_tsv = preprocess.atac_peak_annotation_tsv
		Array[Array[File]] cellranger_arc_feature_linkage_bedpe = preprocess.feature_linkage_bedpe

		# Preprocess
		Array[Array[File]] cellbender_report_html = preprocess.report_html
		Array[Array[File]] cellbender_removed_background_counts = preprocess.removed_background_counts
		Array[Array[File]] cellbender_filtered_removed_background_counts = preprocess.filtered_removed_background_counts
		Array[Array[File]] cellbender_cell_barcodes_csv = preprocess.cell_barcodes_csv
		Array[Array[File]] cellbender_graph_pdf = preprocess.graph_pdf
		Array[Array[File]] cellbender_log = preprocess.log
		Array[Array[File]] cellbender_metrics_csv = preprocess.metrics_csv
		Array[Array[File]] cellbender_posterior_probability = preprocess.posterior_probability
		Array[Array[File]] rna_initial_adata_object = preprocess.rna_initial_adata_object
		Array[Array[File]] atac_initial_adata_object = preprocess.atac_initial_adata_object

		#######################################
		################# RNA #################
		#######################################
		# Project cohort analysis outputs
		## List of samples included in the cohort
		Array[File?] rna_project_cohort_sample_list = rna_project_cohort_analysis.cohort_sample_list

		# Merged adata objects, QC plots, filtered adata objects, MMC results, normalized adata objects
		Array[File?] rna_project_merged_adata_object = rna_project_cohort_analysis.merged_adata_object
		Array[File?] rna_project_qc_initial_metadata_csv = rna_project_cohort_analysis.qc_initial_metadata_csv
		Array[Array[File]?] rna_project_qc_plots_png = rna_project_cohort_analysis.qc_plots_png
		Array[File?] rna_project_filtered_adata_object = rna_project_cohort_analysis.filtered_adata_object
		Array[File?] rna_project_mmc_extended_results_json = rna_project_cohort_analysis.mmc_extended_results_json
		Array[File?] rna_project_mmc_results_csv = rna_project_cohort_analysis.mmc_results_csv
		Array[File?] rna_project_mmc_log_txt = rna_project_cohort_analysis.mmc_log_txt
		Array[File?] rna_project_normalized_adata_object = rna_project_cohort_analysis.normalized_adata_object
		Array[File?] rna_project_all_genes_csv = rna_project_cohort_analysis.all_genes_csv
		Array[File?] rna_project_hvg_genes_csv = rna_project_cohort_analysis.hvg_genes_csv
		Array[File?] rna_project_mmc_adata_object = rna_project_cohort_analysis.mmc_adata_object
		Array[File?] rna_project_mmc_results_parquet = rna_project_cohort_analysis.mmc_results_parquet

		# Clustering outputs
		Array[File?] rna_project_integrated_adata_object = rna_project_cohort_analysis.integrated_adata_object
		Array[File?] rna_project_scvi_model_tar_gz = rna_project_cohort_analysis.scvi_model_tar_gz
		Array[File?] rna_project_labeled_cells_adata_object = rna_project_cohort_analysis.labeled_cells_adata_object
		Array[File?] rna_project_scanvi_model_tar_gz = rna_project_cohort_analysis.scanvi_model_tar_gz
		Array[File?] rna_project_scanvi_cell_types_parquet = rna_project_cohort_analysis.scanvi_cell_types_parquet
		Array[File?] rna_project_umap_clustered_adata_object = rna_project_cohort_analysis.umap_clustered_adata_object

		# PCA and Harmony integrated adata objects and artifact metrics
		Array[File?] rna_project_final_adata_object = rna_project_cohort_analysis.final_adata_object
		Array[File?] rna_project_final_metadata_csv = rna_project_cohort_analysis.final_metadata_csv
		Array[File?] rna_project_scib_report_results_csv = rna_project_cohort_analysis.scib_report_results_csv
		Array[File?] rna_project_scib_report_results_svg = rna_project_cohort_analysis.scib_report_results_svg

		# Groups and features plots
		Array[File?] rna_project_groups_umap_plot_png = rna_project_cohort_analysis.groups_umap_plot_png
		Array[File?] rna_project_features_umap_plot_png = rna_project_cohort_analysis.features_umap_plot_png

		Array[Array[File]?] rna_preprocess_manifests = rna_project_cohort_analysis.preprocess_manifest_tsvs
		Array[Array[File]?] rna_project_manifests = rna_project_cohort_analysis.cohort_analysis_manifest_tsvs

		# Cross-team cohort analysis outputs
		## List of samples included in the cohort
		File? rna_cohort_sample_list = rna_cross_team_cohort_analysis.cohort_sample_list

		# Merged adata objects, QC plots, filtered adata objects, MMC results, normalized adata objects
		File? rna_cohort_merged_adata_object = rna_cross_team_cohort_analysis.merged_adata_object
		File? rna_cohort_qc_initial_metadata_csv = rna_cross_team_cohort_analysis.qc_initial_metadata_csv
		Array[File]? rna_cohort_qc_plots_png = rna_cross_team_cohort_analysis.qc_plots_png
		File? rna_cohort_filtered_adata_object = rna_cross_team_cohort_analysis.filtered_adata_object
		File? rna_cohort_mmc_extended_results_json = rna_cross_team_cohort_analysis.mmc_extended_results_json
		File? rna_cohort_mmc_results_csv = rna_cross_team_cohort_analysis.mmc_results_csv
		File? rna_cohort_mmc_log_txt = rna_cross_team_cohort_analysis.mmc_log_txt
		File? rna_cohort_normalized_adata_object = rna_cross_team_cohort_analysis.normalized_adata_object
		File? rna_cohort_all_genes_csv = rna_cross_team_cohort_analysis.all_genes_csv
		File? rna_cohort_hvg_genes_csv = rna_cross_team_cohort_analysis.hvg_genes_csv
		File? rna_cohort_mmc_adata_object = rna_cross_team_cohort_analysis.mmc_adata_object
		File? rna_cohort_mmc_results_parquet = rna_cross_team_cohort_analysis.mmc_results_parquet

		# Clustering outputs
		File? rna_cohort_integrated_adata_object = rna_cross_team_cohort_analysis.integrated_adata_object
		File? rna_cohort_scvi_model_tar_gz = rna_cross_team_cohort_analysis.scvi_model_tar_gz
		File? rna_cohort_labeled_cells_adata_object = rna_cross_team_cohort_analysis.labeled_cells_adata_object
		File? rna_cohort_scanvi_model_tar_gz = rna_cross_team_cohort_analysis.scanvi_model_tar_gz
		File? rna_cohort_scanvi_cell_types_parquet = rna_cross_team_cohort_analysis.scanvi_cell_types_parquet
		File? rna_cohort_umap_clustered_adata_object = rna_cross_team_cohort_analysis.umap_clustered_adata_object

		# PCA and Harmony integrated adata objects and artifact metrics
		File? rna_cohort_final_adata_object = rna_cross_team_cohort_analysis.final_adata_object
		File? rna_cohort_final_metadata_csv = rna_cross_team_cohort_analysis.final_metadata_csv
		File? rna_cohort_scib_report_results_csv = rna_cross_team_cohort_analysis.scib_report_results_csv
		File? rna_cohort_scib_report_results_svg = rna_cross_team_cohort_analysis.scib_report_results_svg

		# Groups and features plots
		File? rna_cohort_groups_umap_plot_png = rna_cross_team_cohort_analysis.groups_umap_plot_png
		File? rna_cohort_features_umap_plot_png = rna_cross_team_cohort_analysis.features_umap_plot_png

		Array[File]? rna_cohort_manifests = rna_cross_team_cohort_analysis.cohort_analysis_manifest_tsvs

		#######################################
		################# ATAC ################
		#######################################
		# Project cohort analysis outputs
		## List of samples included in the cohort
		Array[File?] atac_project_cohort_sample_list = atac_project_cohort_analysis.cohort_sample_list

		# Merged adata objects, QC plots, processed bins adata object
		Array[File?] atac_project_merged_adata_object = atac_project_cohort_analysis.merged_adata_object
		Array[File?] atac_project_qc_initial_metadata_csv = atac_project_cohort_analysis.qc_initial_metadata_csv
		Array[Array[File]?] atac_project_qc_plots_png = atac_project_cohort_analysis.qc_plots_png
		Array[File?] atac_project_processed_bins_adata_object = atac_project_cohort_analysis.processed_bins_adata_object

		# Sc integration outputs
		Array[File?] atac_project_harmony_integrated_adata_object = atac_project_cohort_analysis.harmony_integrated_adata_object
		Array[File?] atac_project_harmony_clustered_adata_object = atac_project_cohort_analysis.harmony_clustered_adata_object
		Array[File?] atac_project_harmony_clustered_umap_png = atac_project_cohort_analysis.harmony_clustered_umap_png
		Array[File?] atac_project_harmony_merged_peaks_adata_object = atac_project_cohort_analysis.harmony_merged_peaks_adata_object
		Array[File?] atac_project_harmony_merged_peaks_csv = atac_project_cohort_analysis.harmony_merged_peaks_csv
		Array[File?] atac_project_harmony_peaks_matrix_adata_object = atac_project_cohort_analysis.harmony_peaks_matrix_adata_object

		Array[File?] atac_project_peakvi_integrated_adata_object = atac_project_cohort_analysis.peakvi_integrated_adata_object
		Array[File?] atac_project_peakvi_model_tar_gz = atac_project_cohort_analysis.peakvi_model_tar_gz
		Array[File?] atac_project_peakvi_clustered_adata_object = atac_project_cohort_analysis.peakvi_clustered_adata_object
		Array[File?] atac_project_peakvi_clustered_umap_png = atac_project_cohort_analysis.peakvi_clustered_umap_png
		Array[File?] atac_project_peakvi_merged_peaks_adata_object = atac_project_cohort_analysis.peakvi_merged_peaks_adata_object
		Array[File?] atac_project_peakvi_merged_peaks_csv = atac_project_cohort_analysis.peakvi_merged_peaks_csv
		Array[File?] atac_project_peakvi_peaks_matrix_adata_object = atac_project_cohort_analysis.peakvi_peaks_matrix_adata_object

		# Benchmark sc integration tools
		Array[File?] atac_project_scib_report_results_csv = atac_project_cohort_analysis.scib_report_results_csv
		Array[File?] atac_project_scib_report_results_svg = atac_project_cohort_analysis.scib_report_results_svg

		# Gene matrix adata object
		Array[File?] atac_project_gene_matrix_adata_object = atac_project_cohort_analysis.gene_matrix_adata_object
		Array[File?] atac_project_processed_gene_matrix_adata_object = atac_project_cohort_analysis.processed_gene_matrix_adata_object
		Array[File?] atac_project_all_genes_csv = atac_project_cohort_analysis.all_genes_csv
		Array[File?] atac_project_hvg_genes_csv = atac_project_cohort_analysis.hvg_genes_csv
		Array[File?] atac_project_imputed_gene_matrix_adata_object = atac_project_cohort_analysis.imputed_gene_matrix_adata_object

		# MMC from sc RNA-seq pipeline
		Array[File?] atac_project_mmc_extended_results_json = atac_project_cohort_analysis.mmc_extended_results_json
		Array[File?] atac_project_mmc_results_csv = atac_project_cohort_analysis.mmc_results_csv
		Array[File?] atac_project_mmc_log_txt = atac_project_cohort_analysis.mmc_log_txt
		Array[File?] atac_project_mmc_adata_object = atac_project_cohort_analysis.mmc_adata_object
		Array[File?] atac_project_mmc_results_parquet = atac_project_cohort_analysis.mmc_results_parquet

		# Differential chromatin analysis
		Array[File?] atac_project_celltype_merged_peaks_adata_object = atac_project_cohort_analysis.celltype_merged_peaks_adata_object
		Array[File?] atac_project_celltype_merged_peaks_csv = atac_project_cohort_analysis.celltype_merged_peaks_csv
		Array[File?] atac_project_celltype_peaks_matrix_adata_object = atac_project_cohort_analysis.celltype_peaks_matrix_adata_object
		Array[File?] atac_project_motifs_parquet = atac_project_cohort_analysis.motifs_parquet

		# Groups and features plots
		Array[File?] atac_project_groups_umap_plot_png = atac_project_cohort_analysis.groups_umap_plot_png
		Array[File?] atac_project_features_umap_plot_png = atac_project_cohort_analysis.features_umap_plot_png

		# Final artifacts
		Array[File?] atac_project_final_adata_object = atac_project_cohort_analysis.final_adata_object
		Array[File?] atac_project_final_metadata_csv = atac_project_cohort_analysis.final_metadata_csv

		Array[Array[File]?] atac_preprocess_manifests = atac_project_cohort_analysis.preprocess_manifest_tsvs
		Array[Array[File]?] atac_project_manifests = atac_project_cohort_analysis.cohort_analysis_manifest_tsvs

		# Cross-team cohort analysis outputs
		## List of samples included in the cohort
		File? atac_cohort_sample_list = atac_cross_team_cohort_analysis.cohort_sample_list

		# Merged adata objects, QC plots, processed bins adata object
		File? atac_cohort_merged_adata_object = atac_cross_team_cohort_analysis.merged_adata_object
		File? atac_cohort_qc_initial_metadata_csv = atac_cross_team_cohort_analysis.qc_initial_metadata_csv
		Array[File]? atac_cohort_qc_plots_png = atac_cross_team_cohort_analysis.qc_plots_png
		File? atac_cohort_processed_bins_adata_object = atac_cross_team_cohort_analysis.processed_bins_adata_object

		# Sc integration outputs
		File? atac_cohort_harmony_integrated_adata_object = atac_cross_team_cohort_analysis.harmony_integrated_adata_object
		File? atac_cohort_harmony_clustered_adata_object = atac_cross_team_cohort_analysis.harmony_clustered_adata_object
		File? atac_cohort_harmony_clustered_umap_png = atac_cross_team_cohort_analysis.harmony_clustered_umap_png
		File? atac_cohort_harmony_merged_peaks_adata_object = atac_cross_team_cohort_analysis.harmony_merged_peaks_adata_object
		File? atac_cohort_harmony_merged_peaks_csv = atac_cross_team_cohort_analysis.harmony_merged_peaks_csv
		File? atac_cohort_harmony_peaks_matrix_adata_object = atac_cross_team_cohort_analysis.harmony_peaks_matrix_adata_object

		File? atac_cohort_peakvi_integrated_adata_object = atac_cross_team_cohort_analysis.peakvi_integrated_adata_object
		File? atac_cohort_peakvi_model_tar_gz = atac_cross_team_cohort_analysis.peakvi_model_tar_gz
		File? atac_cohort_peakvi_clustered_adata_object = atac_cross_team_cohort_analysis.peakvi_clustered_adata_object
		File? atac_cohort_peakvi_clustered_umap_png = atac_cross_team_cohort_analysis.peakvi_clustered_umap_png
		File? atac_cohort_peakvi_merged_peaks_adata_object = atac_cross_team_cohort_analysis.peakvi_merged_peaks_adata_object
		File? atac_cohort_peakvi_merged_peaks_csv = atac_cross_team_cohort_analysis.peakvi_merged_peaks_csv
		File? atac_cohort_peakvi_peaks_matrix_adata_object = atac_cross_team_cohort_analysis.peakvi_peaks_matrix_adata_object

		# Benchmark sc integration tools
		File? atac_cohort_scib_report_results_csv = atac_cross_team_cohort_analysis.scib_report_results_csv
		File? atac_cohort_scib_report_results_svg = atac_cross_team_cohort_analysis.scib_report_results_svg

		# Gene matrix adata object
		File? atac_cohort_gene_matrix_adata_object = atac_cross_team_cohort_analysis.gene_matrix_adata_object
		File? atac_cohort_processed_gene_matrix_adata_object = atac_cross_team_cohort_analysis.processed_gene_matrix_adata_object
		File? atac_cohort_all_genes_csv = atac_cross_team_cohort_analysis.all_genes_csv
		File? atac_cohort_hvg_genes_csv = atac_cross_team_cohort_analysis.hvg_genes_csv
		File? atac_cohort_imputed_gene_matrix_adata_object = atac_cross_team_cohort_analysis.imputed_gene_matrix_adata_object

		# MMC from sc RNA-seq pipeline
		File? atac_cohort_mmc_extended_results_json = atac_cross_team_cohort_analysis.mmc_extended_results_json
		File? atac_cohort_mmc_results_csv = atac_cross_team_cohort_analysis.mmc_results_csv
		File? atac_cohort_mmc_log_txt = atac_cross_team_cohort_analysis.mmc_log_txt
		File? atac_cohort_mmc_adata_object = atac_cross_team_cohort_analysis.mmc_adata_object
		File? atac_cohort_mmc_results_parquet = atac_cross_team_cohort_analysis.mmc_results_parquet

		# Differential chromatin analysis
		File? atac_cohort_celltype_merged_peaks_adata_object = atac_cross_team_cohort_analysis.celltype_merged_peaks_adata_object
		File? atac_cohort_celltype_merged_peaks_csv = atac_cross_team_cohort_analysis.celltype_merged_peaks_csv
		File? atac_cohort_celltype_peaks_matrix_adata_object = atac_cross_team_cohort_analysis.celltype_peaks_matrix_adata_object
		File? atac_cohort_motifs_parquet = atac_cross_team_cohort_analysis.motifs_parquet

		# Groups and features plots
		File? atac_cohort_groups_umap_plot_png = atac_cross_team_cohort_analysis.groups_umap_plot_png
		File? atac_cohort_features_umap_plot_png = atac_cross_team_cohort_analysis.features_umap_plot_png

		# Final artifacts
		File? atac_cohort_final_adata_object = atac_cross_team_cohort_analysis.final_adata_object
		File? atac_cohort_final_metadata_csv = atac_cross_team_cohort_analysis.final_metadata_csv

		Array[File]? atac_cohort_manifests = atac_cross_team_cohort_analysis.cohort_analysis_manifest_tsvs

		###########################################
		################# MULTIOME ################
		###########################################
		# Project cohort analysis outputs
		Array[File?] multiome_project_filtered_mudata_object = multiome_project_cohort_analysis.filtered_mudata_object
		Array[File?] multiome_project_multivi_integrated_mudata_object = multiome_project_cohort_analysis.multivi_integrated_mudata_object
		Array[File?] multiome_project_multivi_model_tar_gz = multiome_project_cohort_analysis.multivi_model_tar_gz
		Array[File?] multiome_project_multivi_clustered_mudata_object = multiome_project_cohort_analysis.multivi_clustered_mudata_object
		Array[File?] multiome_project_multivi_clustered_umap_png = multiome_project_cohort_analysis.multivi_clustered_umap_png
		Array[File?] multiome_project_annotated_mudata_object = multiome_project_cohort_analysis.annotated_mudata_object
		Array[File?] multiome_project_atac_annotated_adata_object = multiome_project_cohort_analysis.atac_annotated_adata_object
		Array[File?] multiome_project_feature_linked_mudata_object = multiome_project_cohort_analysis.feature_linked_mudata_object
		Array[File?] multiome_project_feature_linkage_csv = multiome_project_cohort_analysis.feature_linkage_csv

		Array[Array[File]?] multiome_preprocess_manifests = multiome_project_cohort_analysis.preprocess_manifest_tsvs
		Array[Array[File]?] multiome_project_manifests = multiome_project_cohort_analysis.cohort_analysis_manifest_tsvs

		# Cross-team cohort analysis outputs
		File? multiome_cohort_filtered_mudata_object = multiome_cross_team_cohort_analysis.filtered_mudata_object
		File? multiome_cohort_multivi_integrated_mudata_object = multiome_cross_team_cohort_analysis.multivi_integrated_mudata_object
		File? multiome_cohort_multivi_model_tar_gz = multiome_cross_team_cohort_analysis.multivi_model_tar_gz
		File? multiome_cohort_multivi_clustered_mudata_object = multiome_cross_team_cohort_analysis.multivi_clustered_mudata_object
		File? multiome_cohort_multivi_clustered_umap_png = multiome_cross_team_cohort_analysis.multivi_clustered_umap_png
		File? multiome_cohort_annotated_mudata_object = multiome_cross_team_cohort_analysis.annotated_mudata_object
		File? multiome_cohort_atac_annotated_adata_object = multiome_cross_team_cohort_analysis.atac_annotated_adata_object
		File? multiome_cohort_feature_linked_mudata_object = multiome_cross_team_cohort_analysis.feature_linked_mudata_object
		File? multiome_cohort_feature_linkage_csv = multiome_cross_team_cohort_analysis.feature_linkage_csv

		Array[File]? multiome_cohort_manifests = multiome_cross_team_cohort_analysis.cohort_analysis_manifest_tsvs
	}

	meta {
		description: "Harmonized human postmortem-derived brain sequencing (PMDBS) brain sc/sn Multiome (RNA + ATAC) workflow."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files during cross-team cohort analysis."}
		projects: {help: "The project ID, set of samples and their associated reads and metadata, output bucket locations, sc data type, and whether or not to run project-level cohort analysis."}
		cellranger_arc_reference_data: {help: "CellRanger ARC reference data; see https://www.10xgenomics.com/support/software/cell-ranger-arc/downloads."}
		cellbender_fpr: {help: "Cellbender false positive rate. [0.0]"}
		pct_counts_mt_max: {help: "Maximum percentage of mitochondrial gene counts allowed per cell. [10]"}
		doublet_score_max: {help: "Maximum doublet detection score threshold. [0.2]"}
		total_counts_limits: {help: "Minimum and maximum total UMI (unique molecular identifier) counts per cell. [100, 100000]"}
		n_genes_by_counts_limits: {help: "Minimum and maximum number of genes detected per cell (genes with at least one count). [100, 10000]"}
		allen_brain_mmc_precomputed_stats_h5: {help: "A precomputed statistics file from the Allen Brain Cell Atlas containing reference statistics (the average gene expression profile per cell type cluster and cell type taxonomy)."}
		norm_target_sum: {help: "The total count value that each cell will be normalized to. [10000]"}
		n_top_genes: {help: "Number of HVG genes to keep. [3000]"}
		n_comps: {help: "Number of principal components to compute. [30]"}
		scvi_latent_key: {help: "Latent key to save the scVI latent to. ['X_scVI']"}
		scanvi_latent_key: {help: "Latent key to save the scANVI latent to. ['X_scANVI']"}
		scanvi_predictions_key: {help: "scANVI cell type predictions column name. ['C_scANVI']"}
		batch_key: {help: "Key in AnnData object for batch information. ['batch_id']"}
		n_neighbors: {help: "The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. [15]"}
		leiden_res: {help: "Leiden resolutions which are the parameter values controlling the coarseness of the clustering. [0.2, 0.5, 1.0]"}
		rna_groups: {help: "RNA groups to produce umap plots for. ['sample', 'batch', 'cell_type', 'leiden_res_0.05', 'leiden_res_0.10', 'leiden_res_0.20', 'leiden_res_0.40']"}
		rna_features: {help: "RNA features to produce umap plots for. ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rb', 'doublet_score', 'S_score', 'G2M_score']"}
		peakvi_latent_key: {help: "Latent key to save the peakVI latent to. ['X_peakVI']"}
		peakvi_max_epochs: {help: "The maximum number of full passes through the training data during PeakVI model training. If the model converges early, training will halt before this limit is reached. [300]"}
		peakvi_n_hidden: {help: "Number of nodes per hidden layer for PeakVI (i.e., the number of neurons per fully-connected layer between the input features and the latent space). [128]"}
		peakvi_batch_size: {help: "Training batch size for PeakVI. Controls how many cells are processed per training step. [128]"}
		atac_groups: {help: "ATAC groups to produce umap plots for. ['sample', 'batch', 'team', 'dataset', 'batch_id', 'leiden']"}
		atac_features: {help: "ATAC features to produce umap plots for. ['n_fragment', 'tsse', 'frac_dup', 'frac_mito', 'doublet_score', 'doublet_probability']"}
		min_overlap_frac: {help: "Minimum fraction of the smaller modality that must be shared; below this the run fails as a likely barcode-key mismatch. [0.5]"}
		multivi_latent_key: {help: "Latent key to save the MultiVI latent to. ['X_multiVI']"}
		multivi_max_epochs: {help: "The maximum number of full passes through the training data during MultiVI model training. If the model converges early, training will halt before this limit is reached. [300]"}
		multivi_n_hidden: {help: "Number of nodes per hidden layer for MultiVI (i.e., the number of neurons per fully-connected layer between the input features and the latent space). [128]"}
		multivi_batch_size: {help: "Training batch size for MultiVI. Controls how many cells are processed per training step. [128]"}
		peak_gene_only: {help: "Keep only peak-gene lineages (drop peak-peak and gene-gene). [true]"}
		min_linkage_corr: {help: "Minimum absolute feature linkage correlation to retain. [0.1]"}
		min_linkage_sig: {help: "Minimum feature linkage significance to retain. [0.0]"}
		run_cross_team_cohort_analysis: {help: "Whether to run downstream harmonization steps on all samples across projects. If set to false, only preprocessing steps (cellranger and generating the initial adata object(s)) will run for samples. [false]"}
		cohort_raw_data_bucket: {help: "Bucket to upload cross-team cohort intermediate files to."}
		cohort_staging_data_buckets: {help: "Set of buckets to stage cross-team cohort analysis outputs in."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}
