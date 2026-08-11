version 1.0

# Run steps in the cohort analysis

import "../../wf-common/wdl/tasks/write_cohort_sample_list.wdl" as WriteCohortSampleList
import "../../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs

workflow cohort_analysis {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids

		# If provided, these files will be uploaded to the staging bucket alongside other intermediate files made by this workflow
		Array[String] preprocessing_output_file_paths = []

		Array[File] cellranger_arc_feature_linkage_bedpe

		File filtered_rna_adata_object
		File filtered_atac_adata_object

		Float min_overlap_frac

		String batch_key
		String multivi_latent_key
		Int multivi_max_epochs
		Int multivi_n_hidden
		Int multivi_batch_size

		Boolean peak_gene_only
		Float min_linkage_corr
		Float min_linkage_sig

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

	call joint_qc {
		input:
			cohort_id = cohort_id,
			filtered_rna_adata_object = filtered_rna_adata_object,
			filtered_atac_adata_object = filtered_atac_adata_object,
			min_overlap_frac = min_overlap_frac,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call integrate_multivi {
		input:
			cohort_id = cohort_id,
			filtered_mudata_object = joint_qc.filtered_mudata_object, #!FileCoercion
			batch_key = batch_key,
			multivi_latent_key = multivi_latent_key,
			multivi_max_epochs = multivi_max_epochs,
			multivi_n_hidden = multivi_n_hidden,
			multivi_batch_size = multivi_batch_size,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call parse_feature_linkage {
		input:
			cohort_id = cohort_id,
			cellranger_arc_feature_linkage_bedpe = cellranger_arc_feature_linkage_bedpe,
			multivi_integrated_mudata_object = integrate_multivi.multivi_integrated_mudata_object, #!FileCoercion
			peak_gene_only = peak_gene_only,
			min_linkage_corr = min_linkage_corr,
			min_linkage_sig = min_linkage_sig,
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

	Array[String] cohort_analysis_final_output_paths = flatten([
		[
			joint_qc.filtered_mudata_object
		],
		[
			integrate_multivi.multivi_integrated_mudata_object,
			integrate_multivi.multivi_model_tar_gz
		],
		[
			parse_feature_linkage.feature_linked_mudata_object,
			parse_feature_linkage.feature_linkage_csv
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

		# 
		File filtered_mudata_object = joint_qc.filtered_mudata_object #!FileCoercion
		File multivi_integrated_mudata_object = integrate_multivi.multivi_integrated_mudata_object #!FileCoercion
		File multivi_model_tar_gz = integrate_multivi.multivi_model_tar_gz #!FileCoercion
		File feature_linked_mudata_object = parse_feature_linkage.feature_linked_mudata_object #!FileCoercion
		File feature_linkage_csv = parse_feature_linkage.feature_linkage_csv #!FileCoercion

		Array[File] preprocess_manifest_tsvs = upload_preprocess_files.manifests #!FileCoercion
		Array[File] cohort_analysis_manifest_tsvs = upload_cohort_analysis_files.manifests #!FileCoercion
	}

	meta {
		description: "Runs QC filtering on combined Multiome data (RNA + ATAC) for a cohort."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		project_sample_ids: {help: "Associated team ID, sample ID, and dataset DOI URL; used to generate a sample list."}
		preprocessing_output_file_paths: {help: "Selected preprocessed output files to upload to the staging bucket alongside selected cohort analysis output files."}
		filtered_rna_adata_object: {help: "Cohort-level filtered RNA AnnData ({cohort}.filtered.h5ad from the RNA `filter` task)"}
		filtered_atac_adata_object: {help: "Cohort-level filtered ATAC AnnData ({cohort}.merged_filtered.h5ad from the ATAC `merge_and_qc` task)"}
		min_overlap_frac: {help: "Minimum fraction of the smaller modality that must be shared; below this the run fails as a likely barcode-key mismatch [0.5]"}
		batch_key: {help: "Key in AnnData object for batch information. ['batch_id']"}
		multivi_latent_key: {help: "Latent key to save the MultiVI latent to. ['X_multiVI']"}
		multivi_max_epochs: {help: "The maximum number of full passes through the training data during MultiVI model training. If the model converges early, training will halt before this limit is reached. [300]"}
		multivi_n_hidden: {help: "Number of nodes per hidden layer for MultiVI (i.e., the number of neurons per fully-connected layer between the input features and the latent space). [128]"}
		multivi_batch_size: {help: "Training batch size for MultiVI. Controls how many cells are processed per training step. [128]"}
		peak_gene_only: {help: "Keep only peak-gene linages (drop peak-peak and gene-gene). [true]"}
		min_linkage_corr: {help: "Minimum absolute feature linkage correlation to retain. [0.1]"}
		min_linkage_sig: {help: "Minimum feature linkage significance to retain. [0.0]"}
		workflow_name: {help: "Workflow name; stored in the file-level manifest and final manifest with all saved files."}
		workflow_version: {help: "Workflow version; stored in the file-level manifest and final manifest with all saved files."}
		workflow_release: {help: "GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		crn_release_version: {help: "CRN Cloud release version; used to organize outputs and for the CRN Cloud release."}
		run_timestamp: {help: "UTC timestamp; stored in the file-level manifest and final manifest with all saved files."}
		raw_data_path_prefix: {help: "Raw data bucket path prefix; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis`)."}
		staging_data_buckets: {help: "Array of staging data buckets to upload intermediate files to (i.e., DEV or UAT buckets depending on internal QC status)."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}

task joint_qc {
	input {
		String cohort_id

		File filtered_rna_adata_object
		File filtered_atac_adata_object
		Float min_overlap_frac

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size([filtered_rna_adata_object, filtered_atac_adata_object], "GB") * 8 + 40)
	Int disk_size = ceil(size([filtered_rna_adata_object, filtered_atac_adata_object], "GB") * 4 + 40)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		joint_qc \
			--rna-adata-input ~{filtered_rna_adata_object} \
			--atac-adata-input ~{filtered_atac_adata_object} \
			--min-overlap-frac ~{min_overlap_frac} \
			--mudata-output "~{cohort_id}.filtered.h5mu"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.filtered.h5mu"
	>>>

	output {
		String filtered_mudata_object = "~{raw_data_path}/~{cohort_id}.filtered.h5mu"
	}

	runtime {
		docker: "~{container_registry}/sc_multiome_tools:1.0.0"
		cpu: 16
		cpuPlatform: "Intel Cascade Lake"
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		bootDiskSizeGb: 40
		zones: zones
	}

	meta {
		description: "Intersect filtered RNA and ATAC barcodes into a paired Multiome MuData."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		filtered_rna_adata_object: {help: "Cohort-level filtered RNA AnnData ({cohort}.filtered.h5ad from the RNA `filter` task)."}
		filtered_atac_adata_object: {help: "Cohort-level filtered ATAC AnnData ({cohort}.merged_filtered.h5ad from the ATAC `merge_and_qc` task)"}
		min_overlap_frac: {help: "Minimum fraction of the smaller modality that must be shared; below this the run fails as a likely barcode-key mismatch. [0.5]"}
		raw_data_path: {help: "Raw data bucket path for outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}

task integrate_multivi {
	input {
		String cohort_id
		File filtered_mudata_object

		String batch_key
		String multivi_latent_key
		Int multivi_max_epochs
		Int multivi_n_hidden
		Int multivi_batch_size

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	# N1 custom machine type must be compatible with gpuType
	Int calc_mem_gb = ceil(size(filtered_mudata_object, "GB") * 8 + 50)
	Int mem_gb = if calc_mem_gb > 312 then 312 else calc_mem_gb
	Int disk_size = ceil(size(filtered_mudata_object, "GB") * 4 + 50)

	command <<<
		set -euo pipefail

		nvidia-smi

		mkdir multivi_dir

		/usr/bin/time -v \
		integrate_multivi \
			--mudata-input ~{filtered_mudata_object} \
			--batch-key ~{batch_key} \
			--latent-key ~{multivi_latent_key} \
			--max-epochs ~{multivi_max_epochs} \
			--n-hidden ~{multivi_n_hidden} \
			--batch-size ~{multivi_batch_size} \
			--mudata-output ~{cohort_id}.multivi_integrated.h5mu \
			--output-multivi-dir "~{cohort_id}_multivi_model"

		# Model name cannot be changed because scvi models serialization expects a path containing a model.pt object
		tar -czvf "~{cohort_id}.multivi_model.tar.gz" "~{cohort_id}_multivi_model"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.multivi_integrated.h5mu" \
			-o "~{cohort_id}.multivi_model.tar.gz"
	>>>

	output {
		String multivi_integrated_mudata_object = "~{raw_data_path}/~{cohort_id}.multivi_integrated.h5mu"
		String multivi_model_tar_gz = "~{raw_data_path}/~{cohort_id}.multivi_model.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/scvi_tools:1.0.0"
		cpu: 8
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		bootDiskSizeGb: 30
		zones: zones
		gpuType: "nvidia-tesla-t4"
		gpuCount: 1
	}

	meta {
		description: "Trains a MultiVI variational autoencoder on paired and unpaired sc RNAseq and sc ATACseq data to learn a joint, batch-corrected latent representation, and exports the trained MultiVI model."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		filtered_mudata_object: {help: "Filtered paired Multiome MuData object."}
		batch_key: {help: "Key in AnnData object for batch information. ['batch_id']"}
		multivi_latent_key: {help: "Latent key to save the MultiVI latent to. ['X_multiVI']"}
		multivi_max_epochs: {help: "The maximum number of full passes through the training data during MultiVI model training. If the model converges early, training will halt before this limit is reached. [300]"}
		multivi_n_hidden: {help: "Number of nodes per hidden layer for MultiVI (i.e., the number of neurons per fully-connected layer between the input features and the latent space). [128]"}
		multivi_batch_size: {help: "Training batch size for MultiVI. Controls how many cells are processed per training step. [128]"}
		raw_data_path: {help: "Raw data bucket path for outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}

task parse_feature_linkage {
	input {
		String cohort_id

		Array[File] cellranger_arc_feature_linkage_bedpe

		File multivi_integrated_mudata_object

		Boolean peak_gene_only
		Float min_linkage_corr
		Float min_linkage_sig

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	String peak_gene_only_flag = if peak_gene_only then "--peak-gene-only" else ""

	Int mem_gb = ceil((size(multivi_integrated_mudata_object, "GB") + size(flatten([cellranger_arc_feature_linkage_bedpe]), "GB")) * 2 + 20)
	Int disk_size = ceil((size(multivi_integrated_mudata_object, "GB") + size(flatten([cellranger_arc_feature_linkage_bedpe]), "GB")) * 4 + 50)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		parse_feature_linkage \
			--mudata-input ~{multivi_integrated_mudata_object} \
			--bedpe-inputs ~{sep=" " cellranger_arc_feature_linkage_bedpe} \ #TODO
			~{peak_gene_only_flag} \
			--corr-threshold ~{min_linkage_corr} \
			--significance-threshold ~{min_linkage_sig} \
			--mudata-output "~{cohort_id}.feature_linked.h5mu" \
			--links-output "~{cohort_id}.feature_linkage_per_sample.csv"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.feature_linked.h5mu" \
			-o "~{cohort_id}.feature_linkage_per_sample.csv"
	>>>

	output {
		String feature_linked_mudata_object = "~{raw_data_path}/~{cohort_id}.feature_linked.h5mu"
		String feature_linkage_csv = "~{raw_data_path}/~{cohort_id}.feature_linkage_per_sample.csv"
	}

	runtime {
		docker: "~{container_registry}/sc_multiome_tools:1.0.0"
		cpu: 4
		cpuPlatform: "Intel Cascade Lake"
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		bootDiskSizeGb: 40
		zones: zones
	}

	meta {
		description: "Stack per-sample Cell Ranger ARC feature_linkage.bedpe files into the merged Multiome MuData. Peaks are per-sample ARC calls, so links are kept sample-stamped rather than collapsed to a cross-sample consensus."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		cellranger_arc_feature_linkage_bedpe: {help: "Per-sample Cell Ranger ARC feature_linkage.bedpe files, one per sample, from preprocess."}
		multivi_integrated_mudata_object: {help: "Paired Multiome MuData object with the joint, batch-corrected MultiVI latent representation stored in obsm."}
		peak_gene_only: {help: "Keep only peak-gene linages (drop peak-peak and gene-gene). [true]"}
		min_linkage_corr: {help: "Minimum absolute feature linkage correlation to retain. [0.1]"}
		min_linkage_sig: {help: "Minimum feature linkage significance to retain. [0.0]"}
		raw_data_path: {help: "Raw data bucket path for outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}
