version 1.0

# Integrate MultiVI and cluster

workflow multivi_integration {
	input {
		String cohort_id
		File filtered_mudata_object

		String batch_key
		String multivi_latent_key
		Int multivi_max_epochs
		Int multivi_n_hidden
		Int multivi_batch_size

		Int n_neighbors
		Array[Float] leiden_res

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	call integrate_multivi {
		input:
			cohort_id = cohort_id,
			filtered_mudata_object = filtered_mudata_object,
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

	call cluster_multivi {
		input:
			cohort_id = cohort_id,
			multivi_integrated_mudata_object = integrate_multivi.multivi_integrated_mudata_object,
			multivi_latent_key = multivi_latent_key,
			n_neighbors = n_neighbors,
			leiden_res = leiden_res,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		File multivi_integrated_mudata_object = integrate_multivi.multivi_integrated_mudata_object
		File multivi_model_tar_gz = integrate_multivi.multivi_model_tar_gz #!FileCoercion
		File multivi_clustered_mudata_object = cluster_multivi.multivi_clustered_mudata_object
		File multivi_clustered_umap_png = cluster_multivi.multivi_clustered_umap_png #!FileCoercion
	}

	meta {
		description: "Integrates samples by training a PeakVI variational autoencoder for batch correction and clusters cells using Leiden community detection, producing an integrated AnnData object, PeakVI model, and UMAP plot."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		filtered_mudata_object: {help: "Filtered paired Multiome MuData object."}
		batch_key: {help: "Key in AnnData object for batch information. ['batch_id']"}
		multivi_latent_key: {help: "Latent key to save the MultiVI latent to. ['X_multiVI']"}
		multivi_max_epochs: {help: "The maximum number of full passes through the training data during MultiVI model training. If the model converges early, training will halt before this limit is reached. [300]"}
		multivi_n_hidden: {help: "Number of nodes per hidden layer for MultiVI (i.e., the number of neurons per fully-connected layer between the input features and the latent space). [128]"}
		multivi_batch_size: {help: "Training batch size for MultiVI. Controls how many cells are processed per training step. [128]"}
		n_neighbors: {help: "The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. [15]"}
		leiden_res: {help: "Leiden resolutions which are the parameter values controlling the coarseness of the clustering. [0.2, 0.5, 1.0]"}
		raw_data_path: {help: "Raw data bucket path for merged adata and QC plots outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
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
			-o "~{cohort_id}.multivi_model.tar.gz"
	>>>

	output {
		File multivi_integrated_mudata_object = "~{cohort_id}.multivi_integrated.h5mu"
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

task cluster_multivi {
	input {
		String cohort_id
		File multivi_integrated_mudata_object

		String multivi_latent_key
		Int n_neighbors
		Array[Float] leiden_res

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int calc_mem_gb = ceil(size(multivi_integrated_mudata_object, "GB") * 5 + 50)
	Int mem_gb = if calc_mem_gb > 624 then 624 else calc_mem_gb
	Int disk_size = ceil(size(multivi_integrated_mudata_object, "GB") * 4 + 50)

	command <<<
		set -euo pipefail

		/usr/bin/time -v \
		cluster_multivi \
			--mudata-input ~{multivi_integrated_mudata_object} \
			--latent-key ~{multivi_latent_key} \
			--n-neighbors ~{n_neighbors} \
			--leiden-res ~{sep=" " leiden_res} \
			--mudata-output ~{cohort_id}.multivi_clustered.h5mu

		mv "figures/umap.png" "./~{cohort_id}.multivi_umap.png"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.multivi_umap.png"
	>>>

	output {
		File multivi_clustered_mudata_object = "~{cohort_id}.multivi_clustered.h5mu"
		String multivi_clustered_umap_png = "~{raw_data_path}/~{cohort_id}.multivi_umap.png"
	}

	runtime {
		docker: "~{container_registry}/scvi_tools:1.0.0"
		cpu: 4
		cpuPlatform: "Intel Cascade Lake"
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}

	meta {
		description: "Constructs a k-nearest neighbor graph on the MultiVI latent embedding and applies Leiden community detection to cluster cells."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		multivi_integrated_adata_object: {help: "MultiVI-integrated AnnData object."}
		multivi_latent_key: {help: "Latent key to save the MultiVI latent to. ['X_multiVI']"}
		n_neighbors: {help: "The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. [15]"}
		leiden_res: {help: "Leiden resolutions which are the parameter values controlling the coarseness of the clustering. [0.2, 0.5, 1.0]"}
		raw_data_path: {help: "Raw data bucket path for outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in. ['us-central1-c us-central1-f']"}
	}
}
