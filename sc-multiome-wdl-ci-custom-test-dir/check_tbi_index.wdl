version 1.0

# Validate input TBI index files
# Input type: TBI file

task check_tbi_index {
	input {
		File current_run_output
		File validated_output
	}

	Int disk_size = ceil(size([current_run_output, validated_output], "GB") + 10)

	command <<<
		set -euo pipefail

		err() {
			message=$1

			echo -e "[ERROR] $message" >&2
		}

		expected_magic=$(printf '\x54\x42\x49\x01')

		check_tbi() {
			local file="$1"
			local label="$2"
			local magic
			# Read first 4 decompressed bytes; tolerate SIGPIPE from head
			magic=$(gzip -cd "$file" 2>/dev/null | head -c 4 || true)
			if [[ "$magic" != "$expected_magic" ]]; then
				err "$label [$(basename "$file")] is not a valid TBI index file"
				return 1
			fi
			echo "$label [$(basename "$file")] is a valid TBI index file"
		}

		check_tbi "~{validated_output}"   "Validated output"   || exit 1
		check_tbi "~{current_run_output}" "Current run output" || exit 1
	>>>

	output {
	}

	runtime {
		docker: "dnastack/dnastack-wdl-ci-tools:0.1.1"
		cpu: 2
		memory: "3.75 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 1
	}
}
