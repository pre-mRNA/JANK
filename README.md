# Ultra-fast kmer counter
<p align="center">
  <img src="./img/ufkc.png" alt="UFKC logo" width="200" height="200">
</p>
JANK is a high-performance kmer counter developed in Rust, designed to efficiently process both FASTA and FASTQ files using parallel processing.

## Features
- Compatible with FASTA and FASTQ input files
- Supports gzipped files and automatic file-type identification
- In-built parallel processing for optimal performance on multi-core systems
- Offers the option to count reverse complement k-mers


## Dependencies
- Rust programming language (version >= 1.54.0)
- clap (version >= 3.0.0)
- flate2 (version >= 1.0.20)
- rayon (version >= 1.5.1)
**Note: Dependencies will be automatically installed during compilation.**

## Installation
```
# Clone the repository
git clone https://github.com/pre-mRNA/UFKC.git

# Compile the code
cd UFKC && cargo build --release

# Add jank to your path
export PATH="$PATH:$(pwd)/target/release"
```
The compiled binary will be located in the target/release directory, and added to your current PATH.

## Usage
```
UFKC [FLAGS] [OPTIONS] <input_file>
```
<input_file>: The input FASTA or FASTQ file (required)
### Flags
-r, --count-reverse-complement: Count reverse complement k-mers as well
### Options
-k, --kmer-size <k_value>: Size of k-mers to count (default: 9)
-c, --cpu-count <cpu_count>: Number of CPU cores to use (default: 0, which means the program will use all available cores)

```
UFKC -k 7 -c 4 -r input.fasta
```
This command will count 7-mers in the input.fasta file, using 4 CPU cores, and also counting reverse complement k-mers.
