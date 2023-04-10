use std::fs::File;
use std::io::{BufRead, BufReader, Write, stderr};
use std::collections::HashMap;
use clap::{Arg, App, value_t};
use rayon::prelude::*;
use std::sync::Mutex;
use flate2::bufread::GzDecoder;

fn main() {
    let matches = App::new("Just ANother Kmer counter")
    .arg(Arg::with_name("input_file")
    .help("Input FASTA or FASTQ file")
    .required(true)
    .index(1))
    .arg(Arg::with_name("k_value")
    .short("k")
    .long("kmer-size")
    .help("Size of k-mers to count")
    .takes_value(true)
    .default_value("9"))
    .arg(Arg::with_name("cpu_count")
    .short("c")
    .long("cpu-count")
    .help("Number of CPU cores to use")
    .takes_value(true)
    .default_value("0"))
    .arg(Arg::with_name("count_reverse_complement")
    .short("r")
    .long("count-reverse-complement")
    .help("Count reverse complement k-mers as well")
    .takes_value(false))
    .get_matches();

    let input_file = matches.value_of("input_file").unwrap();
    let kmer_size: usize = value_t!(matches.value_of("k_value"), usize).unwrap_or_else(|e| e.exit());
    let cpu_count: usize = value_t!(matches.value_of("cpu_count"), usize).unwrap_or_else(|e| e.exit());
    let count_reverse_complement = matches.is_present("count_reverse_complement");

    if cpu_count > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(cpu_count).build_global().unwrap();
    }

    let sequences = read_sequences(input_file);
    let kmer_counts = Mutex::new(HashMap::new());

    let total_length: usize = sequences.par_iter()
    .map(|seq| {
        if seq.len() >= kmer_size {
            for idx in 0..=seq.len() - kmer_size {
                let kmer = String::from_utf8_lossy(&seq[idx..idx + kmer_size]).to_string();
                let mut kmer_counts = kmer_counts.lock().unwrap();
                kmer_counts.entry(kmer.clone()).and_modify(|c| *c += 1).or_insert(1);
                if count_reverse_complement {
                    let reverse_complement = reverse_complement(&kmer);
                    kmer_counts.entry(reverse_complement).and_modify(|c| *c += 1).or_insert(1);
                }
            }
            seq.len()
        } else {
            writeln!(stderr(), "Discarding short sequence with length: {} (less than k-mer size: {})", seq.len(), kmer_size).unwrap();
            0
        }
    })
    .reduce(|| 0, |a, b| a + b);

    let total_kmers: usize = kmer_counts.lock().unwrap().values().map(|v| *v).sum();

    writeln!(stderr(), "Total length of input sequences: {}", total_length).unwrap();
    writeln!(stderr(), "Total number of k-mers counted: {}", total_kmers).unwrap();
    writeln!(stderr(), "Using {} CPU cores", rayon::current_num_threads()).unwrap();

    println!("# Command: kmer_counter {} -k {}", input_file, kmer_size);
    for (kmer, count) in kmer_counts.lock().unwrap().iter() {
        println!("{} {}", kmer, count);
    }
}

fn read_sequences(filename: &str) -> Vec<Vec<u8>> {
    let file = File::open(filename).expect("Error opening file");
    let is_gzipped = filename.to_lowercase().ends_with(".gz");

    let reader: Box<dyn BufRead> = if is_gzipped {
        Box::new(BufReader::new(GzDecoder::new(BufReader::new(file))))
    } else {
        Box::new(BufReader::new(file))
    };

    let is_fasta = filename.to_lowercase().contains(".fasta") || filename.to_lowercase().contains(".fa");
    let is_fastq = filename.to_lowercase().contains(".fastq") || filename.to_lowercase().contains(".fq");

    let mut sequences = Vec::new();
    let mut current_seq = Vec::new();
    let mut in_sequence = false;
    let mut ignore_lines = 0;

    eprintln!("Reading sequences...");
    for line in reader.lines() {
        let line = line.unwrap();
        if ignore_lines > 0 {
            ignore_lines -= 1;
            continue;
        }

        if line.starts_with('>') && is_fasta || line.starts_with('@') && is_fastq {
            if !current_seq.is_empty() {
                sequences.push(current_seq);
                current_seq = Vec::new();
            }
            in_sequence = true;
        } else if in_sequence && is_fastq && line.starts_with('+') {
            ignore_lines = 1;
            in_sequence = false;
        } else if in_sequence {
            current_seq.extend(line.bytes());
        }
    }

    if !current_seq.is_empty() {
        sequences.push(current_seq);
    }

    eprintln!("Finished reading sequences. Sequences count: {}", sequences.len());

    sequences
}

fn reverse_complement(kmer: &str) -> String {
    let mut rev_comp = String::with_capacity(kmer.len());
    for base in kmer.chars().rev() {
        let complement = match base {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => panic!("Invalid base in k-mer: {}", base),
        };
        rev_comp.push(complement);
    }
    rev_comp
}
