mod fast_fishers_exact_test;

use fast_fishers_exact_test::fishers_exact_test;
use linked_hash_map::LinkedHashMap;
use statrs::distribution::{Binomial, ContinuousCDF, DiscreteCDF, Normal};
use std::collections::HashMap;
use std::fmt;
use std::io::{BufRead, BufReader};
use std::process::{Command, Stdio};
use std::str;

#[derive(Debug, Copy, Clone)]
struct StrandCount {
    fwd: u64,
    rev: u64,
}

#[derive(Debug, Copy, Clone)]
enum Strand {
    Both,
    Forward,
    Reverse,
}

impl StrandCount {
    fn total(&self) -> u64 {
        self.fwd + self.rev
    }

    fn of_strand(&self, strand: Strand) -> u64 {
        match strand {
            Strand::Both => self.total(),
            Strand::Forward => self.fwd,
            Strand::Reverse => self.rev,
        }
    }

    fn new() -> StrandCount {
        StrandCount { fwd: 0, rev: 0 }
    }
}

static EMPTY_STRANDCOUNTS: StrandCount = StrandCount { fwd: 0, rev: 0 };

// Stores counts extracted from a pileup string of a single locus and sample
#[derive(Debug)]
struct PileupCount {
    snv_count: LinkedHashMap<char, StrandCount>,
    ins_count: LinkedHashMap<String, StrandCount>,
    del_count: LinkedHashMap<String, StrandCount>,
    depth: StrandCount,
}

fn count_pileup(pileup: &[u8]) -> PileupCount {
    // Record SNVs, indels and total depth
    let mut snv_count = LinkedHashMap::new();
    let mut ins_count = LinkedHashMap::new();
    let mut del_count = LinkedHashMap::new();
    let mut depth = StrandCount::new();

    // Iterate over the pileup string and parse symbols to get counts
    // Strand information is not recorded as we rely on samtools merging overlapping reads instead
    let mut i = 0;
    let n = pileup.len();
    while i < n {
        let symbol = pileup[i] as char;
        match symbol {
            // Reference base
            '.' | ',' => {
                let count = snv_count.entry('.').or_insert(StrandCount::new());
                if symbol == '.' {
                    count.fwd += 1;
                    depth.fwd += 1;
                } else {
                    count.rev += 1;
                    depth.rev += 1;
                }
                i += 1;
            }
            // Mismatching base
            'A' | 'C' | 'G' | 'T' | 'N' | 'a' | 'c' | 'g' | 't' | 'n' => {
                let count = snv_count
                    .entry(symbol.to_ascii_uppercase())
                    .or_insert(StrandCount::new());
                if symbol.is_ascii_uppercase() {
                    count.fwd += 1;
                    depth.fwd += 1;
                } else {
                    count.rev += 1;
                    depth.rev += 1;
                }
                i += 1;
            }
            // Insertion or deletion
            '+' | '-' => {
                let mut n_digits = 1;
                while pileup[i + n_digits + 1] >= b'0' && pileup[i + n_digits + 1] <= b'9' {
                    n_digits += 1;
                }
                let indel_len = &pileup[i + 1..i + 1 + n_digits];
                let indel_len: usize = str::from_utf8(indel_len).unwrap().parse().unwrap();
                let indel_start = i + 1 + n_digits;
                let indel_seq = &pileup[indel_start..indel_start + indel_len];
                let pos_strand = indel_seq[0].is_ascii_uppercase();
                let indel_seq = str::from_utf8(indel_seq).unwrap().to_ascii_uppercase();
                let count = match symbol {
                    '+' => ins_count
                        .entry(String::from(indel_seq))
                        .or_insert(StrandCount::new()),
                    '-' => del_count
                        .entry(String::from(indel_seq))
                        .or_insert(StrandCount::new()),
                    _ => panic!(),
                };
                if pos_strand {
                    count.fwd += 1;
                } else {
                    count.rev += 1;
                }
                i += 1 + n_digits + indel_len;
            }
            // Deleted base & reference skip symbols
            '*' | '>' => {
                depth.fwd += 1;
                i += 1;
            }
            '#' | '<' => {
                depth.rev += 1;
                i += 1;
            }
            // Skip start-of-read event including the mapping quality symbol
            '^' => i += 2,
            // Skip end-of-read symbol
            '$' => i += 1,
            // Unknown symbol
            _ => panic!("Unknown symbol in pileup"),
        }
    }

    // Empty locus returns "*" instead of empty string, temp hack
    if pileup == b"*" {
        depth = StrandCount::new();
    }

    // Return counts from this pileup string
    PileupCount {
        snv_count,
        ins_count,
        del_count,
        depth,
    }
}

// A locus in the genome
#[derive(Debug, Clone)]
struct Locus {
    chrom: String,
    coord: u64,
    ref_base: char,
}

// Process one line of output from samtools mpileup
fn process_mpileup_line(line: &str) -> (Locus, Vec<PileupCount>) {
    // Split fields and individual samples from the mpileup output line
    let mut pileup_counts = Vec::new();
    let fields: Vec<&str> = line.split('\t').collect();
    let pileups = &fields[3..];
    let chrom = String::from(fields[0]);
    let coord = fields[1].parse().unwrap();
    let ref_base = fields[2].chars().next().unwrap();
    let locus = Locus {
        chrom,
        coord,
        ref_base,
    };
    eprintln!(
        "\nLocus: {} {} {}",
        locus.chrom, locus.coord, locus.ref_base
    );
    // Count SNVs and indels from each sample separately
    for i in (0..pileups.len()).step_by(3) {
        let depth: u64 = pileups[i].parse().unwrap();
        let pileup = pileups[i + 1];
        let pileup_count = count_pileup(pileup.as_bytes());
        eprintln!(
            "SNV counts {:?}. Ins counts {:?}. Del counts {:?}",
            pileup_count.snv_count, pileup_count.ins_count, pileup_count.del_count
        );
        assert_eq!(depth, pileup_count.depth.total());
        pileup_counts.push(pileup_count)
    }
    (locus, pileup_counts)
}

// Possible reference mismatchs at a locus
#[derive(Debug, Clone)]
enum AltAllele {
    SNV(char),
    Ins(String),
    Del(String),
}

impl fmt::Display for AltAllele {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AltAllele::SNV(x) => write!(f, "{}", x),
            AltAllele::Ins(x) => write!(f, "+{}", x),
            AltAllele::Del(x) => write!(f, "-{}", x),
        }
    }
}

// Chooses the most common ref mismatch as the alt allele (no multi-allelics)
fn choose_alt_alleles(pileup_counts: &[PileupCount], min_count: u64) -> Vec<Option<AltAllele>> {
    let mut alt_alleles = Vec::new();
    for pileup_count in pileup_counts {
        let mut top_count = 0;
        let mut top_allele = None;
        for (snv, count) in &pileup_count.snv_count {
            if *snv == '.' || count.total() < min_count {
                continue;
            }
            if count.total() > top_count {
                top_allele = Some(AltAllele::SNV(*snv));
                top_count = count.total();
            }
        }
        for (ins, count) in &pileup_count.ins_count {
            if count.total() > top_count {
                top_allele = Some(AltAllele::Ins(String::from(ins)));
                top_count = count.total();
            }
        }
        for (del, count) in &pileup_count.del_count {
            if count.total() > top_count {
                top_allele = Some(AltAllele::Del(String::from(del)));
                top_count = count.total();
            }
        }
        alt_alleles.push(top_allele);
    }
    alt_alleles
}

// Get total counts of a single variant from pileup counts
fn lookup_alt_total_count(pileup_count: &PileupCount, alt_allele: &AltAllele) -> u64 {
    match alt_allele {
        AltAllele::SNV(x) => pileup_count
            .snv_count
            .get(x)
            .unwrap_or(&EMPTY_STRANDCOUNTS)
            .total(),
        AltAllele::Ins(x) => pileup_count
            .ins_count
            .get(x)
            .unwrap_or(&EMPTY_STRANDCOUNTS)
            .total(),
        AltAllele::Del(x) => pileup_count
            .del_count
            .get(x)
            .unwrap_or(&EMPTY_STRANDCOUNTS)
            .total(),
    }
}

// Get total counts of a single variant from pileup counts
fn lookup_alt_strand_count(pileup_count: &PileupCount, alt_allele: &AltAllele) -> StrandCount {
    let strand_counts = match alt_allele {
        AltAllele::SNV(x) => pileup_count.snv_count.get(x).unwrap_or(&EMPTY_STRANDCOUNTS),
        AltAllele::Ins(x) => pileup_count.ins_count.get(x).unwrap_or(&EMPTY_STRANDCOUNTS),
        AltAllele::Del(x) => pileup_count.del_count.get(x).unwrap_or(&EMPTY_STRANDCOUNTS),
    };
    *strand_counts
}

// Calculate Fisher's exact test p-values for each sample at this locus
fn calc_fisher_pvals(
    pileup_counts: &[PileupCount],
    alt_alleles: &[Option<AltAllele>],
    strand: Strand,
) -> (Vec<f64>, Vec<[u64; 4]>) {
    let mut fisher_pvals = Vec::new();
    let mut count_tables = Vec::new();
    let n = pileup_counts.len();
    assert_eq!(n, alt_alleles.len());
    for i in 0..n {
        if let Some(alt_allele) = &alt_alleles[i] {
            let alt_count_this =
                lookup_alt_strand_count(&pileup_counts[i], alt_allele).of_strand(strand);
            let oth_count_this = pileup_counts[i].depth.of_strand(strand) - alt_count_this;
            let mut alt_count_rest = 0;
            let mut oth_count_rest = 0;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let alt_count_that =
                    lookup_alt_strand_count(&pileup_counts[j], alt_allele).of_strand(strand);
                alt_count_rest += alt_count_that;
                let oth_count_that = pileup_counts[j].depth.of_strand(strand) - alt_count_that;
                oth_count_rest += oth_count_that;
            }

            let p = fishers_exact_test(
                alt_count_this as i64,
                oth_count_this as i64,
                alt_count_rest as i64,
                oth_count_rest as i64,
            );
            fisher_pvals.push(p);
            let count_table = [
                alt_count_this,
                oth_count_this,
                alt_count_rest,
                oth_count_rest,
            ];
            count_tables.push(count_table);
        } else {
            fisher_pvals.push(1.0);
            count_tables.push([0, 0, 0, 0]);
        }
    }
    (fisher_pvals, count_tables)
}

fn odds_ratio(x: [u64; 4]) -> f64 {
    (x[0] as f64 * x[3] as f64) / (x[2] as f64 * x[1] as f64)
}

// Calculate strand bias metrics for each sample at this locus
fn calc_strand_bias(
    pileup_counts: &[PileupCount],
    alt_alleles: &[Option<AltAllele>],
) -> (Vec<f64>, Vec<f64>, Vec<[u64; 4]>) {
    let mut strand_pvals = Vec::new();
    let mut strand_odds_ratios = Vec::new();
    let mut count_tables = Vec::new();
    let n = pileup_counts.len();
    assert_eq!(n, alt_alleles.len());
    for i in 0..n {
        if let Some(alt_allele) = &alt_alleles[i] {
            let StrandCount {
                fwd: alt_fwd,
                rev: alt_rev,
            } = lookup_alt_strand_count(&pileup_counts[i], alt_allele);
            let StrandCount {
                fwd: depth_fwd,
                rev: depth_rev,
            } = pileup_counts[i].depth;
            let oth_fwd = depth_fwd - alt_fwd;
            let oth_rev = depth_rev - alt_rev;
            let p = fishers_exact_test(
                alt_fwd as i64,
                alt_rev as i64,
                oth_fwd as i64,
                oth_rev as i64,
            );
            strand_pvals.push(p);
            let count_table = [alt_fwd, alt_rev, oth_fwd, oth_rev];
            strand_odds_ratios.push(odds_ratio(count_table));
            count_tables.push(count_table);
        } else {
            strand_pvals.push(1.0);
            strand_odds_ratios.push(1.0);
            count_tables.push([0, 0, 0, 0]);
        }
    }
    (strand_pvals, strand_odds_ratios, count_tables)
}

fn median(l: &[f64]) -> f64 {
    let mut l = l.to_vec();
    l.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = l.len() / 2;
    if l.len() % 2 == 1 {
        l[mid]
    } else {
        (l[mid] + l[mid + 1]) / 2.0
    }
}

fn mad(l: &[f64]) -> f64 {
    let m = median(&l);
    let abs_devs: Vec<f64> = l.iter().map(|x| (x - m).abs()).collect();
    median(&abs_devs)
}

// Calculate normal distribution p-vals for each sample at this locus
fn calc_normal_pvals(pileup_counts: &[PileupCount], alt_alleles: &[Option<AltAllele>]) -> Vec<f64> {
    let mut normal_pvals = Vec::new();
    let n = pileup_counts.len();
    assert_eq!(n, alt_alleles.len());
    for i in 0..n {
        if let Some(alt_allele) = &alt_alleles[i] {
            let mut log_allelic_fracs = Vec::new();
            for j in 0..n {
                let alt_count = lookup_alt_total_count(&pileup_counts[j], alt_allele);
                let depth = pileup_counts[i].depth.total();
                let allelic_frac = if depth == 0 {
                    0.0
                } else {
                    alt_count as f64 / depth as f64
                };
                let log_allelic_frac = (allelic_frac + 1e-6).log2();
                log_allelic_fracs.push(log_allelic_frac);
            }
            let this_log_allelic_frac = log_allelic_fracs.remove(i);
            let median_log_af = median(&log_allelic_fracs);
            let mad_log_af = mad(&log_allelic_fracs);
            if mad_log_af == 0.0 {
                if this_log_allelic_frac > median_log_af {
                    normal_pvals.push(0.0);
                } else {
                    normal_pvals.push(1.0);
                }
            } else {
                let dist = Normal::new(median_log_af, mad_log_af * 1.4826).unwrap();
                let normal_pvals_est = 1.0 - dist.cdf(this_log_allelic_frac);
                normal_pvals.push(normal_pvals_est);
            }
        } else {
            normal_pvals.push(1.0);
        }
    }
    normal_pvals
}

// Iterator struct for consuming samtools mpileup's output
struct SamtoolsMpileup {
    stdout_iterator: Box<dyn Iterator<Item = std::io::Result<String>>>,
}

impl SamtoolsMpileup {
    fn new(ref_file: &str, bed_file: &str, bam_files: &[&str]) -> SamtoolsMpileup {
        let mut samtools_args: Vec<&str> =
            "mpileup -R -B -Q 60 -q 60 -d 999999".split(' ').collect();
        samtools_args.extend(vec!["-f", ref_file, "-l", bed_file]);
        samtools_args.extend(bam_files);
        let samtools_proc = Command::new("samtools")
            .args(samtools_args)
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        let samtools_reader = BufReader::new(samtools_proc.stdout.unwrap());
        let stdout_iterator = Box::new(samtools_reader.lines());
        SamtoolsMpileup { stdout_iterator }
    }
}

impl Iterator for SamtoolsMpileup {
    type Item = (Locus, Vec<PileupCount>);

    fn next(&mut self) -> Option<Self::Item> {
        match self.stdout_iterator.next() {
            Some(line) => Some(process_mpileup_line(&line.unwrap())),
            None => None,
        }
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
enum NoiseContext {
    SNV(char, char),
    Ins(char),
    Del(char),
}

fn init_noise_counts(
    n_samples: usize,
) -> (Vec<HashMap<NoiseContext, u64>>, Vec<HashMap<char, u64>>) {
    let mut noise_counts: Vec<HashMap<NoiseContext, u64>> = Vec::new();
    let mut cum_depths: Vec<HashMap<char, u64>> = Vec::new();
    for _ in 0..n_samples {
        noise_counts.push(HashMap::new());
        cum_depths.push(HashMap::new());
    }
    (noise_counts, cum_depths)
}

fn update_noise_counts(
    locus: &Locus,
    pileup_counts: &[PileupCount],
    noise_counts: &mut [HashMap<NoiseContext, u64>],
    cum_depths: &mut [HashMap<char, u64>],
) {
    let ref_base = locus.ref_base;
    let n_samples = pileup_counts.len();
    for i in 0..n_samples {
        // Get SNV counts from this locus
        let pileup_count = &pileup_counts[i];
        let noise_count = &mut noise_counts[i];
        for alt_base in ['A', 'C', 'G', 'T'] {
            if alt_base == ref_base {
                continue;
            }
            let pileup_snv_n = pileup_count
                .snv_count
                .get(&alt_base)
                .unwrap_or(&EMPTY_STRANDCOUNTS);
            let noise_n = noise_count
                .entry(NoiseContext::SNV(ref_base, alt_base))
                .or_insert(0);
            *noise_n += pileup_snv_n.total();
        }
        // Get insertion counts from this locus
        let mut pileup_ins_n = 0;
        for (_ins, count) in &pileup_count.ins_count {
            pileup_ins_n += count.total();
        }
        let ins_noise_n = noise_count.entry(NoiseContext::Ins(ref_base)).or_insert(0);
        *ins_noise_n += pileup_ins_n;
        // Get deletion counts from this locus
        let mut pileup_del_n = 0;
        for (_del, count) in &pileup_count.del_count {
            pileup_del_n += count.total();
        }
        let del_noise_n = noise_count.entry(NoiseContext::Del(ref_base)).or_insert(0);
        *del_noise_n += pileup_del_n;
        // Get total depth of this location
        let cum_depth = cum_depths[i].entry(ref_base).or_insert(0);
        *cum_depth += pileup_counts[i].depth.total();
    }
}

fn calc_noise_estimates(
    noise_counts: &[HashMap<NoiseContext, u64>],
    cum_depths: &[HashMap<char, u64>],
) -> Vec<HashMap<NoiseContext, f64>> {
    let mut noise_estimates = Vec::new();
    for i in 0..noise_counts.len() {
        let mut noise_estimate = HashMap::new();
        for (noise_context, noise_n) in &noise_counts[i] {
            let ref_base = match noise_context {
                NoiseContext::SNV(x, _) => x,
                NoiseContext::Ins(x) => x,
                NoiseContext::Del(x) => x,
            };
            let cum_depth = cum_depths[i].get(&ref_base).unwrap_or(&0);
            let frac = if *cum_depth == 0 {
                0.0
            } else {
                *noise_n as f64 / *cum_depth as f64
            };
            noise_estimate.insert(*noise_context, frac);
        }
        noise_estimates.push(noise_estimate);
    }
    noise_estimates
}

#[derive(Debug, Clone)]
struct VariantCall {
    sample_i: usize,
    locus: Locus,
    alt_allele: AltAllele,
    allelic_frac: f64,
    depth: u64,
    count_table: [u64; 4],
    fisher_p: f64,
    normal_p: f64,
    noise_p: f64,
    noise_ratio: f64,
    fisher_p_fwd: f64,
    fisher_p_rev: f64,
    strand_count_table: [u64; 4],
    strand_bias_fisher_p: f64,
    strand_odds_ratio: f64,
    //    filtered: bool
}

fn calc_noise_pval(
    variant_call: &VariantCall,
    noise_estimates: &HashMap<NoiseContext, f64>,
) -> (f64, f64) {
    let ref_allele = variant_call.locus.ref_base;
    let noise_freq = match variant_call.alt_allele {
        AltAllele::SNV(alt_allele) => noise_estimates
            .get(&NoiseContext::SNV(ref_allele, alt_allele))
            .unwrap_or(&0.0),
        AltAllele::Ins(_) => noise_estimates
            .get(&NoiseContext::Ins(ref_allele))
            .unwrap_or(&0.0),
        AltAllele::Del(_) => noise_estimates
            .get(&NoiseContext::Del(ref_allele))
            .unwrap_or(&0.0),
    };
    let alt_count = variant_call.count_table[0];
    let depth = variant_call.depth;
    let dist = Binomial::new(*noise_freq, depth).unwrap();
    let noise_p = 1.0 - dist.cdf(alt_count - 1);
    let allelic_frac = alt_count as f64 / depth as f64;
    let noise_ratio = (allelic_frac / *noise_freq).log10();
    (noise_p, noise_ratio)
}

fn passes_filters(call: &VariantCall, n: usize) -> bool {
    if call.fisher_p > 0.05 / n as f64 {
        false
    } else if call.noise_p > 0.05 / n as f64 {
        false
    } else if call.allelic_frac >= 1.0 / 3.0 {
        false
    } else {
        true
    }
}

fn print_header() {
    print!("sample\t");
    print!("chrom\t");
    print!("coord\t");
    print!("ref\t");
    print!("alt\t");
    print!("allelic_frac\t");
    print!("depth\t");
    print!("n_alt_this\t");
    print!("n_non_alt_this\t");
    print!("n_alt_others\t");
    print!("n_non_alt_others\t");
    print!("fisher_p\t");
    print!("normal_p\t");
    print!("noise_p\t");
    print!("noise_log10_ratio\t");
    print!("n_alt_fwd\t");
    print!("n_alt_rev\t");
    print!("n_oth_fwd\t");
    print!("n_oth_rev\t");
    print!("fisher_p_fwd\t");
    print!("fisher_p_rev\t");
    print!("strand_bias_fisher_p\t");
    print!("strand_bias_or\n");
}

fn print_variant_call(call: &VariantCall, sample_name: &str) {
    print!("{}\t", sample_name);
    print!("{}\t", call.locus.chrom);
    print!("{}\t", call.locus.coord);
    print!("{}\t", call.locus.ref_base);
    print!("{}\t", call.alt_allele);
    print!("{}\t", call.allelic_frac);
    print!("{}\t", call.depth);
    print!("{}\t", call.count_table[0]);
    print!("{}\t", call.count_table[1]);
    print!("{}\t", call.count_table[2]);
    print!("{}\t", call.count_table[3]);
    print!("{}\t", call.fisher_p);
    print!("{}\t", call.normal_p);
    print!("{}\t", call.noise_p);
    print!("{}\t", call.noise_ratio);
    print!("{}\t", call.strand_count_table[0]);
    print!("{}\t", call.strand_count_table[1]);
    print!("{}\t", call.strand_count_table[2]);
    print!("{}\t", call.strand_count_table[3]);
    print!("{}\t", call.fisher_p_fwd);
    print!("{}\t", call.fisher_p_rev);
    print!("{}\t", call.strand_bias_fisher_p);
    print!("{}\n", call.strand_odds_ratio);
}

pub fn call_variants(bams: &[&str]) {
    let mpileup = SamtoolsMpileup::new("ref.fasta", "target.bed", bams);
    let n_samples = bams.len();
    let mut variant_calls = Vec::new();
    let (mut noise_counts, mut cum_depths) = init_noise_counts(n_samples);
    for (locus, pileup_counts) in mpileup {
        update_noise_counts(&locus, &pileup_counts, &mut noise_counts, &mut cum_depths);
        let alt_alleles = choose_alt_alleles(&pileup_counts, 5);
        eprintln!("Alts: {:?}", alt_alleles);
        let (fisher_pvalues, count_tables) =
            calc_fisher_pvals(&pileup_counts, &alt_alleles, Strand::Both);
        let (fisher_pvalues_fwd, _) =
            calc_fisher_pvals(&pileup_counts, &alt_alleles, Strand::Forward);
        let (fisher_pvalues_rev, _) =
            calc_fisher_pvals(&pileup_counts, &alt_alleles, Strand::Reverse);
        let normal_pvalues = calc_normal_pvals(&pileup_counts, &alt_alleles);
        let (strand_pvalues, strand_ors, strand_tables) =
            calc_strand_bias(&pileup_counts, &alt_alleles);
        for i in 0..n_samples {
            if let Some(alt_allele) = &alt_alleles[i] {
                let depth = pileup_counts[i].depth.total();
                let allelic_frac = count_tables[i][0] as f64 / depth as f64;
                let variant_call = VariantCall {
                    sample_i: i,
                    locus: locus.clone(),
                    alt_allele: alt_allele.clone(),
                    allelic_frac: allelic_frac,
                    depth: depth,
                    count_table: count_tables[i],
                    fisher_p: fisher_pvalues[i],
                    normal_p: normal_pvalues[i],
                    noise_p: 1.0,
                    noise_ratio: 0.0,
                    fisher_p_fwd: fisher_pvalues_fwd[i],
                    fisher_p_rev: fisher_pvalues_rev[i],
                    strand_count_table: strand_tables[i],
                    strand_bias_fisher_p: strand_pvalues[i],
                    strand_odds_ratio: strand_ors[i],
                };
                variant_calls.push(variant_call);
            }
        }
    }
    let noise_estimates = calc_noise_estimates(&noise_counts, &cum_depths);
    print_header();
    for variant_call in &variant_calls {
        let sample_i = variant_call.sample_i;
        let (noise_p, noise_ratio) = calc_noise_pval(variant_call, &noise_estimates[sample_i]);
        let variant_call = VariantCall {
            noise_p: noise_p,
            noise_ratio: noise_ratio,
            ..variant_call.clone()
        };
        if passes_filters(&variant_call, variant_calls.len()) {
            print_variant_call(&variant_call, &bams[sample_i]);
        }
    }
}
