use edlib_rs::edlibrs::{edlibAlignRs, EdlibAlignConfigRs};
use noodles_vcf::{
    self as vcf, record::alternate_bases::allele, record::genotypes::sample::value::Genotype,
    record::info::field, record::Filters,
};
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
pub enum Svtype {
    Ins,
    Del,
    Dup,
    Inv,
    Snp,
    Unk,
    //Repl should be one
}

impl FromStr for Svtype {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "DEL" => Ok(Svtype::Del),
            "INS" => Ok(Svtype::Ins),
            "DUP" => Ok(Svtype::Dup),
            "INV" => Ok(Svtype::Inv),
            "SNP" => Ok(Svtype::Snp),
            _ => Ok(Svtype::Unk),
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum Gt {
    Ref,
    Het,
    Hom,
    Non,
    Unk,
    //Hemi should be a thing
}

impl Gt {
    pub fn new(gt: Genotype) -> Self {
        let g1 = match gt.first().unwrap().position() {
            Some(g) => g.to_string().chars().next().unwrap(),
            None => '.',
        };
        // This isn't going to work if its hemi?
        let g2 = match gt.last().unwrap().position() {
            Some(g) => g.to_string().chars().next().unwrap(),
            None => '.',
        };
        if g1 == g2 {
            match g1 {
                '1' => Gt::Hom,
                '0' => Gt::Ref,
                '.' => Gt::Non,
                _ => Gt::Unk,
            }
        } else if (g1 == '1') || (g2 == '1') {
            Gt::Het
        } else {
            Gt::Unk
        }
    }
}

pub fn entry_boundaries(entry: &vcf::Record, ins_inflate: bool) -> (usize, usize) {
    let mut start = usize::from(entry.position()) - 1;
    let mut end = usize::from(entry.end().expect("No Variant End"));
    if ins_inflate & (entry_variant_type(entry) == Svtype::Ins) {
        let size = entry_size(entry);
        start -= size / 2;
        end += size / 2;
    }
    (start, end)
}

pub fn entry_size(entry: &vcf::Record) -> usize {
    let svlen = entry
        .info()
        .get(&field::Key::from_str("SVLEN").expect("No SVLEN INFO"));
    if let Some(Some(field::Value::Integer(svlen))) = svlen {
        return svlen.unsigned_abs() as usize;
    } else if let Some(Some(field::Value::Array(field::value::Array::Integer(svlen)))) = svlen {
        return svlen[0].expect("Bad SVLEN").unsigned_abs() as usize;
    }

    let r_len = entry.reference_bases().len();
    let a_len = match entry.alternate_bases().first() {
        Some(allele::Allele::Bases(alt)) => alt.len(),
        Some(allele::Allele::Symbol(_alt)) => {
            let (start, end) = entry_boundaries(entry, false);
            start.abs_diff(end) + 1 // I don't understand why I have to add 1.
        }
        _ => 0,
    };

    if r_len == a_len {
        if r_len == 1 {
            return 0;
        } else {
            return r_len;
        }
    }

    r_len.abs_diff(a_len)
}

pub fn entry_variant_type(entry: &vcf::Record) -> Svtype {
    match entry
        .info()
        .get(&field::Key::from_str("SVTYPE").expect("Unable to make key"))
    {
        // INFO/SVTYPE
        Some(Some(field::Value::String(svtype))) => svtype.parse().expect("Bad SVTYPE"),
        Some(Some(field::Value::Array(field::value::Array::String(svtype)))) => {
            svtype[0].clone().expect("Bad SVTYPE").parse().unwrap()
        }
        // Direct from REF/ALT
        _ => match entry.alternate_bases().first() {
            Some(allele::Allele::Bases(alt)) => {
                let asz = alt.len();
                let rsz = entry.reference_bases().len();
                if asz > rsz {
                    Svtype::Ins
                } else if asz < rsz {
                    Svtype::Del
                } else if asz == 1 {
                    Svtype::Snp
                } else {
                    Svtype::Unk
                }
            }
            Some(allele::Allele::Symbol(alt)) => {
                Svtype::from_str(&alt.to_string()).expect("Bad Symbolic Alt")
            }
            _ => Svtype::Unk,
        },
    }
}

pub fn entry_distance(entry_a: &vcf::Record, entry_b: &vcf::Record) -> (isize, isize) {
    let (astart, aend) = entry_boundaries(entry_a, false);
    let (bstart, bend) = entry_boundaries(entry_b, false);
    (
        (astart as isize - bstart as isize),
        (aend as isize - bend as isize),
    )
}

pub fn entry_same_variant_type(
    entry_a: &vcf::Record,
    entry_b: &vcf::Record,
    dup_to_ins: bool,
) -> bool {
    let mut a_type = entry_variant_type(entry_a);
    let mut b_type = entry_variant_type(entry_b);
    if dup_to_ins & (a_type == Svtype::Dup) {
        a_type = Svtype::Ins;
    }

    if dup_to_ins & (b_type == Svtype::Dup) {
        b_type = Svtype::Ins;
    }
    a_type == b_type
}

pub fn entry_reciprocal_overlap(entry_a: &vcf::Record, entry_b: &vcf::Record) -> f32 {
    let (astart, aend) = entry_boundaries(entry_a, true);
    let (bstart, bend) = entry_boundaries(entry_b, true);
    reciprocal_overlap(astart, aend, bstart, bend)
}

pub fn reciprocal_overlap(astart: usize, aend: usize, bstart: usize, bend: usize) -> f32 {
    let ovl_start = std::cmp::max(astart, bstart);
    let ovl_end = std::cmp::min(aend, bend);
    if ovl_start < ovl_end {
        return (ovl_end as isize - ovl_start as isize) as f32
            / std::cmp::max(aend - astart, bend - bstart) as f32;
    }
    0.0
}

pub fn entry_size_similarity(entry_a: &vcf::Record, entry_b: &vcf::Record) -> (f32, isize) {
    sizesim(entry_size(entry_a), entry_size(entry_b))
}

pub fn sizesim(size_a: usize, size_b: usize) -> (f32, isize) {
    if ((size_a == 0) || (size_b == 0)) && size_a == size_b {
        return (1.0, 0);
    }
    let pct = std::cmp::max(std::cmp::min(size_a, size_b), 1) as f32
        / std::cmp::max(std::cmp::max(size_a, size_b), 1) as f32;
    let diff = size_a as isize - size_b as isize;
    (pct, diff)
}

pub fn overlap_percent(astart: usize, aend: usize, bstart: usize, bend: usize) -> f32 {
    if (astart >= bstart) & (aend <= bend) {
        return 1.0;
    }
    let ovl_start = std::cmp::max(astart, bstart);
    let ovl_end = std::cmp::min(aend, bend);
    if ovl_start < ovl_end {
        return (ovl_end - ovl_start) as f32 / (aend - astart) as f32;
    }
    0.0
}

pub fn overlaps(s1: usize, e1: usize, s2: usize, e2: usize) -> bool {
    std::cmp::max(s1, s2) < std::cmp::min(e1, e2)
}

pub fn entry_gt_comp(
    entry_a: &vcf::Record,
    entry_b: &vcf::Record,
    sample_a: usize,
    sample_b: usize,
) -> bool {
    let gt_a = Gt::new(
        entry_a
            .genotypes()
            .get_index(sample_a)
            .expect("Bad sample index")
            .genotype()
            .expect("Unable to parse genotype")
            .unwrap(),
    );
    let gt_b = Gt::new(
        entry_b
            .genotypes()
            .get_index(sample_b)
            .expect("Bad sample index")
            .genotype()
            .expect("Unable to parse genotype")
            .unwrap(),
    );
    gt_a == gt_b
}

pub fn entry_is_present(entry: &vcf::Record, sample: usize) -> bool {
    let gt = Gt::new(
        entry
            .genotypes()
            .get_index(sample)
            .expect("Bad sample index")
            .genotype()
            .expect("Unable to parse genotype")
            .unwrap(),
    );
    return (gt == Gt::Het) || (gt == Gt::Hom);
}

pub fn entry_is_filtered(entry: &vcf::Record) -> bool {
    // Filter is None or PASS not in filter
    match entry.filters() {
        Some(map) => *map != Filters::Pass,
        None => false,
    }
}

pub fn seqsim(seq_a: &String, seq_b: &String) -> f32 {
    let align_res = edlibAlignRs(
        seq_a.as_bytes(),
        seq_b.as_bytes(),
        &EdlibAlignConfigRs::default(),
    );
    let totlen: f32 = (seq_a.len() + seq_b.len()) as f32;
    (totlen - align_res.editDistance as f32) / totlen
}

pub fn entry_seq_similarity(entry_a: &vcf::Record, entry_b: &vcf::Record) -> f32 {
    let mut a_seq = match entry_variant_type(entry_a) {
        Svtype::Del => entry_a.reference_bases().to_string(),
        _ => entry_a
            .alternate_bases()
            .first()
            .expect("Should only be called on sequence resolved alts")
            .to_string(),
    };
    let mut b_seq = match entry_variant_type(entry_b) {
        Svtype::Del => entry_b.reference_bases().to_string(),
        _ => entry_b
            .alternate_bases()
            .first()
            .expect("Should only be called on sequence resolved alts")
            .to_string(),
    };

    let (mut st_dist, mut ed_dist) = entry_distance(&entry_a, &entry_b);
    if (st_dist == 0) || (ed_dist == 0) {
        return seqsim(&a_seq, &b_seq);
    }
    if st_dist < 0 {
        st_dist *= -1
    } else {
        std::mem::swap(&mut a_seq, &mut b_seq);
    }

    seqsim(&a_seq, &b_seq)
        .max(unroll_compare(&a_seq, &b_seq, st_dist as usize, true))
        .max(unroll_compare(&b_seq, &a_seq, st_dist as usize, false))
}

pub fn unroll_compare(a_seq: &String, b_seq: &String, p: usize, up: bool) -> f32 {
    let b_len = b_seq.len() as usize;
    let f = p % b_len;
    let position = (b_len - f) as usize; // I'm worried about signs here
    if position >= b_len {
        return 0.0 // should never be called unless Symbolic alts are present, in which case we
                   // can't compare
    }
    // If up, a_seq is upstream of b_seq
    let rolled = match up {
        true => format!("{}{}", &b_seq[position..], &b_seq[..position]),
        false => format!("{}{}", &b_seq[..position], &b_seq[position..])
    };
    seqsim(&a_seq, &rolled)
}

/* TODO
Only if we're doing vcf2df (probably not...?)
def entry_to_hash(entry, hasher=hashlib.sha1):
def entry_to_key(entry, prefix="", bounds=False):

Deprecating reference context comparison
def create_pos_haplotype(a1, a2, ref, min_len=0):
def entry_shared_ref_context(entryA, entryB, ref, use_ref_seq=False, min_len=0):
 */
