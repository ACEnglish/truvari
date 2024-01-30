use noodles_vcf::record::genotypes::sample::value::Genotype;
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
