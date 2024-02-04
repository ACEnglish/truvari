use crate::types::Gt;
use crate::comparisons;
use noodles_vcf::{self as vcf};
use std::cmp::Ordering;

#[derive(Debug)]
pub struct MatchResult {
    pub base: Option<vcf::Record>,
    pub comp: Option<vcf::Record>,
    pub base_gt: Option<Gt>,
    pub base_gt_count: u8,
    pub comp_gt: Option<Gt>,
    pub comp_gt_count: u8,
    pub state: bool,
    pub seqsim: Option<f32>,
    pub sizesim: Option<f32>,
    pub ovlpct: Option<f32>,
    pub sizediff: Option<i64>,
    pub st_dist: Option<i64>,
    pub ed_dist: Option<i64>,
    pub gt_match: Option<bool>,
    pub multi: Option<bool>,
    pub score: Option<f32>,
    pub matid: Option<String>,
}

impl MatchResult {
    pub fn new() -> Self {
        MatchResult {
            base: None,
            comp: None,
            base_gt: None,
            base_gt_count: 0,
            comp_gt: None,
            comp_gt_count: 0,
            state: false,
            seqsim: None,
            sizesim: None,
            ovlpct: None,
            sizediff: None,
            st_dist: None,
            ed_dist: None,
            gt_match: None,
            multi: None,
            score: None,
            matid: None,
        }
    }

    pub fn calc_score(mut self) {
        self.score = match (self.seqsim, self.sizesim, self.ovlpct) {
            (Some(seq), Some(size), Some(ovl)) => Some((seq + size + ovl) / 3.0 * 100.0),
            _ => None,
        }
    }
}

impl PartialOrd for MatchResult {
    fn partial_cmp(&self, other: &MatchResult) -> Option<Ordering> {
        Some(other.cmp(self))
    }
}

impl Ord for MatchResult {
    fn cmp(&self, other: &MatchResult) -> Ordering {
        if self.state != other.state {
            if self.state < other.state {
                Ordering::Greater
            } else {
                Ordering::Less
            }
        } else {
            match (self.score, other.score) {
                (Some(sscore), Some(oscore)) => {
                    if sscore < oscore {
                        Ordering::Greater
                    } else if sscore > oscore {
                        Ordering::Less
                    } else {
                        Ordering::Equal
                    }
                }
                _ => Ordering::Equal, // Not good?
            }
        }
    }
}

impl PartialEq for MatchResult {
    fn eq(&self, other: &MatchResult) -> bool {
        (self.state == other.state) & (self.score == other.score)
    }
}

impl Eq for MatchResult {}

#[derive(Debug)]
pub struct MatchParams {
    pub reference: Option<String>,
    pub refdist: usize ,
    pub pctseq: f32,
    pub pctsize: f32,
    pub pctovl: f32,
    pub typeignore: bool,
    pub chunksize: usize,
    pub bSample: usize,
    pub cSample: usize,
    pub dup_to_ins: bool,
    pub sizemin: usize,
    pub sizefilt: usize,
    pub sizemax: usize,
    pub passonly: bool,
    pub no_ref: char,
    // pub pick // Picker Object
    pub ignore_monref: bool,
    pub check_multi: bool,
    pub check_monref: bool,
}

impl Default for MatchParams {
    fn default() -> MatchParams {
        MatchParams {
            reference: None,
            refdist: 500,
            pctseq:  0.70,
            pctsize:  0.70,
            pctovl:  0.0,
            typeignore:  false,
            chunksize:  1000,
            bSample:  0,
            cSample:  0,
            dup_to_ins:  false,
            sizemin:  50,
            sizefilt:  30,
            sizemax:  50000,
            passonly:  false,
            no_ref:  'o', // maybe should be some enum
            // params.pick:  'single',
            ignore_monref: true,
            check_multi:  true,
            check_monref:  true,
        }
    }
}

#[derive(Debug)]
pub struct Matcher {
    pub params: MatchParams

}

impl Matcher {
    pub fn new(m_params: MatchParams) -> Self {
        Matcher { params: m_params }
    }

    pub fn filter_call(&self, entry: &vcf::Record, base: bool) -> bool {
        if self.params.check_monref & (entry.alternate_bases().len() == 0) {
            return true;
        }

        if self.params.check_multi & (entry.alternate_bases().len() > 1) {
            //panic and exit
            return true;
        }

        if self.params.passonly & comparisons::entry_is_filtered(&entry) {
            return true;
        }

        let size = comparisons::entry_size(&entry);
        if (size > self.params.sizemax) || (base & (size < self.params.sizemin)) || (!base & (size < self.params.sizefilt)) {
            return true;
        }
        let (samp, prefix) = if base { (self.params.bSample, 'b') } else { (self.params.cSample, 'c') };
        if (self.params.no_ref == prefix) || (self.params.no_ref == 'a') { // self.params.pick == 'ac'
            return true;
        }

        false
    }
}
    //check_monref and alts is None
    //# check_multi and 
/* TODO
class Matcher():
    make_match_params with defaults...
    make_match_params_from_args with the thing
    def filter_call(self, entry, base=False):
    def build_match(self, base, comp, matid=None, skip_gt=False, short_circuit=False):
def file_zipper(*start_files):
def chunker(matcher, *files):
 */
