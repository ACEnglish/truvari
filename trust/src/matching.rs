use crate::types::Gt;
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

/* TODO
pub struct Matcher {
}
class Matcher():
    make_match_params with defaults...
    make_match_params_from_args with the thing
    def filter_call(self, entry, base=False):
    def build_match(self, base, comp, matid=None, skip_gt=False, short_circuit=False):
def file_zipper(*start_files):
def chunker(matcher, *files):
 */
