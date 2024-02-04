use crate::comparisons;
use crate::types::Gt;
use noodles_vcf::{self as vcf};
use std::cmp::Ordering;

#[derive(Debug)]
pub struct MatchResult<'a> {
    pub base: Option<&'a vcf::Record>,
    pub comp: Option<&'a vcf::Record>,
    pub base_gt: Option<Gt>,
    pub base_gt_count: u8,
    pub comp_gt: Option<Gt>,
    pub comp_gt_count: u8,
    pub state: bool,
    pub seqsim: Option<f32>,
    pub sizesim: Option<f32>,
    pub ovlpct: Option<f32>,
    pub sizediff: Option<isize>,
    pub st_dist: Option<isize>,
    pub ed_dist: Option<isize>,
    pub gt_match: Option<bool>,
    pub multi: Option<bool>,
    pub score: Option<f32>,
    pub matid: Option<String>,
}

impl Default for MatchResult<'_> {
    fn default() -> MatchResult<'static> {
        MatchResult {
            base: None,
            comp: None,
            base_gt: None,
            base_gt_count: 0,
            comp_gt: None,
            comp_gt_count: 0,
            state: true,
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
}

impl MatchResult<'_> {
    pub fn calc_score(&mut self) {
        self.score = match (self.seqsim, self.sizesim, self.ovlpct) {
            (Some(seq), Some(size), Some(ovl)) => Some((seq + size + ovl) / 3.0 * 100.0),
            _ => None,
        }
    }
}

impl PartialOrd for MatchResult<'_> {
    fn partial_cmp(&self, other: &MatchResult) -> Option<Ordering> {
        Some(other.cmp(self))
    }
}

impl Ord for MatchResult<'_> {
    fn cmp(&self, other: &MatchResult) -> Ordering {
        if self.state != other.state {
            if !self.state & other.state {
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

impl PartialEq for MatchResult<'_> {
    fn eq(&self, other: &MatchResult) -> bool {
        (self.state == other.state) & (self.score == other.score)
    }
}

impl Eq for MatchResult<'_> {}

#[derive(Debug)]
pub struct Matcher {
    pub reference: Option<String>,
    pub refdist: usize,
    pub pctseq: f32,
    pub pctsize: f32,
    pub pctovl: f32,
    pub typeignore: bool,
    pub chunksize: usize,
    pub b_sample: usize,
    pub c_sample: usize,
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

impl Default for Matcher {
    fn default() -> Matcher {
        Matcher {
            reference: None,
            refdist: 500,
            pctseq: 0.70,
            pctsize: 0.70,
            pctovl: 0.0,
            typeignore: false,
            chunksize: 1000,
            b_sample: 0,
            c_sample: 0,
            dup_to_ins: false,
            sizemin: 50,
            sizefilt: 30,
            sizemax: 50000,
            passonly: false,
            no_ref: 'o', // maybe should be some enum
            // params.pick:  'single',
            ignore_monref: true,
            check_multi: true,
            check_monref: true,
        }
    }
}

impl Matcher {
    pub fn filter_call(&self, entry: &vcf::Record, base: bool) -> bool {
        if self.check_monref & (entry.alternate_bases().len() == 0) {
            return true;
        }

        if self.check_multi & (entry.alternate_bases().len() > 1) {
            //panic and exit
            return true;
        }

        if self.passonly & comparisons::entry_is_filtered(entry) {
            return true;
        }

        let size = comparisons::entry_size(entry);
        if (size > self.sizemax)
            || (base & (size < self.sizemin))
            || (!base & (size < self.sizefilt))
        {
            return true;
        }
        let (samp, prefix) = if base {
            (self.b_sample, 'b')
        } else {
            (self.c_sample, 'c')
        };
        if (self.no_ref == prefix) || (self.no_ref == 'a') {
            // self.pick == 'ac'
            return true;
        }

        false
    }

    pub fn build_match<'a>(
        &'a self,
        base: &'a vcf::Record,
        comp: &'a vcf::Record,
        matid: Option<String>,
        skip_gt: bool,
        short_circuit: bool,
    ) -> MatchResult {
        let mut ret = MatchResult {
            base: Some(base),
            comp: Some(comp),
            matid,
            ..Default::default()
        };

        if !self.typeignore
            & !comparisons::entry_same_variant_type(base, comp, self.dup_to_ins)
        {
            ret.state = false;
            if short_circuit {
                return ret;
            }
        }

        // Don't want to call this twice, but also might want to add in insertion inflation...
        // Revisit after the alpha
        let (bstart, bend) = comparisons::entry_boundaries(base, false);
        let (cstart, cend) = comparisons::entry_boundaries(comp, false);

        if !comparisons::overlaps(
            bstart - self.refdist,
            bend + self.refdist,
            cstart,
            cend,
        ) {
            ret.state = false;
            if short_circuit {
                return ret;
            }
        }

        let (szsim, szdiff) = comparisons::entry_size_similarity(base, comp);
        ret.sizesim = Some(szsim);
        ret.sizediff = Some(szdiff);
        if ret.sizesim < Some(self.pctsize) {
            ret.state = false;
            if short_circuit {
                return ret;
            }
        }

        if !skip_gt {
            ret.base_gt = Some(Gt::new(base, self.b_sample));
            //ret.base_gt_count = 1;
            ret.comp_gt = Some(Gt::new(comp, self.c_sample));
            //ret.comp_gt_count = 1;
            // ret.gt_match = ret.base_gt_count.abs_diff(ret.comp_gt_count)
        }

        ret.ovlpct = Some(comparisons::reciprocal_overlap(bstart, bend, cstart, cend));
        if ret.ovlpct < Some(self.pctovl) {
            ret.state = false;
            if short_circuit {
                return ret;
            }
        }
       
        if self.pctseq > 0.0 {
            ret.seqsim = Some(comparisons::entry_seq_similarity(base, comp));
            if ret.seqsim < Some(self.pctseq) {
                ret.state = false;
                if short_circuit {
                    return ret;
                }
            }
        }

        ret.st_dist = Some(bstart as isize - cstart as isize);
        ret.ed_dist = Some(bend as isize - cend as isize);
        ret.calc_score();

        ret
    }
}
