use noodles_vcf::{self as vcf};
mod comparisons;
mod matching;
mod types;

fn main() {
    let mut reader = vcf::reader::Builder::default()
        .build_from_path("sample.vcf.gz")
        .expect("Unable to parse vcf");
    let header = reader.read_header().expect("Unable to parse header");
    let mut mmatch1 = matching::MatchResult {
        ..Default::default()
    };
    mmatch1.base_gt_count = 4;
    mmatch1.score = Some(4.0);
    let mut mmatch2 = matching::MatchResult {
        ..Default::default()
    };
    mmatch2.score = Some(5.0);
    let mut mmatch3 = matching::MatchResult {
        ..Default::default()
    };
    mmatch3.state = true;
    let mut parts = vec![mmatch1, mmatch2, mmatch3];
    parts.sort();
    parts.reverse();
    for i in parts {
        println!("{:?}", i);
    }

    let mat = matching::Matcher { 
        ..Default::default()
    };
    let mut up_record = vcf::Record::default();
    reader
        .read_record(&header, &mut up_record)
        .expect("Unable to parse record");
    for result in reader.records(&header) {
        let dn_record = result.expect("Unable to parse record");
        // Need to put a sequence resolved guard on seq_similarity
        if mat.filter_call(&dn_record, false) {
            //println!("filtering {:?}", dn_record);
            continue;
        }
        println!(
            "size={:?} type={:?} bound={:?} dist={:?} ovl={:?} szsim={:?} gtcmp={:?} pres={:?} filt={:?} sim={:?}",
            comparisons::entry_size(&dn_record),
            comparisons::entry_variant_type(&dn_record),
            comparisons::entry_boundaries(&dn_record, false),
            comparisons::entry_distance(&up_record, &dn_record),
            comparisons::entry_reciprocal_overlap(&up_record, &dn_record),
            comparisons::entry_size_similarity(&up_record, &dn_record),
            comparisons::entry_gt_comp(&up_record, &dn_record, 0, 0),
            comparisons::entry_is_present(&dn_record, 0),
            comparisons::entry_is_filtered(&dn_record),
            comparisons::entry_seq_similarity(&up_record, &dn_record),
        );
        let x = mat.build_match(&up_record, &dn_record, Some("".to_string()), false, false);
        println!("{:?}", x);
        up_record = dn_record;
    }
}
