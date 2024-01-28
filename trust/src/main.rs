use noodles_vcf::{self as vcf};
mod comparisons;

fn main() {
    let mut reader = vcf::reader::Builder::default()
        .build_from_path("sample.vcf.gz")
        .expect("Unable to parse vcf");
    let header = reader.read_header().expect("Unable to parse header");

    let mut up_record = vcf::Record::default();
    reader
        .read_record(&header, &mut up_record)
        .expect("Unable to parse record");

    for result in reader.records(&header) {
        let dn_record = result.expect("Unable to parse record");
        println!(
            "size={:?} type={:?} bound={:?} dist={:?} ovl={:?} szsim={:?} gtcmp={:?} pres={:?} filt={:?}",
            comparisons::entry_size(&dn_record),
            comparisons::entry_variant_type(&dn_record),
            comparisons::entry_boundaries(&dn_record, false),
            comparisons::entry_distance(&up_record, &dn_record),
            comparisons::entry_reciprocal_overlap(&up_record, &dn_record),
            comparisons::entry_size_similarity(&up_record, &dn_record),
            comparisons::entry_gt_comp(&up_record, &dn_record, 0, 0),
            comparisons::entry_is_present(&dn_record, 0),
            comparisons::entry_is_filtered(&dn_record),
        );
        up_record = dn_record;
    }
}
