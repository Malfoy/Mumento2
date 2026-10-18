# Human chromosome 19 benchmark inputs

These are real haplotype-resolved human assemblies from the Human Pangenome Reference Consortium (HPRC). Each input contains the chromosome 19 sequence from one assembled haplotype. Haplotype labels 1 and 2 distinguish chromosome copies and do not universally imply paternal or maternal origin.

The frozen 512-input collection contains 376 haplotypes from official HPRC Release 2 and 136 from public R3 Verkko v2.3.2 hybrid assembly submissions. It represents 276 individuals: 236 contribute two haplotypes and 40 contribute one. The collection contains 31,558,712,967 bases; input lengths range from 58,700,692 to 65,739,657 bases, averaging 61,638,111 bases.

Release 2 metadata comes from the [official HPRC assembly index](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/assemblies_release2_v1.0.index.csv). Assembly FASTA records were obtained from the public human-pangenomics S3 bucket. For example, the first input is [HG00097 haplotype 1](https://human-pangenomics.s3.us-west-2.amazonaws.com/working/HPRC/HG00097/assemblies/release2/HG00097_hap1_hprc_r2_v1.0.1.fa.gz). Chromosome 19 was identified using UCSC chromosome aliases and the matching FASTA index record.

The added R3 inputs come from the public `submissions/bed38f3f-cac4-4703-9852-181bc362cea1--R3_verkko-v2.3.2_hybrid_assembly/` S3 submission, using its `verkko-hi-c` assemblies. For example, [HG02497 haplotype 1](https://human-pangenomics.s3.us-west-2.amazonaws.com/submissions/bed38f3f-cac4-4703-9852-181bc362cea1--R3_verkko-v2.3.2_hybrid_assembly/HG02497/verkko-hi-c/HG02497.assembly.haplotype1.fasta). These are public assembly submissions rather than an assertion that they constitute a final R3 release. Chromosome identity comes from each assembly's supplied chromosome assignment file. Record headers, lengths, object sizes and ETags were checked during extraction.

Eligibility required a unique chromosome-19 assignment and a record length of 40–90 Mb; this does not establish that every sequence is telomere-to-telomere. The original Release 2 prefix was preserved: eligible haplotype 1 records in sample-ID order, then eligible haplotype 2 records. R3 haplotypes 1 and 2 from new individuals were appended in sample-ID order. Smaller collections are nested prefixes. Selection is based on availability and chromosome assignment, rather than balanced population sampling.

| Haplotypes | Individuals | Total sequence (Gb, decimal) |
|---:|---:|---:|
| 2 | 2 | 0.1234 |
| 4 | 4 | 0.2458 |
| 8 | 8 | 0.4899 |
| 16 | 16 | 0.9820 |
| 32 | 32 | 1.9625 |
| 64 | 64 | 3.9328 |
| 128 | 128 | 7.8723 |
| 256 | 195 | 15.7821 |
| 512 | 276 | 31.5587 |

The [public input manifest](human-chr19-inputs-20261018.csv) lists all 512 inputs in benchmark order, including assembly URLs and checksums. The local frozen manifest records sample ID, haplotype, assembly accession/name, source URL, chromosome record, sequence length, object ETag and checksums. Both tools receive the same collection in the same order. Uniqueness is evaluated within these chromosome-19 inputs; whole-genome uniqueness and performance have not been measured.
