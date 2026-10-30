# Supplementary File cohort_inventory. Public sources and analysis units

## Overview and scope

Supplementary File cohort_inventory documents the public microRNA datasets, screening inventory, analysis units, source annotations, and data citations.
Recorded person or specimen identifiers define cross-validation grouping but do not establish undocumented donor relationships across accessions.
The benchmark corpus contains 34 public analysis units comprising 23,192 fitted assay profiles (12,636 cases and 10,556 controls) grouped into 22,938 biological groups (12,509 case groups and 10,429 control groups).
Exclusion of 396 near-duplicate profiles from held-out scoring (394 in the deduplicated 3D-Gene serum unit and 2 in PRJEB30542) leaves 22,796 evaluated profiles (12,634 cases and 10,162 controls).
Thirty-two benchmark units form the primary active-cancer cell-free and extracellular vesicle synthesis (22,714 fitted and 22,318 evaluated profiles across 22,580 biological groups).
Two benchmark units are outside that 32-unit synthesis, namely GSE145176 (acute lymphoblastic leukaemia survivorship, 240 profiles across 120 people) and GSE31568 (whole-blood cellular RNA, 238 profiles across 238 people).
GSE268524 is a separate pancreatobiliary application outside the 34-unit benchmark corpus in which 1,403 technical arrays (668 cancer and 735 control arrays) collapse by median into 477 biological specimens (227 cases and 250 controls).
The deduplicated 3D-Gene serum unit audits 10 source accessions of which 8 contribute canonical rows (GSE211692 with 12,659, GSE106817 with 3,698, GSE113740 with 1,238, GSE113486 with 840, GSE164174 with 288, GSE137140 with 239, GSE112264 with 126, and GSE139031 with 55).
The remaining two audited 3D-Gene accessions (GSE110317 and GSE122497) contribute zero canonical rows after cross-accession deduplication and detector-consistency quarantine.

## Column definitions for cohort_inventory.tsv and cohort_inventory_notes.tsv

cohort_inventory.tsv summarises the 34 benchmark analysis units and the standalone application unit GSE268524.
Analysis unit gives the unit identifier and manuscript reference numbers.
Clinical comparison states the target malignancy or clinical contrast and marks the two corpus-only units outside the 32-unit synthesis.
Specimen / assay states the biofluid fraction and measurement technology.
Profiles and Case / control profiles report the fitted assay profile total and its case and control profile split for each unit.
Evaluated profiles and Evaluated case / control profiles report the held-out evaluated profile total and its case and control split after the 396 near-duplicate exclusions.
Groups and Cases / controls report the biological group total and the biological group case and control counts kept together during grouped cross-validation.
Control composition summarises the clinical composition of the non-cancer comparator arm.
The trailing columns verbatim_table_row, source_path, and source_line record manuscript table provenance.
cohort_inventory_notes.tsv provides the table footnote text in verbatim_text with provenance columns source_path and source_line.

## Column definitions for cohort_inventory_dataset_characteristics.tsv

cohort_inventory_dataset_characteristics.tsv provides unit-level technical and provenance metadata for the 34 benchmark units and GSE268524.
Public analysis unit and Role in this study identify each unit and distinguish the 32 primary synthesis units, the 2 corpus-only benchmark units, and the standalone application.
Cancer comparison, Specimen type, and Measurement platform describe the clinical contrast, biofluid matrix, and assay platform.
Assay profiles, Cancer profiles, and Control profiles give fitted profile counts.
Evaluated assay profiles, Evaluated cancer profiles, and Evaluated control profiles give held-out evaluated profile counts after near-duplicate exclusions.
Independent biological groups, Cancer groups, and Control groups give the grouped donor or specimen counts used for outer fold partitioning.
Measured miRNA features gives the number of microRNA features in the analysis matrix.
Source accessions, Source-study citations, and Number of source accessions list the repository accessions, linked publications, and audited versus contributing accession counts.
Repeated-measure handling and Repository data references describe replicate or multi-platform grouping and reference numbers, followed by provenance columns evidence_file, archive_member, and source_row.

## Column definitions for cohort_inventory_public_dataset_screening.tsv

cohort_inventory_public_dataset_screening.tsv documents all 117 screened candidate records.
Public record, Screening disposition, Scientific screening category, and Analysis unit if included identify each candidate record, its screening decision, its disposition category, and its target analysis unit.
Material used, Material not used, and Scientific reason state which deposited components entered the benchmark, which components were excluded, and the scientific rationale.
Cancer comparison, Measurement platform, Repository platform identifier, and Specimen type describe the clinical contrast, assay technology, platform accession, and biofluid or tissue matrix.
Cancer profiles, Control profiles, and Other profiles report the deposited sample counts by depositor annotation for malignant cancer, healthy or non-cancer control, and other diagnoses or non-target sample types.
Smaller case/control arm equals the minimum of deposited Cancer profiles and deposited Control profiles for every row.
Analysed cancer profiles and Analysed control profiles report the fitted case and control profile counts contributed by each record to its analysis unit, or zero for excluded records.
Deposit to analysis difference reason states in one sentence why deposited and analysed counts differ when subsetting, comparator pooling, deduplication, or array collapse occurred.
Processed expression available, Raw data available, Source paper available, PubMed identifier, and DOI record data and publication availability.
Expression scale, Logarithm base, Required score direction, and Cross-record independence handling record measurement scale and cross-accession grouping rules.
Repository data references, evidence_file, archive_member, and source_row record data citations and archive provenance.

## Column definitions for citations, curation, and reference crosswalk tables

cohort_inventory_accession_citations.tsv links each repository accession (accession) to its analysis unit (analysis_units), role (inventory_scope), dataset citation (data_reference_numbers and data_citation_verbatim), publication citation (source_publication_reference_numbers, source_publication_verbatim, and source_publication_identifiers_verbatim), and provenance fields (data_citation_state, source_publication_state, citation_source, identifier_source, and identifier_members).
cohort_inventory_source_curation.tsv records structured curation entries by analysis unit (unit), source accession (accession), curated attribute name (field), curated description (value), evidence location (source), supporting excerpt (quote), confidence rating (confidence), verification status (second_reading), and archive provenance (evidence_file, archive_member, and source_row).
cohort_inventory_reference_list.tsv lists each numbered citation in reference_number and verbatim_reference with its source file in source_path.
cohort_inventory_reference_number_crosswalk.tsv maps each analysis_unit and its source_accessions between archived_inventory_reference_numbers and specified_manuscript_reference_numbers with provenance columns archived_inventory_source and manuscript_reference_source.
