Mock Data Processing for Platanthera project
- _A. thaliana_ reads to be stripped from mock
- _P. zijinensis_ reads to be generated _in silico_ using **NGSNGS** 
- Dilute mock data with _P. zijinensis_ reads to benchmark metagenomics analysis based on host-contamination.
  - Dilute to 0, 50, and 95% host-contamination (maybe 80% also)


*SCRIPTS*
- `abundance_info`: In an effort to understand the abundance metric, this script sums both the value and the abs(log(val)) for each file.
- `filter_thaliana`: Gets the reads that map to taxid `3702` in the `reads_mapping.tsv.gz` and removes them from the fastq file.
- `get_mock`: Downloads 10 datasets from the _A. thaliana_ rhizosphere challenge data from **CAMI2** challenge.
- `unpack_mock`: Unpacks the downloaded challenge data 


