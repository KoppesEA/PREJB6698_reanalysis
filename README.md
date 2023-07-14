# PREJB6698_reanalysis
## Summary:
This project is a (re)-analysis of a reduced-representation bisulfite sequencing (RRBS) dataset from DNMT1-TET-OFF mESCs generated by my PhD thesis mentor J.R. Chaillet (U.Pitt) and  initially analyzed by S. McGraw and J. Trasler (McGill University).
1.https://pubmed.ncbi.nlm.nih.gov/25578964/ 
2.https://www.ebi.ac.uk/ena/browser/view/PRJEB6698

## Data Mining Aquisition:
The McGraw et al. RRBS data is found on the European Nucleotide Archive (ENA) Record PRJEB6698.
1. Download read file TSV report directly from ENA
2. Use Awk one-liner to extract R1 and R2 FTP locations for 10 samples
```
cat filereport_read_run_PRJEB6698.tsv | awk -F"\t" '{print $7}' | awk -F";" -v OFS="\n" 'NR>1 {print $1, $2}' > PRJEB6698_ACC_List.txt
```
3. Using the wget bash script `PRJEB6698_ERR_wgetdownload.bash` download corresponding _R1 and _R2 fastq records, check out log to make sure no errors, should add checksum feature
5. Download read experiment XML
6. Use '.bash' or the below command to extract sample metadata
```
cat ena_PRJEB6698_read_experiment.xml | \
grep -o -E ">SAMEA\d{7}<|>ERR\d{6}<|refname\=\"[dN].*>" | \
tr -d '>|<' | tr '\n' ' ' | \
awk -v OFS="\t" ' {print $3, $2, $1, "\n", $6, $5, $4, "\n", $9, $8, $7, "\n", $12, $11, $10, "\n", $15, $14, $13, "\n", $18, $17, $16, "\n", $21, $20, $19, "\n", $24, $23, $22, "\n", $27, $26, $25, "\n", $30, $29, $28}' | \
sed 's/^\t//' | sed 's/refname=\"//' | sed 's/\"//' > PRJEB6698_metadata.tsv
```

## Bismark Genome Preperation:
1. Download the GRC38 Mouse Reference Genome to a Reference folder using the commands below:
```
wget -O CAST_EiJ_v1.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus_casteij/dna/Mus_musculus_casteij.CAST_EiJ_v1.dna.toplevel.fa.gz
```
2. Prepare bisulfite converted reference genome with CT and GA stranded C-->T deamination transition
```
module load bismark/0.20.0
bismark_genome_preparation ./
```

## Fastq Trimming and QC
1. Implement ``
