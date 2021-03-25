#!/bin/bash
set -euo 

# requires:
#  - text files with lists of accessions for all 3 papers/datasets in accession_lists
#  - downloaded hashes for the two ENA datasets
#  - sra tools installed and configured
#  - aspera connect installed and configured

python scripts/get_ena_checksums.py accession_lists/gouliouris_2018_uk_accessions.txt > checksums/gouliouris_2018_uk_ena
mkdir -p gouliouris_2018_uk
cd gouliouris_2018_uk
perl ../scripts/sra_download.pl --ascp ../accession_lists/gouliouris_2018_uk_accessions.txt
find . -name '*.gz'  -exec md5sum {} \; >> ../checksums/gouliouris_2018_uk_downloaded
cd ..
python scripts/compare_checksums.py checksums/gouliouris_2018_uk_ena checksum_check/gouliouris_2018_uk_downloaded

python scripts/get_ena_checksums.py accession_lists/raven_2016_uk_hospital_accessions.txt > checksums/raven_2016_uk_hospital_ena
mkdir -p raven_2016_uk_hospital
cd raven_2016_uk_hospital
perl ../scripts/sra_download.pl --ascp ../accession_lists/raven_2016_uk_hospital_accessions.txt
find . -name '*.gz'  -exec md5sum {} \; >> ../checksums/raven_2016_uk_hospital_downloaded
cd ..
python scripts/compare_checksums.py checksums/raven_2016_uk_hospital_ena checksums/raven_2016_uk_hospital_downloaded

mkdir -p zaheer_2020_alberta_sra
cd zaheer_2020_alberta_sra
prefetch --progress --verify yes --resume yes --output-directory . --option-file ../accession_lists/zaheer_2020_alberta_accessions.txt
vdb-validate ./*
mkdir -p ../zaheer_2020_alberta
for acc in *
    do 
        mkdir -p ../ab_dataset_reads/"$acc"
        fasterq-dump -p -O ../ab_dataset_reads/"$acc" --threads 4 "$acc"
        pigz -p 4 ../ab_dataset_reads/"$acc"/*.fastq
done
rm -rf zaheer_2020_alberta_sra
