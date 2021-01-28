#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import requests
import numpy as np
import bs4
import re
from collections import Counter


# # UK Enterococcus Paper Metadata
# 
# ## Add sample and project accessions to UK metadata
# 
# Starting from table S2 in https://mbio.asm.org/content/9/6/e01780-18 which contains `run_accessions` for deposited raw sequencing data in ENA, use the ENA API to automatically add `study_accession` and `sample_accession`
# 

# In[2]:


# table S2 from https://mbio.asm.org/content/9/6/e01780-18
uk_metadata = pd.read_csv('datasets/inline-supplementary-material-3.csv', sep='\t', skiprows=1)
uk_metadata = uk_metadata.rename(columns={'Accession number': 'Run_accession', 
                                          'Isolate ID': 'Isolate_name',
                                          'ST': 'Sequence_Type'})

# break into 500 accession chunks to more efficiently query the API
api_responses = ""

run_accessions = uk_metadata.loc[uk_metadata['Run_accession'].str.startswith('ERR'), 'Run_accession'].values
for run_accession_chunk in [run_accessions[i:i + 500] for i in range(0, len(run_accessions), 500)]:
    run_accession_chunk = ",".join(run_accession_chunk)
    api_response = requests.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{run_accession_chunk}")    
    api_responses += api_response.text

# parse combined xml records into a soup for easier traversal 
api_response_soup = bs4.BeautifulSoup(api_responses, 'lxml')


def add_study_and_sample_metadata(row, api_response_soup):
    """
    Use the ENA API to get the study and sample accessions
    """
    run_accession = row['Run_accession']
        
    run_soup = api_response_soup.find(accession=run_accession)
    if run_soup:
        for xref in run_soup.find_all('xref_link'):
            if xref.db.text == 'ENA-STUDY':
                row['Study_accession'] = xref.id.text.strip()
            elif xref.db.text == 'ENA-SAMPLE':
                row['Sample_accession'] = xref.id.text.strip()
    else:
        row['Study_accession'] = np.nan
        row['Sample_accession'] = np.nan
    return row

uk_metadata = uk_metadata.apply(lambda x: add_study_and_sample_metadata(x, api_response_soup), axis=1)


# Then let's tidy up the UK metadata sheet by renaming some fields and adding some extra metadata to make our life easier merging in with the alberta data later. 

# In[3]:


# tidy up UK metadata for merging
uk_metadata = uk_metadata.rename(columns={'Origin': 'Origin',
                                          'BAPS group': 'BAPS_group',
                                          'Location': 'Location',
                                          'Ampicillin resistance': 'Ampicillin',
                                          'Vancomycin resistance': 'Vancomycin'})

uk_metadata.loc[uk_metadata['Removed in deduplication'] == 'Removed', 'Metadata_status'] = 'Removed for deduplication in original paper (10.1128/mBio.01780-18)'
uk_metadata = uk_metadata.drop('Removed in deduplication', axis=1)

# all in this paper are E. faecium
uk_metadata['Species'] = 'Enterococcus faecium'
uk_metadata['Country/Province'] = 'United Kingdom'

# drop reference strains
uk_metadata = uk_metadata.loc[uk_metadata['Origin'] != 'Reference strain']


# ## Add read data for UK reads to the UK metadata
# 
# Then parse the read data for all the reads (or accessions) given to me by HS and archived on `deivos.research.cs.dal.ca`

# In[4]:


read_data = pd.read_csv('datasets/Enterococcus_reads.tsv', sep='\t')
read_data = read_data.rename(columns={'Sample_name': "Strain_name"})
read_data = read_data.replace({'E_faecalis': 'Enterococcus faecalis',
                               'E_faecium': 'Enterococcus faecium'})

# split the datasets for now
uk_read_data = read_data[read_data['Strain_name'].str.startswith('ERR')]
uk_read_data['Run_accession'] = uk_read_data['Strain_name']
alberta_read_data = read_data[~read_data['Strain_name'].str.startswith('ERR')]


# And merge the UK metadata and UK read data using the `Run_accession` ensuring injective mapping.
# 
# Add explanation about why readsets are missing, download protocol initially used by HS:
# 1. get all accessions
# 2. try to download
# 3. keep the subset that successfully downloaded as a "random sample with similar habitat distribution as all Alberta E. faecium metadata"

# In[5]:


uk_metadata_and_reads = pd.merge(uk_metadata, 
                                 uk_read_data, 
                                 on='Run_accession', 
                                 validate='one_to_one', 
                                 suffixes=['_metadata', '_reads'],
                                 how='outer')
uk_metadata_and_reads.loc[uk_metadata_and_reads['Read_1'].isna(), 'Read_status'] = "Not successfully downloaded in initial attempt by HS"

# drop extraneous columns
uk_metadata_and_reads = uk_metadata_and_reads.drop(['Strain_name', 'Species_reads'], axis=1)

# make species information from metadata the official species for the UK data as this seems reliable for this dataset
uk_metadata_and_reads = uk_metadata_and_reads.rename(columns={'Species_metadata': 'Species'})

# re-order columns into natural groupings
uk_metadata_and_reads = uk_metadata_and_reads[
                      ['Study_accession', 'Sample_accession', 
                       'Run_accession', 'Isolate_name', 
                       'Metadata_status',
                       'Species', 'BAPS_group', 'Sequence_Type',
                       'Country/Province', 'Origin', 'Location',
                      'Ampicillin', 'Vancomycin', 'Read_status', 
                      'Read_1', "Read_2"]]


# # AAFC Enterococcus Paper Metadata
# 
# ## Link NCBI metadata to metadata published with paper
# 
# https://www.nature.com/articles/s41598-020-61002-5
# 
# Raw sequencing data was not deposited for this paper so only a `BioProject` accession (`PRJNA604849`) i.e., `study_accession` above was provided, this study on NCBI contains `sample_accessions` and the raw assembly files.
# 
# Therefore, download the metadata from this study and tidy up column names to make merging easier

# In[6]:


with open('datasets/PRJNA604849_full_biosample_list.xml') as fh:
    PRJNA604849_xml = bs4.BeautifulSoup(fh, 'lxml')

parsed_xml_data = {}
for biosample in PRJNA604849_xml.find_all('biosample'):
    biosample_data = {}
    biosample_data['Sample_name'] = biosample.find(db_label='Sample name').text
    biosample_data['Species'] = biosample.find('organismname').text
    biosample_data['Strain_name'] = biosample.find(attribute_name="strain").text
    biosample_data['Study_accession'] = biosample.find(target="bioproject")['label']
    parsed_xml_data[biosample['accession']] = biosample_data
    
alberta_ncbi = pd.DataFrame(parsed_xml_data).T.reset_index().rename(columns={'index': 'Sample_accession'})

# tidy up as we don't use sample_name for matching because strain_name was a closer match to the isolate name
# in the paper metadata
alberta_ncbi = alberta_ncbi.drop('Sample_name', axis=1)


# Now we need to parse and tidy up the metadata from https://www.nature.com/articles/s41598-020-61002-5  `41598_2020_61002_MOESM2_ESM.csv` as this seems to have some disconnections with the metadata via NCBI e.g., species assignments we will use the NCBI data when in doubt. Given that NCBI confirm species assignments this seems a prudent choice.

# In[7]:


alberta_metadata = pd.read_csv('datasets/41598_2020_61002_MOESM2_ESM.csv', sep='\t')

# drop the 2 inexplicably duplicated rows 
alberta_metadata = alberta_metadata.drop_duplicates()

# get rid of the identical rows apart from _ vs - 
alberta_metadata = alberta_metadata[alberta_metadata['ISOLATE'] != 'SWEntR-0393']
# identical row apart from incorrect species (when compared to NCBI assembly)
alberta_metadata = alberta_metadata[alberta_metadata['ISOLATE'] != 'SWEntR 0262']

# tidy column names to help with merging later
alberta_metadata = alberta_metadata.rename(columns={'ISOLATION SOURCE': 'Origin',
                                                    'SPECIFIC LOCATION': 'Location',
                                                    'VNCO': 'Vancomycin',
                                                    'AMPI': 'Ampicillin',
                                                    'TEIC': 'Teicoplanin',
                                                    'DOXY': 'Doxycycline',
                                                    'ERTH': 'Erythromycin',
                                                    'GENT': 'Gentamicin',
                                                    'LNZD': 'Linezolid',
                                                    'LVFL': 'Levofloxacin',
                                                    'QUIN': 'Quinolone',
                                                    'STEP': 'Streptomycin',
                                                    'NTRO': 'Nitrofurantoin',
                                                    'TGC': 'Tigecycline',
                                                    'SPECIATION': 'Species',
                                                    'SOURCE CODE': 'Source_code',
                                                    'ISOLATE': 'Isolate_name_paper'})

# get rid of useless columns and add extra column to help with merging UK and AB metadata
alberta_metadata = alberta_metadata.drop(['Unnamed: 18', "Resistance count", 'LOCATION'], axis=1)
alberta_metadata['Country/Province'] = 'Canada/Alberta'

# Remove trailing spaces that were left in the paper metadata
alberta_metadata['Species'] = alberta_metadata['Species'].str.strip()

# To ensure we merge the correct identifiers we are going to use the Source_code
# information AS WELL as the isolate_name, therefore let's create a new identifier
# out of the source code and isolate_name (and then remove spaces/underscores/hyphens etc)
def combine_source_and_isolate_name(row):
    """
    Try to combine isolate name with the source code field
    if it isn't already prefixed by the source code information
    """
    # Spaces, hyphens and dashes are a major source of disconnect so just remove them and try to map
    metadata_isolate_name = row['Isolate_name_paper'].replace('_', '').replace('-', '').replace(' ', '')
    if metadata_isolate_name.startswith(row['Source_code']):
        return metadata_isolate_name
    else:
        return row['Source_code'] +  metadata_isolate_name


alberta_metadata['Source_code_and_isolate_name'] = alberta_metadata.apply(                                                            combine_source_and_isolate_name, axis=1)


# Now we want to merge in the accession data from the bioproject.
# 
# The NCBI biosamples have both a distinct `Strain_name` and a distinct `Sample_name` whereas the paper metadata just has an `Isolate_name_paper` (original `ISOLATE`). 
# 
# `Isolate_name_paper` seems to correspond to either one of these NCBI identifiers with no apparent pattern. However, `Strain_name` does match the paper metadata `Isolate_name_paper` more often and seems to be generally closer.
# 
# Therefore, we are going to try and identify mappings between the NCBI `Strain_name` and the paper metadata `Isolate_name_paper`.  Then we are going to treat the NCBI `Strain_name` as the true `Isolate_name`.
# 
# Fortunately, it is usually fairly obvious what the correct mapping is as they usually only differ in hyphens/dashes/spacing/number of leading 0's. 
# The `Strain_name` in NCBI often contain part of what is the `Source_code` in the paper metadata (i.e., source), therefore I combined `Source_code` with `Isolate_name_paper` as `Source_code_and_isolate_name` and then searched for mappings using that.   
# 
# This means I can match `Strain_name` to `Source_code_and_isolate_name` by finding the `Source_code_and_isolate_name` that shares the longest suffix with each `Strain_name` in NCBI.
# 
# For extra security to prevent mis-assignments, I also added the following filter conditions: 
# 
# 1. All the numerical portions of `Strain_name` and `Source_code_and_isolate_name` had to match to be valid
# 2. The relationship between the two sets of names had to be one-to-one i.e., each `Strain_name` was assigned to one and only one `Source_code_and_isolate_name`.
# 
# Any remaining `Strain_name` in the NCBI data that wasn't assigned a `Source_code_and_isolate_name` can then be manually reviewed.

# In[8]:


def attempt_to_reconcile_name_sets(set1, set2):
    """
    Try and find the closest match between strings in set1 to set2
    First find the closest matching string by longest shared suffix
    Then as a safety move extract all digits and confirm they match between
    the biosample identifier and the closest paper identifier
    
    returns: dictionary mapping likely paper identifiers from set1 to set2
    """  
    closest_match_strings_by_suffix_length = []
    # for each name in set1 clean it up then and revese it (to make the suffix a prefix)
    for set1_name in set1:
        set1_name_clean = set1_name.replace('-', '').replace('_', '').lower()[::-1]
        
        # compare to each cleaned name in set2
        distances_between_set1_name_and_all_set2 = []
        longest_suffix = {'suffix_length': 0, 'set2_name': '', 'set1_name': set1_name}
        for set2_name in set2:
            set2_name_clean = set2_name.replace('-', '').replace('_', '').lower()[::-1]
            
            # and recover the string that has the longest shared suffix 
            suffix_length = 0
            for set1_name_char, set2_name_char in zip(set1_name_clean, set2_name_clean):
                if set1_name_char == set2_name_char:
                    suffix_length += 1
                else:
                    break
                    
            if suffix_length > longest_suffix['suffix_length']:
                longest_suffix['suffix_length'] = suffix_length
                longest_suffix['set2_name'] = set2_name

            closest_match_strings_by_suffix_length.append(longest_suffix)
    
    # then check that the best matchs ALSO share all the same numerical components
    closest_match_strings_by_suffix_size_and_numerical = []
    for closest_pairing in closest_match_strings_by_suffix_length:
        
        # handle mapping long date names (failing due to numbers in the date portion)
        if closest_pairing['set1_name'].startswith('ES-'):
            set1_numerical = "".join(re.findall(r'\d+', closest_pairing['set1_name'].split('-')[-1]))
        else:
            set1_numerical = "".join(re.findall(r'\d+', closest_pairing['set1_name']))

        if closest_pairing['set2_name'].startswith('ES-'):
            set2_numerical =  "".join(re.findall(r'\d+', closest_pairing['set2_name'].split('-')[-1]))
        else:
            set2_numerical =  "".join(re.findall(r'\d+', closest_pairing['set2_name']))
        
        if set1_numerical == set2_numerical:
            closest_match_strings_by_suffix_size_and_numerical.append((closest_pairing['set1_name'], 
                                                                       closest_pairing['set2_name']))
        else:
            # check if there is an issue with different numbers of leading 0s
            # as this is a common issue
            if set1_numerical.lstrip('0') == set2_numerical.lstrip('0'):
                 closest_match_strings_by_suffix_size_and_numerical.append((closest_pairing['set1_name'], 
                                                                            closest_pairing['set2_name']))
    
    # create a dictionary from the closest hits and extract any in set1 without a mapping to set2
    closest_match_strings_by_suffix_size_and_numerical = dict(closest_match_strings_by_suffix_size_and_numerical)
    
    # check for any non-one to one mappings in the dicionary i.e., different set1_name -> the same set2_name
    # delete and manually resolve 
    set2_name_counter = Counter(closest_match_strings_by_suffix_size_and_numerical.values())
    duplicates = []
    closest_match_strings_by_suffix_size_and_numerical_no_duplicates = {}
    for set1_name, set2_name in closest_match_strings_by_suffix_size_and_numerical.items():
        # i.e. drop all those who don't have one-to-one mapping
        if set2_name_counter[set2_name] == 1:
            closest_match_strings_by_suffix_size_and_numerical_no_duplicates[set1_name] = set2_name
    
    set1_names_without_matches = set(set1) - set(closest_match_strings_by_suffix_size_and_numerical_no_duplicates.keys())
    return closest_match_strings_by_suffix_size_and_numerical_no_duplicates, set1_names_without_matches

alberta_ncbi_to_metadata_match, alberta_ncbi_without_metadata_match = attempt_to_reconcile_name_sets(alberta_ncbi['Strain_name'].values, 
                                                                                                    alberta_metadata['Source_code_and_isolate_name'].values)


# Which NCBI strain_names didn't get a metadata match:

# In[9]:


alberta_ncbi_without_metadata_match


# Now we can try and manually fix these and identify cases of missing accessions or unresolvably ambiguous assignments (i.e., >1 perfectly valid appearing mapping from `Strain_name` to `Source_code_and_isolate_name` was possible).

# In[10]:


missing = {'HC_NS0026',  # ambiguous 
           'HC_NS0078', # ambiguous
           'HC_NS0854', # missing
           'HC_NS1042', # missing
           'HC_NS1090', # missing
           'HC_NS1104', # missing
           'HC_VRE0078', # missing
           'HC_SS0026', #missing
           'WW_0060M', # missing
           'WW_0089I'} # missing

manual_fixes = {'CB_0150': 'CBEntR0150',
                 'CB_0182': 'CBEntR0182',
                 'CB_0383': 'CBSWEntR0383',
                 'ES-C-ST002-07DEC15-0142B': 'WW0142B',
                 'FC_0142B': 'FC0142B',
                 'HC_NS0150': 'NSSNS0150',
                 'HC_NS0238': 'NSSNS0238',
                 'HC_NS0383': 'NSSNS0383',
                 'HC_NS210': 'NSSNS0210',
                 'HC_SS0002': 'SS0002',
                 'HC_SS0025': 'SS0025',
                 'SW_0002': 'NWSEnt0002',
                 'SW_0025': 'NWSEnt0025',
                 'SW_0182': 'NWSEnt0182',
                 'SW_0238': 'NWSEnt0238'}

alberta_ncbi_to_metadata_match.update(manual_fixes)


# Unfortunately, this still leaves us with 10 isolates that are in the deposited NCBI data but don't seem to have any supplied metadata

# In[11]:


unresolved_ncbi = alberta_ncbi[~alberta_ncbi['Strain_name'].isin(alberta_ncbi_to_metadata_match.keys())].sort_values('Strain_name')['Strain_name'].values
# add a status to the NCBI data indicating missing
alberta_ncbi.loc[alberta_ncbi['Strain_name'].isin(unresolved_ncbi), 'Metadata_status'] = 'No metadata in original paper (10.1038/s41598-020-61002-5)'


# Using our mapping from `Strain_name` to `Source_code_and_isolate_name` let's add a `Source_code_and_isolate_name` to the NCBI metadata sheet and merge using that (ensuring valid one-to-one merging).  
# 
# We are doing a `left` merge because we only want to keep the metadata from the paper for isolates that correspond to the NCBI bioproject.

# In[12]:


# translate the incorrect isolate_names in the paper metadata to the correct ones deposited in NCBI
alberta_ncbi['Source_code_and_isolate_name'] = alberta_ncbi['Strain_name'].apply(lambda x: alberta_ncbi_to_metadata_match[x] if x in alberta_ncbi_to_metadata_match else f"UNMATCHED: {x}")

# add the 10 unresolved identifiers to the paper metadata before merging
alberta_metadata = pd.concat([alberta_metadata, pd.DataFrame({'Source_code_and_isolate_name': unresolved_ncbi})])

alberta_merged = pd.merge(alberta_ncbi, alberta_metadata, 
                          validate='one_to_one', how='left', 
                          on='Source_code_and_isolate_name', suffixes=['_ncbi', '_paper'])


# Finally, we want to tidy up the alberta metadata to make our ultimate merging with the UK data easier

# In[13]:


# rename NCBI species information as official species information
# and change the NCBI `Strain_name` to `Isolate_name` to match better with other datasets (even though what
# NCBI calls things should remain the master)
alberta_merged = alberta_merged.rename(columns={'Species_ncbi': 'Species',
                                                'Strain_name': 'Isolate_name'})

# drop the superflous and likely erroneous species information from paper metadata sheet
alberta_merged = alberta_merged.drop('Species_paper', axis=1)

# reoroder the columns into logical groupings
alberta_merged = alberta_merged[['Study_accession', 'Sample_accession', 'Isolate_name', 'Isolate_name_paper', 'Metadata_status',
                                 'Species', 'Country/Province', "Origin", "Location", "Source_code", 'Ampicillin', 'Vancomycin', 'Teicoplanin',
                                 'Doxycycline', 'Erythromycin', 'Nitrofurantoin',
                                 'Gentamicin', 'Linezolid', 'Levofloxacin', 'Quinolone', 'Streptomycin',
                                 'Tigecycline', 'Source_code_and_isolate_name']]


# ## Add read data to Alberta metadata
# 
# Despite reconciling this at various stages the constant changes to metadata has meant this pairing has broken again, I'm applying the same approach as used for merging in the paper metadata to NCBI again. 
# 
# Specifically I need to link the "Strain_name" in `alberta_read_data` to the `Source_code_and_isolate_name` in the `alberta_merged` metadata dataframe
# 
# Let's apply the same approach as last time: suffix + numerical matching and guaranteeing injectivity

# In[14]:


metadata_to_read_data,     metadata_to_read_data_mismatches = attempt_to_reconcile_name_sets(alberta_merged['Source_code_and_isolate_name'].dropna().values,
                                                                      alberta_read_data['Strain_name'].values)


# Now let's fix the mismatches that we can manually

# In[15]:


metadata_to_read_data_mismatches


# In[16]:


missing_reads = {'BPC1383', # missing
                 'BPH112E1', # missing
                 'BPH532', # missing
                 'CB0104A', # missing
                 'CB0139J', # missing
                 'CBEnt0179', # missing
                 'FC0029A', # ambiguous
                 'FC0064E', # missing
                 'FC0142B', #ambiguous
                 'FC0194A', # missing
                 'FC0512A', # missing
                 'FC0555C', # missing
                 'FC0606I', # missing
                 'FC0616B', # missing
                 'FC0670I', # missing
                 'FC0834B', # missing 
                 'FC0845J', # missing
                 'NSSNS0176', # missing
                 'NSSNS0210', # missing
                 'NSSNS0554', # missing
                 'NSSNS0564', # missing
                 'NWSSWEnt0978', # missing
                 'NWSSWEnt1013', # missing
                 'NWSSWEnt1077', # missing
                 'NWSSWEnt1143', # missing
                 'NWSSWEntR0331', # missing   
                 'SS0018', # missing
                 'SS0030', # missing
                 'VRE0030', # missing 
                 'VRE0067', # missing
                 'VRE0069', # missing
                 'VRE0084', # missing
                 'VRE0090', # missing
                 'WW0039J', # ambiguous
                 'WW0050I', # missing
                 'WW0050M', # missing
                 'WW0081E', # missing
                 'WW0124I', # missing
                 'WW0137C', # missing
                 'WW0141J', # missing
                 'WW0141M', # missing
                 'WW0142A', # missing
                 'WW0142B', # missing
                 'WW0142C', # missing 
                 'WW0142E', # missing
                 'WW0142I', # missing
                 'WW0146B'} # missing

                 
manual_fix_reads = {'NSSNS0331': 'NS0331',
                    'VRE0018': 'VRE-0018',
                    'WW0029A': 'ES-M-ST001-21JUL14-0029A'}


metadata_to_read_data.update(manual_fix_reads)


# In[17]:


alberta_merged.loc[~alberta_merged['Source_code_and_isolate_name'].isin(metadata_to_read_data.keys()), 'Read_status'] = "Read data missing in dataset from HS"
alberta_merged['Isolate_name_reads'] = alberta_merged['Source_code_and_isolate_name'].apply(lambda x: metadata_to_read_data[x] if x in metadata_to_read_data else f"UNMATCHED: {x}")


# In[18]:


alberta_read_data = alberta_read_data.rename(columns={'Strain_name': 'Isolate_name_reads'})

alberta_data_all = pd.merge(alberta_merged, alberta_read_data, 
                          validate='one_to_one', how='left', 
                          on='Isolate_name_reads', suffixes=['_metadata', '_reads'])


# Final tidy up before merging of UK and AAFC data

# In[19]:


# finally get rid of our franken colunm with the source code and the species names from the reads
alberta_data_all = alberta_data_all.drop(['Source_code_and_isolate_name', 'Species_reads'], axis=1)

# rename Species_metadata to Species again
alberta_data_all = alberta_data_all.rename(columns={'Species_metadata': 'Species'})

# finally reorder the columns into logical groupings to better understand the data (even though we have to do this again
# after the big final merge)
alberta_data_all = alberta_data_all[['Study_accession', 'Sample_accession', 'Isolate_name',
                                     'Isolate_name_paper', 'Isolate_name_reads', 
                                     'Metadata_status', 'Species', 'Country/Province',
                                     'Origin', 'Location', 'Source_code', 'Ampicillin', 'Vancomycin',
                                     'Teicoplanin', 'Doxycycline', 'Erythromycin', 'Nitrofurantoin',
                                     'Gentamicin', 'Linezolid', 'Levofloxacin', 'Quinolone', 'Streptomycin',
                                     'Tigecycline', 'Read_status', 'Read_1', 'Read_2']]


# # Merging the two datasets
# 
# Combine the two datasets and tidy the names

# In[20]:


all_data = pd.concat([alberta_data_all, uk_metadata_and_reads])

all_data = all_data[['Study_accession', 'Sample_accession', 'Run_accession', 'Isolate_name',
                     'Isolate_name_paper', 'Isolate_name_reads', 
                     'Metadata_status', 'Species', 'BAPS_group','Sequence_Type', 'Country/Province', 
                     'Origin', 'Location', 'Source_code', 
                     'Ampicillin', 'Vancomycin', 'Teicoplanin', 'Doxycycline',
                     'Erythromycin', 'Nitrofurantoin', 'Gentamicin', 'Linezolid',
                     'Levofloxacin', 'Quinolone', 'Streptomycin', 'Tigecycline',
                     'Read_status', 'Read_1', 'Read_2']]

all_data.to_csv('all_combined_enterococcus_metadata.tsv', index=False, sep='\t')

