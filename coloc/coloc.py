#!/usr/bin/env python3
import requests
import re
from functools import lru_cache
from pysam import VariantFile
import numpy as np
import pandas as pd
from collections import Counter

# url_server = 'https://rest.ensembl.org'
url_server = 'https://grch37.rest.ensembl.org'

annotation_vcf_file_path = '/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/vep-annotation/all_studies.sampleids_cleaned_to_lowercase.bcf.gz.annot.vcf.gz'

# from pysam import TabixFile

# tbx = TabixFile('/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/vep-annotation/all_studies.sampleids_cleaned_to_lowercase.bcf.gz.annot.head.chrompos.gz')
# for header in tbx.header:
#     pass
# header = header.decode().split('\t')

# rows = []
# for row in tbx.fetch("1", 727655-1, 727655+1):
#     rows.append(row.split('\t'))

re_location = re.compile(r'(?P<chrom>\w+):(?P<start>\d+)-(?P<end>\d+)')

DTYPES = {
    'AFR_AF':  'float64',
    'AMR_AF':  'float64',
    'CADD_PHRED':  'float64',
    'CADD_RAW':  'float64',
    'CDS_position':  'float64',
    'DISTANCE':  'float64',
    'EAS_AF':  'float64',
    'EUR_AF':  'float64',
    'HGNC_ID':  'str',
    'INTRON':  'str',
    'PHENO':  'str',
    'Protein_position':  'float64',
    'SAS_AF':  'float64',
    'STRAND':  'float64',
    'cDNA_position':  'float64',
    'gnomAD_AF':  'float64',
    'gnomAD_AFR_AF':  'float64',
    'gnomAD_AMR_AF':  'float64',
    'gnomAD_ASJ_AF':  'float64',
    'gnomAD_EAS_AF':  'float64',
    'gnomAD_FIN_AF':  'float64',
    'gnomAD_NFE_AF':  'float64',
    'gnomAD_OTH_AF':  'float64',
    'gnomAD_SAS_AF':  'float64',
    'TYPED':  'bool',
    'AN':  'float64',
    'RefPanelAF':  'float64',
    'AC':  'float64',
    'chrom':  'str',
    'start':  'float64',
    'stop':  'float64',
}

# from http://genetics.bwh.harvard.edu/pph/pph_help_text.html#PredictionBasis
POLYPHEN_PREDICTIONS_LO_HI = {
    'benign': 0,
    'unknown': 1,
    'possibly damaging': 2,
    'probably damaging': 3,
}

SIFT_PREDICTIONS_LO_HI = {
    'tolerated': 0,
    'deleterious': 1,
}

def extract_chrom_start_end(location):
    m = re_location.match(location)
    return m.groupdict()


def url_get(url):
#     print(url)
    r = requests.get(url, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        raise Exception()
    return r.json()


def get_variants_in_chrom_start_end(chrom, start, end, minimum_maf=None, **kwargs):
#     print(chrom, start, end, minimum_maf)
    url_endpoint = '/overlap/region/human/'
    url_region = f'{chrom}:{start}-{end}'
    url_options = '?feature=variation;content-type=application/json'
    url = ''.join([url_server, url_endpoint, url_region, url_options])
    variants = url_get(url)

#     if minimum_maf:
#         variants = [x for x in variants if variant_has_greater_equal_MAF(x['id'], minimum_maf)]

    return variants


def parse_children(child):
    results = []
    if 'children' in child:
        for children in child['children']:
            results += parse_children(children)
    if 'sequence' in child:
        c = child['sequence']
        if 'name' in c:
            name = c['name']
            location = c['location']
            info = extract_chrom_start_end(location)
            assert len(c['id']) ==  1
            info['id'] = c['id'][0]['accession']
            info['name'] = c['name']
            results.append(info)
    return results


def get_region_info_for_gene(gene):
    """
    there is a chance this will return > 1 regions for 1 gene - process both in the pipeline
    """
    url_endpoint = '/genetree/member/symbol/homo_sapiens/'
    url_options = '?prune_species=homo_sapiens;content-type=application/json'
    url = ''.join([url_server, url_endpoint, gene, url_options])
    data = url_get(url)
    
    """
    try:
    
    except e:
        filter on name to include gene name ie "IL10RB-001"
    """
#     return data
    tree = data['tree']
    if 'sequence' in tree:
        location = tree['sequence']['location']
        info = extract_chrom_start_end(location)
        info['id'] = tree['id']['accession']
    else:
        info = parse_children(tree)
        info = [x for x in info if gene in x['name']]
        assert len(info) == 1
        info = info[0]
    return info


def get_eqtl_info_for_stable_id(stable_id, p_value=2):
    # p_value removes rubbish
    url_endpoint = '/eqtl/id/homo_sapiens/'
    url_options = '?statistic=p-value;content-type=application/json'
    url = ''.join([url_server, url_endpoint, stable_id, url_options])
    data = url_get(url)
    return data


def variant_has_greater_equal_MAF(rsid, minimum_maf=0.05, population='EUR'):
    variation_frequency = get_variation_frequency(rsid, population)
    minor_allele = variation_frequency['minor_allele']
    if minor_allele is None:
        print(variation_frequency)
        raise Exception('abc')
    print('minor_allele', minor_allele)
    passed = [x for x in variation_frequency['populations'] if 
                x['population'] == f'1000GENOMES:{population}' and
                x['allele'] == minor_allele and 
                float(x['frequency']) >= minimum_maf]
    if len(passed) == 1:
        return True
    if len(passed) > 1:
        raise Exception('More than one minor allele found: ' + str(passed))
    return False


def get_variation_frequency(region):
    pass
#     gnomAD_NFE_AF


def sort_eqtl_info_by_minus_log10_p_value(eqtl_info):
    return sorted(eqtl_info, key=lambda x:float(x['minus_log10_p_value']), reverse=True)


def get_genotype_accumulation_from_variant_file(variants, variant_file_path):
    """
    Discard samples
    
    eventually use a set of files
    
        seq_region_name
        start
        end
    """
    variant_file = VariantFile(variant_file_path)

    n_samples = get_number_of_samples(variant_file_path)
    n_variants = len(variants)
    GENOTYPES = np.zeros((n_samples, 2*n_variants), dtype=np.float16)

    for variant_index, variant in enumerate(variants):
        print('variant_index', variant_index)
        chrom = variant['seq_region_name']
        start = variant['start'] - 1
        end = variant['end']

        records = list(variant_file.fetch(chrom, start, end))

        n_records = 0
        for record in records:
            n_records += 1
            if n_records > 1:
                raise Exception(f'More than one record found for {chrom}:{start}-{end}')
            print(chrom, start, end)
            genotypes = [s['GT'] for s in record.samples.values()]
            genotypes = np.array(genotypes, dtype=np.float16)
            GENOTYPES += genotypes

    #         print('len(np.where(genotypes > 0)[0])', len(np.where(genotypes > 0)[0]))
    #         print('len(np.where(GENOTYPES > 0)[0])', len(np.where(GENOTYPES > 0)[0]))
    
    return GENOTYPES, variant_file.header.samples


def get_genotypes_from_variant_file(variants, variant_file_path):

    variant_file = VariantFile(variant_file_path)
    
    n_samples = get_number_of_samples(variant_file_path)
    n_variants = len(variants)
    GENOTYPES = np.zeros((n_samples, n_variants, 2), dtype=np.float16)

    for variant_index, variant in enumerate(variants):
        chrom = variant['seq_region_name']
        start = variant['start'] - 1
        end = variant['end']

        records = list(variant_file.fetch(chrom, start, end))

        n_records = 0
        for record in records:
            n_records += 1
            if n_records > 1:
                raise Exception(f'More than one record found for {chrom}:{start}-{end}')
            genotypes = [s['GT'] for s in record.samples.values()]
            genotypes = np.array(genotypes, dtype=np.float16)
            GENOTYPES[:,variant_index] = genotypes
            
    return GENOTYPES, variant_file.header.samples


def test_classify_sample_genotypes():
    hom = np.array([
           [1, 1],
           [0, 1],
           [1, 0],
           [1, 0],
           [0, 1],
           [1, 0],
           [1, 0],
           [np.nan, 1],
           [np.nan, np.nan],
           [0, 0]])

    chet = np.array([
           [0, 1],
           [0, 1],
           [1, 0],
           [1, 0],
           [0, 1],
           [1, 0],
           [1, 0],
           [np.nan, 1],
           [np.nan, np.nan],
           [0, 0]])

    het = np.array([
           [0, 1],
           [0, 1],
           [0, 0],
           [0, 0],
           [0, 1],
           [0, 0],
           [0, 0],
           [np.nan, 0],
           [np.nan, np.nan],
           [0, 0]])

    ref = np.array([
           [0, 0],
           [0, 0],
           [0, 0],
           [0, 0],
           [0, 0],
           [0, 0],
           [0, 0],
           [np.nan, 0],
           [np.nan, np.nan],
           [0, 0]])

    classify_sample_genotypes(hom)
    classify_sample_genotypes(chet)
    classify_sample_genotypes(het)
    classify_sample_genotypes(ref)


def classify_sample_genotypes(genotypes):
    homozygous = np.array([1.,1.])
    found_homozygous = (homozygous == genotypes).all(axis=1)
    n_homozygous = sum(found_homozygous)
    if n_homozygous > 0:
        return 'homozygous', n_homozygous

    genotypes_nans_to_zero = genotypes.copy()
    nans = np.isnan(genotypes_nans_to_zero)
    genotypes_nans_to_zero[nans] = 0
    
    n_alt = genotypes_nans_to_zero.sum(axis=0)

    found_compound_heterozygous = (n_alt > 0).all()
    if found_compound_heterozygous:
        return 'compound_heterozygous', n_alt.sum()

    found_heterozygous = (n_alt > 0).any()
    if found_heterozygous:
        return 'heterozygous', max(n_alt)
    
    reference = np.array([0.,0.])
    found_reference = (reference == genotypes_nans_to_zero).all(axis=1)
    if found_reference.all():
        n_nan_rows = nans.any(axis=1).sum()
        n_reference = found_reference.sum() - n_nan_rows
        return 'reference', n_reference
    
    
def classify_all_genotypes(GENOTYPES, samples):
    L = []
    for genotypes, sample in zip(GENOTYPES, samples):
        genotype_class, n = classify_sample_genotypes(genotypes)
        L.append({
            'genotype_class': genotype_class,
            'n': n,
            'sanger_sample_id': sample
        })
    df = pd.DataFrame(L)
    return df
    
#     c = Counter([x[0] for x in classes])
#     print(c)
#     return classes

# # @lru_cache(maxsize=None)
# def get_compound_hets_from_variant_file(variants, variant_file_path):
#     GENOTYPES, samples = get_genotype_accumulation_from_variant_file(variants, variant_file_path)    
#     i = (GENOTYPES[:,0] > 0) & (GENOTYPES[:,1] > 0)
#     # return number of hits
#     return GENOTYPES, samples, i

# def get_hom_from_variant_file(variants, variant_file_path):
#     """
#     For any variant, if 1 is hom (1|1) return True
#     """
    
# def get_het_from_variant_file(variants, variant_file_path):
#     """
#     For any variant, if 1 is het (0|1) or (1|0) return True
#     (either column is >= 1)
#     """


# def get_compound_homs_from_variant_file(variants, variant_file_path):
#     GENOTYPES, samples = get_genotype_accumulation_from_variant_file(variants, variant_file_path)
#     i = (GENOTYPES[:,0] > 0) & (GENOTYPES[:,1] == 0)
#     i |= (GENOTYPES[:,0] == 0) & (GENOTYPES[:,1] > 0)
#     return GENOTYPES, samples, i


def get_number_of_samples(variant_file_path):
    query_variant_file = VariantFile(variant_file_path)
    return len(query_variant_file.header.samples)


def get_most_consequence_atleast_as_bad():
    pass
    
    
def get_consequence_types():
    url_endpoint = '/info/variation/consequence_types?'
    url = ''.join([url_server, url_endpoint])
    consequence_types = url_get(url)
    consequence_types = [x['SO_term'] for x in consequence_types]
    return consequence_types


def are_input_consequences_valid(input_consequence_types):
    input_consequence_types = set(input_consequence_types)
    valid_consequence_types = set(get_consequence_types())
    
    if not valid_consequence_types.issuperset(input_consequence_types):
        invalid = input_consequence_types - valid_consequence_types
        invalid = ''.join(sorted([f'\t{x}\n' for x in invalid]))
        valid = ''.join(sorted([f'\t{x}\n' for x in valid_consequence_types]))
        print('The following provided consequences types are not valid:')
        print(invalid)
        print('Please select from the following:')
        print(valid)
        return False
    
    return True

def filter_annotation_data_maf(annotation_data, minimum_maf=None, maximum_maf=None):
    """
    The following are available form the annotation data:
        AFR_AF
        AMR_AF
        EAS_AF
        EUR_AF
        SAS_AF
        gnomAD_AF
        gnomAD_AFR_AF
        gnomAD_AMR_AF
        gnomAD_ASJ_AF
        gnomAD_EAS_AF
        gnomAD_FIN_AF
        gnomAD_NFE_AF
        gnomAD_OTH_AF
        gnomAD_SAS_AF
        RefPanelAF

    
        1. gnomAD_AF
        if not present use
        2. RefPanelAF
        3. return the number not found
        UPDATE: use only RefPanelAF
    """
#     af_column = 'RefPanelAF'
#     annotation_data = annotation_data.copy()
#     annotation_data['maf_used_for_filtering'] = annotation_data['gnomAD_AF']
#     null_index =  annotation_data['maf_used_for_filtering'].isnull()
#     annotation_data.loc[null_index, 'maf_used_for_filtering'] = annotation_data.loc[null_index, 'RefPanelAF']
#     not_null_index =  annotation_data['maf_used_for_filtering'].notnull()
#     assert not_null_index.all()
#     index = annotation_data['maf_used_for_filtering'] > minimum_maf
#     return annotation_data[index]

    if minimum_maf is None and maximum_maf is None:
          raise Exception('One of minimum_maf or maximum_maf must be set.')
    
    maf_column = 'RefPanelAF'
    maf = annotation_data[maf_column]
    assert maf.notnull().all()
    index = pd.Series([True]*len(annotation_data))
    if minimum_maf is not None:
        index &= maf >= minimum_maf
    if maximum_maf is not None:
        index &= maf <= maximum_maf
    return annotation_data[index]


# @lru_cache(maxsize=None)
def get_annotation_data(chrom, start, end, **kwargs):
    """
    input
    consequences = ['5_prime_UTR', 'missense']
    
    return
    consequenceA poly fiff cad
    consequenceB poly fiff cad
    """
    
    chrom = str(chrom)
    start = int(start)
    end = int(end)
    variant_file = VariantFile(annotation_vcf_file_path)
    t = variant_file.header.info['CSQ']
    description = t.description.replace('Consequence annotations from Ensembl VEP. Format: ', '')
    CSQ_fields = description.split('|')
    f = variant_file.fetch(chrom, start, end)

    dfs = []
    for site in f:
        CSQs = []
        for CSQ in site.info['CSQ']:
            CSQ_split = CSQ.split('|')
            CSQs.append({k:v for k,v in zip(CSQ_fields, CSQ_split)})
        df = pd.DataFrame(CSQs)

        df['TYPED'] = site.info.get('TYPED')
        df['AN'] = site.info['AN']

        RefPanelAF = site.info['RefPanelAF']
        assert len(RefPanelAF) ==1
        df['RefPanelAF'] = RefPanelAF[0]
        AC = site.info['AC']
        assert len(AC) == 1
        df['AC'] = AC[0]
        df['chrom'] = site.chrom
        df['start'] = site.start
        df['stop'] = site.stop
        df['id'] = site.id
        dfs.append(df)

    df = pd.concat(dfs)
    df = df.reset_index(drop=True)
    df = df.replace('', pd.np.nan)
    df = df.astype(DTYPES)

    return df


def filter_annotation_data_on_consequences(df, consequences_without_severity_measures, 
                                           consequences_with_severity_measures):
    f = lambda row: filter_annotation_data_row_on_consequences(
            row,
            consequences_without_severity_measures, 
            consequences_with_severity_measures
    )
    index = df.apply(f, axis=1)
    return df[index]


def filter_annotation_data_row_on_consequences(row, consequences_without_severity_measures,
                                               consequences_with_severity_measures):
    """
    TODO
    ----
    Currently we allow multiple Consequences ie
        synonymous variant&NMD transcript variant
    But it is not clear how to treat them at the moment.
    """
    observed_consequences = row['Consequence']
    if '&' not in observed_consequences:
        observed_consequences = set([observed_consequences])
    else:
        observed_consequences = set('&'.split(observed_consequences))
    
    # consequences_without_severity_measures
    wanted_consequences = set(consequences_without_severity_measures)
    if wanted_consequences & observed_consequences:
        return True
    
    # consequences_with_severity_measures
    wanted_consequences = set(consequences_with_severity_measures['consequences'])
    if wanted_consequences & observed_consequences:
        SIFT = row['SIFT']
        SIFT_encoded = SIFT_PREDICTIONS_LO_HI.get(SIFT, pd.np.nan)
        CADD_PHRED = row['CADD_PHRED']
        PolyPhen = row['PolyPhen']
        PolyPhen_encoded = POLYPHEN_PREDICTIONS_LO_HI.get(PolyPhen, pd.np.nan)
        
        _ = consequences_with_severity_measures['severity_measures']['CADD_PHRED']
        CADD_PHRED_test = CADD_PHRED >= _
        
        _ = consequences_with_severity_measures['severity_measures_encoded']['POLYPHEN']
        PolyPhen_test = PolyPhen_encoded >= _
        
        _ = consequences_with_severity_measures['severity_measures_encoded']['SIFT']
        SIFT_test = SIFT_encoded >= _

        tests = [CADD_PHRED_test, PolyPhen_test, SIFT_test]

        conjnction = consequences_with_severity_measures['filters_conjunction']
        if conjnction ==  'or':
            outcome = any(tests)
        elif conjnction ==  'and':
            outcome = all(tests)
        else:
            raise Exception(f'Unrecognized conjuction: {conjnction}')
        if outcome:
            return True
    
    return False


def encode_severity_measures(severity_measures):
    severity_measures['severity_measures_encoded'] = {
        'SIFT': SIFT_PREDICTIONS_LO_HI[severity_measures['severity_measures']['SIFT']],
        'POLYPHEN': POLYPHEN_PREDICTIONS_LO_HI[severity_measures['severity_measures']['POLYPHEN']],
    }
    return severity_measures


def get_info_score_from_variant_file(chrom, start, end, variant_file, **kwargs):
    """
    return the INFO score to be appended to the annotation data dataframe
    """
    records = list(variant_file.fetch(chrom, start, end))
    n_records = 0
    for record in records:
        n_records += 1
        if n_records > 1:
            raise Exception(f'More than one record found for {chrom}:{start}-{end}')
        info_score = record.info['INFO']
    return info_score


def append_info_scores_from_variant_file(df, variant_file):
    columns = ['id', 'chrom', 'start', 'stop']
    variants = df[columns].drop_duplicates('id').copy()
    f = lambda x: get_info_score_from_variant_file(x['chrom'], x['start'], x['stop'], variant_file)
    variants['INFO_score'] = variants.apply(f, axis=1)
    return pd.merge(df, variants, on=columns, how='left')


def filter_annotation_data_on_INFO_score(df, minimimum_INFO_score):
    return df[df['INFO_score'] >= minimimum_INFO_score]


def get_variants_from_df(df):
    variants = df.drop_duplicates(['chrom', 'start', 'stop'])
    L = []
    for _, variant in variants.iterrows():
        L.append({
            'seq_region_name': variant['chrom'],
            'start': variant['start'],
            'end': variant['stop']
        })
    return L


def main():
    gene_name = 'NOD2' #'IL10RA' #or 'IL10R'
    consequences_without_severity_measures = [
        'frameshift_variant',
        'stop_gained',
        'start_lost',
        'transcript_ablation',
        'feature_truncation',
        'incomplete_terminal_codon_variant',
        'inframe_insertion',
        'regulatory_region_variant',
        '3_prime_UTR_variant',
        'stop_gained',
        'stop_retained_variant',
        'regulatory_region_amplification',
        'protein_altering_variant',
        'transcript_amplification',
        'transcript_ablation',
        'splice_acceptor_variant',
        'coding_sequence_variant',
        'inframe_deletion',
        'feature_elongation',
        'NMD_transcript_variant',
        'regulatory_region_ablation',
        'splice_region_variant',
        'stop_lost',
        'splice_donor_variant',
        'start_lost',
        'TFBS_ablation',
        'TFBS_amplification',
        'frameshift_variant',
    ]

    consequences_with_severity_measures = {
        'consequences': {
            'missense_variant',
        },
        'filters_conjunction': 'or', # in the case of 'and' NA are counted as False and should not be returned, or
        'severity_measures': {
                'SIFT': 'deleterious',
                'POLYPHEN': 'possibly damaging',
                'CADD_PHRED': 10,
        }
    }
    maximum_maf = 0.1
    
    variant_file_path = '/lustre/scratch119/realdata/mdt2/teams/anderson/users/dr9/new_imputation_dan/all_studies.sampleids_cleaned_to_lowercase.bcf.gz' 

    consequences_with_severity_measures = encode_severity_measures(consequences_with_severity_measures)
    consequences = consequences_with_severity_measures['consequences']
    assert are_input_consequences_valid(consequences)

    region_info = get_region_info_for_gene(gene_name)
    print(region_info)
    annotation_data = get_annotation_data(**region_info)

    maf_filtered_annotation_data = filter_annotation_data_maf(annotation_data, maximum_maf=maximum_maf)
    print(len(maf_filtered_annotation_data))

    print(min(annotation_data['RefPanelAF']))
    print(max(annotation_data['RefPanelAF']))
    print(min(maf_filtered_annotation_data['RefPanelAF']))
    print(max(maf_filtered_annotation_data['RefPanelAF']))

    consequence_maf_filtered_annotation_data = filter_annotation_data_on_consequences(
        maf_filtered_annotation_data,
        consequences_without_severity_measures,
        consequences_with_severity_measures
    )

    variants = get_variants_from_df(consequence_maf_filtered_annotation_data)
    GENOTYPES, samples = get_genotypes_from_variant_file(variants, variant_file_path)
    classes = classify_all_genotypes(GENOTYPES, samples)

    consequence_maf_filtered_annotation_data.to_csv('consequences.tsv', index=False, sep='\t')
    classes.to_csv('classes.tsv', index=False, sep='\t')

if __name__ == '__main__':
    main()
