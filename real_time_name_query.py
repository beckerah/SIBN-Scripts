"""
This script reproduces the functionality of the Real-Time Name Query web app
found at https://www.globalgeno.me/gaps/live. It is intended to be run as 
a script at the command line. The correct usage is:

python real_time_name_query.py <names file path> <NCBI API key>
e.g. python real_time_name_query.py querynames.txt 01234567890ABCDEFGHIK

Names files should be text files formatted with one name on each line.

A NCBI API key can be obtained by creating an account on the NCBI website.

**Note:
Species-level queries are not available for GGBN, and any GGBN results obtained for
species names will not be accurate.
"""

import argparse, requests
import pandas as pd
    
def inCOL(tname):
    """Query Catalog of Life for taxonomic name and return boolean."""

    api_url = 'http://www.catalogueoflife.org/col/webservice'
    params = {
        'format': 'json',
        'response': 'terse',
        'name': tname
    }
    col = requests.get(api_url, params).json()
    numResults = int(col["total_number_of_results"])
    return numResults > 0
    
def inGGBN(tname):
    """Query GGBN for taxonomic name and return boolean."""

    api_url = 'http://data.ggbn.org/ggbn_portal/api/search'
    params = {
        'getCounts': True,
        'name': tname + "*"
    }
    ggbn = requests.get(api_url, params).json()
    numResults = int(ggbn["nbSamples"])
    return numResults > 0
    
def inGenBankAny(tname, ncbikey):
    """Query GenBank Eutils for any records, return count."""

    api_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    params = {
        'db': 'nuccore',
        'retmode': 'json',
        'rettype': 'count',
        'api_key': ncbikey,
        'term': (tname + '[Organism]')
    }
    result = requests.get(api_url, params).json()
    try:
        return int(result['esearchresult']['count'])
    except KeyError:
        return 0

def inGenBankBarcode(tname, ncbikey):
    """Query GenBank Eutils for barcode records and return counts 
    for coi, rbcl, matk, its"""

    api_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    record_counts = {'its': 0}
    
    # COI, rbcL, matK 
    for gene in ['coi', 'rbcl', 'matk']:
        params = {
            'db': 'nuccore',
            'retmode': 'json',
            'rettype': 'count',
            'api_key': ncbikey,
            'term': (
                tname + "[Organism] AND "+ 
                gene + "[gene] " +
                "AND 500:99999999[slen] " + 
                "AND src country[prop] " + 
                "AND src specimen voucher[prop] " +
                "AND src pcr primers[prop]"
            )
        }
        result = requests.get(api_url, params).json()
        try:
            record_counts[gene] = int(result['esearchresult']['count'])
        except KeyError:
            record_counts[gene] = 0
    
    # ITS
    its_params = {
        'db': 'nuccore',
        'retmode': 'json',
        'rettype': 'count',
        'api_key': ncbikey,
        'term': (
            tname + '[Organism] '
            'AND internal transcribed spacer[All Fields] ' + 
            'AND src country[prop] ' + 
            'AND src specimen voucher[prop] ' + 
            'AND src pcr primers[prop]'
        )
    }
    its_result = requests.get(api_url, its_params).json()
    try:
        record_counts['its'] = int(its_result['esearchresult']['count'])
    except KeyError:
        record_counts['its'] = 0

    return record_counts  

def processTaxNames(namelist, ncbikey):
    """Retrieve results from CoL, GGBN, and GenBank for each name in a supplied
    list of taxonomic names. Return results as a list of dictionaries."""

    allResults = []
    for t in namelist:
        result = {}
        result["name"] = t
        result["col"] = inCOL(t)
        result["ggbn"] = inGGBN(t)
        result["genbank_any"] = inGenBankAny(t, ncbikey)
        if result['genbank_any'] > 0:
            barcodes = inGenBankBarcode(t, ncbikey)
            result['genbank_barcodes_coi'] = barcodes['coi']
            result['genbank_barcodes_rbcl'] = barcodes['rbcl']
            result['genbank_barcodes_matk'] = barcodes['matk']
            result['genbank_barcodes_its'] = barcodes['its']
        else:
            result['genbank_barcodes_coi'] = 0
            result['genbank_barcodes_rbcl'] = 0
            result['genbank_barcodes_matk'] = 0
            result['genbank_barcodes_its'] = 0
        allResults.append(result)
    return allResults
    
if __name__ == '__main__':

    # Retrieve name file and source options supplied as argument
    parser = argparse.ArgumentParser()
    parser.add_argument('namesfile', help = '.txt file containing taxonomic names to query', type = str)
    parser.add_argument('ncbikey', help = 'API Key supplied by NCBI', type = str)
    args = parser.parse_args()

    # Read name file from supplied parameters
    with open(args.namesfile, 'r') as f:
        taxnames = [x.strip() for x in f.readlines()]

    # Retrieve results from name queries
    searchresults = processTaxNames(taxnames, args.ncbikey)

    # Output results as a CSV
    df = pd.DataFrame(searchresults)
    df = df[[
        'name', 
        'col', 
        'ggbn', 
        'genbank_any', 
        'genbank_barcodes_coi', 
        'genbank_barcodes_rbcl',
        'genbank_barcodes_matk',
        'genbank_barcodes_its',
    ]]
    df.to_csv('real_time_name_query_results.csv', index=False)