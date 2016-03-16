#!/usr/bin/env python
import argparse
import vcf  
import sys
import random

def passed(record):
    """ Did this variant pass all of its filters? """
    return record.FILTER is None or len(record.FILTER) == 0 or record.FILTER == "PASS"

def populate_dict(filename, caller, onlypasses=True):
    """
    Read each variant in filename, creating a dictionary
    with key (chrom, pos, ref, alt) and value VAF
    """
    caller_dict = {}
    if filename is None:
        return {}
    vcf_reader = vcf.Reader(open(filename, 'r'))
    for variant in vcf_reader:
        if not passed(variant) and onlypasses:
            continue
        usecallers = variant.INFO['Callers'] if 'Callers' in variant.INFO else [caller]
        for alt in variant.ALT:
            caller_dict[(variant.CHROM, variant.POS, variant.REF, str(alt))] = usecallers
    return caller_dict

def compare_against_merged(calldict, mergedict, ncalls=1000):
    ncalls = min(len(calldict), ncalls)
    variants = random.sample(calldict.keys(), ncalls)
    for variant in variants:
        if not variant in mergedict or calldict[variant][0] not in mergedict[variant]:
            print("Error: "+calldict[variant][0]+" not in merged dict for variant "+str(variant))
            print(mergedict[variant])
            return False
    return True

def compare_merged_against_callsets(mergedict, calldicts, ncalls=1000):
    ncalls = min(len(mergedict), ncalls)
    variants = random.sample(mergedict.keys(), ncalls)
    allcallers = calldicts.keys()

    for variant in variants:
        for caller in allcallers:
            if caller in mergedict[variant] and not variant in calldicts[caller]:
                print("Error: "+variant+" found in merged set for callers "+mergedict[variant]+
                       "not found in callset for "+caller+".")
                return False
            if not caller in mergedict[variant] and variant in calldicts[caller]:
                print("Error: "+variant+" not found in merged set for callers "+mergedict[variant]+
                       "but found in callset for "+caller+".")
                return False
    return True

def main():
    parser = argparse.ArgumentParser(description='Annotate merged vcf with VAF information where available')
    parser.add_argument('mergedvcf', type=str, help="Merged VCF file")
    parser.add_argument('-b', '--broad', type=str, help="Broad file")
    parser.add_argument('-d', '--dkfz', type=str, help="DKFZ file")
    parser.add_argument('-s', '--sanger', type=str, help="Sanger file")
    parser.add_argument('-m', '--muse', type=str, help="Muse file")
    parser.add_argument('-n', '--nvariants', type=int, default=1000, help="Number of variants to check per test")
    args = parser.parse_args()

    ncallers = 0
    for filename in [args.broad, args.dkfz, args.sanger]:
        if filename is not None:
            ncallers += 1

    dicts = {"broad":populate_dict(args.broad, "broad"),
             "dkfz":populate_dict(args.dkfz, "dkfz"),
             "sanger":populate_dict(args.sanger, "sanger"),
             "muse":populate_dict(args.muse, "muse")}

    mergeddict = populate_dict(args.mergedvcf, "merged", onlypasses=False)

    # for each call set, make sure all of a sample of its calls
    # are correctly tagged in the merged set
    success = True
    for caller in ["broad", "dkfz", "sanger"]:
        success = success and compare_against_merged(dicts[caller], mergeddict, args.nvariants)

    if not success:
        return 1

    # for merged call set, make sure a sample of its calls 
    # reflect the individual callsets
    success = success and compare_merged_against_callsets(mergeddict, dicts, args.nvariants)
    if not success:
        return 1

    print("PASSED: "+str(args.nvariants*(ncallers+1))+" checks.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
