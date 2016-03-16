#!/usr/bin/env python
import argparse
import vcf  
import sys
import numpy

def passed(record):
    """ Did this variant pass all of its filters? """
    return record.FILTER is None or len(record.FILTER) == 0 or record.FILTER == "PASS"

def round_three(num):
    scale=1000
    return int(num*scale)*1./scale

def get_vaf(record, broad=False, sanger=False, muse=False, dkfz=False, SNV=True):
    """ Returns the VAF corresponding to the record and the caller.
        Indels not yet implemented."""
    if SNV:
        if broad:
            return float(record.INFO['tumor_f'][0])
        if dkfz:
            return record.INFO['AF'][1]
        if sanger:
            return record.samples[1]['PM']
        if muse:
            return 1.0*record.samples[0]['AD'][1]/record.samples[0]['DP']
        return -100.
    else:  #indel
        return -100.

def populate_dict(filename, *args, **kwargs):
    """
    Read each variant in filename, creating a dictionary
    with key (chrom, pos, ref, alt) and value VAF
    """
    vaf_dict = {}
    if filename is None:
        return {}
    vcf_reader = vcf.Reader(open(filename, 'r'))
    for variant in vcf_reader:
        if not passed(variant):
            continue
        vaf = get_vaf(variant, *args, **kwargs)
        for alt in variant.ALT:
            vaf_dict[(variant.CHROM, variant.POS, variant.REF, str(alt))] = vaf
    return vaf_dict

def main():
    parser = argparse.ArgumentParser(description='Annotate merged vcf with VAF information where available')
    parser.add_argument('mergedvcf', type=argparse.FileType('r'), default=sys.stdin, help="Merged VCF file")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)")
    parser.add_argument('-b', '--broad', type=str, help="Broad file")
    parser.add_argument('-d', '--dkfz', type=str, help="DKFZ file")
    parser.add_argument('-s', '--sanger', type=str, help="Sanger file")
    parser.add_argument('-m', '--muse', type=str, help="Muse file")
    parser.add_argument('-i', '--indel', action='store_true',
                        help="Variant type == indel (default:snv_mnv)")
    args = parser.parse_args()
    snvs = not args.indel

    dicts = [populate_dict(args.broad, broad=True, SNV=snvs),
             populate_dict(args.dkfz, dkfz=True, SNV=snvs),
             populate_dict(args.sanger, sanger=True, SNV=snvs),
             populate_dict(args.muse, muse=True, SNV=snvs)]

    vcf_reader = vcf.Reader(args.mergedvcf)
    vcf_writer = vcf.Writer(args.output, vcf_reader)
    for variant in vcf_reader:
        key = variant.CHROM, variant.POS, variant.REF, str(variant.ALT[0])
        vafs = [vaf_dict[key]
                for vaf_dict in dicts
                if key in vaf_dict]
        roundvafs = [round_three(vaf) for vaf in vafs]
        if len(vafs) > 0:
            variant.INFO['VAFs'] = roundvafs
            variant.INFO['medianVAF'] = round_three(numpy.median(vafs))
        vcf_writer.write_record(variant)

    return 0

if __name__ == "__main__":
    sys.exit(main())
