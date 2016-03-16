#!/bin/bash

readonly SAMPLE=$1

if [ -z ${SAMPLE} ] 
then
    >&2 echo "Usage: $0 sample-id"
    >&2 echo "  combines the merged VCF file with annotations from all callers"
    >&2 echo "invocation: $0 $1 $2"
    exit
fi

readonly merged=results/snv_mnv/${SAMPLE}.merged.somatic.snv_mnv.vcf

if [ ! -f $merged ]
then
    >&2 echo "Invalid sample ${SAMPLE}"
    exit
fi

filenamelist=$( ./analysis/filenames-from-sample.sh $SAMPLE snv_mnv )
read -r -a filenames <<< "${filenamelist}"

module load python/2.7.2
module load gcc/4.8.1 openblas python-packages/2

mkdir -p annotated/snv_mnv

./analysis/annotate_vaf.py $merged \
    -b ${filenames[0]} \
    -d ${filenames[1]} \
    -s <( zcat ${filenames[2]} | sed '/=$/d' ) \
    -m ${filenames[3]} \
    | sed -e '/#CHROM/i\
##INFO=<ID=VAFs,Number=.,Type=Float,Description="VAFs identified by callers, in order of callers in Callers record">\
##INFO=<ID=medianVAF,Number=1,Type=Float,Description="median VAF identified by callers">\
##FILTER=<ID=LOWSUPPORT,Description="Not called by enough callers in ensemble">' \
    > annotated/snv_mnv/${SAMPLE}.merged.somatic.snv_mnv.vcf
