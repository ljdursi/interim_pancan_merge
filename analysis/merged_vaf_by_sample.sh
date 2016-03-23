#!/bin/bash

readonly SAMPLE=$1
readonly VARIANT=$2

readonly USAGE="Usage: $0 sample-id snv_mnv|indel|indel_normed"

if [ -z ${SAMPLE} ] || [ -z ${VARIANT} ]
then
    >&2 echo $USAGE
    >&2 echo "  combines the merged VCF file with annotations from all callers"
    >&2 echo "invocation: $0 $1 $2"
    exit
fi

if [ ${VARIANT} != "snv_mnv" ] && [ ${VARIANT} != "indel" ] && [ ${VARIANT} != "indel_normed" ]
then
    >&2 echo $USAGE
    >&2 echo " Invalid VARIANT type ${VARIANT}"
    >&2 echo "invocation: $0 $1 $2"
    exit
fi

readonly merged=results/${VARIANT}/${SAMPLE}.merged.somatic.${VARIANT}.vcf

if [ ! -f $merged ]
then
    >&2 echo $USAGE
    >&2 echo "Invalid sample ${SAMPLE}"
    exit
fi

filenamelist=$( ./analysis/filenames-from-sample.sh $SAMPLE ${VARIANT} )
read -r -a filenames <<< "${filenamelist}"

module load python/2.7.2
module load gcc/4.8.1 openblas python-packages/2

mkdir -p annotated/${VARIANT}

function addheaders {
    sed -e '/#CHROM/i\
##INFO=<ID=VAFs,Number=.,Type=Float,Description="Variant Allele Fractions identified by callers, in order of callers in Callers record">\
##INFO=<ID=medianVAF,Number=1,Type=Float,Description="median VAF identified by callers">\
##FILTER=<ID=LOWSUPPORT,Description="Not called by enough callers in ensemble">' 
}

if [ ${VARIANT} == "snv_mnv"]
then
    ./analysis/annotate_vaf.py $merged \
        -b ${filenames[0]} \
        -d ${filenames[1]} \
        -s <( zcat ${filenames[2]} | sed '/=$/d' ) \
        -m ${filenames[3]} \
        | addheaders \
            > annotated/${VARIANT}/${SAMPLE}.merged.somatic.${VARIANT}.vcf
else
    ./analysis/annotate_vaf.py $merged \
        -b ${filenames[0]} \
        -d ${filenames[1]} \
        -s <( zcat ${filenames[2]} | sed '/=$/d' ) \
        -i \
        | addheaders \
            > annotated/${VARIANT}/${SAMPLE}.merged.somatic.${VARIANT}.vcf
fi
