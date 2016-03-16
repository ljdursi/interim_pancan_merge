#!/bin/bash

readonly BASEDIR=/oicr/data/pancanxfer
readonly SANGER_DIR=${BASEDIR}/Sanger_workflow_variants/batch07/
readonly DKFZ_DIR=${BASEDIR}/DKFZ_EMBL_workflow_variants/
readonly MUSE_DIR=${BASEDIR}/Broad_workflow_variants/batch02/muse_vcf/
readonly BROAD_DIR=${BASEDIR}/Broad_workflow_variants/batch02/broad_vcf/

readonly SAMPLE=$1
readonly VARIANT=$2

if [ -z ${SAMPLE} ] || [ -z ${VARIANT} ]
then
    >&2 echo "Usage: $0 sample-id snv_mnv|indel"
    >&2 echo "  returns the sanger, dkfz, muse, and broad filenames"
    >&2 echo "  corresponding to the variant calls on that sample"
    >&2 echo "invocation: $0 $1 $2"
    >&2 exit
fi

if [ ${VARIANT} != "snv_mnv" ] && [ ${VARIANT} != "indel" ]
then
    >&2 echo "Error: variant $VARIANT unrecognized"
    >&2 echo "       must be one of: snv_mnv indel"
    >&2 echo "invocation: $0 $1 $2"
    exit
fi

sangerfile=$( find $SANGER_DIR -type f -name "${SAMPLE}*somatic.${VARIANT}.vcf.gz" )
dkfzfile=$( find $DKFZ_DIR -type f -name "${SAMPLE}*somatic.${VARIANT}.vcf.gz" )
broadfile=$( find $BROAD_DIR -type f -name "${SAMPLE}*somatic.${VARIANT}.vcf.gz" )

if [ ${VARIANT} == "snv_mnv" ]
then
    musefile=$( find $MUSE_DIR -type f -name "${SAMPLE}*somatic.${VARIANT}.vcf.gz" )
fi

echo $broadfile $dkfzfile $sangerfile $musefile
