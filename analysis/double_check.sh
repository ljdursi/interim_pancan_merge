#!/bin/bash

readonly MERGEDFILE=$1
readonly VARIANT=$2

readonly USAGE="Usage: $0 mergedfile snv_mnv|indel|indel_normed
  checks the merged VCF against the input callset"

if [ ! -f ${MERGEDFILE} ]  || [ -z ${VARIANT} ]
then
    >&2 echo "${USAGE}"
    >&2 echo "invocation: $0 $1 $2"
    exit
fi

if [ ${VARIANT} != "snv_mnv" ] && [ ${VARIANT} != "indel" ] && [ ${VARIANT} != "indel_normed" ]
then
    >&2 echo "${USAGE}"
    >&2 echo "invocation: $0 $1 $2"
    >&2 echo "invalid variant type: ${VARIANT}"
    exit
fi

function passvariants {
    awk '($1 !~ /^#/) && ($7 == "." || $7 == "PASS") { split($5, alts, ","); for (i in alts) {printf "%s_%s_%s_%s\n",$1,$2,$4,alts[i];}}' | sort
}

function allvariants {
    awk '($1 !~ /^#/) {split($5, alts, ","); for (i in alts) {printf "%s_%s_%s_%s\n",$1,$2,$4,alts[i];}}' | sort 
}

function compare {
    local mergedfilename=$1
    local callername=$2
    local callerfile=$3

    # I'd prefer not to have the uniq in here, but sanger indel caller will sometimes double-report an indel variant
    n_caller_in_merged=$( grep $callername $mergedfilename | wc -l )
    n_caller_in_individual=$( zcat $callerfile | passvariants | uniq | wc -l )
    n_in_common=$( join <( grep $callername $mergedfilename | allvariants ) \
                        <( zcat $callerfile | passvariants | uniq )  \
                     | wc -l )

    if [ $n_caller_in_merged -ne $n_caller_in_individual ] || [ $n_in_common -ne $n_caller_in_merged ]
    then
        >&2 echo "Warning: inconsistency between $callerfile ($callername) and $mergedfilename "
        >&2 echo "$n_caller_in_merged $n_caller_in_individual $n_in_common"
        exit 1
    fi
}

readonly SAMPLE=$( basename $MERGEDFILE | cut -f 1 -d . )
readonly filenamelist=$( ./analysis/filenames-from-sample.sh $SAMPLE $VARIANT )
IFS=' ' read -r -a filenames <<< "${filenamelist}"

compare $MERGEDFILE broad ${filenames[0]}
compare $MERGEDFILE dkfz ${filenames[1]}
compare $MERGEDFILE sanger ${filenames[2]}

echo "${SAMPLE}: PASSED"
