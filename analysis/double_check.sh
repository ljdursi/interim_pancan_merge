#!/bin/bash

readonly MERGEDFILE=$1

if [ ! -f ${SAMPLE} ] 
then
    >&2 echo "Usage: $0 mergedfile"
    >&2 echo "  checks the merged VCF against the input callsets"
    >&2 echo "invocation: $0 $1"
    exit
fi

function passvariants {
    awk '($1 !~ /^#/) && ($7 == "." || $7 == "PASS") { split($5, alts, ","); for (alt in alts) {printf "%s_%s_%s_%s\n",$1,$2,$4,alt;}}' | sort
}

function allvariants {
    awk '($1 !~ /^#/) {split($5, alts, ","); for (alt in alts) {printf "%s_%s_%s_%s\n",$1,$2,$4,alt;}}' | sort
}

function compare {
    local mergedfilename=$1
    local callername=$2
    local callerfile=$3

    n_caller_in_merged=$( grep $callername $mergedfilename | wc -l )
    n_caller_in_individual=$( zcat $callerfile | passvariants | wc -l )
    n_in_common=$( join <( grep $callername $mergedfilename | allvariants ) \
                        <( zcat $callerfile | passvariants )  \
                     | wc -l )

    if [ $n_caller_in_merged -ne $n_caller_in_individual ] || [ $n_in_common -ne $n_caller_in_merged ]
    then
        >&2 echo "Warning: inconsistency between $callerfile ($callername) and $mergedfilename "
        >&2 echo "$n_caller_in_merged $n_caller_in_individual $n_in_common"
        exit 1
    fi
}

readonly SAMPLE=$( basename $MERGEDFILE | cut -f 1 -d . )
readonly filenamelist=$( ./analysis/filenames-from-sample.sh $SAMPLE snv_mnv )
IFS=' ' read -r -a filenames <<< "${filenamelist}"

compare $MERGEDFILE broad ${filenames[0]}
compare $MERGEDFILE dkfz ${filenames[1]}
compare $MERGEDFILE sanger ${filenames[2]}

echo "${SAMPLE}: PASSED"
