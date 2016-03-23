# Interim Consensus Callset

A set of routines for generating and verifying interim merged callsets,
using a simple "2 of m" approach.

Relies on http://github.com/ljdursi/mergevcf

To use:
```bash
 $ ./gen_filelist      # map donor ids to files; depends on location of files on local FS
 $ ./normalize_indels  # bcftools norm the indel files, generate local copies
 $ ./merge             # merge the vcfs and label by number of callers
 $ ./annotate          # annotate with VAFs from callers
 $ ./test_merge        # trust, but verify
```
