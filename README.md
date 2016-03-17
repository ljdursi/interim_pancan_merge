# Interim Consensus Callset

A set of routines for generating and verifying interim merged callsets,
using a simple "2 of m" approach.

Relies on http://github.com/ljdursi/mergevcf

To use:
```
 $ ./gen_filelist      # map donor ids to files; depends on location of files on local FS
 $ ./merge             # merge the vcfs and label by number of callers
 $ ./annotate_snvs     # annotate with VAFs from callers
 $ ./annotate_indels   # annotate with VAFs from callers
 $ ./test_merge        # trust, but verify
```
