To run
------
./sample-sieve.sh <path/to/yaml>

Example YAML file
-----------------
gene_name: NOD2
maximum_maf: 0.1
consequences_without_severity_measures: 
  - frameshift_variant
  - stop_gained
consequences_with_severity_measures: 
  consequences: 
    - missense_variant
  filters_conjunction: or
  severity_measures: 
    CADD_PHRED: 10
    POLYPHEN: possibly damaging
    SIFT: deleterious
out_directory: /lustre/scratch115/teams/barrett/users/dr9/sample-sieve-output


YAML file issues
----------------
For issues with the YAML format please copy/paste your file into http://www.yamllint.com/


Notes
-----
source activate ./sample-sieve-env
source deactivate ./sample-sieve-env
