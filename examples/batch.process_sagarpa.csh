## To run this script:
## - replace $PROJECT with the path to your clone of this GitHub repo
## - replace $SAGARPA_ROOT with the path to the top-level directory where SAGARPA files are stored
##
## Also, you must run this script from $PROJECT/tools/SAGARPA or add $PROJECT/tools/SAGARPA to your $PATH environment variable.

remove_nondelim_commas_in_csv.pl $SAGARPA_ROOT/orig/Agt_cierre_1980_2012.csv > $SAGARPA_ROOT/no_extra_commas/Agt_cierre_1980_2012.csv
process_sagarpa_files.pl $SAGARPA_ROOT/no_extra_commas/Agt_cierre_1980_2012.csv > $SAGARPA_ROOT/summ_by_mun/summ_by_mun.csv
