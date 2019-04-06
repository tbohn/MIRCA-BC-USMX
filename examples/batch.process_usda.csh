## To run this script:
## - replace $PROJECT with the path to your clone of this GitHub repo
## - replace $USDA_ROOT with the path to the top-level directory where USDA files are stored
##
## Also, you must run this script from $PROJECT/tools/USDA or add $PROJECT/tools/USDA to your $PATH environment variable.

wrap_parse_usda_pdfs.pl $USDA_ROOT/orig table 1992 $PROJECT/data/county_codes.csv $USDA_ROOT/reformat >& log.wrap_parse_usda_pdfs.pl.txt
wrap_parse_usda_pdfs.pl $USDA_ROOT/orig table 2002 $PROJECT/data/county_codes.csv $USDA_ROOT/reformat >& log.wrap_parse_usda_pdfs.pl.txt
wrap_parse_usda_pdfs.pl $USDA_ROOT/orig table 2012 $PROJECT/data/county_codes.csv $USDA_ROOT/reformat >& log.wrap_parse_usda_pdfs.pl.txt
collate_usda_county_data.pl $USDA_ROOT/reformat 1992 $USDA_ROOT/summ >& log.collate_usda_county_data.pl.txt
collate_usda_county_data.pl $USDA_ROOT/reformat 2002 $USDA_ROOT/summ >& log.collate_usda_county_data.pl.txt
collate_usda_county_data.pl $USDA_ROOT/reformat 2012 $USDA_ROOT/summ >& log.collate_usda_county_data.pl.txt
