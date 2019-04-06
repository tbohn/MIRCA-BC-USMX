#!/usr/bin/perl

$indir = shift;
$prefix = shift;
$year = shift;
$county_table = shift;
$outdir = shift;

$tmpdir = $indir;
$tmpdir =~ s/orig/asc/g;
`mkdir -p $tmpdir`;
`mkdir -p $outdir`;

opendir(DIR,$indir) or die "$0: ERROR: cannot open dir $indir for reading\n";
@files = grep /$year/, grep /^$prefix/, readdir(DIR);
closedir(DIR);

foreach $file (sort(@files)) {

  $cmd = "pdftotext -layout $indir/$file asc/$file.asc";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

  if ($file =~ /table.(\d\d)\./) {
    $statenum = $1;
  }
  $outfile = $file;
  $outfile =~ s/pdf/csv/;
  $cmd = "parse_usda_irrig_table.pl $tmpdir/$file.asc $statenum $year $county_table > $outdir/$outfile";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}
