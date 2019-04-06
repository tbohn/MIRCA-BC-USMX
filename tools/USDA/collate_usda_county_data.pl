#!/usr/bin/perl

$indir = shift;
$year = shift;
$outdir = shift;

`mkdir -p $outdir`;

opendir(DIR,$indir) or die "$0: ERROR: cannot open dir $indir for reading\n";
@files = grep /$year/, grep /^table/, readdir(DIR);
close(DIR);

$outfile = "summ.$year.csv";
$first = 1;
foreach $file (sort(@files)) {
  if ($first) {
    $cmd = "cat $indir/$file > $outfile";
    $first = 0;
  }
  else {
    $cmd = "grep -v \"^#\" $indir/$file >> $outdir/$outfile";
  }
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
}
