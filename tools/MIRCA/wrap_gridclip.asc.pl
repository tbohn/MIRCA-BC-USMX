#!/usr/bin/perl

$indir = shift;
$prefix = shift;
$south = shift;
$north = shift;
$west = shift;
$east = shift;
$outdir = shift;
$verbose = shift;

`mkdir -p $outdir`;

opendir(DIR,$indir) or die "$0: ERROR: cannot open dir $indir for reading\n";
@files = grep /^$prefix/, readdir(DIR);
closedir(DIR);

foreach $file (sort(@files)) {

  $cmd = "gridclip $indir/$file float $west $east $south $north 5 4 $outdir/$file";
  if ($verbose) { print "$cmd\n"; }
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}


