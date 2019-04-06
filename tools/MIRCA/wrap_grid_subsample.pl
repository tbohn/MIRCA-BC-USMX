#!/usr/bin/perl

$indir = shift;
$prefix = shift;
$new_cellsize = shift;
$outdir = shift;

`mkdir -p $outdir`;

opendir(DIR,$indir) or die "$0: ERROR: cannot open dir $indir for reading\n";
@files = grep /^$prefix/, readdir(DIR);
closedir(DIR);

foreach $file (sort(@files)) {

  $cmd = "grid_subsample $indir/$file float $new_cellsize nn 0 float 5 4 $outdir/$file";
  print "$cmd\n";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";

}


