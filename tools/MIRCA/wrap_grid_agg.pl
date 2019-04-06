#!/usr/bin/perl

$indir = shift;
$prefix = shift;
$type = shift;
$res_ratio = shift;
$coordprec = shift;
$dataprec = shift;
$outdir = shift;

`mkdir -p $outdir`;

opendir(DIR,$indir) or die "$0: ERROR: cannot open dir $indir for reading\n";
@files = grep /^$prefix/, readdir(DIR);
closedir(DIR);

foreach $file (sort(@files)) {
  $cmd = "grid_agg $indir/$file $type $res_ratio $coordprec $dataprec $outdir/$file";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
}
