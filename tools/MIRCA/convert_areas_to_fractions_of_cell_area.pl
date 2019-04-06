#!/usr/bin/perl

$indir = shift;
$prefix = shift;
$outdir = shift;

`mkdir -p $outdir`;

opendir(DIR,$indir) or die "$0: ERROR: cannot open dir $indir for reading\n";
@files = grep /^$prefix/, readdir(DIR);
closedir(DIR);

$first = 1;
foreach $file (sort(@files)) {
  if ($first) {
    # Compute grid cell areas in km2
    $cmd = "create_cell_area_field.pl $indir/$file > $outdir/area_field.asc";
    (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
    # Convert to hectares
    $cmd = "grid_math $outdir/area_field.asc float 100 const \"*\" -1 float 0 16 2 $outdir/area_field.ha.asc";
    (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
    $first = 0;
  }
  $cmd = "grid_math $indir/$file float $outdir/area_field.ha.asc float / -1 float 0 16 4 $outdir/$file";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
}
