#!/usr/bin/perl

$indir = shift;
$outdir = shift;

`mkdir -p $outdir`;

foreach $irrtype ("irrigated", "rainfed") {

  opendir(DIR,$indir) or die "$0: ERROR: cannot open dir $indir for reading\n";
  @files = grep /001.asc/, grep /$irrtype/, readdir(DIR);
  closedir(DIR);

  $outpfx = "crop_all_$irrtype";
  for ($m=1; $m<=12; $m++) {
    $month = sprintf "%03d", $m;
    $outfile = "$outpfx\_$month.asc";
    $first = 1;
    foreach $file (sort(@files)) {
      $infile = $file;
      $infile =~ s/001.asc/$month.asc/;
      if ($first) {
        $cmd = "cp $indir/$file $outdir/$outfile.tmp";
        (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
        $first = 0;
      }
      else {
        $cmd = "grid_math $indir/$file float $outdir/$outfile.tmp float + -1 float 0 16 4 $outdir/$outfile.tmp2";
        (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
        $cmd = "mv $outdir/$outfile.tmp2 $outdir/$outfile.tmp";
        (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
      }
    }
    $cmd = "mv $outdir/$outfile.tmp $outdir/$outfile";
    (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
  }

}

for ($m=1; $m<=12; $m++) {
  $month = sprintf "%03d", $m;
  $file1 = "crop_all_irrigated_$month.asc";
  $file2 = "crop_all_rainfed_$month.asc";
  $outfile = "crop_all_all_$month.asc";
  $cmd = "grid_math $outdir/$file1 float $outdir/$file2 float + -1 float 0 16 4 $outdir/$outfile";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
}
