#!/usr/bin/perl

$indir = shift;
$outdir = shift;

`mkdir -p $outdir`;

%outpfx_map = (
  "crop_all_all_" => "fplant_of_cell_",
  "crop_all_irrigated_" => "firr_of_cell_",
);

foreach $inpfx (sort(keys(%outpfx_map))) {
  opendir(DIR,$indir) or die "$0: ERROR: cannot open dir $indir for reading\n";
  @files = grep /^$inpfx/, readdir(DIR);
  closedir(DIR);

  foreach $file (@files) {
    $outfile = $file;
    $outfile =~ s/$inpfx/$outpfx_map{$inpfx}/g;
    $cmd = "cp $indir/$file $outdir/$outfile";
    (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
  }
}

for ($m=0; $m<12; $m++) {
  $month = sprintf "%03d", $m+1;
  $cmd = "grid_math $outdir/firr_of_cell_$month.asc float $outdir/fplant_of_cell_$month.asc float / -1 float 0 5 4 $outdir/firr_of_fplant_$month.asc";
  (system($cmd)==0) or die "$0: ERROR: $cmd failed\n";
}
