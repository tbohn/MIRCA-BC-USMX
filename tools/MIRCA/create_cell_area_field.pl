#!/usr/bin/perl
use Math::Trig;

$gridfile = shift;

open(FILE,$gridfile) or die "$0: ERROR: cannot open file $gridfile for reading\n";
for ($count=0; $count<6; $count++) {
  $line = <FILE>;
  chomp $line;
  @fields = split /\s+/, $line;
  if ($fields[0] =~ /ncols/i) {
    $ncols = $fields[1];
  }
  elsif ($fields[0] =~ /nrows/i) {
    $nrows = $fields[1];
  }
  elsif ($fields[0] =~ /xllcorner/i) {
    $xllcorner = $fields[1];
  }
  elsif ($fields[0] =~ /yllcorner/i) {
    $yllcorner = $fields[1];
  }
  elsif ($fields[0] =~ /cellsize/i) {
    $cellsize = $fields[1];
  }
  elsif ($fields[0] =~ /nodata/i) {
    $nodata = $fields[1];
  }
  print "$line\n";
}
close(FILE);

$R = 6371.009; # km
$nsteps = 100;
$dlat = $cellsize / $nsteps;
for ($row=$nrows-1; $row>=0; $row--) {
  $lat_south = $yllcorner + $row * $cellsize;
  $lat_center = $lat_south + 0.5 * $cellsize;
  $area = 0;
  for ($y=0; $y<$nsteps; $y++) {
    $lat = $lat_south + ($y+0.5) * $dlat;
    $width = 2 * pi * ($R * cos(deg2rad($lat))) * ($cellsize / 360);
    $height = 2 * pi * $R * ($cellsize / $nsteps / 360);
    $area += $width * $height;
  }
  for ($col=0; $col<$ncols; $col++) {
    printf "%.4f", $area;
    if ($col < $ncols-1) {
      print " ";
    }
    else {
      print "\n";
    }
  }
}
