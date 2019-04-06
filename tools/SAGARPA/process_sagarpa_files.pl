#!/usr/bin/perl

$file = shift;

@irrtypes = ("irr","noirr");

open(FILE,$file) or die "$0: ERROR: cannot open file $file for reading\n";
foreach (<FILE>) {
  if (!/^\d/) { next; }
  chomp;
  @fields = split /,/;
  $year = $fields[0];
  $stateID = $fields[1];
  $state = $fields[2];
  $munID = $fields[7];
  $mun = $fields[8];
  if ($fields[15] == 1) {
    $season = "w";
  }
  elsif ($fields[15] == 2) {
    $season = "s";
  }
  elsif ($fields[15] == 3) {
    $season = "p";
  }
  $irrig = $fields[17]; # 0 = rainfed, 1 = irrigated
  if ($fields[17] == 1) {
    $irrig = "irr";
  }
  elsif ($fields[17] == 2) {
    $irrig = "noirr";
  }
  $Aplanted = $fields[19]*0.01; # ha to km2
  $Aharvest = $fields[20]*0.01; # ha to km2
  if ($munID > 0) {
    $state_mun_key = sprintf "%02d%03d", $stateID, $munID;
    $AplantedTotal{$state_mun_key}{$year}{$irrig}{$season} += $Aplanted;
    $AharvestTotal{$state_mun_key}{$year}{$irrig}{$season} += $Aharvest;
  }
  $state_mun_key = sprintf "%02d000", $stateID;
  $AplantedTotal{$state_mun_key}{$year}{$irrig}{$season} += $Aplanted;
  $AharvestTotal{$state_mun_key}{$year}{$irrig}{$season} += $Aharvest;
  $AplantedTotal{"entire"}{$year}{$irrig}{$season} += $Aplanted;
  $AharvestTotal{"entire"}{$year}{$irrig}{$season} += $Aharvest;
}
close(FILE);

foreach $state (sort(keys(%AplantedTotal))) {
  foreach $year (sort(keys(%{$AplantedTotal{$state}}))) {
    print "$state,$year,";
    foreach $irrig (@irrtypes) {
      print "$irrig,";
      foreach $season ("p","w","s") {
        if ($AplantedTotal{$state}{$year}{$irrig}{$season} eq "") {
          $AplantedTotal{$state}{$year}{$irrig}{$season} = 0;
        }
        if ($AharvestTotal{$state}{$year}{$irrig}{$season} eq "") {
          $AharvestTotal{$state}{$year}{$irrig}{$season} = 0;
        }
        print "$season,$AplantedTotal{$state}{$year}{$irrig}{$season},$AharvestTotal{$state}{$year}{$irrig}{$season},";
      }
      $tmp1 = $AplantedTotal{$state}{$year}{$irrig}{"p"} + $AplantedTotal{$state}{$year}{$irrig}{"w"};
      $tmp2 = $AharvestTotal{$state}{$year}{$irrig}{"p"} + $AharvestTotal{$state}{$year}{$irrig}{"w"};
      print "wp,$tmp1,$tmp2,";
      $tmp1 = $AplantedTotal{$state}{$year}{$irrig}{"p"} + $AplantedTotal{$state}{$year}{$irrig}{"s"};
      $tmp2 = $AharvestTotal{$state}{$year}{$irrig}{"p"} + $AharvestTotal{$state}{$year}{$irrig}{"s"};
      print "sp,$tmp1,$tmp2,";
      $tmp1 = $AplantedTotal{$state}{$year}{$irrig}{"p"} + $AplantedTotal{$state}{$year}{$irrig}{"s"} + $AplantedTotal{$state}{$year}{$irrig}{"w"};
      $tmp2 = $AharvestTotal{$state}{$year}{$irrig}{"p"} + $AharvestTotal{$state}{$year}{$irrig}{"s"} + $AharvestTotal{$state}{$year}{$irrig}{"w"};
      print "wsp,$tmp1,$tmp2,";
    }
    print "\n";
  }
}
