#!/usr/bin/perl

$file = shift;
$statenum = shift; # number from 1 to 50
$year = shift;
$county_table = shift;

$factor = 0.004046856422; # conv acres to km2

# US state statecodes in the pop stats tables
$c = 1;
for ($i=0; $i<50; $i++) {
  $key = sprintf "%02d", $i+1;
  $statecodes{$key} = sprintf "%02d", $c;
  if ($i == 1 || $i == 4 || $i == 7 || $i == 9 || $i == 37 || $i == 45) {
    $c++;
  }
  $c++;
}
$statenum = sprintf "%02d", $statenum;
$this_statecode = $statecodes{$statenum};

open(FILE,$county_table) or die "$0: ERROR: cannot open file $county_table for reading\n";
foreach (<FILE>) {
  if (/^#/) { next; }
  chomp;
  @fields = split /,/;
  if ($fields[2] ne "") {
    $county_name_usda{$fields[0]} = $fields[2];
    $county_name_census{$fields[0]} = $fields[1];
    if ($fields[0] =~ /(\d\d?)\d\d\d/) {
      $statecode_tmp = sprintf "%02d", $1;
    }
    $county_code{$statecode_tmp}{$fields[2]} = $fields[0];
  }
}
close(FILE);

open(FILE,$file) or die "$0: ERROR: cannot open file $file for reading\n";
$newpage = 0;
$count = 0;
foreach (<FILE>) {
#print "count $count\n";
#print;
  chomp;
  if (($year == 1992 && /Table 8/) || ($year != 1992 && /Table 10/)) {
    $newpage = 1;
    $count = 0;
  }
  if (($year == 1992 && $count == 5) || ($year != 1992 && $count == 3)) {
#print "$_\n";
    s/(\w\.?) (\w)/$1_$2/g;
#print "$_\n";
    s/[ ]+/,/g;
#print "$_\n";
    s/^,//;
#print "$_\n";
    if ($year != 1992) {
      s/Item,//;
    }
#print "$_\n";
    @tmparray = split /,/;
    $any_valid = 0;
    for ($i=0; $i<@tmparray; $i++) {
      $code = $county_code{$statecodes{$statenum}}{$tmparray[$i]};
      if ($code eq "") {
        $county_tmp = $tmparray[$i];
        $county_tmp =~ s/_//g;
        $code = $county_code{$statecodes{$statenum}}{$county_tmp};
        if ($code ne "") { $tmparray[$i] = $county_tmp; }
      }
      if ($code ne "") { $any_valid = 1; }
    }
    if ($any_valid) {
      foreach $tmp (@tmparray) {
        push @counties, $tmp;
      }
      $ncounties_on_this_page = @tmparray;
    }
    else {
      print "$0: ERROR: no valid county names found\n";
    }
  }
  else {
#print "$_\n";
    s/(\d),(\d)/$1$2/g;
    s/(\d) ([^\s^\-])/$1$2/g;
#print "$_\n";
    s/[ ]+/,/g;
#print "$_\n";
    if (($year == 1992 && $count == 14) || ($year == 2002 && $count == 10) || ($year == 2012 && $count == 12)) {
#print "parsing harvest\n";
      s/^.*acres,,$year,(--,)?//;
#print "$_\n";
      s/\(D\)/0/g;
      s/,-/,0/g;
#print "$_\n";
      s/\,[^\w^\d]/,/g;
#print "$_\n";
      @tmparray = split /,/;
      $ntmp = @tmparray;
      if ($ntmp != $ncounties_on_this_page) {
        print "$0: ERROR: ndata ($ntmp) != ncounties ($ncounties_on_this_page)\n";
      }
      foreach $tmp (@tmparray) {
        $tmp =~ s/[^\w^\d]//g;
        if ($tmp eq "") { $tmp = 0; }
        push @harvest_crop, $tmp;
      }
    }
    elsif (($year == 1992 && $count == 24) || ($year == 2002 && $count == 18) || ($year == 2012 && $count == 20)) {
#print "parsing pasture\n";
      s/^.*acres,,$year,(--,)?//;
#print "$_\n";
      s/\(D\)/0/g;
      s/,-/,0/g;
#print "$_\n";
      s/\,[^\w^\d]/,/g;
#print "$_\n";
      @tmparray = split /,/;
      $ntmp = @tmparray;
      if ($ntmp != $ncounties_on_this_page) {
        print "$0: ERROR: ndata ($ntmp) != ncounties ($ncounties_on_this_page)\n";
      }
      foreach $tmp (@tmparray) {
        $tmp =~ s/[^\w^\d]//g;
        if ($tmp eq "") { $tmp = 0; }
        push @harvest_pasture, $tmp;
      }
    }
    elsif (($year == 1992 && $count == 27) || ($year == 2002 && $count == 20) || ($year == 2012 && $count == 23)) {
#print "parsing irrigated\n";
      s/^.*acres,,$year,(--,)?//;
#print "$_\n";
      s/\(D\)/0/g;
      s/,-/,0/g;
#print "$_\n";
      s/\,[^\w^\d]/,/g;
#print "$_\n";
      @tmparray = split /,/;
      $ntmp = @tmparray;
      if ($ntmp != $ncounties_on_this_page) {
        print "$0: ERROR: ndata ($ntmp) != ncounties ($ncounties_on_this_page)\n";
      }
      foreach $tmp (@tmparray) {
#print "before $tmp\n";
        $tmp =~ s/[^\w^\d]//g;
        if ($tmp eq "") { $tmp = 0; }
#print "after $tmp\n";
        push @irrigated, $tmp;
      }
    }
    elsif (($year == 1992 && $count == 31) || ($year == 2002 && $count == 24) || ($year == 2012 && $count == 27)) {
      s/^.*acres,,$year,(--,)?//;
#print "$_\n";
      s/\(D\)/0/g;
      s/,-/,0/g;
#print "$_\n";
      s/\,[^\w^\d]/,/g;
#print "$_\n";
      @tmparray = split /,/;
      $ntmp = @tmparray;
      if ($ntmp != $ncounties_on_this_page) {
        print "$0: ERROR: ndata ($ntmp) != ncounties ($ncounties_on_this_page)\n";
      }
      foreach $tmp (@tmparray) {
        $tmp =~ s/[^\w^\d]//g;
        if ($tmp eq "") { $tmp = 0; }
        push @irrig_harvest_crop, $tmp;
      }
    }
    elsif (($year == 1992 && $count == 36) || ($year == 2002 && $count == 28) || ($year == 2012 && $count == 31)) {
      s/^.*acres,,$year,(--,)?//;
#print "$_\n";
      s/\(D\)/0/g;
      s/,-/,0/g;
#print "$_\n";
      s/\,[^\w^\d]/,/g;
#print "$_\n";
      @tmparray = split /,/;
      $ntmp = @tmparray;
      if ($ntmp != $ncounties_on_this_page) {
        print "$0: ERROR: ndata ($ntmp) != ncounties ($ncounties_on_this_page)\n";
      }
      foreach $tmp (@tmparray) {
        $tmp =~ s/[^\w^\d]//g;
        if ($tmp eq "") { $tmp = 0; }
        push @irrig_pasture, $tmp;
      }
    }
  }
  $count++;
}
close(FILE);

print "#Code,State,County,AHarvest,AIrr\n";
$nCounties = @counties;
# skip first element, which is statewide total
for ($i=1; $i<$nCounties; $i++) {
  $code = sprintf "%05d", $county_code{$statecodes{$statenum}}{$counties[$i]};
  if ($code eq "") {
    print "$0: ERROR: no county code found for county $counties[$i]\n";
  }
  $AHarvest = $harvest_crop[$i] + $harvest_pasture[$i];
  $AHarvest *= $factor;
  $AIrr = $irrigated[$i];
  $AIrr *= $factor;
  $AHarvest = sprintf "%.2f", $AHarvest;
  $AIrr = sprintf "%.2f", $AIrr;
  print "$code,$counties[0],$counties[$i],$AHarvest,$AIrr\n";
}
