#!/usr/bin/perl

$file = shift;

open(FILE,$file) or die "$0: ERROR: cannot open file $file for reading\n";
foreach (<FILE>) {
  chomp;
  s/\s+$//;
#print "$_\n";
  @fields = split /,/;
  $ncols = @fields;
#print "ncols $ncols\n";
  @newfields = ();
  $tmp = "";
  $col = 0;
  $read = 1;
  while ($col<$ncols) {
#print "col $col fields[col] $fields[$col]\n";
    if ($fields[$col] !~ /^"/) {
      $tmp = $fields[$col];
#print "tmp $tmp read $read\n";
    }
    else {
      while ($read) {
        $tmp = $tmp . $fields[$col];
        if ($fields[$col] =~ /"$/) {
          $read = 0;
        }
#print "tmp $tmp read $read\n";
        $col++;
      }
      $col--;
    }
    $tmp =~ s/"//g;
    push @newfields, $tmp;
#print "newfields @newfields\n";
    $tmp = "";
    $read = 1;
    $col++;
  }
  $line = join ",", @newfields;
  print "$line\n";
}
close(FILE);
