#!/usr/bin/perl

use English;


my $infilename=shift;
my $jobid=shift;
my $ppn=shift;

$jobid="$infilename" if length($jobid)<1;
chomp $infilename;
my($infile, $outfile);
open($infile, "<", "$infilename") or die "cannot open $infilename";
my $nodelistname=".nodelist.$jobid";
open($outfile, ">","$nodelistname") or die "cannot open $nodelistname";
print $outfile "group main +shell ssh\n";
print "creating nodelist file $nodelistname\n";
my $countppn=1;
my $prevnode="";
my $numnodes=0;
foreach my $node (<$infile>)
{
  chomp $node;
  print $outfile "host $node\n";
  if($prevnode eq $node)
  {
      $countppn++;
  }
  else
  {
      $countppn=1;
      $numnodes++;
  }
  $prevnode=$node;
}
close $infile;
close $outfile;

my $procspernode=1;
if(length($ppn)>0)
{
    $procspernode=$countppn/$ppn;
}
else
{
    $ppn=$countppn;
}
$nodelistname.=".smp";
open($infile, "sort -u $infilename |") or die "cannot open $infilename";
open($outfile, ">","$nodelistname") or die "cannot open $nodelistname";
print $outfile "group main +shell ssh\n";
print "creating nodelist file $nodelistname\n";
foreach my $node (<$infile>)
{
  chomp $node;
  for(my $i=0; $i< $procspernode; $i++)
  {
      print $outfile "host $node +ncpus $ppn\n";
  }
}
close $infile;
close $outfile;

