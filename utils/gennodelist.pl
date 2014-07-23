#!/usr/bin/perl

use English;


my $infilename=shift;
my $jobid=shift;
$jobid="$infilename" if length($jobid)<1;
chomp $infilename;
my $infile, $outfile;
open($infile, "<", "$infilename") or die "cannot open $infilename";
my $nodelistname=".nodelist.$jobid";
open($outfile, ">","$nodelistname") or die "cannot open $nodelistname";
print $outfile "group main +shell ssh\n";
print "creating nodelist file $nodelistname\n";
foreach my $node (<$infile>)
{
    chomp $node;
    print $outfile "host $node\n";
}
close $infile;
close $outfile;

