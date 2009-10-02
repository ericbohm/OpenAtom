#!/usr/bin/perl

use English;
use strict;
#compare the results of 2 executions to N places.

my $precision=shift;
my $answerfile=shift;
my $executionfile=shift;
my $columnskip=0;
if($ARGV>=0)
{
    $columnskip=shift;
}
else
{
    $columnskip=0;
}

my $ansfile=open(ANSFILE,'<', $answerfile) || die "cannot open $answerfile";
my $exefile=open(EXEFILE,'<', $executionfile) || die "cannot open $executionfile";

while(my $ans =<ANSFILE>)
{
    my $exe=<EXEFILE>;
    my @ans=split(/\s+/,$ans);
    my @exe=split(/\s+/,$exe);
    my @ans_p;
    my @exe_p;
    if($columnskip>0)
    {
	@ans=@ans[$columnskip..$#ans];
	@exe=@exe[$columnskip..$#exe];
    }
    foreach my $avalue (@ans)
    {
	push(@ans_p,trunc_to_precision($precision,$avalue));
    }
    foreach my $evalue (@exe)
    {
	push(@exe_p,trunc_to_precision($precision,$evalue));
    }
    for(my $elem=0;$elem<=$#ans_p;$elem++)
    {
#	print "comparing $ans_p[$elem] : $exe_p[$elem]\n";
	die "mismatch on line $INPUT_LINE_NUMBER $ans_p[$elem] != $exe_p[$elem]\n" if ($ans_p[$elem] != $exe_p[$elem]);
    }
}

sub trunc_to_precision
{
    my($precision,$value)=@_;
    my $format="%.$precision"."g";
    return sprintf($format,$value);
}
