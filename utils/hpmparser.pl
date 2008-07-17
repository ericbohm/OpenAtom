#!/usr/bin/perl

use English;
use strict;

# HPM files a hpm_data.N where N is the node number

# by default each node counts chip 0 1, can be conf to count 2,3

# so methodologically if you want a per chip amount you need 2 runs to get
# both sets.  If you can assume reasonably uniform distribution you can use 1
# run and multiply by 2.

#first cut a this script assumes latter for simplicity
my @pernode;
my %sumpernode;

#HPM author suggests FDIV is 13 ops, based on the compiler doing fancy tricks
#but UPC is what the chip was actually asked to do so that seems irrelevant
#unless you're trying to reconcile a manually computed FLOPS with counted FLOPS

# UPC authors state FMA_2 is 2 flops
# UPC authors state FMA_4 is 4 flops
my %FP_KEYS=( 'BGP_PU0_FPU_ADD_SUB_1' => 1,
		'BGP_PU0_FPU_MULT_1' => 1,
		'BGP_PU0_FPU_FMA_2' => 2,
		'BGP_PU0_FPU_DIV_1' => 1,
		'BGP_PU0_FPU_OTHER_NON_STORAGE_OPS' => 1,
		'BGP_PU0_FPU_ADD_SUB_2' => 1,
		'BGP_PU0_FPU_MULT_2' => 1,
		'BGP_PU0_FPU_FMA_4' => 4,
		'BGP_PU0_FPU_DUAL_PIPE_OTHER_NON_STORAGE_OPS' => 1,
		'BGP_PU1_FPU_ADD_SUB_1' => 1,
		'BGP_PU1_FPU_MULT_1' => 1, 
		'BGP_PU1_FPU_FMA_2' => 2,
		'BGP_PU1_FPU_DIV_1' => 1,
		'BGP_PU1_FPU_OTHER_NON_STORAGE_OPS' => 1,
		'BGP_PU1_FPU_ADD_SUB_2' => 1,
		'BGP_PU1_FPU_MULT_2' => 1,
		'BGP_PU1_FPU_FMA_4' => 4,
		'BGP_PU1_FPU_DUAL_PIPE_OTHER_NON_STORAGE_OPS' => 1
		);

foreach my $infile (@ARGV)
{
    my $file= open(FILE,"<$infile") || die "cannot open $infile";
    my $node=0;
    if($infile =~ /.*\.(\d+)/)
    {
	$node= $1 ;
    }
    my %thisnode;
    while(<FILE>)
    {
	# extract per node data
	chomp;
	my @column=split;
#	print "columns:$column[0] $column[1] $column[2] \n";
	next if $column[0] !~ /\d+/;
	$thisnode{$column[2]}=$column[1];
	$sumpernode{$column[2]}+=$column[1];
    }
    $pernode[$node]=\%thisnode;
    # keep running totals

}
my $rawflops=0;
my $scaledflops=0;
foreach my $fpreg (keys %FP_KEYS)
{
    my $scaled = $sumpernode{$fpreg}*$FP_KEYS{$fpreg};
    print "$fpreg total $sumpernode{$fpreg} scaled $scaled by $FP_KEYS{$fpreg}\n";
    $rawflops+= $sumpernode{$fpreg};
    $scaledflops+= $scaled;
}
my $allchipsrflops=$rawflops*2;
my $allchipssflops=$scaledflops*2;
print "raw flops $rawflops scaled $scaledflops\n";
print "raw flops by 2 $allchipsrflops scaled $allchipssflops\n";
