#!/usr/bin/perl

use English;

die "args $#ARGV usage: beadmaker.pl numberofbeads spreadcoordbinary relativepathtodatabasedir molecularsystemdirectory molsysteminput"  if ($#ARGV != 4);
my $numbeads=shift;
my $spreadx=shift;
my $databasedir=shift;
my $molsystem=shift;
my $molsysteminput=shift;

`mkdir BEADS_$numbeads`;
chdir("BEADS_$numbeads");
`ln -s ../$databasedir`;
`cp -ar ../$molsystem .`;
chdir("$molsystem");
my $dir="ATOM_COORDS_IN";
#write spreadcoord input file
open( my $spf, ">", "spread.input") || die "cannot open spread.input for writing\n";
print $spf "$numbeads 1475145.0 initial\n";
print $spf "10.0 300.0\n";
print $spf "$dir/Bead.0_Temper.0/water.coords_initial $dir water.coords_initial\n";
close $spf;

#make ATOM directories
for(my $i=1; $i<$numbeads;$i++)
{

#input
    my $target=$dir.'/Bead.'.$i.'_Temper.0';
    `mkdir $target`;
#output
    $target=ATOM_COORDS_OUT.'/Bead.'.$i.'_Temper.0';
    `mkdir $target`;
    `touch $target/ChkDirExistOAtom`;
}


#make STATE directories
for(my $i=1; $i<$numbeads;$i++)
{

#input
    my $target='STATES/Spin.0_Kpt.0_Bead.'.$i.'_Temper.0';
    `mkdir $target`;
    `(cd $target; cp ../Spin.0_Kpt.0_Bead.0_Temper.0/* .)`;
#output
    $target='STATES_OUT/Spin.0_Kpt.0_Bead.'.$i.'_Temper.0';
    `mkdir $target`;
    `touch $target/ChkDirExistOAtom`;
}

my $beadsminus1=$numbeads-1;
system("$spreadx spread.input $beadsminus1");
my $outputfile=$molsysteminput.".bead";
`cp $molsysteminput $outputfile`;
open($conf,">>",$outputfile) ||die "cannot open $outputfile for appending";
print $conf '~sim_pimd_def['."\n".'\path_int_beads{'.$numbeads.'}'."\n".']'."\n";



