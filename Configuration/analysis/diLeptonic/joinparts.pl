#!/usr/bin/perl

use strict;
use warnings;

open my $FH, '<joined.txt' or die $!;
my $joined = <$FH>;
if (-e $joined) { exit }
$joined =~ s/.root/_partly.root/;

my ($f) = glob("j_*.sh");
$f =~ s/^j_//;
$f =~ s/\.sh$//;

my @toJoin;
for my $job (1..100) {
    if (-e "$f-$job.root") {
        push @toJoin, "$f-$job.root";
    } else {
        if (-e "$f-$job.root.part.1") {
            push @toJoin, "$f-$job.root.part.*";
        }
    }
}

if (@toJoin) {
    print "hadd -f $joined @toJoin\n";
    system("hadd -f $joined @toJoin\n");
}
