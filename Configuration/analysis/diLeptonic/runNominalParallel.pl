#!/usr/bin/perl

use strict;
use warnings;

open my $in, '<', 'selectionList.txt' or die $!;
chomp(my @files = grep { /Nominal/ } grep { /\w/ } grep { !/^\s*#/ } <$in>);

while (@files) {
    my $file = shift @files;
    system("./load_Analysis -f $file &");
}

