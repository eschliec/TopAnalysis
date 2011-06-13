#!/usr/local/bin/perl -w

# script to run over a list of input files (use of wildcard is allowed) and
# to grep different cuts. These cuts are expected to be named AND ordered as:
#
#   (1) muonEtaCut
#   (2) muonPtCut
#   (3) jetsEtaCut
#   (4) jetsPtCut
#
# It sums num1 ... num8 before and after the cut was applied and puts the
# result out to screen.

# get list of input files
my @myfiles = @ARGV;

# initialize error counters
my $err1 = 0;

# initialize event counters
my $num1 = 0;
my $num2 = 0;
my $num3 = 0;
my $num4 = 0;
my $num5 = 0;
my $num6 = 0;
my $num7 = 0;
my $num8 = 0;

# loop over each file and accumulate the event numbers
foreach my $file (@myfiles) {

    # print current file name
    print "... processing file $file \n";

    # get the grep output from the shell
    my @beforeMuonEtaCut = `grep -B 1 -A 8 muonEtaCut $file | grep 'Events Before Cut'`;
    my @beforeMuonPtCut  = `grep -B 1 -A 8 muonPtCut  $file | grep 'Events Before Cut'`;
    my @beforeJetsEtaCut = `grep -B 1 -A 8 jetsEtaCut $file | grep 'Events Before Cut'`;
    my @beforeJetsPtCut  = `grep -B 1 -A 8 jetsPtCut  $file | grep 'Events Before Cut'`;
    my @afterMuonEtaCut  = `grep -B 1 -A 8 muonEtaCut $file | grep 'Events After  Cut'`;
    my @afterMuonPtCut   = `grep -B 1 -A 8 muonPtCut  $file | grep 'Events After  Cut'`;
    my @afterJetsEtaCut  = `grep -B 1 -A 8 jetsEtaCut $file | grep 'Events After  Cut'`;
    my @afterJetsPtCut   = `grep -B 1 -A 8 jetsPtCut  $file | grep 'Events After  Cut'`;

    # check the pattern and extract the numbers, if succesful
    # accumulate else voice error message ; first occurrence
    # semi leptonic second dileptonic preselection
    my $failed = 0;
    if($beforeMuonEtaCut[0] =~ /Events\sBefore\sCut\s\(Weighted\):\s+(\d+)\s\(\s+\d+\w/) {
      $num1 = $num1 + $1;
    } else {
      $failed = 1;
      print ">> ERROR: could not get first number for muonEtaCut from $file <<\n";
    };
    if($beforeMuonPtCut[0]  =~ /Events\sBefore\sCut\s\(Weighted\):\s+(\d+)\s\(\s+\d+\w/) {
      $num3 = $num3 + $1;
    } else {
      $failed = 1;
      print ">> ERROR: could not get first number for muonPtCut  from $file <<\n";
    };
    if($beforeJetsEtaCut[0] =~ /\sEvents\sBefore\sCut\s\(Weighted\):\s+(\d+)\s\(\s+\d+\w/) {
      $num5 = $num5 + $1;
    } else {
      $failed = 1;
      print ">> ERROR: could not get first number for jetsEtaCut from $file <<\n";
    };
    if($beforeJetsPtCut[0]  =~ /Events\sBefore\sCut\s\(Weighted\):\s+(\d+)\s\(\s+\d+\w/) {
      $num7 = $num7 + $1;
    } else {
      $failed = 1;
      print ">> ERROR: could not get first number for jetsPtCut  from $file <<\n";
    };

    if($afterMuonEtaCut[0]  =~ /Events\sAfter\s+Cut\s\(Weighted\):\s+(\d+)\s\(\s+\d+\w/) {
      $num2 = $num2 + $1;
    } else {
      $failed = 1;
      print ">> ERROR: could not get last  number for muonEtaCut from $file <<\n";
    };
    if($afterMuonPtCut[0]   =~ /Events\sAfter\s+Cut\s\(Weighted\):\s+(\d+)\s\(\s+\d+\w/) {
      $num4 = $num4 + $1;
    } else {
      $failed = 1;
      print ">> ERROR: Could not get last  number for muontPtCut  from $file <<\n";
    };
    if($afterJetsEtaCut[0]  =~ /Events\sAfter\s+Cut\s\(Weighted\):\s+(\d+)\s\(\s+\d+\w/) {
      $num6 = $num6 + $1;  
    } else {
      $failed = 1;
      print ">> ERROR: Could not get last  number for jetsEtaCut from $file <<\n";
    };
    if($afterJetsPtCut[0]   =~ /Events\sAfter\s+Cut\s\(Weighted\):\s+(\d+)\s\(\s+\d+\w/) { 
      $num8 = $num8 + $1;  
    } else {
      $failed = 1;
      print ">> ERROR: Could not get last  number for jetsPtCut  from $file <<\n";
    };
    if($failed == 1){ $err1 = $err1 + 1 }
  }; # end loop over files

# output record
print "\n";
print "----------------------------------------------";
print "\n  Accumulated errors    :\n";
print "----------------------------------------------";
print "\n";
print "  N/A record: \t $err1 \n";
print "----------------------------------------------";
print "\n  Accumulated statistics:\n";
print "----------------------------------------------";
print "\n            \t before: \t after: \n";
print "  muonEtaCut: \t $num1 \t $num2 \n";
print "  muonPtCut : \t $num3 \t $num4 \n";
print "  jetsEtaCut: \t $num5 \t $num6 \n";
print "  jetsPtCut : \t $num7 \t $num8 \n";
print "----------------------------------------------\n";

exit;
