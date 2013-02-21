#!/usr/bin/perl -w

# % cat Plots_temp/PDF_CENTRAL/combined/HypToppTLaTeX.txt
# Variable: HypToppT   Channel: Dilepton Combined
# BinCenter & LowXbinEdge  &  HighXbinEdge  &   DiffXSec  &  StatError(\%)  & SystError(\%)  & TotalError(\%) \\
# \hline
# $40$ & $0$ to $80$   &  0.0051167  &   1.33 &    0 &    1.33 \\
# $105$ & $80$ to $130$   &  0.0058073  &   1.5 &    0 &    1.5 \\
# $165$ & $130$ to $200$   &  0.0030112  &   1.45 &    0 &    1.45 \\
# $250$ & $200$ to $300$   &  0.00084305  &   1.76 &    0 &    1.76 \\
# $350$ & $300$ to $400$   &  0.00016227  &   3.4 &    0 &    3.4 \\

use strict;
use warnings;
use List::Util qw(max);
use File::Path qw(mkpath);

use Data::Dumper;

my $outputPath = "UnfoldingResults/PDF_";

sub pdfSum {
    my ($quantity, $channel) = @_;
    my $nominal = UnfoldedResult->new("Plots/PDF_CENTRAL/$channel/${quantity}LaTeX.txt");
    #print Dumper $result;
    my (@up, @down);
    for my $var_no (1..22) { 
        push @up, UnfoldedResult->new("Plots/PDF_${var_no}_UP/$channel/${quantity}LaTeX.txt");
        push @down, UnfoldedResult->new("Plots/PDF_${var_no}_DOWN/$channel/${quantity}LaTeX.txt");
    }
    
    my $filename = "$outputPath/$channel/${quantity}Results.txt";
    open my $OUTFH, '>', $filename or die $!;
    print "Writing to $filename\n";
    
    for my $bin (0..@{$nominal->{bins}}-1) {
        my $binNom = $nominal->{bins}->[$bin];
        my $nominalXsec = $binNom->{-xsec};
        my ($upSUM2, $downSUM2) = (0,0);
        for my $no (0..21) {
            my $upMnom = $up[$no]{bins}[$bin]{-xsec} - $nominalXsec;
            my $nomMdown = $nominalXsec - $down[$no]{bins}[$bin]{-xsec};
            $upSUM2 += max(0, $upMnom, $nomMdown)**2;
            $downSUM2 += max(0, -$upMnom, -$nomMdown)**2;
        }
        my $relUncUp = sqrt($upSUM2)/$nominalXsec;
        my $relUncDown = sqrt($downSUM2)/$nominalXsec;
#         my $unc = 0.5 * ($relUncUp + $relUncDown);
        my $unc = max($relUncUp, $relUncDown);
        #my $line = "XAxisbinCenters[bin]: $binNom->{-center} bin: $binNom->{-from} to $binNom->{-to} SystematicRelError: $unc (+$relUncUp -$relUncDown)\n";
        my $line = "XAxisbinCenters[bin]: $binNom->{-center} bin: $binNom->{-from} to $binNom->{-to} SystematicRelError: $unc\n";
        print $line;
        print $OUTFH $line;
    }
    
}

for my $channel qw(ee emu mumu combined) {
    mkpath("$outputPath/$channel");
    for my $plot qw(
    HypLeptonBjetMass
    HypBJetpTLead
    HypBJetpTNLead
    HypBJetEtaLead
    HypBJetEtaNLead
    HypTopRapidityLead
    HypTopRapidityNLead
    HypToppTLead
    HypToppTNLead
    HypLeptonpTLead
    HypLeptonpTNLead
    HypLeptonEtaLead
    HypLeptonEtaNLead
    HypLeptonEta
    HypLeptonpT
    HypBJetEta
    HypBJetpT
    HypTopRapidity
    HypToppT
    HypTTBarRapidity
    HypTTBarpT
    HypTTBarMass
    HypLLBarpT
    HypLLBarMass
    HypNeutrinopT
    ) 
    {
        print "Unc for $plot in channel $channel:\n";
        eval{ pdfSum($plot, $channel); };
        print "\n";
    }
}

package UnfoldedResult;

sub new {
    my ($class, $filename) = @_;
    my $self = {};
    $self->{filename} = $filename;
    bless $self, $class;
    $self->{bins} = $self->readFile();
    return $self;
}

sub readFile {
    my $self = shift;
    open my $f, '<', $self->{filename} or die $!;
    scalar <$f> for (1..3); #skip first 3 lines
    my @bins;
    while (my $line = <$f>) {
        $line =~ s/\$//g;
        my @cols = split / +/, $line;
        push @bins, {-center => $cols[0], 
                     -from => $cols[2], -to => $cols[4], 
                     -xsec => $cols[6], -xsecstatunc => $cols[8]};
    }
    return \@bins;
}



