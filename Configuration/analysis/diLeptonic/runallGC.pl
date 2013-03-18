#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use File::Spec;
use File::Path;
use POSIX;

my %arg;
getopts('d:c:r:sj:m:f:R:gh', \%arg);

unless ( $arg{d} && $arg{c} && -f $arg{c}) {
    print <<'USAGE';
Syntax:
 $ runallGC.pl [-g] -d directoryName -c configFile.py [-r regexp] [-s]
        [-j NJobs] [-x factor] [-m maxEventsPerJob] [-f files.txt]

Run runallGC.pl to run over all data samples given in files.txt using the
configuration file configFile.py - and use GridControl to submit jobs

A directory on the sonar user space named "GridControl_workingDir" will
be created to contain working directories and output.

Set environment variable HN_USER to your hypernews name!

-g      run on the grid

-r regexp
If a regexp is given only jobs matching regexp will be submitted (case sensitive).

-R regexp
If a regexp is given only jobs NOT matching the regexp will be submitted (case sensitive).

-s                    really submit jobs using GridControl
-j NJobs              submit NJobs for each dataset. 
-m maxEventsPerJob    Overwrites default value from files.txt

-o output            output directory when output is downloaded from dcache

-f files.txt          (default value: files.txt)
files.txt must contain lines of
cff_file root_output number_of_jobs
# comments with # and empty lines are allowed
Variables can be added via ${VAR} and will be replaced with environment variable $VAR

USAGE
    exit;
}

my $source = $arg{f} ? $arg{f} : 'files.txt';
open my $IN, '<', $source or die "Cannot open input file: $!\n";
my @createDirs;
my @runGCs;
my @checkGCs;
my @getGCs;
my @killGCs;
my @forHadd;
my @forJson;

my $hypernewsName = $arg{h} || $ENV{HN_USER};
die "Who are you? Set your hypernews name (HN_USER).\n" unless $hypernewsName;

# prepare workingdir

my $globalGcWorkingdir="/scratch/hh/dust/naf/cms/user/$ENV{USER}/GridControl_workingDir";

unless(-e $globalGcWorkingdir){
    print "creating global grid-control working dir at $globalGcWorkingdir...\n";
    mkdir($globalGcWorkingdir);
}
unless(-e "${globalGcWorkingdir}/grid-control"){
    print "checking out stable grid-control version...\n";
    system("svn co -q https://ekptrac.physik.uni-karlsruhe.de/public/grid-control/tags/stable/grid-control ${globalGcWorkingdir}/grid-control");
    print "\n\n";
}

my $workDirWithTime="$arg{d}_" . strftime("%FT%H-%M-%S",localtime);

if ($arg{g}) {
    push @createDirs, "srmmkdir srm://dcache-se-cms.desy.de:8443/pnfs/desy.de/cms/tier2/store/user/$hypernewsName/$workDirWithTime\n";
}


while(my $line = <$IN>) {
    chomp $line;
    next if $line =~ /^\s*#/; #skip comments
    next if exists $arg{'R'} && $line =~ /$arg{'R'}/;
    next unless $line =~ /\w/; #skip empty lines
    next unless (!exists $arg{'r'} || $line =~ /$arg{'r'}/);
    
    my ($eventsPerJob, $dataset, $outputFile, $options, $jsonFile) =
        map { s!\${(\w+)}!$ENV{$1}!g; $_ }
        split / {2,}/, $line;
        
    next unless $options;
    
    my $path = "${globalGcWorkingdir}/$workDirWithTime/GC-$outputFile";
    mkpath($path);
    my $jobdir = "output-$outputFile"; #. strftime("%FT%H-%M-%S",localtime);
    
    my $t = getGCtemplate();

    my $json = File::Spec->rel2abs($jsonFile);

    for ($t) {
        s/##LOCAL_OR_GRID##/$arg{g} ? "grid" : "local"/eg;
        s/##WORKDIR##/$path/g;
        s/##JOBS##/$arg{j} ? $arg{j} : ''/eg;
        s/##CMSSW_BASE##/$ENV{CMSSW_BASE}/g;
        s/##CMSSWConfigFile##/File::Spec->rel2abs($arg{c})/eg;
        s/##EVENTS_PER_JOB##/$eventsPerJob/g;
        s/##DATASET##/$dataset/g;
        s/##OUTPUTFILE##/$outputFile/g;
        s/##HN_USER##/$hypernewsName/g;
        s/##SEOUTDIR##/$path/g;
        s/##OPTIONS##/$options/g;
	s/##FILTERLUMI##/$jsonFile ? '' : ';'/eg
        s/##JSON##/$json/g;
        s/##JOBDIR##/$jobdir/g;
        s/##SE_DCACHE##/$arg{g} ? '' : ';'/eg;
        s/##SE_SCRATCH##/$arg{g} ? ';' : ''/eg;        
    }
    my $gcConfig = "${globalGcWorkingdir}/$workDirWithTime/GC-$outputFile.conf";
    open my $OUT, '>', $gcConfig or die $!;
    print $OUT $t;
    close $OUT;
    
    if ($arg{g}) {
        push @createDirs, "srmmkdir srm://dcache-se-cms.desy.de:8443/pnfs/desy.de/cms/tier2/store/user/$hypernewsName/$workDirWithTime/$jobdir\n";
	push @getGCs, ".${globalGcWorkingdir}/grid-control/scripts/downloadFromSE.py -m -o ${globalGcWorkingdir}/$workDirWithTime/$jobdir $gcConfig\n";
    } else {
        push @createDirs, "mkdir -p ${globalGcWorkingdir}/$workDirWithTime/$jobdir\n";
    }
    push @runGCs, "echo '\\n\\n---------############----------\\n\\nstarting $gcConfig'\n.${globalGcWorkingdir}/grid-control/go.py -i $gcConfig\n";
    push @checkGCs, "echo '\\n\\n---------############----------\\n\\nchecking $gcConfig'\n.${globalGcWorkingdir}/grid-control/go.py -r -R $gcConfig\n";
    push @killGCs, "echo '\\n\\n---------############----------\\n\\nkilling $gcConfig'\n.${globalGcWorkingdir}/grid-control/go.py -d ALL $gcConfig\n";
    push @forHadd, "hadd ${globalGcWorkingdir}/$workDirWithTime/${outputFile}.root ${globalGcWorkingdir}/$workDirWithTime/$jobdir/*.root\n";
    push @forJson, ".${globalGcWorkingdir}/grid-control/scripts/lumiInfo.py -j $gcConfig\n";

}

open my $OUT, '>', "${globalGcWorkingdir}/$workDirWithTime/createDirsAndRun.sh" or die $!;
print $OUT $_ for @createDirs;
print $OUT "\n" x 3;
print $OUT $_ for @runGCs;
close $OUT;
chmod 0755, "${globalGcWorkingdir}/$workDirWithTime/createDirsAndRun.sh";
print "run .${globalGcWorkingdir}/$workDirWithTime/createDirsAndRun.sh to create dirs and submit jobs\n";

open my $OUTC, '>', "${globalGcWorkingdir}/$workDirWithTime/checkAllJobs.sh" or die $!;
print $OUTC $_ for @checkGCs;
close $OUTC;
chmod 0755, "${globalGcWorkingdir}/$workDirWithTime/checkAllJobs.sh";
print "run .${globalGcWorkingdir}/$workDirWithTime/checkAllJobs.sh to check dirs and resubmit jobs\n";

if ($arg{g}) {
    open my $OUTG, '>', "${globalGcWorkingdir}/$workDirWithTime/downloadAllFromSE.sh" or die $!;
    print $OUTG $_ for @getGCs;
    close $OUTG;
    chmod 0755, "${globalGcWorkingdir}/$workDirWithTime/downloadAllFromSE.sh";
    print "run .${globalGcWorkingdir}/$workDirWithTime/downloadAllFromSE.sh to download from dcache to your sonar user space\n";
}

open my $OUTD, '>', "${globalGcWorkingdir}/$workDirWithTime/killAllJobs.sh" or die $!;
print $OUTD $_ for @killGCs;
close $OUTD;
chmod 0755, "${globalGcWorkingdir}/$workDirWithTime/killAllJobs.sh";
print "run .${globalGcWorkingdir}/$workDirWithTime/killAllJobs.sh to kill all jobs\n";

open my $OUTR, '>', "${globalGcWorkingdir}/$workDirWithTime/haddAllRoot.sh" or die $!;
print $OUTR $_ for @forHadd;
close $OUTR;
chmod 0755, "${globalGcWorkingdir}/$workDirWithTime/haddAllRoot.sh";
print "run .${globalGcWorkingdir}/$workDirWithTime/haddAllRoot.sh to hadd all output\n";

open my $OUTL, '>', "${globalGcWorkingdir}/$workDirWithTime/getJsons.sh" or die $!;
print $OUTL $_ for @forJson;
close $OUTL;
chmod 0755, "${globalGcWorkingdir}/$workDirWithTime/getJsons.sh";
print "run .${globalGcWorkingdir}/$workDirWithTime/getJsons.sh to get all processed Jsons\n";


sub getGCtemplate
{
    return <<'ENDE_TEMPLATE';
[global]
module  = CMSSW
backend = ##LOCAL_OR_GRID## ; (local or grid) ; switch is enough
workdir = ##WORKDIR## ; (when running locally also the sandbox will be placed here ; 

[local]
wms = SGE

[SGE]
site = hh 

[grid]
sites     =  desy.de physik.rwth-aachen.de -colorado.edu -hep.wisc.edu -hep.ucl.ac.uk -in2p3.fr

[jobs]
##JOBS##
wall time    = 40:59:00 ;queue dingens // thats the 12 hour cue!!!
; in queue     = $maxjobsinqueue ; 1000 -1
memory       = 3700
monitor jobs = dashboard
max retry    = 1


[CMSSW]
project area       = ##CMSSW_BASE##
se runtime         = True
config file        = ##CMSSWConfigFile## ; 
prepare config     = True
arguments          = ##OPTIONS## outputFile=##OUTPUTFILE##  ; 

dataset splitter   = EventBoundarySplitter
events per job     = ##EVENTS_PER_JOB##   ; I think better than to restrict total number of jobs


dataset            = ##DATASET## ; /MuEG/Run2012A-recover-06Aug2012-v1/AOD
##FILTERLUMI##lumi filter        = ##JSON## ;

[storage]
se output files    = *.root ;
##SE_DCACHE##se path            = srm://dcache-se-cms.desy.de:8443/pnfs/desy.de/cms/tier2/store/user/##HN_USER##/##JOBDIR##
##SE_SCRATCH##se path           = dir://##SEOUTDIR##/##JOBDIR##





; $workdir must be created before to avoid script from stopping
; srmmkdir srm://dcache-se-cms.desy.de:8443/pnfs/desy.de/cms/tier2/store/user/${HN_USER}/##JOBDIR##
; 
; ./grid-control/go.py -i [gc_configfile]  ; inits working dir, runs the jobs, then quits
; ./grid-control/scripts/downloadFromSE.py -ml --out=$outputdir [gc_configfile]  ; loops in background, checks if files are there and if ok, moves them to outputdir
; ./grid-control/go.py [gc_configfile] ; checks and resubmits
ENDE_TEMPLATE

}

