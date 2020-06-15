#!/usr/bin/perl -w

use IO::File;
use strict;

opendir(DIR, ".");
my @fastqc_count_files = grep(/\.fastqc\.count/, readdir(DIR));
closedir(DIR);

opendir(DIR, ".");
my @flagstat_count_files = grep(/\.flagstat\.count/, readdir(DIR));
closedir(DIR);

my $in_fh = new IO::File;
my %counts;
foreach my $file (@fastqc_count_files)
{
    $file =~ /(.+)\.fastqc\.count/;
    my $name = $1;
    
    $in_fh->open($file, "r") or die "Could not open $file\n$!";
    my $line = <$in_fh>;
    chomp $line;
    $in_fh->close();

    my ($read1, $read2) = split(/\s+/, $line);
    $counts{$name}{'fastqc'}{'R1'} = $read1;
    $counts{$name}{'fastqc'}{'R2'} = $read2;
}

foreach my $file (@flagstat_count_files)
{
    $file =~ /(.+)\.flagstat\.count/;
    my $name = $1;
    
    $in_fh->open($file, "r") or die "Could not open $file\n$!";
    my $line = <$in_fh>;
    chomp $line;
    $in_fh->close();

    my ($read1, $read2) = split(/\s+/, $line);
    $counts{$name}{'flagstat'}{'R1'} = $read1;
    $counts{$name}{'flagstat'}{'R2'} = $read2;
}

print "Sample\tFlagstat_R1\tFlagstat_R2\tFastQC_R1\tFastQC_R2\tDiff_R1\tDiff_R2\n";
foreach my $name (sort keys %counts)
{
    printf "$name\t%d\t%d\t%d\t%d\t%d\t%d\n",
    $counts{$name}{'flagstat'}{'R1'},
    $counts{$name}{'flagstat'}{'R2'},
    $counts{$name}{'fastqc'}{'R1'},
    $counts{$name}{'fastqc'}{'R2'},
    $counts{$name}{'flagstat'}{'R1'} - $counts{$name}{'fastqc'}{'R1'},
    $counts{$name}{'flagstat'}{'R2'} - $counts{$name}{'fastqc'}{'R2'};
}
