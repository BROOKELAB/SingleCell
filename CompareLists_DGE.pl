#!/usr/bin/perl

# ************
# Author: J. Cristobal Vera
# email: jcvera@illinois.edu
##compare multiple lists of DGE results (e.g. Seurat_StatusPB2.tsv) and return the intersection based on a minimum representation parameter
##(e.g. result must occure in two out of three lists)

#use strict;
use Cwd;
use File::Spec;

# defaults, initializations, and constants
my $help = "\n\nCompareLists_DGE.\n".
            "\t-d  Option: Input directory.  Required.\n".
            "\t-o  Option: Output file. Default=STDOUT.\n".
            "\t-x  Option: Minimum intersect. Default=[Nlist].\n".
            "\n************\nAuthor: J. Cristobal Vera, email: jvera8888\@gmail.com\n"; 
my $usage = "\nCompareLists_DGE.pl -d [Input Dir] -o [Output File] -x [Min Intersect]\n";
my $outfh = 'STDOUT';
my ($indir);
my $fdr_min = 0.01;
my $log2FC_min = 0.5;
my $min_intersect = 0;
my $i = my $n = 0;
my %lists;


#process command line custom script tags
my %cmds = ReturnCmds(@ARGV);
die "\n$help\n$usage\n" if ($cmds{'h'});
if ($cmds{'o'}) {
  $outfh = 'OUT';
  open ($outfh, ">$cmds{'o'}") or die "Cannot open $cmds{'o'}: $!\n";
}
$indir = $cmds{'d'} if ($cmds{'d'});
$min_intersect = $cmds{'x'} if ($cmds{'x'});

##make absolute paths
die "\nError: please specify input directory!\n\n" if (!$indir);
$indir = File::Spec->rel2abs($indir);

##get all
opendir(INDIR,$indir) or die "Can't open directory: $indir: $!\n";
my @inputfiles = grep {m/\.tsv$/ && -f "$indir/$_"} readdir(INDIR);
close(INDIR);
@inputfiles = sort @inputfiles;
print STDERR "\nVCF files found:\n";
foreach my $h (@inputfiles){
  $n++;
  print STDERR "\t$n: $h\n";
}
$n = 0;


##set 'max' min intersect (i.e. size of file list)
if (!$min_intersect){
  $min_intersect = scalar @inputfiles;
}

my @listnames;
foreach my $inputfile (@inputfiles){
  my $listname = $inputfile;
  $listname =~ s/^([^_]+)_([^_]+)_.+/$1\_$2/;
  push @listnames,$listname;
  open (IN, "<$indir/$inputfile") or die "Cannot open $indir/$inputfile: $!\n";
  $n += 1;
  while (my $line = <IN>) {
    $i += 1;
    chomp $line;
    my @line = split /\t/,$line;
    my ($id,$fdr,$log2FC);
    ###get specific DGE tool results (so far: NBID & Seurat)
    ($id,$fdr,$log2FC) = ($line[0],$line[5],$line[2]) if ($listname =~ m/mast/i);
    ($id,$fdr,$log2FC) = ($line[0],$line[1],$line[7]) if ($listname =~ m/nbid/i);
    ($id,$fdr,$log2FC) = ($line[0],$line[3],$line[4]) if ($listname =~ m/combined/i);
    #die "\nError: file $inputfile does not appear to be a DGE file: $indir/$inputfile\n\n" if ($i == 1 and ($id ne 'ID');
    next if ($i == 1);
    next if ($line eq '');
    $lists{$id} .= ";$fdr:$log2FC:$listname" if (exists $lists{$id} and $fdr <= $fdr_min and abs($log2FC) >= $log2FC_min);
    $lists{$id} = "$fdr:$log2FC:$listname" if (!exists $lists{$id} and $fdr <= $fdr_min and abs($log2FC) >= $log2FC_min);
  }
  close(IN);
}
print STDERR "\nTotal DGE lists parsed: $n\n\n";

##print headers
print $outfh "ID";
foreach my $listname (@listnames){
  print $outfh "\tFDR_$listname\tlog2FC_$listname";
}
print $outfh "\n";
##print out IDs that pass min intersect check
my $x = my $y = 0;
foreach my $id (sort keys %lists){
  $y++;
  my @lists = split /;/,$lists{$id};
  if (scalar @lists >= $min_intersect){
    $x++;
    print $outfh "$id";
    foreach my $listname (@listnames){
      my $list = ReturnList($listname,@lists);
      if ($list){
        my ($fdr,$log2fc,$ln) = split /:/,$list;
        print $outfh "\t$fdr\t$log2fc";
      }
      else{
        print $outfh "\t\t";
      }
    }
    print $outfh "\n";
  }
}
print STDERR "\nTotal unique gene IDs: $y\n";
print STDERR "Pass-filter unique gene IDs: $x\n\n";

##return a list from a list of lists based on value match
sub ReturnList{
  my ($value,@lists) = @_;
  foreach my $list (@lists){
    return $list if ($list =~ m/$value/);
  }
  return 0;
}
sub ReturnCmds{
  my (@cmds) = @_;
  my ($opt);
  my %cmds;
  foreach my $cmd (@cmds){
    if (!$opt and $cmd =~ m/^-([a-zA-Z])/) {
      $opt = $1;
    }
    elsif ($opt and $cmd =~ m/^-([a-zA-Z])/){
      $cmds{$opt} = 1;
      $opt = $1;
    }
    elsif ($opt){
      $cmds{$opt} = $cmd;
      $opt = '';
    }
  }
  $cmds{$opt} = 1 if ($opt);
  return %cmds;
}
