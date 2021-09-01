#!/usr/bin/perl
use strict;
use warnings;


my $aa_file = $ARGV[0];
my $ntfile = $ARGV[1];

my $outfile =$ARGV[2];
my $nt_outfile = $ARGV[3];


open (my $ntFH, '<', $ntfile)
or die "Could not open $ntfile: $!\n";	

if ($aa_file =~ m/(.*)\.ali\.sth$/){
open (my $FH, '<', "$aa_file");

#aa output
	my $eog = $1;
	
	#my $outfile = $aa_file;
	#$outfile =~s/.ali.sth/.cut.fas/g;
	
	open (my $outfh, '>>', $outfile)
    or die "Could not write to $outfile: $!\n";

	#my $nt_outfile = $ntfile;
	#$nt_outfile =~ s/.nt_ali.fas/.nt_cut.fas/g;
	
    open (my $nt_outfh, '>>', $nt_outfile)
	or die "Could not write to $nt_outfile: $!\n";

my @names;
my %seq;
my %ntseq;
my $rf = '';

		while (my $line = <$FH>){
			chomp $line;
			if ($line =~ m/^([A-Z]\S+) +(.*)$/){
				my $name = $1;
				if (exists $seq{$name}){
					$seq{$name} .= $2;
				}
				else {
					push @names, $name;
					$seq{$name} = $2;
				}
			}
			if ($line =~ m/^#\=GC RF +(.*)$/){
				$rf .= $1;
			}
		}
		my $start;
		my $end;
		if ($rf =~ m/^(\.*)x.*x(\.*)$/){
			$start = length($1);
			$end = length($2);
		}

		my @ntlines = <$ntFH>;
		for (my $x = 0; $x < @ntlines; ++$x){
			if ($ntlines[$x] =~ m/^>(.*)$/){
				my $taxon = $1;
				$ntseq{$taxon} = '';
				my $y = 1;
				while (1){
					if ($x+$y == @ntlines or $ntlines[$x+$y] =~ m/^>/){
						last;
					}
					else {
						chomp $ntlines[$x+$y];
						$ntseq{$taxon} .= $ntlines[$x+$y];
						++$y;
					}
				}
			}
		}
		foreach my $name (@names){
			my @pos = split //, $seq{$name};
			my $length = @pos;
			my $newseq = '';
			for (my $x = $start; $x < $length - $end; ++$x){
				if ($pos[$x] eq '.'){
					$newseq .= '-';
				}
				else {
					$newseq .= $pos[$x];
				}
			}
			print {$outfh} '>', $name, "\n";
			print {$outfh} $newseq, "\n";
			my @nt_pos = split //, $ntseq{$name};
			my $nt_newseq = '';
			for (my $x = 3 * $start; $x < 3 * ($length - $end); ++$x){
				$nt_newseq .= $nt_pos[$x];
			}
			print {$nt_outfh} '>', $name, "\n";
			print {$nt_outfh} $nt_newseq, "\n";
		}
	}

