#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;

my $inputfile = '';

my @scriptpath = split /\//, $0;
my $scriptname = pop @scriptpath;
my $scriptdir  = join '/', @scriptpath;

my $GetOpt = GetOptions( 'fasta=s'  => \$inputfile );

my %me2x5 = &makescorematrix("$scriptdir/splicemodels/me2x5");
my %seq   = &makesequencematrix("$scriptdir/splicemodels/splice5sequences");

my %bgd   = (
	'A' => 0.27,
	'C' => 0.23,
	'G' => 0.23,
	'T' => 0.27 );

my $SequenceObj = Bio::SeqIO->new(
	'-file'   => "$inputfile",
	'-format' => 'fasta' );

SCORE: while (my $SiteObj = $SequenceObj->next_seq()) {

		my $str = $SiteObj->seq();
		my $id  = $SiteObj->id();
		$str    = uc($str);

		unless ($str =~ /[ACGT]{9}/) {
			print STDERR "Invalid character in $id $str, skipping!\n";
			next SCORE;
		}

		my $score = '0';
		$score    = &log2(&scoreconsensus($str)*$me2x5{$seq{&getrest($str)}});

		print "$id\t$score\n";
}

sub makesequencematrix {

	my $file = shift;

	my %matrix = ();
	my $n      = 0;

	open(SCOREF, $file) || die "Can't open $file!\n";
	while(<SCOREF>) { 
		chomp;
		$_=~ s/\s//;
		$matrix{$_} = $n;
		$n++;
	}
	close(SCOREF);

	return %matrix;
}

sub makescorematrix {

	my $file = shift;

	my %matrix = ();
	my $n      = 0;

	open(SCOREF, $file) || die "Can't open $file!\n";
	while(<SCOREF>) {
		chomp;
		$_=~ s/\s//;
		$matrix{$n} = $_;
		$n++;
	}
	close(SCOREF);

	return %matrix;
}

sub getrest {

	my $seq = shift;

	my @seqa = split(//,uc($seq));

	return $seqa[0].$seqa[1].$seqa[2].$seqa[5].$seqa[6].$seqa[7].$seqa[8];
}

sub scoreconsensus {

	my $seq  = shift;

	my @seqa = split(//,uc($seq));
	my %bgd  = ();

	$bgd{'A'} = 0.27;
	$bgd{'C'} = 0.23;
	$bgd{'G'} = 0.23;
	$bgd{'T'} = 0.27;

	my %cons1 = ();

	$cons1{'A'} = 0.004;
	$cons1{'C'} = 0.0032;
	$cons1{'G'} = 0.9896;
	$cons1{'T'} = 0.0032;

	my %cons2 = ();

	$cons2{'A'} = 0.0034;
	$cons2{'C'} = 0.0039;
	$cons2{'G'} = 0.0042;
	$cons2{'T'} = 0.9884;

	my $addscore = $cons1{$seqa[3]}*$cons2{$seqa[4]}/($bgd{$seqa[3]}*$bgd{$seqa[4]}); 

	return $addscore;
}

sub log2 {

	my ($val) = @_;
	return log($val)/log(2);
}
