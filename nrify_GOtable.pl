#!/usr/bin/perl

$usage= "

nrify_GOtable.pl:

removes duplicate entries for a gene from gene<tab>semicolon-separated GOterms table
concatenates nonredundant categories for each gene

Misha Matz July 2013, matz\@utexas.edu

";

$inp=shift or die $usage;

open IN, $inp or die "cannot open input $inp\n";

my %gos={};
my $gene="";
my $goline="";

while(<IN>){
	chomp;
	($gene,$goline)=split('\t',$_);
	if (!$gos{$gene}) {
		$gos{$gene}=$goline;
		next;
	}
	my @goo=split(';',$goline);
	foreach $g (@goo){
		if ($gos{$gene}=~/$g/){next;}
		$gos{$gene}=$gos{$gene}.";".$g;
	}
}

foreach $g (keys %gos){
	if ($g=~/HASH/){next;}
	print $g,"\t",$gos{$g},"\n";
}