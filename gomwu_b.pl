#!/usr/bin/env perl

my $usage= "

gomwu_b.pl  (v. Feb 2015):

This is the second script in the GO database slimming and reformatting procedure,
called automatically by goStats.R 

See README_GO_MWU.txt file for details.

Mikhail Matz, UT Austin; matz@utexas.edu

";

my $gen2go=shift;
my $measure=shift;
my $div=shift or die "$usage\nNot enough arguments for gomwu_b.pl\n";

my $clfile="cl_dissim0_".$div."_".$gen2go; 
open CLF, $clfile or die "cannot locate primary clustering file $clfile\n";
my %clgo={};
my $go;
my $cl;

<CLF>;
while(<CLF>){
	($go, $cl)=split(/,/,$_);
	push @{$clgo{$cl}},$go;
}
close CLF;

unlink $clfile;
unlink "dissim0_".$div."_".$gen2go; 

opendir THISDIR, ".";
my @donealready=grep /$gen2go/, readdir THISDIR;
my $dones=" "."@donealready"." ";
#print "DONE: $dones\n";

my $inname2=$measure.".".$div.".tmp";
open TAB, $inname2 or die "go_nrify: cannot open input table $inname2\n";
<TAB>;
my %level={};
my %desc={};
my %value={};
my $des;
my $go; 
my $l;
my $gn;
my $val;
my @gos=();
my %genes={};
my @gcount=();
my %gci={};
my %goi={};

while (<TAB>){
	chomp;
	($des,$go,$l,$val,$gn)=split(/\t/,$_);
	$value{$gn}=$val;
	$desc{$go}=$des;
	push @{$genes{$go}},$gn;
	push @gcount, $gn unless ($gci{$gn}==1) ;
	$gci{$gn}=1;
#	push @genes,$gn;
	unless ($goi{$go}==1){
		$desc{$go}=$des;
		$level{$go}=$l;
		push @gos, $go;
		$goi{$go}=1;
	}
}
close TAB;
unlink $inname2;

#--------------------


my @nrgos=();
my %nrlev={};
my @nrgos=();
my %nrgenes={};
my %nrdesc={};
my $gcount=0;
my $cl;
my $go;
my $gene;

foreach $cl (keys %clgo){
	my $largest=${$clgo{$cl}}[0];	
	my $maxgenes=$#{$genes{$largest}}+1;
	my $maxlevel=$level{$largest};
	$gcount+=$#{$genes{$largest}}+1;
	my @nrgens=@{$genes{$largest}};
	if ($#{$clgo{$cl}}>0) {
		foreach $go (@{$clgo{$cl}}) {
			if ($maxgenes<($#{$genes{$go}}+1)) {
				$maxgenes=($#{$genes{$go}}+1);
				$largest=$go;
			}
			elsif ($maxgenes==($#{$genes{$go}}+1)) {
				if ($maxlevel<$level{$go}){
					$maxlevel=$level{$go};
					$largest=$go;
				}
			}
			foreach $gene (@{$genes{$go}}){
				push @nrgens, $gene unless(" @nrgens "=~/ $gene /);
			}
		}			
	}
	my $goos=join(";",@{$clgo{$cl}});
	push @nrgos, $goos;
	$nrdesc{$goos}=$desc{$largest};
	$nrlev{$goos}=$maxlevel;
	@{$nrgenes{$goos}}=@nrgens;	
	$gcount+=$#nrgens+1;
}

print $#nrgos+1," non-redundant GO categories of good size\n-------------\n"; 
$outname=$div."_".$measure;
open OUT, ">$outname" or die "gomwu_b: cannot create output $outname\n";
print {OUT} "name\tterm\tlev\tseq\tvalue\n";

foreach $go (@nrgos) {
	foreach $gene (@{$nrgenes{$go}}){
		print {OUT} "$nrdesc{$go}\t$go\t$nrlev{$go}\t$gene\t$value{$gene}\n";
	}
}
close OUT;
my %level={};
my %desc={};
my %value={};
my $des;
my $go; 
my $l;
my $gn;
my $val;
my @gos=();
my %genes={};
my @gcount=();
my %nrlev={};
my @nrgos=();
my %nrgenes={};
my %nrdesc={};
my $gcount=0;
my %dnds;

####################
# building dissimilarity matrix

my $inname4="dissim_".$div."_".$gen2go;
my $inname3=$div."_".$measure;

#if($dones!~/ $inname4 /) { 

	use List::Util qw[min max];
	open TAB, $inname3 or die "go_cluster: cannot open input table $inname3\n";
	<TAB>;

	my $des;
	my $go; 
	my $l;
	my $gn;
	my $val;
	my @gos=();
	my %genes={};
	my %gosi={};

	print"\nSecondary clustering:\ncalculating similarities....\n";
	while (<TAB>){
		chomp;
		($des,$go,$l,$gn,$val)=split(/\t/,$_);
		push @{$genes{$go}},$gn;
		unless ($gosi{$go}==1 ){
			push @gos, $go;
			$gosi{$go}=1;
		}
	}

	my %dissim={};
	for ($g1=0;$g1<$#gos;$g1++){
		my $go=@gos[$g1];
	#if ($go eq "unknown") { print "unknown as go\n";}
		for ($g2=$g1+1;$g2<=$#gos;$g2++){
			my $go2=@gos[$g2];		
			if ($go2 eq "unknown") {
	#print "$go against $go2\n";
				$dissim{$go,$go2}=$dissim{$go2,$go}=1;
				next;
			}
			my %seen={}; my $count=0;
			foreach $g (@{$genes{$go}},@{$genes{$go2}}){
				unless($seen{$g}==1 ){
					$count++;
					$seen{$g}=1;
				}
			}
			my $shared=$#{$genes{$go}}+1+$#{$genes{$go2}}+1-$count;
			my $ref=min($#{$genes{$go}}+1,$#{$genes{$go2}}+1);
			$dissim{$go,$go2}=$dissim{$go2,$go}=sprintf("%.3f",1-$shared/$ref);
		}
	}

	open OUT, ">$inname4" or die "gomwu_b: cannot create output $inname4\n";
	print {OUT} join("\t",@gos),"\n";

	foreach $go (@gos) {
		$dissim{$go,$go}=0;
		foreach $go2 (@gos){
			print {OUT} "$dissim{$go,$go2}";
			print {OUT} "\t" unless ($go2 eq $gos[$#gos]);
		}
		print {OUT} "\n";
	}
#}
