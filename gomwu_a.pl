#!/usr/bin/env perl

my $usage= "

gomwu_a.pl  (v. Feb 2015):

This is the fist script in the GO database slimming and reformatting procedure,
called automatically by goStats.R 

See README_GO_MWU.txt file for details.

Mikhail Matz, UT Austin; matz@utexas.edu

";

print "@ARGV";

my $onto=$ARGV[0] or die $usage;
my $gen2go=$ARGV[1] or die $usage;
my $measure=$ARGV[2] or die $usage;
my $div=$ARGV[3] or die $usage;
my $altern="t";
if ("@ARGV"=~/alternative=(\w)/) { $altern=$1; }
my $toomany=0.1;
if ("@ARGV"=~/largest=(0\.\d+)/) { $toomany=$1; }
my $mingenes=5;
if ("@ARGV"=~/smallest=(\d+)/) { $mingenes=$1; }
my $cutHeight=0.25;
if ("@ARGV"=~/cutHeight=(\S+)/) { $cutHeight=$1; }
my $padj="BH";
if ("@ARGV"=~/p.adjust=(\w+)/) { $padj=$1; }


print "

Run parameters:

largest GO category as fraction of all genes (largest)  : $toomany
         smallest GO category as # of genes (smallest)  : $mingenes
                clustering threshold (clusterCutHeight) : $cutHeight

";

my $division;
if ($div eq "BP") { $division="biological_process";}
elsif ($div eq "MF") { $division="molecular_function";}
elsif ($div eq "CC") { $division="cellular_component";}
else { die "unrecognized division: $div\n";}

my $inname2=$measure.".".$div.".tmp";
my $inname3=$div."_".$measure;
my $inname31="dissim0_".$div."_".$gen2go;
my $inname4="dissim_".$div."_".$gen2go;

my @donealready=();

open GO, $onto or die "cannot open ontology $onto\n";
open CORAL, $gen2go or die "cannot open genes2go table $gen2go\n";
open DNDS, $measure or die "cannot open measure table $measure\n";
my $outname=$inname2;
open VOOL, ">$outname" or die "cannot create output file $outname\n";

#################################################
# reading GO hierarchy, only the $divterm division
print "-----------------\nretrieving GO hierarchy, reformatting data...\n\n";
my %parents;
my %name;
my $term;
my %goodterms;
#my $divterm;
#my $division="molecular_function";

$bads=0;

while (<GO>){
	if ($_=~/^id:\s(GO\S+)/){
		$term=$1;
	}
	elsif ($_=~/^namespace:\s(\S+)/){
			if ($1 ne $division) {
				$term="";
			}
			else {
				$goodterms{$term}="OK" unless (!$term);
			}
	}

	elsif ($_=~/^is_a:\s(GO\S+)\s/) {
		push @{$parents{$term}}, $1 unless (!$term);
		push @allparents, $1;
	}
	elsif ($_=~/^name:\s(.+)/) {
		$name{$term}=$1;
		chomp($name{$term}) unless (!$term);
		if ($name{$term} eq $division) { $divterm=$term;}
	}
}
$goodterms{$term}="OK" unless (!$term);

# calculating hierarchy levels, putting all parents together

my %golevel;
$golevel{$divterm}=0;

foreach $term (keys %parents) {
	$extra=$#{$parents{$term}}+1;
	for ($in=$#{$parents{$term}};$extra;$in+=$extra){
		$golevel{$term}++;
		$tstart=$in-$extra+1;
		$extra=0;
		for ($tt=$tstart; $tt<=$in;$tt++){ 
			$t=${$parents{$term}}[$tt];
			foreach $t0 (@{$parents{$t}}) {
				next if ("@{$parents{$term}}"=~/$t0/ || !$goodterms{$t0});
				push @{$parents{$term}}, $t0;
				$extra++;
			}
		}	
	}
}
###################################################

# reading measures
my $l=0;
my %dnds;
my $seq;

while (<DNDS>){
	if (!$l){
		$l++;
		next;
	}
	chomp;
	($seq,$ns)=split(/,/, $_);
	if ($seq=~/SEQ/) { $seq.="_s";}
	$dnds{$seq}=$ns;
}

#reformatting

$bads=0;
$goods=0;

my $seq;
my $goline;
my @terms;
my @orphans;
my @nolevel;

print {VOOL} "id\tterm\tlev\tvalue\tseq\n";

my $count;
while (<CORAL>){
	chomp;
	($seq,$goline)=split(/\t/,$_);
	if ($dnds{$seq}!~/\d/) {
		push @orphans, $seq;
		next;
	}
	if ($goline=~/unknown/){
		print {VOOL} "\"$goline\"\t$goline\t5\t$dnds{$seq}\t$seq\n";
		next;
	}
	@terms=split(/;/,$goline);
	$nt0=$#terms+1;
	$count++;
#print "-------------------\n$count | $seq terms: $nt0\n";
	my @collect;
	foreach $term (@terms) {
		if (!$goodterms{$term}) {
				$bads++;
#print "$term is not $division$1\n";
				next;
		}
#print "$term\n";
		$goods++;
		if (!$golevel{$term}){
			push @nolevel,$term;
			$golevel{$term}=-1;
		}
		push @collect,($term,@{$parents{$term}});
	}
	
$ncoll=$#collect+1;
#print "collected terms: $ncoll\n";

	my @nrcollect;
	foreach $term (@collect) { 
		push @nrcollect, $term unless ("@nrcollect"=~/$term/);
	}

$ncoll=$#nrcollect+1;
#print "       nr terms: $ncoll\n";
	
	foreach $term (@nrcollect) { 
		print {VOOL} "\"$name{$term}\"\t$term\t$golevel{$term}\t$dnds{$seq}\t$seq\n";
	}
}
#print "terms matching division: $goods\nterms not matching division: $bads\n";
$nor=$#orphans+1;
$nnl=$#nolevel+1;
print "-------------\ngo_reformat:\nGenes with GO annotations, but not listed in measure table: $nor\n";
print "\nTerms without defined level (old ontology?..): $nnl\n-------------\n";
my %parents;
my %name;
my $term;
my %goodterms;
close VOOL;
#########################
# selecting good sized categroies, collapsing redundant ones
#if($dones!~/ $inname3 /) { 
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
#unlink $inname2;

$gc=$#gcount+1;
$goc=$#gos+1;
print "-------------\ngo_nrify:\n$goc categories, $gc genes; size range $mingenes-",$gc*$toomany,"\n";

# excluding too large and  too small categories
my %gonego={};
my $toobroad=0;
my $toosmall=0;
for ($g1=0;$g1<=$#gos;$g1++){
	$go=@gos[$g1];
	my $golen=$#{$genes{$go}}+1;
	if ($golen > $gc*$toomany) { 
		unless ($go eq "unknown") {
			$gonego{$go}=1;
			$toobroad++;
		}
#warn "$go: too broad ($golen genes)\n";
	}
	elsif ($golen < $mingenes) { 
		$gonego{$go}=1 ;
		$toosmall++;
#warn "\t\t$go: too narrow ($golen genes) @{$genes{$go}}\n";
	}
}
print "\t$toobroad too broad\n\t$toosmall too small\n\t", $goc-$toobroad-$toosmall," remaining

removing redundancy:

calculating GO term similarities based on shared genes...\n";

my @goodgo=();
foreach $go (@gos) { push @goodgo, $go unless ($gonego{$go}==1); }
@gos=@goodgo;

################################

	#warn "comparing categories...\n"; 
#my $clfile="cl_".$inname31;
#if($dones!~/ $clfile /) { 

	use List::Util qw[min max];
	for ($g1=0;$g1<=$#gos;$g1++){
		my $go=@gos[$g1];
		next if ($gonego{$go}==1);
	#warn "----------------\n$go  $desc{$go}  level $level{$go}\n";
		my $goos=$go;
		my $lev=$level{$go};
		my $dsc=$desc{$go};
		for ($g2=$g1+1;$g2<=$#gos;$g2++){
			my $go2=@gos[$g2];	
			next if ($gonego{$go2}==1);
			next if ($ggi{$go2}==1);
			my %seen={}; 
			my $count=0;
			my @combo=();
			if ($lump<=1) {
				foreach $g (@{$genes{$go}},@{$genes{$go2}}){
					unless($seen{$g}==1 ){
						$count++;
						$seen{$g}=1;
						push @combo, $g;
					}
				}
				my $shared=$#{$genes{$go}}+1+$#{$genes{$go2}}+1-$count;
				$overlap{$go,$go2}=min($shared/($#{$genes{$go}}+1),$shared/($#{$genes{$go2}}+1));			
				$overlap{$go2,$go}=min($shared/($#{$genes{$go}}+1),$shared/($#{$genes{$go2}}+1));			
			}
		}
	}

	open OUT, ">$inname31" or die "gomwu_b: cannot create output $inname31\n";
	
	print {OUT} join("\t",@gos),"\n";

	foreach $go (@gos) {
		$overlap{$go,$go}=1;
		foreach $go2 (@gos){
			print {OUT} sprintf("%.3f",1-$overlap{$go,$go2});;
			print {OUT} "\t" unless ($go2 eq $gos[$#gos]);
		}
		print {OUT} "\n";
	}
	close OUT;
#}

#print "calling clusteringGOs.R script ....\n";
#	my $err=`Rscript clusteringGOs.R $inname31 $cutHeight `;
#	print $err;





