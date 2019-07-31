#!/usr/bin/perl -w
=head1 Name

  homogene_identification.pl

=head1 Version

  Author: Wang Lu, double_giraffe@qq.com
  Version:1.0,  Date: 2017-09-26

=head1 Usage

  perl homogene_identification.pl <input.pep> <ref.pep> <outdir>

=cut

use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);

unless(@ARGV==3){&help();};

my $software="$Bin/software.txt";

#software
my %config;
open IN,$software;
while(<IN>){
	chomp;
	if($_=~/^(\S+?)\=(\S+)/){
		$config{$1}=$2;
	}
}
close IN;
my $pfam=$config{"pfam"};
my $blast=$config{"blast"};
my $perl=$config{"perl"};


my $input=$ARGV[0];
my $ref=$ARGV[1];
my $outdir=$ARGV[2];

#get absolute path
$input=abs_path($input);
$ref=abs_path($ref);
$outdir=abs_path($outdir);

print "Referance protein sequences:$ref\nSubject protein sequences:$input\n";

#make blast database
if ( -e "$ref\.phr"){
	print "Blast database is found\n";
}else{
	print "There is not blast database\nMaking blast database\n";
	system("$blast/makeblastdb -in $ref -dbtype prot -out $ref");
}

#blastp
print "Blastp\n";
system("$blast/blastp -query $input -out $outdir/blastp.m8 -db $ref -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs\" -evalue 1e-3 -num_threads 2 ");
my %hash0;
open IN,"$outdir/blastp.m8";
while(<IN>){
	chomp;
	my @a=split/\t/;
	if($a[2]>=20 && $a[12]>=40){
		if(exists $hash0{$a[0]}){
			if($hash0{$a[0]}<=$a[2]){$hash0{$a[0]}=$a[2];}
		}else{$hash0{$a[0]}=$a[2];}
	}
}
close IN;

#max identity
my $max=0;
foreach(keys %hash0){
	if($hash0{$_}>=$max){
		$max=$hash0{$_};
	}
}

#filt blastp
print "Filting blastp\n";
open IN,"$input";
open OUT,">$outdir/blastp.m8.filt.fa";
$/='>';
<IN>;
my %haxi;
while(<IN>){
	chomp;
	my $id=$1 if($_=~/^(\S+)/);
	$_=~s/^.+?\n//;
	$_=~s/\s//g;
	$haxi{$id}=$_;
	if(exists $hash0{$id}){
		print OUT ">$id\n$_\n";
	}
}
close IN;
close OUT;
$/="\n";

#ref pfam
if ( -e "$outdir/ref.pfam"){
	print "File ref.pfam already exists\n";
}else{
	print "Running pfam for referance sequences\n";
	system("$perl $pfam/pfam_scan.pl -fasta $ref -dir $pfam/library -outfile $outdir/ref.pfam -cpu 2");
}


#input pfam
print "Running pfam for Subject sequences\n";
if ( -e "$outdir/blastp.m8.filt.fa.pfam"){
	`rm $outdir/blastp.m8.filt.fa.pfam`;
}
system("$perl $pfam/pfam_scan.pl -fasta $outdir/blastp.m8.filt.fa -dir $pfam/library -outfile $outdir/blastp.m8.filt.fa.pfam -cpu 2");

#filt pfam
print "Filting subject pfam\n";
my %hash;
my $n=`grep '>' $ref | wc -l`;
chomp $n;
print "There are $n referance sequences\n";

open IN,"$outdir/ref.pfam";
while(<IN>){
	next if($_=~/^\#/ || $_=~/^\s/);
	chomp;
	my @a=split /\s+/,$_;
	if(exists $hash{$a[6]}){
		$hash{$a[6]}++;
	}else{
		$hash{$a[6]}=1;
	}
}
close IN;
my @domain;
foreach(keys %hash){
	if($hash{$_}>=$n*0.9){push @domain,$_}
}

my %hash1;
open IN,"$outdir/blastp.m8.filt.fa.pfam";
while(<IN>){
	next if($_=~/^\#/ || $_=~/^\s/);
	chomp;
	my @a=split /\s+/,$_;
	if(exists $hash1{$a[0]}){
		$hash1{$a[0]}="$hash1{$a[0]}\t$a[6]"
	}else{$hash1{$a[0]}=$a[6]}
}
close IN;
my $tmp=0;
open OUT,">$outdir/blastp.m8.filt.fa.pfam.filt";
foreach(keys %hash1){
	my @a=split /\t/,$hash1{$_};
	foreach my $key(@domain){
		for(my $i=0;$i<=$#a;$i++){
			if($a[$i] eq $key){
				$tmp=1;last;
			}
		}
		last;
		$tmp=0
	}
	if($tmp==1 && $hash0{$_}>=0.8*$max){
		print OUT ">$_\n$haxi{$_}\n";
	}
}
close IN;
close OUT;

system("ln -s $outdir/blastp.m8.filt.fa.pfam.filt $outdir/final.fa");
print "The final result is $outdir/final.fa\n";

sub help{
        system("pod2text $0");
	exit();
}

