#!/bin/usr/perl
use strict;
use warnings;
####$type=1; merege $type=2 getseq
my ($type, @parameter) = @ARGV;

if($type==1){
	my($IN, $transLen, $start, $end, $break)=@parameter;
	$break=int($break);
	if(not defined $IN){
		print STDERR 'please provide the input file\n';
		exit;
	}

	open IN,$transLen;
	my %transLen = ();
	while(<IN>){
		chomp;
		my @line = split /\t/;
		#print $line[0]."\t".$line[1]."\n";
		$transLen{$line[0]}=$line[1];
	}

	open IN,$IN|| die 'can not open the input file\n';
	my %trans=();
	while(<IN>){
		chomp;
		my ($part1,$part2, $trans, $sindex, $eindex) = ();
		my ($position, @line) = split /\t/;
		
		if($position=~/^ENS/){
			($trans, $sindex, $eindex) = split /_/, $position;
		}elsif($position=~/^N/){
			($part1,$part2, $sindex, $eindex) = split /_/, $position;
			$trans = $part1."_".$part2;	
		}
		#print $position."\t".$sindex."\t".$eindex."\n";
		if( int($sindex) < $start ){
			$trans{$trans}{"0\t".($eindex+$end-1)} =1 ;	
				
		}elsif( int($eindex) + $end > $transLen{$trans}){
			$trans{$trans}{($sindex-$start)."\t".($transLen{$trans}-1)}=1;
				
		}else{
			$trans{$trans}{($sindex-$start)."\t".($eindex+$end-1)}=1;	
		}	
	
	}

	open OUT,">window.tmp.bed";
	foreach my $key1 (keys %trans){
		foreach my $key2 (keys %{$trans{$key1}}){
			print OUT $key1."\t".$key2."\n";
		}
	}

	my $tmpbed="./window.tmp.bed";
	my $sortbed="./window.sort.bed";
	my $mergebed="./window.merge.bed";
	my $mergeBed="/software/biosoft/software/MeRIP-PF/tools/BEDTools-Version-2.16.2/bin/mergeBed";
	system("sort -k1,1 -k2,2n $tmpbed > $sortbed" );
	system("$mergeBed -d 0  -i $sortbed > $mergebed");

	print $break;
	open IN,"./window.merge.bed";
	open OUT,">./window.used.bed";
	while(<IN>){
		chomp;
		my @tmp = split /\t/;
		my $peak = int($tmp[2]-$tmp[1]);
		if( $peak < $break){
			print  OUT  $_."\n";
		}else{
			my $mode = $peak%$break;
			my $num = $peak/$break;
			my $intNum = ();
			if($mode == 0){
				$intNum  =  $num;
			}elsif($mode > $break/2){
				$intNum  =  $num+1;
			}elsif($mode <= $break/2){
				$intNum  =  $num;
			}
			#print $intNum."\n";
		
			if($intNum == 1){
				print OUT $_."\n";
			}else{
				foreach my $index ( 1..$intNum){
					if($index == 1){
						print OUT $tmp[0]."\t".$tmp[1]."\t".($tmp[1]+$break)."\n";
					}elsif($index < $intNum && $index>1){
						#print $index;
						print OUT $tmp[0]."\t".($tmp[1]+($index-1)*$break+1)."\t".($tmp[1]+($index)*$break)."\n";
					}elsif($index == $intNum){
						print OUT $tmp[0]."\t".($tmp[1]+($intNum-1)*$break+1)."\t".($tmp[2]-1)."\n";
					}
				}	
			}
	
		}
	}

}else{
	my($input, $tranFa)=@parameter;
	my %trans=();
	open IN,$input;
	while(<IN>){
		chomp;
		my @tmp = split /\t/;
		$trans{$tmp[0]}{$tmp[1]."\t".$tmp[2]}++;
	}
	close IN;
	
	open IN,$tranFa||die "can't open the input file";
	$/=">";
	while(<IN>){
		chomp;
		$_=~/(^\w.*?)\s+/;
		my $tr=$1;
		my @tmp = split /\n/;
		if(exists $trans{$tr}){
			foreach my $index (keys %{$trans{$tr}}){
				my ($start, $end) = split /\t/,  $index; 
				my $seq = substr(join('',@tmp[1..$#tmp]),$start, $end-$start+1);
				my $start2 = $start+1;
				my $end2 = $end+1;			
				$seq=~tr/Tt/Uu/;
				print ">$tr:$start2-$end2\n";
				print $seq."\n";	
			}
		}
	}


}





































































