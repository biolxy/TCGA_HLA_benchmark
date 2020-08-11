use strict;
use File::Basename;
use FindBin '$Bin';

=head1 sample

        perl MHC_ld2pretype_score.pl ld_file version(hg19) > pretype_file

=head1 done

=cut
die `pod2text $0` if(@ARGV < 1);
#########################################################################################

my$version = $ARGV[1]||"hg19";
my($sample,$gene) = (split /[\._]/,basename$ARGV[0])[0,-2];
my$path = dirname$ARGV[0];
my$pos_mhc = "$Bin/MHC.$version.database";
my$mhc_type = "$Bin/Data_type_$version/$gene.type";
my($pos_start,$pos_end,%data,%hash);

open IN,"< $pos_mhc" or die "$!\n";
while(<IN>){
	chomp;
	my($exon,$pos1,$pos2) = (split)[0,2,3];
	if($exon eq $gene){$pos_start=$pos1;$pos_end=$pos2;}
}

open IN,"< $mhc_type" or die "$!\n";
while(<IN>){
	chomp;my@ln=split;
	my@indel;foreach my$i(3..$#ln){push @indel,$ln[$i];}
	my@snp = split /-/,$ln[2];
	@{$data{$ln[0]}{"snp"}}=@snp;@{$data{$ln[0]}{"indel"}}=@indel;
}
#########################################################################################

my$cluster_core;my$most_reads=0;my%novel;
open IN,"< $ARGV[0]" or die "$!->$ARGV[0]\n";
while(my$ln1=<IN>){
	my$ln2=<IN>;chomp($ln1,$ln2);
	my@ln=split /\s+/,$ln1;my@reads=split /\s+/,$ln2;my@hap=split /\*/,$ln[0];
	if($ln1=~/Novel_mutation/){$novel{$ln[1]}=$ln[2]."\t".$ln[3];next;}

	my%tem;@tem{@reads}=();@reads=keys %tem;
	next if(@ln<11);
	my$match=0;my@type;

	foreach my$kes(sort keys %data){
		my@snp_data;my@indel_data;
		my@snp = @{$data{$kes}{"snp"}};my@indel = @{$data{$kes}{"indel"}};
		foreach my$snp(@snp){
			my$pos=(split /:/,$snp)[0];
			if($pos>=$hap[0] && $pos<=$hap[1]){push @snp_data,$snp;}
		}
		foreach my$indel(@indel){
			my$pos=(split /-/,$indel)[0];
			if($pos>=$hap[0] && $pos<=$hap[1]){push @indel_data,$indel;}
		}
		my$snp_data=join("-",@snp_data)||"None";
		my$indel_data=join(":",@indel_data)||"None";
		if($snp_data eq $hap[2] && $indel_data eq $hap[3]){push @type,$kes;$match++;}
	}

	if($match){
		my@dep;
		my$cov=sprintf("%.2f",($hap[1]-$hap[0]+1)/($pos_end-$pos_start+1));
		my$mean;foreach my$n(1..$#ln){$mean += $ln[$n];}$mean = $mean/(@ln-1);
		my$x2;foreach my$n(1..$#ln){$x2 += ($ln[$n]-$mean)**2/$mean;push @dep,$ln[$n];}
		my$p_value=$x2/log(@dep-1);if($p_value<=60){$p_value=60;}
		my$hap=join("*",@hap);
		my$reads_total=0;
		foreach my$kes(@reads){
			my($reads,$num)=split /_/,$kes;
			$reads_total+=$num;
		}
		my$total_score=($cov**2)*$reads_total/log($p_value);

		$hash{$hap}{"score"}=$total_score;
		@{$hash{$hap}{"reads"}}=@reads;
		@{$hash{$hap}{"type"}}=@type;
		if($total_score>$most_reads){$most_reads=$total_score;$cluster_core=$hap;}
	}
}
unless($cluster_core){die "No reads cover $sample-$gene gene\n";}

my%check_dep;
foreach my$hap(sort keys %hash){
	my$reads_hap=0;
	foreach my$kes(@{$hash{$hap}{"reads"}}){
		my($reads,$num)=split /_/,$kes;
		$reads_hap+=$num;
	}

	my$score=$hash{$hap}{"score"}/$hash{$cluster_core}{"score"};
	my@type=@{$hash{$hap}{"type"}};
	my@reads=@{$hash{$hap}{"reads"}};
	my$mean_dep=$reads_hap*90/($pos_end-$pos_start+1);
#	my$score_trans=log($mean_dep)/log(30);if($score_trans>=1){$score_trans=1;}
#	$score=$score_trans*$score;

	foreach my$kes(keys %novel){
		my($times,$snp)=split /\s+/,$novel{$kes};
		if($hap=~/$snp/){push @{$check_dep{$kes}},$mean_dep;}
	}

	print "$hap\t$score\t@type\n@reads\n";
}

foreach my$kes(sort keys %check_dep){
	my$best_dep=0;
	foreach my$dep(@{$check_dep{$kes}}){if($best_dep<$dep){$best_dep=$dep;}}
	my($times,$snp)=split /\s+/,$novel{$kes};
	if($times/$best_dep>0.8){print "Novel_mutation\t$kes\t$snp\n";}
}

