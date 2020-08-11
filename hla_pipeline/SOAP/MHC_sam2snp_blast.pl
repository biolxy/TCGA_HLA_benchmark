use strict;
use File::Basename;
use FindBin '$Bin';

=head1 sample

        perl MHC_sam2snp_blast.pl sam_sort_file version(hg19) > reads_snp&indel_file

=head1 done

=cut
die `pod2text $0` if(@ARGV < 1);
############################################################################################################

my$version = $ARGV[1]||"hg19";
my($sample,$gene) = (split /[\._]/,basename$ARGV[0])[0,-2];
my$path = dirname$ARGV[0];
my$pos_mhc = "$Bin/MHC.$version.database";
my$mhc_type = "$Bin/Data_type_$version/$gene.type";
my($pos_start,$pos_end,%data,%reads);

open IN,"< $pos_mhc" or die "$!->$pos_mhc\n";
while(<IN>){
	chomp;
	my($exon,$pos1,$pos2) = (split)[0,2,3];
	if($exon eq $gene){$pos_start=$pos1;$pos_end=$pos2;}
}

open IN,"< $mhc_type" or die "$!->$mhc_type\n";
while(<IN>){
	chomp;my@ln = split;
	my@indel;foreach my$i(3..$#ln){push @indel,$ln[$i];}
	my@snp = split /-/,$ln[2];
	@{$data{$ln[0]}{"snp"}}=@snp;
	@{$data{$ln[0]}{"indel"}}=@indel;
	$data{$ln[0]}{"base"}=$ln[1];
}
##############################################################################################################

open IN,"< $ARGV[0]" or die "$!->$ARGV[0]\n";open OUT,"> $path/$sample\_$gene.fasta";
while(<IN>){
	chomp;my@ln = split;
	my($index,$len) = check($ln[10],"#");

	if($index == 0){$ln[9] = substr($ln[9],$len);}
	if($index > 0){$ln[9] = substr($ln[9],0,$index);}
	if(length$ln[9]<60){next;}

	print OUT ">$ln[0]\n$ln[9]\n";$reads{$ln[0]}=$ln[9];
}

if(($pos_end-$pos_start)>50){
	open IN,"< $path/$sample\_single_unmap.sam";
	while(<IN>){
		chomp;my@ln = split;next if($ln[9]=~/N/);
		my($index,$len) = check($ln[10],"#");

		if($index == 0){$ln[9] = substr($ln[9],$len);}
		if($index > 0){$ln[9] = substr($ln[9],0,$index);}
		if(length$ln[9]<60){next;}

		print OUT ">$ln[0]\n$ln[9]\n";$reads{$ln[0]}=$ln[9];
	}
}

if($gene=~/DRB12/ || $gene=~/DRB13/ || $gene=~/DQB12/ || $gene=~/DQB13/){
	open IN,"< $path/$sample\_pair_unmap.sam";
	while(<IN>){
		chomp;my@ln = split;next if($ln[9]=~/N/);
		my($index,$len) = check($ln[10],"#");

		if($index == 0){$ln[9] = substr($ln[9],$len);}
		if($index > 0){$ln[9] = substr($ln[9],0,$index);}
		if(length$ln[9]<60){next;}

		print OUT ">$ln[0]\n$ln[9]\n";$reads{$ln[0]}=$ln[9];
	}
}
#########################################################################################################

my%hash;my%hash_mat;my%hash_snp;
`$Bin/blastall -i $path/$sample\_$gene.fasta -d $Bin/Data_type_$version/$gene/$gene.fa -o $path/$sample\_$gene.blast -p blastn -e 10 -v 1 -b 1 -F F -m 8`;
open IN,"< $path/$sample\_$gene.blast" or die "$!->blast\n";
while(<IN>){
	chomp;my@ln=split;my$check_mis=0;
	if($ln[5]!=0 || $ln[4]/$ln[3]>0.02 || $ln[4]>1 || $ln[3]<15){next;}
	if($ln[8]>$ln[9]){my$reverse_temp=$ln[8];$ln[8]=$ln[9];$ln[9]=$reverse_temp;}
	if($ln[8]==1 || $ln[9]==length$data{$ln[1]}{"base"} || $ln[3]>=length($reads{$ln[0]})){
		if($ln[4]!=0){$check_mis=1;}

		my($index,$len) = check($data{$ln[1]}{"base"},"*");
		my($pos_rela_1,$pos_rela_2);
		if($index == 0){$pos_rela_1 = $pos_start+$ln[8]-1+$len;}
		else{$pos_rela_1 = $pos_start+$ln[8]-1;}
		
		my(@snp_reads,@indel_reads);
		if(@{$data{$ln[1]}{"indel"}}==0){
			$pos_rela_2 = $pos_rela_1+$ln[3]-1;
		}
		
		if(${$data{$ln[1]}{"indel"}}[0]=~/-D-/){
			my($pos_indel,$base_indel) = split /-D-/,${$data{$ln[1]}{"indel"}}[0];
			if($pos_indel<$pos_rela_1){
				$pos_rela_1=$pos_rela_1+length$base_indel;$pos_rela_2=$pos_rela_1+$ln[3]-1;
			}
			if($pos_indel>=$pos_rela_1 && $pos_indel<($pos_rela_1+$ln[3])){
				$pos_rela_2=$pos_rela_1+$ln[3]-1+length$base_indel;
				push @indel_reads,@{$data{$ln[1]}{"indel"}};
			}
			if($pos_indel>=($pos_rela_1+$ln[3])){
				$pos_rela_2=$pos_rela_1+$ln[3]-1;
			}
		}

		if(${$data{$ln[1]}{"indel"}}[0] =~ /-I-/){
			my($pos_indel,$base_indel) = split /-I-/,${$data{$ln[1]}{"indel"}}[0];
			if(($pos_indel+length$base_indel)<$pos_rela_1){
				$pos_rela_1=$pos_rela_1-length$base_indel;$pos_rela_2=$pos_rela_1+$ln[3]-1;
			}
			if($pos_indel<=$pos_rela_1 && ($pos_indel+length$base_indel)>$pos_rela_1){
				$pos_rela_1=$pos_indel;$pos_rela_2=$pos_rela_1+$ln[3]-1-length$base_indel;
				push @indel_reads,@{$data{$ln[1]}{"indel"}};
			}
			if($pos_indel>$pos_rela_1 && ($pos_indel+length$base_indel)<$pos_rela_1+$ln[3]){
				$pos_rela_2=$pos_rela_1+$ln[3]-1-length$base_indel;
				push @indel_reads,@{$data{$ln[1]}{"indel"}};
			}
			if($pos_indel<$pos_rela_1 && ($pos_indel+length$base_indel)>=($pos_rela_1+$ln[3])){
				$pos_rela_2=$pos_indel+length$base_indel;
				push @indel_reads,@{$data{$ln[1]}{"indel"}};
			}
			if($pos_indel>($pos_rela_1+$ln[3])){
				$pos_rela_2=$pos_rela_1+$ln[3]-1;
			}
			if($pos_rela_1<$pos_indel && ($pos_rela_1+$ln[3])<($pos_indel+length$base_indel)){
				$pos_rela_1=$pos_indel;$pos_rela_2=$pos_indel+1;
				push @indel_reads,@{$data{$ln[1]}{"indel"}};
			}
		}

		foreach my$snp(@{$data{$ln[1]}{"snp"}}){
			my$pos=(split /:/,$snp)[0];
			if($pos>=$pos_rela_1 && $pos<=$pos_rela_2){push @snp_reads,$snp;}
		}

		my@snp_extra;
		if($check_mis==1){
			my$ref=$data{$ln[1]}{"base"};my$reads=$reads{$ln[0]};
			my$ali_ref=substr($ref,$ln[8],($ln[3]-1));
			my$ali_reads=substr($reads,$ln[6],$ln[3]);
			my@ali_ref=split //,$ali_ref;
			my@ali_reads=split //,$ali_reads;
			
			foreach my$n(0..$#ali_ref){
				if($ali_ref[$n] ne $ali_reads[$n]){push @snp_extra,($n+$pos_rela_1+1).":".$ali_reads[$n];}
			}
			if(@snp_extra>1){next;}
			$hash_snp{$snp_extra[0]}++;
		}

		my$snp_reads=join("-",@snp_reads)||"None";
		my$indel_reads=join(":",@indel_reads)||"None";

		push @{$hash{$pos_rela_1."*".$pos_rela_2."*".$snp_reads."*".$indel_reads}},$ln[0];
		if($check_mis==1){
			$hash_mat{$pos_rela_1."*".$pos_rela_2."*".$snp_reads."*".$indel_reads}=$snp_extra[0];
		}
	}
}

foreach my$kes(sort keys %hash){
	my($pos1,$pos2)=split /\*/,$kes;
	if($pos1>$pos2 || ($pos2-$pos1<15)){next;}
	my%tem;@tem{@{$hash{$kes}}}=();@{$hash{$kes}}=sort keys %tem;
	print "$kes\t@{$hash{$kes}}\n";
}

foreach my$kes(sort keys %hash_mat){
	my($pos1,$pos2,$snp,$indel)=split /\*/,$kes;
	if($pos1>$pos2 || ($pos2-$pos1<15)){next;}
	if($hash_snp{$hash_mat{$kes}}<3){next;}
	
	my$pos_novel=(split /:/,$hash_mat{$kes})[0];
	my@snp=split /-/,$snp;my@indel=split /:/,$indel;
	my$near_len=100;my$near_snp;
	foreach my$snp(@snp){
		my($pos,$base)=split /:/,$snp;
		if(abs($pos-$pos_novel) < $near_len){$near_len=abs($pos-$pos_novel);$near_snp=$snp;}
	}
	foreach my$indel(@indel){
		my($pos,$base)=split /-/,$indel;
		if(abs($pos-$pos_novel) < $near_len){$near_len=abs($pos-$pos_novel);$near_snp=$indel;}
	}
	$near_snp||="None";
	if($near_snp eq "None"){next;}

	print "Novel_mutation\t$hash_mat{$kes}\t$hash_snp{$hash_mat{$kes}}\t$near_snp\n";
}
`rm $path/$sample\_$gene.fasta $path/$sample\_$gene.blast`;
#########################################################################################

sub check{
	my$index = index($_[0],$_[1]);
	my@base = split //,$_[0];my$len=0;
	foreach my$i($index..$#base){if($base[$i] eq $_[1]){$len++;}}
	return($index,$len);
}

