use strict;
use File::Basename;

=head1 sample

	perl MHC_snp2ld_assemble.pl snp_file version(hg19) > ld_file

=head1 done

=cut
die `pod2text $0` if(@ARGV < 1);

my(%hash,%been_done);
open IN,"< $ARGV[0]" or die "$!->$ARGV[0]\n";
while(<IN>){
	chomp;my($snp,@index)=split;
	if(/Novel_mutation/){print "$_\n\n";next}
	else{$hash{$snp}=@index;}
}
######################################################################################################

while(1){
	my($start_hap,%start_hap);my%sub=%hash;
	foreach my$snp(sort keys %sub){
		if(not exists $been_done{$snp}){$start_hap=$snp;$been_done{$start_hap}++;last;}
	}
	unless($start_hap){last;}
	$start_hap{$start_hap}++;

	while(1){
		my@hap_test;
		foreach my$hap(sort keys %start_hap){
			my($hap_pos1,$hap_pos2,$hap_snp,$hap_indel)=split /\*/,$hap;
			$hap_snp=~s/None//g;$hap_indel=~s/None//g;
			my@hap_snp=split /-/,$hap_snp;my@hap_indel=split /:/,$hap_indel;
			my@hap_reads;

			foreach my$sub(sort {$a<=>$b} keys %sub){
				my($sub_pos1,$sub_pos2,$sub_snp,$sub_indel)=split /\*/,$sub;
				$sub_snp=~s/None//g;$sub_indel=~s/None//g;
				my@sub_snp=split /-/,$sub_snp;my@sub_indel=split /:/,$sub_indel;
				
				my%pos;my@pos_com;my@pos_len;
				foreach my$pos($hap_pos1..$hap_pos2){$pos{$pos}++;}
				foreach my$pos($sub_pos1..$sub_pos2){$pos{$pos}++;}
				foreach my$pos(sort {$a<=>$b} keys %pos){if($pos{$pos}==2){push @pos_com,$pos;}}
				@pos_len=sort keys %pos;
				next if(@pos_com<10);
				
				my(@hap_snp_pos,@hap_indel_pos,@sub_snp_pos,@sub_indel_pos);
				foreach my$snp(@hap_snp){
					my$pos=(split /:/,$snp)[0];
					if($pos{$pos}==2){push @hap_snp_pos,$snp;}
				}
				foreach my$indel(@hap_indel){
					my$pos=(split /-/,$indel)[0];
					if($pos{$pos}==2){push @hap_indel_pos,$indel;}
				}
				foreach my$snp(@sub_snp){
					my$pos=(split /:/,$snp)[0];
					if($pos{$pos}==2){push @sub_snp_pos,$snp;}
				}
				foreach my$indel(@sub_indel){
					my$pos=(split /-/,$indel)[0];
					if($pos{$pos}==2){push @sub_indel_pos,$indel;}
				}
	
				my$hap_snp_pos=join("-",@hap_snp_pos)||"None";
				my$hap_indel_pos=join(":",@hap_indel_pos)||"None";
				my$sub_snp_pos=join("-",@sub_snp_pos)||"None";
				my$sub_indel_pos=join(":",@sub_indel_pos)||"None";
				if($hap_snp_pos eq $sub_snp_pos && $hap_indel_pos eq $sub_indel_pos){
					push @hap_reads,$sub;$been_done{$sub}++;
				}
			}
			push @hap_test,@hap_reads;

			my%temp_hap;
			foreach my$hap_reads(@hap_reads){
				my($hap_reads_pos1,$hap_reads_pos2,$hap_reads_snp,$hap_reads_indel)=split /\*/,$hap_reads;
				$hap_reads_snp=~s/None//g;$hap_reads_indel=~s/None//g;
				my@hap_reads_snp=split /-/,$hap_reads_snp;my@hap_reads_indel=split /:/,$hap_reads_indel;
				my$hap_match=0;
	
				foreach my$temp(sort keys %temp_hap){
					my($temp_pos1,$temp_pos2,$temp_snp,$temp_indel)=split /\*/,$temp;
					$temp_snp=~s/None//g;$temp_indel=~s/None//g;
					my@temp_snp=split /-/,$temp_snp;my@temp_indel=split /:/,$temp_indel;
	
					my%pos_temp;
					foreach my$pos($hap_reads_pos1..$hap_reads_pos2){$pos_temp{$pos}++;}
					foreach my$pos($temp_pos1..$temp_pos2){$pos_temp{$pos}++;}
					
					my(@hap_snp_com,@hap_indel_com,@temp_snp_com,@temp_indel_com);
					foreach my$snp(@hap_reads_snp){
						my$pos=(split /:/,$snp)[0];
						if($pos_temp{$pos}==2){push @hap_snp_com,$snp;}
					}
					foreach my$indel(@hap_reads_indel){
						my$pos=(split /-/,$indel)[0];
						if($pos_temp{$pos}==2){push @hap_indel_com,$indel;}
					}
					foreach my$snp(@temp_snp){
						my$pos=(split /:/,$snp)[0];
						if($pos_temp{$pos}==2){push @temp_snp_com,$snp;}
					}
					foreach my$indel(@temp_indel){
						my$pos=(split /-/,$indel)[0];
						if($pos_temp{$pos}==2){push @temp_indel_com,$indel;}
					}
	
					my$hap_snp_com=join("-",@hap_snp_com)||"None";
					my$hap_indel_com=join(":",@hap_indel_com)||"None";
					my$temp_snp_com=join("-",@temp_snp_com)||"None";
					my$temp_indel_com=join(":",@temp_indel_com)||"None";
				
					if($hap_snp_com eq $temp_snp_com && $hap_indel_com eq $temp_indel_com){
						push my@snp_all,@hap_reads_snp,@temp_snp;push my@indel_all,@hap_reads_indel,@temp_indel;
						my%tem;@tem{@snp_all}=();@snp_all=sort keys %tem;
						my%tem;@tem{@indel_all}=();@indel_all=sort keys %tem;

						my@pos_temp=sort keys %pos_temp;
						my$snp_all=join("-",@snp_all)||"None";my$indel_all=join(":",@indel_all)||"None";
						my$new_key=$pos_temp[0]."*".$pos_temp[-1]."*".$snp_all."*".$indel_all;
						$hap_match++;
						if($new_key ne $temp){$temp_hap{$new_key}++;}
					}
				}	
				if($hap_match==0){$temp_hap{$hap_reads}++;}
			}

			foreach my$hap_retain(sort keys %temp_hap){
				my($hap_retain_pos1,$hap_retain_pos2,$hap_retain_snp,$hap_retain_indel)=split /\*/,$hap_retain;
				$hap_retain_snp=~s/None//g;$hap_retain_indel=~s/None//g;
				my@hap_retain_snp=split /-/,$hap_retain_snp;my@hap_retain_indel=split /:/,$hap_retain_indel;
				
				my%pos_retain;my@pos_retain;
				foreach my$pos($hap_pos1..$hap_pos2){$pos_retain{$pos}++;}
				foreach my$pos($hap_retain_pos1..$hap_retain_pos2){$pos_retain{$pos}++;}
				my@pos_retain=sort {$a<=>$b} keys %pos_retain;
	
				push my@hap_snp_total,@hap_snp,@hap_retain_snp;
				push my@hap_indel_total,@hap_indel,@hap_retain_indel;
				my%tem;@tem{@hap_snp_total}=();@hap_snp_total=sort keys %tem;
				my%tem;@tem{@hap_indel_total}=();@hap_indel_total=sort keys %tem;
				my$hap_snp_total=join("-",@hap_snp_total)||"None";
				my$hap_indel_total=join(":",@hap_indel_total)||"None";
	
				my$new_hap_retain=$pos_retain[0]."*".$pos_retain[-1]."*".$hap_snp_total."*".$hap_indel_total;
				$start_hap{$new_hap_retain}++;
			}
		}
		
		foreach my$hap_test(@hap_test){if(exists $sub{$hap_test}){delete($sub{$hap_test});}}
		if(@hap_test==0){last;}
	}

	my@start_hap_done=rmdup(\%start_hap);
	foreach	my$kes(@start_hap_done){
		my@total_reads;my%total_pos;my@total_dep;
		my($kes_pos1,$kes_pos2,$kes_snp,$kes_indel)=split /\*/,$kes;
		my@kes_snp=split /-/,$kes_snp;my@kes_indel=split /:/,$kes_indel;

		foreach my$start_snpt(sort keys %hash){
			my($start_pos1,$start_pos2,$start_snp,$start_indel)=split /\*/,$start_snpt;

			if($start_pos1>=$kes_pos1 && $start_pos2<=$kes_pos2){
				my@kes_snp_in;my@kes_indel_in;
				foreach my$snp(@kes_snp){
					my$pos=(split /:/,$snp)[0];
					if($pos>=$start_pos1 && $pos<=$start_pos2){push @kes_snp_in,$snp;}
				}
				foreach my$indel(@kes_indel){
					my$pos=(split /-/,$indel)[0];
					if($pos>=$start_pos1 && $pos<=$start_pos2){push @kes_indel_in,$indel;}
				}
	
				my$kes_snp_in=join("-",@kes_snp_in)||"None";
				my$kes_indel_in=join(":",@kes_indel_in)||"None";
				if($kes_snp_in eq $start_snp && $kes_indel_in eq $start_indel){
					push @total_reads,$start_snpt."_".$hash{$start_snpt};
					foreach my$pos($start_pos1..$start_pos2){$total_pos{$pos}+=$hash{$start_snpt};}
				}
			}
		}

		next if(@total_reads<=1);
		foreach my$pos(sort keys %total_pos){push @total_dep,$total_pos{$pos};}
		print "$kes\t@total_dep\n@total_reads\n";
	}
}
###########################################################################################################

sub rmdup{
	my@hap_retain_done;my%hash_hap=%{$_[0]};
	while(%hash_hap){
		my@type=sort keys %hash_hap;
		my$type_test=shift@type;delete($hash_hap{$type_test});
		my($type_pos1,$type_pos2,$type_snp,$type_indel)=split /\*/,$type_test;
		my@type_snp=split /-/,$type_snp;my@type_indel=split /:/,$type_indel;
		my$type_count=0;

		foreach my$left_type(sort keys %hash_hap){
			my($left_pos1,$left_pos2,$left_snp,$left_indel)=split /\*/,$left_type;
			my@left_snp=split /-/,$left_snp;my@left_indel=split /:/,$left_indel;

			my@left_pos;my%left_pos;
			foreach my$pos($type_pos1..$type_pos2){$left_pos{$pos}++;}
			foreach my$pos($left_pos1..$left_pos2){$left_pos{$pos}++;}

			my(@type_snp_in,@type_indel_in,@left_snp_in,@left_indel_in);
			foreach my$snp(@type_snp){
				my$pos=(split /:/,$snp)[0];
				if($left_pos{$pos}==2){push @type_snp_in,$snp;}
			}
			foreach my$indel(@type_indel){
				my$pos=(split /-/,$indel)[0];
				if($left_pos{$pos}==2){push @type_indel_in,$indel;}
			}
			foreach my$snp(@left_snp){
				my$pos=(split /:/,$snp)[0];
				if($left_pos{$pos}==2){push @left_snp_in,$snp;}
			}
			foreach my$indel(@left_indel){
				my$pos=(split /-/,$indel)[0];
				if($left_pos{$pos}==2){push @left_indel_in,$indel;}
			}

			my$type_snp_in=join("-",@type_snp_in)||"None";
			my$type_indel_in=join(":",@type_indel_in)||"None";
			my$left_snp_in=join("-",@left_snp_in)||"None";
			my$left_indel_in=join(":",@left_indel_in)||"None";
			
			if($type_snp_in eq $left_snp_in && $type_indel_in eq $left_indel_in){
				if($type_pos1>=$left_pos1 && $type_pos2<=$left_pos2){$type_count++;last;}
			}
		}

		foreach my$left_type(@hap_retain_done){
			my($left_pos1,$left_pos2,$left_snp,$left_indel)=split /\*/,$left_type;
			my@left_snp=split /-/,$left_snp;my@left_indel=split /:/,$left_indel;
			
			my@left_pos;my%left_pos;
			foreach my$pos($type_pos1..$type_pos2){$left_pos{$pos}++;}
			foreach my$pos($left_pos1..$left_pos2){$left_pos{$pos}++;}
			
			my(@type_snp_in,@type_indel_in,@left_snp_in,@left_indel_in);
			foreach my$snp(@type_snp){
				my$pos=(split /:/,$snp)[0];
				if($left_pos{$pos}==2){push @type_snp_in,$snp;}
			}
			foreach my$indel(@type_indel){
				my$pos=(split /-/,$indel)[0];
				if($left_pos{$pos}==2){push @type_indel_in,$indel;}
			}
			foreach my$snp(@left_snp){										
				my$pos=(split /:/,$snp)[0];
				if($left_pos{$pos}==2){push @left_snp_in,$snp;}
			}
			foreach my$indel(@left_indel){
				my$pos=(split /-/,$indel)[0];
				if($left_pos{$pos}==2){push @left_indel_in,$indel;}
			}
			
			my$type_snp_in=join("-",@type_snp_in)||"None";
			my$type_indel_in=join(":",@type_indel_in)||"None";
			my$left_snp_in=join("-",@left_snp_in)||"None";
			my$left_indel_in=join(":",@left_indel_in)||"None";
			
			if($type_snp_in eq $left_snp_in && $type_indel_in eq $left_indel_in){
				if($type_pos1>=$left_pos1 && $type_pos2<=$left_pos2){$type_count++;last;}
			}
		}
		if($type_count==0){push @hap_retain_done,$type_test;}
	}

	return(@hap_retain_done);
}

