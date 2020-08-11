use strict;
use File::Basename;
use Cwd 'abs_path';
use FindBin '$Bin';

=head1 sample

	perl MHC_typing_construct.pl indir_txt outdir version 

=head1 done

=cut
die `pod2text $0` if(@ARGV < 1);
###################################################################################################

my$indir = abs_path $ARGV[0];
my$outdir = abs_path $ARGV[1];`mkdir -m 700 -p $outdir` unless(-d "$outdir");
my$version = $ARGV[2]||"hg19";
my$MHC = "$Bin/MHC.$version.database";
my@type = glob "$indir/*txt";
##################################################################################################

foreach my$type(@type){
	my$type_name = (split /\./,basename($type))[0];
	my($position,$symbol,$stand_type);

	`mkdir -p -m 700 $outdir/$type_name` unless(-d "$outdir/$type_name");
	open OUT,"> $outdir/$type_name.type";
	open OUT1,"> $outdir/$type_name/$type_name.fa";
	open MHC,"< $MHC" or die "$!->$MHC\n";

	while(<MHC>){
		chomp;
		my($chr,$pos) = (split /\s+/,$_)[0,2];
		if($chr eq "$type_name"){
			$position = $pos-1;
			$symbol = (split /\s+/,$_)[1];
		}
	}

	if($type_name =~ /^A\d/){$stand_type = "A*03:01:01:01";}
	elsif($type_name =~ /^B\d/){$stand_type = "B*07:02:01";}
	elsif($type_name =~ /^C\d/){$stand_type = "C*07:02:01:03";}
	elsif($type_name =~ /^DQB1\d/){$stand_type = "DQB1*06:02:01";}
	elsif($type_name =~ /^DRB1\d/){$stand_type = "DRB1*15:01:01:01";}
	elsif($type_name =~ /^DMA\d/){$stand_type = "DMA*01:01:01:01";}
	elsif($type_name =~ /^DMB\d/){$stand_type = "DMB*01:03:01:01";}
	elsif($type_name =~ /^DOA\d/){$stand_type = "DOA*01:01:02:01";}
	elsif($type_name =~ /^DOB\d/){$stand_type = "DOB*01:01:01:02";}
	elsif($type_name =~ /^DPA1\d/){$stand_type = "DPA1*01:03:01:01";}
	elsif($type_name =~ /^DPB1\d/){$stand_type = "DPB1*04:01:01:01";}
	elsif($type_name =~ /^DQA1\d/){$stand_type = "DQA1*01:02:01:01";}
	elsif($type_name =~ /^DRA\d/){$stand_type = "DRA*01:02:03";}
	elsif($type_name =~ /^E\d/){$stand_type = "E*01:03:02:01";}
	elsif($type_name =~ /^F\d/){$stand_type = "F*01:03:01:01";}
	elsif($type_name =~ /^G\d/){$stand_type = "G*01:01:01:01";}
	elsif($type_name =~ /^H\d/){$stand_type = "H*02:04";}
	elsif($type_name =~ /^J\d/){$stand_type = "J*01:01:01:01";}
	elsif($type_name =~ /^K\d/){$stand_type = "K*01:01:01:01";}
	elsif($type_name =~ /^L\d/){$stand_type = "L*01:01:01:01";}
	elsif($type_name =~ /^MICA\d/){$stand_type = "MICA*008:04";}
	elsif($type_name =~ /^MICB\d/){$stand_type = "MICB*004:01:01";}
	elsif($type_name =~ /^P\d/){$stand_type = "P*02:01:01:02";}
	elsif($type_name =~ /^TAP1\d/){$stand_type = "TAP1*01:01:01:01";}
	elsif($type_name =~ /^TAP2\d/){$stand_type = "TAP2*01:01:03:01";}
	elsif($type_name =~ /^V\d/){$stand_type = "V*01:01:01:01";}

	my%typing;
	open TXT,"< $type" or die "$!->$type\n";
	while(<TXT>){
		chomp;
		if(/^\s*#/ || /^\s+$/ || /DNA/){next};
		my@temp = split /\s+/,$_;
		my$type_tem = $temp[1];
		for my$kes(2..$#temp){$typing{$type_tem}.=$temp[$kes];}
	}

	foreach my$typing(sort keys %typing){
		if($symbol eq "-"){
			$typing{$typing}=~tr/ATGCatgc/TACGtacg/;
			$typing{$typing}=reverse$typing{$typing};
		}
	}
	my$tep=$typing{$stand_type};

	my@blank;
	while(1){
		unless($tep =~ /\./){last;}
		my$index = index($tep,".");
		if($tep =~ /(\.+)/){$tep = $`.$';}
		push @blank,$index."-".length($1);
	}

	for my$kes(sort keys %typing){
		my$temp = $typing{$kes};
		my(@temp,@Ins,@Del);
		my$temp_pos = $position;

		for my$blank(@blank){
			my($pos,$len) = split /-/,$blank;
			my$check = substr($temp,$pos,$len);
			$temp = substr($temp,0,$pos).substr($temp,($pos+$len));
			$check=~s/\.//g;
			if($check){push @Ins,($temp_pos+$pos)."-I-".$check;}
		}

		while(1){
			unless($temp =~ /\./){last;}
			my$index = index($temp,".");
			if($temp =~ /(\.+)/){
	       		        $temp = substr($temp,0,$index)."N"x(length($1)).substr($temp,($index+length($1)));
	       		}
			my$base = substr($tep,$index,length($1));
			push @Del,($temp_pos+$index)."-D-".$base;
		}

		my@align1 = split //,$tep;my@align2 = split //,$temp;
		for my$num(0..$#align1){
			$temp_pos++;
			if($align1[$num] ne $align2[$num]){
				if($align2[$num] eq "*"){push @temp,"$temp_pos:N";}
				if($align2[$num] =~ /[ATGCatgc]/){push @temp,"$temp_pos:$align2[$num]";}
			}
		}

		my$SNP = join("-",@temp)||"None";
		$typing{$kes} =~ s/\.//g;
		my$len_check = $typing{$kes} =~ s/\*/\*/g;
		if($len_check/length$typing{$kes}>0.8){next;}

		print OUT1 ">$kes\n$typing{$kes}\n";
		print OUT "$kes\t$typing{$kes}\t$SNP\t@Ins\t@Del\n";
	}
	`$Bin/formatdb -i $outdir/$type_name/$type_name.fa -p F -o T`;
}

