use strict;
use lib 'scripts';
use SeqCIndexFile;

#---------------------------------
# This is a simple pipeline wraper that allow the user to
# move from raw 4C sequences to normalized nearcis profile and
# domainograms representing them. It provide little flexibility
# and is basically not doing much more than invoking the mapper
# and domainogram programs. 
#--------------------------------

sub get_opt {
	my($nm, $default) = @_;
	if(!exists($::conf{$nm})) {
		if(defined($default)) {
			return($default);
		}
		die "cannot find required config opt $nm\n";
	}
	return($::conf{$nm});
}

sub read_conf {
	open(CONF, "4cseqpipe.conf");
	while(<CONF>) {
		if($_=~/^#/) {
			next;
		}
		my($nm, $val) = /(.+)\s*=\s*(.+)/;
		$::conf{$nm} = $val;
	}
}

sub parse_argv_conf {
	my($params) = @_;

	my($i);

	for($i = 0; $i <= $#$params; $i++) {
		if($params->[$i] !~/^\-/) {
			die "bad parameters format (should be -x x_val -y y_val) at $i\n";
		}
		my($nm) = substr($params->[$i], 1);
		if(!exists($::optnames{$nm})) {
			die "unrecogznied option(s) $nm\n";
		}
		if($nm =~/=/) {
			die "Bad option name $nm - cannot use \"=\"\n";
		}
		if($i == $#$params || $params->[$i+1] =~/^\-/) {
			$::conf{$nm} = 1;
		} else {
			$::conf{$nm} = $params->[$i+1];
			$i++;
		}
	}
}

sub fastq2raw {
	my($infn, $outfn) = @_;

    if ($infn =~ /\.gz$/) {
        open(FASTQ,"gunzip -c $infn |") or die "cannot open pipe to $infn: $!\n";
    } elsif ($infn =~ /\.bz2$/) {
        open(FASTQ,"bunzip2 -c $infn |") or die "cannot open pipe to $infn: $!\n";
    } elsif ($infn =~ /\.zip$/) {
            open(FASTQ,"unzip -p $infn | ") or die "cannot open pipe to $infn: $!\n";
    } else {
        open(FASTQ, $infn) || die "cannot open input fastq file $infn\n";
    }

	open(OUTF, ">$outfn") || die "cannot open output raw file $outfn when converting from fastq\n";

	my($convert_qual) = get_opt("convert_qual", 0);
	if($convert_qual) {
		print STDERR "will convert quality scores\n";
	}
	
	my($read_length) = get_opt("read_length", 36);

	my($cnt) = 0;
	my($l1, $l2, $l3, $l4);
	while(<FASTQ>) {
		$l1 = $_;
		$l2 = <FASTQ>;
		$l3 = <FASTQ>;
		$l4 = <FASTQ>;

		chop $l2;
		chop $l4;

		$l2 = substr($l2, 0, $read_length);
		$l4 = substr($l4, 0, $read_length);

		if($convert_qual) {
			my($max_i) = length($l4);
			for(my($i) = 0; $i < $max_i; $i++) {
				substr($l4, $i, 1) = chr(ord(substr($l4,$i,1))+31);
			}
		}
		if($cnt % 100000 == 0) {
			print STDERR "..$cnt";
		}
		$cnt++;
		print OUTF "0\t$l4\t$l2\n";
	}
	print STDERR "..done\n";
}

sub map_raw {
	my($ids) = @_;

	my($index_fn) = $::conf{"index"};
	my($groot) = get_opt("trackdb_root");
	my($trackset) = get_opt("trackset_name");
	my($rawdir) = get_opt("rawdir", "rawdata");
	my($binsize) = get_opt("binsize");

	my($id);
	foreach $id (@$ids) {
		print STDERR "will invoke the 4C mapper for experiment id $id - this can take a while...\n";
		system("scripts/4cmap.sh -dataidx=$index_fn -raw_dir=$rawdir -id=$id -groot=$groot -trackset=$trackset -binsize=$binsize -statfile=stats/$id.stat");
	}
}

sub nearcis2norm {
	my($ids) = @_;

	my($fromx);
	my($tox);
	my($cfrom);
	my($cto);
	my($horiz5);
	my($horiz3);

	if(exists($::conf{"horiz5"}) && exists($::conf{"horiz3"})) {
		$horiz5 = get_opt("horiz5");
		$horiz3 = get_opt("horiz3");
	}
	if(exists($::conf{"calc_from"}) && exists($::conf{"calc_to"})) {
		$cfrom = get_opt("calc_from");
		$cto = get_opt("calc_to");
	}
	
	my($chrom) = $::idxf->bait_chromo($ids->[0]);
	if(exists($::conf{"horiz5"})) {
		$fromx = $::idxf->bait_coord($ids->[0]) - $horiz5;
	}
	if(exists($::conf{"horiz3"})) {
		$tox = $::idxf->bait_coord($ids->[0]) + $horiz3;
	}
	if(exists($::conf{"calc_from"})) {
		$fromx = $cfrom;
	}
	if(exists($::conf{"calc_to"})) {
		$tox = $cto;
	}

	my($binsize) = get_opt("binsize");
	my($trackset) = get_opt("trackset_name");
	my($nonunique) = get_opt("nonunique", 0);

	my($expname);
	if(exists($::conf{"expname"})) {
		$expname = $::conf{"expname"};
	} else {
		$expname = join("", @$ids);
	}
	
	my($offset) = 4+int($fromx/$binsize)*4;
	my($buffer_size) = ($tox-$fromx)/4+1;

#extract data from the relevant tracks
	my(@combined);
	my($head);
	my($id);
	my($groot) = get_opt("trackdb_root");
	my($re_pref) = "$groot/tracks/re/";
	foreach $id (@$ids) {
		my($first_cut) = $::idxf->first_cutter_seq($id);
		my($second_cut) = $::idxf->sec_cutter_seq($id);
		my($track_pref) = "$groot/tracks/$trackset/$id";
		if($chrom ne $::idxf->bait_chromo($id)
		|| ( ($fromx ne $::idxf->bait_coord($ids->[0]) - $horiz5) && ($fromx ne $cfrom) )) {
			die "Exp ids $ids->[0] and $id refer to different baits, aborting\n";
		}
		open(COV3, $track_pref."_3cov.track/chr$chrom") || die "Cannot open cov3 track\ $track_pref"."_3cov_"."$chrom\n";
		open(COV5, $track_pref."_5cov.track/chr$chrom") || die "Cannot open cov5 track\n";
		open(MULT3, $track_pref."_3multi.track/chr$chrom") || die "Cannot open multi3 track\n";
		open(MULT5, $track_pref."_5multi.track/chr$chrom") || die "Cannot open multi5 track\n";
		open(FL, "$re_pref/$first_cut"."_fraglen/chr$chrom") || die "Cannot open fragment lenth track for first cutter $first_cut in db $re_pref, you may build an re db using the build_re_db option\n";
		open(FE5, "$re_pref/$first_cut"."_dist2$second_cut"."_5prime/chr$chrom") || die "Cannot open fend length track for first cutter $first_cut and second cutter $second_cut, you can build the cutter sites tracks using the build_re_db option\n";
		open(FE3, "$re_pref/$first_cut"."_dist2$second_cut"."_3prime/chr$chrom") || die "Cannot open multi5 track\n";

#		open(COV3, $track_pref."_3cov/chr$chrom") || die "Cannot open cov3 track\ $track_pref"."_3cov_"."$chrom\n";
#		open(COV5, $track_pref."_5cov/chr$chrom") || die "Cannot open cov5 track\n";
#		open(MULT3, $track_pref."_3multi/chr$chrom") || die "Cannot open multi3 track\n";
#		open(MULT5, $track_pref."_5multi/chr$chrom") || die "Cannot open multi5 track\n";
#		open(FL, "$re_pref/$first_cut"."_fraglen/chr$chrom") || die "Cannot open fragment lenth track for first cutter $first_cut in db $re_pref, you may build an re db using the build_re_db option\n";
#		open(FE5, "$re_pref/$first_cut"."_dist2$second_cut"."_5prime/chr$chrom") || die "Cannot open fend length track for first cutter $first_cut and second cutter $second_cut, you can build the cutter sites tracks using the build_re_db option\n";
#		open(FE3, "$re_pref/$first_cut"."_dist2$second_cut"."_3prime/chr$chrom") || die "Cannot open multi5 track\n";

		my($cov3_data, $cov5_data, $mult3_data, $mult5_data, $fl_data, $fe3_data, $fe5_data);
		seek(COV3, $offset, 0); read(COV3, $cov3_data, $buffer_size+4);
		seek(COV5, $offset, 0); read(COV5, $cov5_data, $buffer_size+4);
		seek(MULT3, $offset, 0); read(MULT3, $mult3_data, $buffer_size+4);
		seek(MULT5, $offset, 0); read(MULT5, $mult5_data, $buffer_size+4);
		seek(FL, $offset, 0); read(FL, $fl_data, $buffer_size+4);
		seek(FE3, $offset, 0); read(FE3, $fe5_data, $buffer_size+4);
		seek(FE5, $offset, 0); read(FE5, $fe3_data, $buffer_size+4);

		for(my($bin) = 0; $bin <= $buffer_size/4; $bin++) {
			$combined[$bin] .= 
				"\t".unpack("f", substr($cov3_data, $bin*4, 4)).
				"\t".unpack("f", substr($cov5_data, $bin*4, 4)).
				"\t".unpack("f", substr($mult3_data, $bin*4, 4)).
				"\t".unpack("f", substr($mult5_data, $bin*4, 4)).
				"\t".unpack("f", substr($fl_data, $bin*4, 4)).
				"\t".unpack("f", substr($fe3_data, $bin*4, 4)).
				"\t".unpack("f", substr($fe5_data, $bin*4, 4));
			$combined[$bin] =~s/nan/NA/g;
		}
		my($hnames) = "\teid_3\teid_5\teid_3m\teid_5m\teid_fl\teid_fe3\teid_fe5";
		$hnames =~s/id/$id/g;
		$head .= $hnames;
	}
	open(EXTOUT, ">tables/$expname.nearcis.cover.txt");
	print EXTOUT "intervID\tchrom\tstart\t$head\n";
	for(my($i) = 0; $i <= $#combined; $i++) {
		print EXTOUT "1\t$chrom\t".($fromx+$i*$binsize)."$combined[$i]\n";
	}
	close EXTOUT;

#normalize (blind+replicates)
	my($id_list) = join(",",@$ids);
	print STDERR "Analyzing coordinates $fromx - $tox, ID(s) $id_list\n";
	system("$::Rscript scripts/norm_and_combine_nearcis.r tables/$expname.nearcis.cover.txt 3000 $fromx $tox 5000 $id_list tables/$expname.nearcis.norm.txt $nonunique");
}

sub nearcis2plot {
	my($ids) = @_;

	my($pfrom) = 0;
	my($pto) = 0;
	if(exists($::conf{"plot_from"}) && exists($::conf{"plot_to"})) {
		my($cfrom) = get_opt("calc_from");
		my($cto) = get_opt("calc_to");
		my($pfrom_coord) = get_opt("plot_from");
		my($pto_coord) = get_opt("plot_to");
		$pfrom = $pfrom_coord - $cfrom;
		$pto = $cto - $pto_coord;
	}

	my($expname);
	if(exists($::conf{"expname"})) {
		$expname = $::conf{"expname"};
	} else {
		$expname = join("", @$ids);
	}
	my($feat_tab) = "";
	if(exists($::conf{"feat_tab"})) {
		$feat_tab = $::conf{"feat_tab"};
	}else{
		$feat_tab = "NA";
	}
	my($stat) = get_opt("stat_type");
	my($res) = get_opt("trend_resolution");
	my($png_fn) = "figures/".get_opt("figure_fn");
	my($adjust_color) = 1;
	if(exists($::conf{"adjust_col"})) {
		$adjust_color = $::conf{"adjust_col"};
	}
	my($color_def_script) = get_opt("color_def_script", "NA");
	my($int_type) = get_opt("interval_type", "default");
	my($int_down) = get_opt("interval_down", "default");
	my($int_up) = get_opt("interval_up", "default");
	my($trunc_pct) = get_opt("truncate_percentile", 2);
	my($trunc_max_n) = get_opt("max_truncation", 4);
	my($ymax) = get_opt("max_height", 1);


	print STDERR "Creating near-cis plots for ID(s) $expname\n";	
	system("$::Rscript scripts/create_nearcis_plots.r tables/$expname.nearcis.norm.txt $res $feat_tab $stat $adjust_color $png_fn $color_def_script $pfrom $pto $int_type $int_down $int_up $trunc_pct $trunc_max_n $ymax");
}

#------------
# build the restriction site tracks (fragment lengths, distance to second cutters)

sub build_re_db {
	
	my($re_seq) = get_opt("first_cutter");
	my(@second_cutters) = split(",", get_opt("second_cutters"));
	my($groot) = get_opt("trackdb_root");
	my($binsize) = get_opt("binsize");

	if(!-e "$groot/trakcs/re/") {
		system("mkdir $groot/tracks/re");
	}

	system("perl scripts/extract_restrict_sites.pl $re_seq $groot/seq 36 tmp/re_db.txt ".join(" ", @second_cutters));
	system("perl scripts/convert_4c_sites.pl tmp/re_db.txt $groot re.$re_seq $binsize");
	system("rm tmp/re_db.txt");

}

#------------------------------------------------------------------------
# code start here
#-----------------------------------------------------------------------
#read configuration from defailt file and argv

#my(%::optnames) = (
our(%optnames) = (
"index",1,
"trackdb_root",1,
"trackset_name",1,
"ids",1,
"binsize",1,
"rawdir",1,
"feat_tab",1,
"Rscript_path",1,
"build_re_db",1,
"color_def_script",1,
"convert_qual",1,
"dopipe",1,
"expname",1,
"fastq2raw",1,
"fastq_fn",1,
"figure_fn",1,
"first_cutter",1,
"horiz3",1,
"horiz5",1,
"map",1,
"nearcis",1,
"read_length",1,
"second_cutters",1,
"stat_type",1,
"trend_resolution",1,
"interval_type",1,
"interval_down",1,
"interval_up",1,
"truncate_percentile",1,
"max_truncation",1,
"max_height",1,
"nonunique",1,
"calc_from",1,
"calc_to",1,
"plot_from",1,
"plot_to",1
);

if($#ARGV < 0) {
	print STDERR "usage perl 4cseqpipe.pl [-pipe] [-fastq2raw] [-map] [-nearcis] [-build_re_db] -ids x,1,y,z..\n";
	print STDERR "Default configuration file in 4cseqpipe.conf\n";
	die "See 4cseq_pipe_manual.pdf in the distribution directory\n";
}

read_conf();
parse_argv_conf(\@ARGV);

open(HIST,">>.4cpipe.history");
print HIST "perl 4cseqpipe.pl ".join(" ", @ARGV)."\n";
close HIST;

$::Rscript = get_opt("Rscript_path");

#-----------
#build restriction site tracks if needed

if(get_opt("build_re_db", 0)) {
	build_re_db();
	exit 0;
}

#---------
#Read index table and determine on which exp ids we should work

if(!exists($::conf{"index"})) {
	die "cannot find index name in config file - aborting\n";
}

$::idxf = SeqCIndexFile->new();
$::idxf->read_file($::conf{"index"});

if(!exists($::conf{"ids"})) {
	die "cannot find -ids option - aborting\n";
}
my(@ids) = split(",", $::conf{"ids"});

my($dopipe) = get_opt("dopipe", 0);

#----------
#Convert fastq files to raw format (one line per sequence)

if(get_opt("fastq2raw",0) || ($dopipe && get_opt("fastq_fn", 0))) {
	if(!exists($::conf{"fastq_fn"})) {
		die "fastq2raw was specified, but -fastq_fn was not determined\n";
	}
	my($inp_fn) = get_opt("fastq_fn");
	my($rawdir) = get_opt("rawdir", "rawdata");
	my($id);
	foreach $id (@ids) {
		if($::idxf->raw_fname($id) ne $::idxf->raw_fname($ids[0])) {
			die "cannot run in fastq2raw mode with ids that points to different raw files, first convert the fastq files to raw files (one by one) and then complete the pipeline. Aborting now.\n";
		}	
	}
	fastq2raw($inp_fn, "$rawdir/".($::idxf->raw_fname($ids[0])));
}

#---------------
#map raw files, generate binary tracks

if(get_opt("map",0) || $dopipe) {
	map_raw(\@ids);
}

#--------------
#normalize nearcis and generate domainogram

if(get_opt("nearcis", 0) || $dopipe) {
	nearcis2norm(\@ids);
	nearcis2plot(\@ids);
}

