#version 1.02

#extract all cutter sites
#count apperances
#report all sites + multi (foreward and reverse)


use strict;
use warnings;

my(%rc) = ("A", "T", "C", "G", "G", "C", "T", "A", "N", "N");

# returns reverse-complement of the argument
sub get_rc {
	my($str, $len) = @_;
	my($rcstr);
	for (my($i) = $len - 1; $i >= 0; $i--) {
		$rcstr .= $rc{substr($str, $i, 1)};
	}
	return($rcstr);
}

die "Usage: $0 FIRST_CUTTER TRACKDB_SEQ_DIR FRAG_LEN OUT_FILE [SEC_CUTTER]*\n" if $#ARGV < 3;

my($site) = $ARGV[0];
my($site_len) = length($site);
my($frag_len) = $ARGV[2];
my($outfile) = $ARGV[3];
my(@second_sites);
my(@second_site_lens);

my(@seq_fns) = <$ARGV[1]/*.seq>;
#my(@seq_fns) = <$ARGV[1]/*.fa>;

# save second cutters in second_sites array
for (my($i) = 4; $i <= $#ARGV; $i++) {
	my($second_site) = $ARGV[$i];
	push(@second_sites, $second_site);
	push(@second_site_lens, length($second_site));
}

my(@hits);

# verify that the site is reverse-symmetric (i.e. the reverse-complement of the site is equal to the site)
$site eq get_rc($site, $site_len) || die "Site must be reverse-symmetric";

print "Extracting $site\n";

my($frag_id) = 0;

my($ser_id) = 0;
my($seq_fn);
foreach $seq_fn (@seq_fns) {
	open(SEQ, $seq_fn) || die "cannot open fasta $seq_fn\n";

	print "will read $seq_fn\n";
	my($chrom) = $seq_fn=~/chr(.+).seq/;
#	my($chrom) = $seq_fn=~/chr(.+).fa/;

	my($seq) = "";
	$seq = uc(<SEQ>);
	close(SEQ);

	print "done reading $seq_fn, got ".length($seq)." bps\n";
	my($max_i) = length($seq)-$site_len;
	
	my($start_coord) = index($seq, $site);
	my($mid_coord) = index($seq, $site, $start_coord + $site_len);
	my($end_coord) = index($seq, $site, $mid_coord + $site_len);

	while ($mid_coord != -1 && $end_coord != -1) {
		# for each site write its distance to the next primary cutter and distances to the next secondary cutter site
		# strand +
		my($frag) = substr($seq, $mid_coord + $site_len, $frag_len);
		if(length($frag) == $frag_len) {
			my($hit) = $ser_id."\t".($frag_id+1)."\t1\t$chrom\t$mid_coord\t$frag";
			$hit .= "\t" . ($end_coord - $mid_coord);  # fragment length
			for (my($i) = 0; $i <= $#second_sites; $i++) {
				my($sec_coord) = index($seq, $second_sites[$i], $mid_coord + $site_len);
				if ($sec_coord == -1) {
					$hit .= "\t-1";
				} else {
					$hit .= "\t" . ($sec_coord - $mid_coord - $site_len + $second_site_lens[$i]); # distance from the end of mid site to the end of next second site
				}
			}
			push(@hits, $hit);
			$ser_id++;
		}
		# same as before but looking for a complementary strand
		# strand -
		my($rcfrag) = get_rc(substr($seq, $mid_coord - $frag_len, $frag_len), $frag_len);
		if(length($rcfrag) == $frag_len) {
			my($hit) = $ser_id."\t$frag_id\t-1\t$chrom\t$mid_coord\t$rcfrag";
			$hit .= "\t" . ($mid_coord - $start_coord);  # fragment length
			for (my($i) = 0; $i <= $#second_sites; $i++) {
				my($sec_coord) = rindex($seq, $second_sites[$i], $mid_coord - $site_len);
				if ($sec_coord == -1) {
					$hit .= "\t-1";
				} else {
					$hit .= "\t" . ($mid_coord - $sec_coord); # distance from the beginning of mid site to the beginning of previous second site
				}
			}
			push(@hits, $hit);
			$ser_id++;
		}

		$start_coord = $mid_coord;
		$mid_coord = $end_coord;
		$end_coord = index($seq, $site, $mid_coord + $site_len);
		$frag_id++;
	}
	$frag_id++;
	
	print "counted $ser_id at chrom $chrom\n";
}

open(TAB, ">$outfile");

# write the header
print TAB "FEND_ID\tFRAG_ID\tSTRAND\tCHROM\tCOORD\tFRAG\tFRAG_LEN";
for (my($i) = 0; $i <= $#second_sites; $i++) {
	print TAB "\tLEN_TO_$second_sites[$i]";
}
print TAB "\n";

# write hits
for(my($i) = 0; $i <= $#hits; $i++) {
	my(@hit) = split("\t", $hits[$i]);
	print TAB "$hit[0]";
	for (my($i) = 1; $i <= $#hit; $i++) {
		print TAB "\t$hit[$i]";
	}
	print TAB "\n";
}

close(TAB);
