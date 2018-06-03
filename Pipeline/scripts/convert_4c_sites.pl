use strict;
use warnings;

use File::Path;
use IO::Handle;
use POSIX;

die "Usage: $0 4C_SITES_FILE GROOT_DIR TRACK_SET.TRACK_NAME BIN_SIZE\n" if ($#ARGV != 3);

my $sites_fname = $ARGV[0];
my $groot = $ARGV[1];
my $trackname = $ARGV[2];
my $binsize = $ARGV[3];

# replace "." with "/"
$trackname =~ s/\./\//g;

# READ CHROM FILE
my $chromsizes_fname = $groot . "/chrom_sizes.txt";
open(CHROM_SIZES, $chromsizes_fname) || die "Failed to open chrom sizes file $chromsizes_fname\n";

my %chromsizes;

while (<CHROM_SIZES>) {
	chomp;
	my(@f) = split("\t", $_);
	die "Invalid format of chrom_sizes.txt file $chromsizes_fname" if ($#f != 1);
	$chromsizes{$f[0]} = $f[1];
}

close CHROM_SIZES;

# READ SITES FILE
open(SITES, $sites_fname) || die "Failed to open 4c sites file $sites_fname\n";

#   Read the header
my @second_cutters;

my(@head) = split("\t", <SITES>);

for (my $i = 0; $i <= $#head; $i++) {
	if ($head[$i] =~ /LEN_TO_([ACGT]+)/) {
		push(@second_cutters, $1);
	}
}

for (my $i = 0; $i < $#second_cutters; $i++) {
	print STDOUT $second_cutters[$i] . "\n";
}

#   Open output files for write
umask 02;

my %outfiles;
my $outfile;

my %last_bins;

foreach my $chrom (keys %chromsizes) {
	my $dir = "$groot/tracks/$trackname" . "_fraglen";
	mkdir $dir;
	my $fname = "$dir/chr$chrom";
    $outfile = IO::Handle->new();
	open($outfile, ">$fname") || die "Failed to open $fname for writing.";
	$outfiles{"fraglen\t$chrom"} = $outfile;

	for my $strand (1, -1) {
		for (my $i = 0; $i <= $#second_cutters; $i++) {
			$dir = "$groot/tracks/$trackname" . "_dist2" . $second_cutters[$i] . "_" . ($strand + 4) . "prime";
			mkdir $dir;
			$fname = "$dir/chr$chrom";
			$outfile = IO::Handle->new();
			open($outfile, ">$fname") || die "Failed to open $fname for writing.";
			$outfiles{$second_cutters[$i] . "\t$strand\t$chrom"} = $outfile;
		}
		$last_bins{"$strand\t$chrom"} = -1;
	}
}

#   write bin size
foreach my $key (keys %outfiles) {
	print { $outfiles{$key} } pack('L', $binsize);
}

my $last_chrom = "";

#   Start conversion
while (<SITES>) {
	chomp;
	my(@f) = split("\t", $_);
	die "Invalid format of 4c sites file $sites_fname" if ($#f != 7 + $#second_cutters);

	my $strand = $f[2];
	my $chrom = $f[3];
	my $coord = $f[4];
	my $fraglen = $f[6];

	if ($last_chrom ne $chrom) {
		print STDERR "converting chrom $chrom to binary\n";
		$last_chrom = $chrom;
	}
	
	next if (!exists($chromsizes{$chrom}));
	
	my $bin = int(($coord + 0.5 * $strand * $fraglen) / $binsize);
	my $last_bin = $last_bins{"$strand\t$chrom"};
	
	next if ($last_bin >= $bin);

	# write to fragment length track
	if ($strand == 1) {
		my $outfile = $outfiles{"fraglen\t$chrom"};
			
		# fill the gap between the previous and the current bin with NaNs
		for (my $j = $last_bin + 1; $j < $bin; $j++) {
			print { $outfile } pack('f', "Nan");
		}
		print { $outfile } pack('f', $fraglen);
	}
	
	# write to "length to second cutter" tracks.
	for (my $i = 0; $i <= $#second_cutters; $i++) {
		my $outfile = $outfiles{$second_cutters[$i] . "\t$strand\t$chrom"};
		
		# fill the gap between the previous and the current bin with NaNs
		for (my $j = $last_bin + 1; $j < $bin; $j++) {
			print { $outfile } pack('f', "Nan");
		}
		if ($f[7 + $i] != -1) {
			print { $outfile } pack('f', $f[7 + $i]);
		} else {
			# at the end of the chromosome the distance to the second cutter can be -1
			print { $outfile } pack('f', "Nan");
		}
	}

	$last_bins{"$strand\t$chrom"} = $bin;
}

close SITES;

# if  chromosome is longer than the last bin, add NaNs at the end of the track
foreach my $chrom (keys %chromsizes) {
	my $chromsize = $chromsizes{$chrom};
	my $totbins = ceil($chromsize / $binsize);

	my $last_bin = $last_bins{"1\t$chrom"};
	my $outfile = $outfiles{"fraglen\t$chrom"};

	# write to "fraglen" track
	for (my $j = $last_bin + 1; $j < $totbins; $j++) {
		print { $outfile } pack('f', "Nan");
	}
	close $outfile;

	for my $strand (1, -1) {
		my $last_bin = $last_bins{"$strand\t$chrom"};

		# write to "length to second cutter" tracks
		for (my $i = 0; $i <= $#second_cutters; $i++) {
			my $outfile = $outfiles{$second_cutters[$i] . "\t$strand\t$chrom"};
			for (my $j = $last_bin + 1; $j < $totbins; $j++) {
				print { $outfile } pack('f', "Nan");
			}
			close $outfile;
		}
	}
}
