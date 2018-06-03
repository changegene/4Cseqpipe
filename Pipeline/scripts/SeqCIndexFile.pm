package SeqCIndexFile;

use strict;

1;
	
sub new($) {
	my($self) = {};
	bless($self, $_[0]);

	my(@expected_header) = ("id", "run", "lane_no", "exp", "primer_seq", "species_name", "first_cutter_name",
							"first_cutter_seq", "sec_cutter_name", "sec_cutter_seq", "linearization_name",
							"linearization_seq", "bait_chromo", "bait_coord", "seq_len", "raw_fname");
	
	$self->{header} = \@expected_header;
	$self->{ids} = [];
	$self->{lines} = {};
	return $self;
}

sub check_format($$$$$$) {
	my($self, $filename, $lineno, $id, $colname, $reg) = @_;
	if (!exists($self->{lines}{"$id\t$colname"}) || $self->{lines}{"$id\t$colname"} !~ $reg) {
		die "index file $filename, line $lineno: invalid format of $colname\n";
	}
}

sub read_file($$) {
	my($self, $filename) = @_;
	open(INDEX, $filename) || die "cannot open index file $filename\n";

	# read the header
	my $line = <INDEX>;
	chop($line);
	my(@header) = split("\t", $line);

	die "invalid number of columns in index file $filename\n" if $#header != $#{$self->{header}};

	for (my($i) = 0; $i <= $#header; $i++) {
		die "invalid column $header[$i] in index file $filename\n" if ($header[$i] ne $self->{header}[$i]);
	}

	# read the rest of the index file
	my($lineno) = 2;
	while (<INDEX>) {
		chop;
		my(@words) = split("\t", $_);
		my($id) = $words[0];
		push(@{$self->{ids}}, $id);
		
		$self->{lines}{$id} = \@words;

		die "index file $filename, line $lineno: invalid number of columns" if ($#words != $#header);
		die "index file $filename: line $lineno: id $id has been already used" if (exists($self->{lines}{"$id\tid"}));
		for (my($i) = 0; $i <= $#words; $i++) {
			$self->{lines}{"$id\t$header[$i]"} = $words[$i];
		}
		
		check_format($self, $filename, $lineno, $id, "id", qr/^\d+$/);
		check_format($self, $filename, $lineno, $id, "lane_no", qr/^\d+$/);
		check_format($self, $filename, $lineno, $id, "primer_seq", qr/^[ACGT]+$/);
		check_format($self, $filename, $lineno, $id, "first_cutter_seq", qr/^[ACGT]+$/);
		check_format($self, $filename, $lineno, $id, "sec_cutter_seq", qr/^([ACGT]+|NA)$/);
		check_format($self, $filename, $lineno, $id, "linearization_seq", qr/^([ACGT]+|NA)$/);
		check_format($self, $filename, $lineno, $id, "bait_coord", qr/^\d+$/);
		check_format($self, $filename, $lineno, $id, "seq_len", qr/^\d+$/);
		
		$lineno++;
	}
	close(INDEX);
}

sub ids($) {
	my($self) = @_;
	return $self->{ids};
}

sub assert_id($$) {
	my($self, $id) = @_;
	die "id $id does not exist in index file\n" unless (exists($self->{lines}{"$id\tid"}));
}

sub run($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\trun"};
}

sub lane_no($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\tlane_no"};
}

sub exp($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\texp"};
}

sub primer_seq($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\tprimer_seq"};
}

sub species_name($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\tspecies_name"};
}

sub first_cutter_name($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\tfirst_cutter_name"};
}

sub first_cutter_seq($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\tfirst_cutter_seq"};
}

sub sec_cutter_name($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\tsec_cutter_name"};
}

sub sec_cutter_seq($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\tsec_cutter_seq"};
}

sub linearization_name($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\tlinearization_name"};
}

sub linearization_seq($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\tlinearization_seq"};
}

sub bait_chromo($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\tbait_chromo"};
}

sub bait_coord($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\tbait_coord"};
}

sub seq_len($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\tseq_len"};
}

sub raw_fname($$) {
	my($self, $id) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\traw_fname"};
}

sub col_names($) {
	my($self) = @_;
	return $self->{header};
}

sub col_value($$$) {
	my($self, $id, $colname) = @_;
	$self->assert_id($id);
	return $self->{lines}{"$id\t$colname"};
}
