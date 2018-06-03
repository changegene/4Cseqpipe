use strict;
use POSIX;

sub median {
	my($win, $up, $down, $int_type, $int_down, $int_up) = @_;
	if($int_type eq "default"){
		$int_type = "quantile";	
	}
	if($int_down eq "default"){
		$int_down = 20;	
	}
	if($int_up eq "default"){
		$int_up = 80;	
	}
	if($#$win < 4) {
		return("NA");
	} else {
		my(@sw) = sort { $a <=> $b } @$win;
		if($int_type eq "quantile"){
			$$up = $sw[POSIX::ceil($int_up/100*$#sw)];
			$$down = $sw[POSIX::floor($int_down/100*$#sw)];
		}
		return($sw[int($#sw/2)]);
	}
}

sub lmean {
	my($win, $up, $down, $int_type, $int_down, $int_up) = @_;

	if($int_type eq "default"){
		$int_type = "standard_deviation";	
		if($int_down eq "default"){
			$int_down = 1;	
		}
		if($int_up eq "default"){
			$int_up = 1;	
		}
	}
	if($int_type eq "standard_deviation"){
		if($int_down eq "default"){
			$int_down = 1;	
		}
		if($int_up eq "default"){
			$int_up = 1;	
		}
	}
	if($int_type eq "quantile"){
		if($int_down eq "default"){
			$int_down = 20;	
		}
		if($int_up eq "default"){
			$int_up = 80;	
		}
	}

	my(@sw) = sort { $a <=> $b } @$win;
	if($#$win > 0) {
		my($tot_l) = 0;
		my($tot_l2) = 0;
		my($n) = 0;
		for(my($i) = 0; $i <= $#$win; $i++) {
			my($v) = $win->[$i];
			if($v >= 0) {
				$tot_l += log($v+1);
				$tot_l2 += log($v+1)*log($v+1);
				$n++;
			}
		}
		if($n > 3) {
			my($e) = $tot_l/$n;
			my($v) = $tot_l2/$n - $e*$e;
			if($v <= 0) {
				$$up = "NA";
				$$down = "NA";
			} else {
				my($std) = sqrt($v);
				$$up = $e+$int_up*$std;
				$$down = $e-$int_down*$std;
				$$up = exp($$up)-1; 
				$$down = exp($$down)-1; 
				if($int_type eq "quantile"){
					$$up = $sw[POSIX::ceil($int_up/100*$#sw)];
					$$down = $sw[POSIX::floor($int_down/100*$#sw)];
				}
			}
			return(exp($tot_l/$n)-1);
		} else {
			$$up = "NA";
			$$down = "NA";
			return("NA");
		}
	} else {
		$$up = "NA";
		$$down = "NA";
		return("NA");
	}
}

sub lin_trunc_mean {
	my($win, $win_x, $trunc_p, $trunc_n, $center, $width, $do_linear, $do_log, $up, $down, $int_type, $int_down, $int_up) = @_;

	if($int_type eq "default"){
		$int_type = "standard_deviation";	
		if($int_down eq "default"){
			$int_down = 1;	
		}
		if($int_up eq "default"){
			$int_up = 1;	
		}
	}
	if($int_type eq "standard_deviation"){
		if($int_down eq "default"){
			$int_down = 1;	
		}
		if($int_up eq "default"){
			$int_up = 1;	
		}
	}
	if($int_type eq "quantile"){
		if($int_down eq "default"){
			$int_down = 20;	
		}
		if($int_up eq "default"){
			$int_up = 80;	
		}
	}

	if($#$win < 6) {
		return("NA");
	}
	my(@sw) = sort { $a <=> $b } @$win;

	my($n) = $#sw + 1;

	my($max_i) = POSIX::floor($n*(1-$trunc_p/100));
	my($min_i) = POSIX::ceil($n*$trunc_p/100);

	if($n-$max_i < 2) {
		$max_i = $n-3;
		
	}
	if($n-$max_i > $trunc_n) {
		$max_i = $n-$trunc_n - 1;
	}
	if($min_i < 2) {
		$min_i = 2;
	}
	if($min_i > $trunc_n) {
		$min_i = $trunc_n;
	}
	my($max) = $sw[$max_i];
	my($min) = $sw[$min_i];

#mean of weighted truncated counts
	my($tot_v) = 0;
	my($tot_v2) = 0;
	my($tot_n) = 0;
	for(my($i) = 0; $i <= $#$win; $i++) {
		my($v) = $win->[$i];
		if($v > $max) {
			$v = $max;
		}
		if($v < $min) {
			$v = $min;
		}
		if($do_log) {
			$v = log($v+1);
		}
		if($do_linear == 1) {
			my($d) = 1-(2*abs($win_x->[$i]-$center))/$width;
			$tot_v += $d*$v;
			$tot_v2 += $d*$v*$v;
			$tot_n += $d;
		} else {
			$tot_v += $v;
			$tot_v2 += $v*$v;
			$tot_n++;
		}
	}
	if($tot_n > 4) {
		my($e) = $tot_v/$tot_n;
		my($v) = $tot_v2/$tot_n - $e*$e;
		if($v <= 0) {
			$$up = "NA";
			$$down = "NA";
		} else {
			my($std) = sqrt($v);
			$$up = $e+$int_up*$std;
			$$down = $e-$int_down*$std;
			if($int_type eq "quantile"){
				$$up = $sw[POSIX::ceil($int_up/100*$#sw)];
				$$down = $sw[POSIX::floor($int_down/100*$#sw)];
			}
		}
		if($do_log == 1) {
			$$up = exp($$up)-1;
			$$down = exp($$down)-1;
			if($int_type eq "quantile"){
				$$up = $sw[POSIX::ceil($int_up/100*$#sw)];
				$$down = $sw[POSIX::floor($int_down/100*$#sw)];
			}
			return(exp($e)-1);
		} else {
			if($int_type eq "quantile"){
				$$up = $sw[POSIX::ceil($int_up/100*$#sw)];
				$$down = $sw[POSIX::floor($int_down/100*$#sw)];
			}
			return($tot_v/$tot_n);
		}
	} else {
		$$up = "NA";
		$$down = "NA";
		return("NA");
	}
}


open(DATA, $ARGV[0]) || die "cannot open data $ARGV[0]\n";

my($stat_type) = $ARGV[1];
my(%stat_types) = (
"median", 1,
"mean", 1,
"linear_mean", 1,
"trunc_mean", 1,
"log_mean", 1,
"linear_trunc_mean", 1, 
"trunc_log_mean", 1
);

if(!exists($stat_types{$stat_type})) {
	print STDERR "unrecognized stat type $stat_type, using median statistics as a default\n";
}
my($suf) = $ARGV[2];
my($int_type) = $ARGV[4];
my($int_down) = $ARGV[5];
my($int_up) = $ARGV[6];
my($trunc_pct) = $ARGV[7];
my($trunc_max_n) = $ARGV[8];

open(SCALES, ">$ARGV[0].$suf.scales");
open(NORM_T, ">$ARGV[0].$suf.normfactor_t");
open(NORM_M, ">$ARGV[0].$suf.normfactor_m");

my(@vals);
my(@coords);
my($n_na) = 0;
my($width) = 8;
my($min_x);
my($max_x);
my($j) = 1;

while(<DATA>) {
	chop;
	my(@f) = split("\t", $_);

	if($j==1){
		$min_x = $f[0]-16;
	}
	if(eof){
		$max_x = $f[0]+16;
	}
	$j++;

	for(my($i) = 1; $i <= $#f; $i++) {
		if($f[$i] ne "NA") {
			push(@vals, $f[$i]);
			push(@coords, $f[0]);
			$n_na++;
		}
	}
}
my(@stats);
my(@stats_up);
my(@stats_down);

print STDERR "obtained ". ($#vals+1). " values\n";

my($win_i) = 0;
for(my($win) = 2000; $win <= 50000; $win+=1000) {
	my($base_i) = 0;
	$stats[$win_i] = [];
	$stats_up[$win_i] = [];
	$stats_down[$win_i] = [];
	my($coord_i) = ($win/(1000))/2;
	for(my($base_x) = $min_x; $base_x <= $max_x - $win; $base_x+=1000) {	
		while($coords[$base_i] < $base_x) {
			$base_i++;
			last if !exists($coords[$base_i]);
		}
		last if !exists($coords[$base_i + 1]);
		my($max_i) = $base_i;
		my(@win);
		my(@win_x);
		my($tot_l) = 0;
	
		while($coords[$max_i] < $base_x+$win) {
			my($v) = $vals[$max_i];
			push(@win, $v);
			push(@win_x, $coords[$max_i]);
			$max_i++; 
			last if ($coords[$max_i] == $coords[$#coords]);
		}
		my($stat);
		my($stat_up,$stat_down);
		if($stat_type eq "log_mean") {
			$stat = lmean(\@win, \$stat_up, \$stat_down, $int_type, $int_down, $int_up);
		} elsif($stat_type eq "linear_trunc_mean") {
			$stat = lin_trunc_mean(\@win, \@win_x, $trunc_pct, $trunc_max_n,  $base_x+$win/2, $win, 1, 0, \$stat_up,\$stat_down, $int_type, $int_down, $int_up);
		} elsif($stat_type eq "trunc_mean") {
			$stat = lin_trunc_mean(\@win, \@win_x, $trunc_pct, $trunc_max_n,  $base_x+$win/2, $win, 0, 0, \$stat_up, \$stat_down, $int_type, $int_down, $int_up);
		} elsif($stat_type eq "linear_mean") {
			$stat = lin_trunc_mean(\@win, \@win_x, 0, 0,  $base_x+$win/2, $win, 1, 0, \$stat_up, \$stat_down, $int_type, $int_down, $int_up);
		} elsif($stat_type eq "mean") {
			$stat = lin_trunc_mean(\@win, \@win_x, 0, 0,  $base_x+$win/2, $win, 0, 0, \$stat_up, \$stat_down, $int_type, $int_down, $int_up);
		} elsif($stat_type eq "trunc_log_mean") {
			$stat = lin_trunc_mean(\@win, \@win_x, $trunc_pct, $trunc_max_n,  $base_x+$win/2, $win, 0, 1, \$stat_up, \$stat_down, $int_type, $int_down, $int_up);
		} else {
			$stat = median(\@win, \$stat_up, \$stat_down, $int_type, $int_down, $int_up);
		}
		$stats[$win_i]->[$coord_i] = $stat;
		$stats_up[$win_i]->[$coord_i] = $stat_up;
		$stats_down[$win_i]->[$coord_i] = $stat_down;
		$coord_i++;
	}
	$win_i++;
}

my($max_i) = $#{$stats[0]}+1;

my($max_med_at_targwin_trend) = 0;
my($max_med_at_targwin_multi) = 0;
my($max_at_win_trend) = ($ARGV[3]+1)/3-1;
my($max_at_win_multi) = 10;

for(my($i) = 0; $i <= $max_i; $i++) {
	if(defined($stats[$max_at_win_trend]->[$i]) && $stats[$max_at_win_trend]->[$i] > $max_med_at_targwin_trend) {
		$max_med_at_targwin_trend = $stats[$max_at_win_trend]->[$i];
	}
}

for(my($i) = 0; $i <= $max_i; $i++) {
	if(defined($stats[$max_at_win_multi]->[$i]) && $stats[$max_at_win_multi]->[$i] > $max_med_at_targwin_multi) {
		$max_med_at_targwin_multi = $stats[$max_at_win_multi]->[$i];
	}
}

my($norm_factor_trend) = 1/$max_med_at_targwin_trend;
my($norm_factor_multi) = 1/$max_med_at_targwin_multi;

print NORM_T "$norm_factor_trend\n";
print NORM_M "$norm_factor_multi\n";

my($max_win_i) = $#stats;
for(my($i) = 0; $i <= $max_i; $i++) {
	print SCALES ($min_x+1000*$i);
	for(my($win_i) = 0; $win_i < $max_win_i; $win_i++) {
		if(defined($stats[$win_i]->[$i])) {
			print SCALES "\t".($norm_factor_trend*($stats[$win_i]->[$i]));
			print SCALES "\t".($norm_factor_trend*($stats_up[$win_i]->[$i]));
			print SCALES "\t".($norm_factor_trend*($stats_down[$win_i]->[$i]));
		} else {
			print SCALES "\tNA\tNA\tNA";
		}
	}
	print SCALES "\n";
}
