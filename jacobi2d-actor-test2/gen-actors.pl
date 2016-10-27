#!/usr/bin/perl

use List::Util qw(max min);
use List::MoreUtils 'pairwise';

my $UNROLL = 1;
$UNROLL = $ARGV[0] if ($#ARGV>=0);

my $NDIM = 3;
my $VAR_TYPE = "F";
my $VAR_NAME = "B";
my @DIM_NAMES = ("t", "i", "j");
my @DEPS = ([-1, 0, 0], [-1, -1, 0], [-1, 1, 0], [-1, 0, -1], [-1, 0, 1]);
#currently holes are also specified; holes ignore the t dimension, since it is always assumed to be -1
my @HOLES = ([0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1]);
my @SKEW_SHIFT = (0, -1, -1);
my $ORDER = &order(@DEPS);

printf("---unroll factor & stencil order\n");
printf($UNROLL . ' ' . $ORDER."\n");

my @sorted = &lexSort(@DEPS);
my @skewed_deps;

foreach $d (@sorted) {
	my @s = pairwise {$a + $b} @{$d}, @SKEW_SHIFT;
	push @skewed_deps, \@s;
}

my @DEP_MAX = (-999,-999,-999);
my @DEP_MIN = (999,999,999);
foreach $d (@sorted) {
	@DEP_MAX = pairwise {max($a, $b)} @{$d}, @DEP_MAX;
	@DEP_MIN = pairwise {min($a, $b)} @{$d}, @DEP_MIN;
}

printf("---deps (max/min)---\n");
printf (&dump_dep(\@DEP_MAX)."\n");
printf (&dump_dep(\@DEP_MIN)."\n");

printf("---holes---\n");
printf (&dump_deps(@HOLES));

printf("---deps (lex sorted)---\n");
&dump_deps(@sorted);
printf("---skewed deps\n");
&dump_deps(@skewed_deps);

my @frontier  = &frontier(\@skewed_deps);
printf("---innermost frontier---\n");
&dump_deps(@frontier);


#add deps after unrolling
foreach $dep (@frontier) {
	for (my $i=1; $i<$UNROLL; $i+=1) {
		my @u = pairwise {$a + $b} @{$dep}, @{[0, 0, $i]};
		push @skewed_deps, \@u;
	}
}

@skewed_deps = &lexSort(@skewed_deps);

printf("---unrolled deps (id)---\n");
#&dump_deps(@skewed_deps);

my %dep_to_fifo;
my $id = 0;
foreach $d (@skewed_deps) {
	$dep_to_fifo{$d} = $id;
	$id+=1;
}

foreach $d (@skewed_deps) {
	my $dstr =&dump_dep($d);
	printf('(' . $dep_to_fifo{$d} .') ' . $dstr . "\n");
}

my @frontier_j = &frontier_unrolled(\@skewed_deps, 'j');
my @frontier_i = &frontier_unrolled(\@skewed_deps, 'i');
my @frontier_ij = &frontier_unrolled(\@skewed_deps, 'ij');

printf("---j frontier (after unrolling)---\n");
&dump_deps(@frontier_j);
printf("---i frontier (after unrolling)---\n");
&dump_deps(@frontier_i);
printf("---ij frontier (after unrolling)---\n");
&dump_deps(@frontier_ij);


my @cases = &enumerateCases();


#Projection along t
#  0<=t<HALO_0 and 0<=i<S1+HALO_1 and HALO_2<=j<S2+HALO_2
# modulo along t, full halo along i

#Projection along i
#  HALO_0<=t<S0+HALO_0 and 0<=i<HALO_1 and 0<=j<S2+HALO_2
# modulo along i, full halo along j

#Projection along j
#  0<=t<S0+HALO_0 and HALO_1<=i<S1+HALO_1 and 0<=j<HALO_2
# modulo along j, full halo along t

#tile iterations scanned
#  HALO_0<=t<S0+HALO_0 and HALO_1<=i<S1+HALO_1 and HALO_2<=j<S2+HALO_2
# unroll factor x>1 will visit every x iteration of the j dimension. Assumes S2 is a multiple of x.


printf("\n");

my %FIFOs_used;

foreach $case (@cases) {
	printf(&caseCondition($case). "\n");
	foreach $dep (@skewed_deps) {
		my $readOp = &classifyDep($dep, $case);
		printf($readOp."\n");
		$FIFOs_used{$readOp} = 1 if ($readOp =~ /READ_FIFO/);
	}

}

printf("---FIFOs---\n");
foreach $f (sort keys %FIFOs_used) {
	printf($f."\n");
}


printf("---read actors---\n");

print &readActor('Pt');
print &readActor('Pi');
print &readActor('Pj');
print &readActor('aux');
printf("\n");
print &inputShuffleActor();
printf("\n");


printf("---write actors---\n");
print &writeActor('Pt');
print &writeActor('Pi');
print &writeActor('Pj');
print &writeActor('aux');
printf("\n");
print &outputShuffleActor();
printf("\n");

sub writeActor() {
	my $face = $_[0];

	my $fname = 'store_'.$VAR_NAME.'_'.$face;

	my ($d0, $l0, $lb0, $ub0, $s0) = &writeActorLoop($face, 0);
	my ($d1, $l1, $lb1, $ub1, $s1) = &writeActorLoop($face, 1);
	my ($d2, $l2, $lb2, $ub2, $s2) = &writeActorLoop($face, 2);

	my $i0 = $DIM_NAMES[$d0];
	my $i1 = $DIM_NAMES[$d1];
	my $i2 = $DIM_NAMES[$d2];

	my $UNROLLm1 = $UNROLL-1;

	my $mainLoop;

	#it is simple when j is the innermost
	if ($d2 == 2) {
		$mainLoop = << "EOM";
	pack_t output = outputs.read();
	$l0 {
		$l1 {
			$l2 {
				//make sure numerator >0
				int index = ($i2 + $UNROLL - HALO_$d2) % $UNROLL;
				base[count] = output.vals[index];
				count++;
				if (index == $UNROLLm1) output = outputs.read();
			}
		}
	}
EOM
	#otherwise need to keep track of multiple packs
	} elsif ($d1 == 2) {
		$mainLoop = << "EOM";
	pack_t input[$s2];
	$l0 {
		$l1 {
			$l2 {
				//make sure numerator >0
				int index = ($i1 + $UNROLL - HALO_$d1) % $UNROLL;
				input[$i2-$lb2].vals[index] = base[count];
				count++;
				if (index == $UNROLLm1) inputs.write(input[$i2]); 
			}
		}
	}
EOM
	#when the loop order doesn't follow lex order is tricky
	} elsif ($d0 == 2) {
		$mainLoop = << "EOM";
	int jSize = ($s0+$UNROLL-1) / $UNROLL; //ceil
	pack_t input[$s1][$s2][jSize];
	int jCount = 0;
	$l0 {
		//make sure numerator >0
		int index = ($i0 + $UNROLL - HALO_$d0) % $UNROLL;
		$l1 {
			$l2 {
				input[$i1-$lb1][$i2-$lb2][jCount].vals[index] = base[count];
				count++;
			}
		}
		if (index == $UNROLLm1) jCount++;
	}

	//must load the entire face to populate the fifo in lex order

	for ($DIM_NAMES[0]=HALO_0; $DIM_NAMES[0]<S0+HALO_0; $DIM_NAMES[0]++) {
		for ($DIM_NAMES[1]=0; $DIM_NAMES[1]<HALO_1; $DIM_NAMES[1]++) {
			for ($DIM_NAMES[2]=0; $DIM_NAMES[2]<jSize; $DIM_NAMES[2]++) {
				inputs.write(input[$DIM_NAMES[0]][$DIM_NAMES[1]][$DIM_NAMES[2]]); 
			}
		}
	}
EOM
	}


	return << "EOM";
void $fname($VAR_TYPE *base, fifo_t &outputs) {
	int $DIM_NAMES[0], $DIM_NAMES[1], $DIM_NAMES[2];
	int count = 0;
	

$mainLoop
}
EOM
}

sub writeActorLoop() {
	my $face = $_[0];
	my $dim = $_[1];

	my %scan_region;
	$scan_region{'aux'} = [[0, '0', 'HALO_0'],    [1, '0', 'HALO_1'],         [2, '0', 'HALO_2']];
	$scan_region{'Pt'}  = [[1, 'HALO_1', 'S1+HALO_1'], [2, 'HALO_2', 'S2+HALO_2'], [0, '0', 'HALO_0']];
	$scan_region{'Pi'}  = [[2, 'HALO_2', 'S2+HALO_2'], [0, 'HALO_0', 'S0+HALO_0'], [1, '0', 'HALO_1']];
	$scan_region{'Pj'}  = [[0, 'HALO_0', 'S0+HALO_0'], [1, 'HALO_1', 'S1+HALO_1'], [2, '0', 'HALO_2']];

	my $sr = ${$scan_region{$face}}[$dim];
	my $it = $DIM_NAMES[$$sr[0]];
	my $loop = "for ($it=$$sr[1]; $it<$$sr[2]; $it++)";

	#this is a symbolic simplification that works for the current SR only
	my $size = $$sr[2];
	$size =~ s/\+$$sr[1]//;

	return ($$sr[0], $loop, $$sr[1], $$sr[2], $size);
}



sub contains() {
	my $array_ref = $_[0];
	my $query = $_[1];

	foreach (@{$array_ref}) {
		return 1 if ($query == $_);
	}

	return 0;
}

sub classifyProjection() {
	my $dep_ref = $_[0];
	my $case = $_[1];

	my @dep = pairwise {$a + $b} @{$dep_ref}, @{$case};

	return 'aux' if ($dep[0] < 0 && $dep[1] < 0 && $dep[2] < 0);
	return 'Pt'  if ($dep[0] < 0 && $dep[2] >= 0);
	return 'Pi'  if ($dep[1] < 0 && $dep[0] >= 0);
	return 'Pj'  if ($dep[2] < 0 && $dep[1] >= 0);

	return 'err';
}

sub classifyDep() {
	my $dep_ref = $_[0];
	my $case = $_[1];

	my $depID = $dep_to_fifo{$dep_ref};
	my $reuse = 'REUSE('.&dump_dep($dep_ref, ', ').')';
	my $prev  = 'PREV('.&dump_dep($dep_ref, ', ').')';
	my $fifo  = 'READ_FIFO(' . $depID . ', ' . &classifyProjection($dep_ref, $case) . ')';

	#t == 0
	if ($$case[0] == 0) {
		#very first iteration -- everything from fifo
		if ($$case[1] == 0 && $$case[2] == 0) {
			return $fifo;
		}

		#first column; frontier j are from FIFO
		if ($$case[1] == 0 && $$case[2] > 0) {
			if (&contains(\@frontier_j, $dep_ref)) {
				return $fifo;
			} else {
				return $reuse;
			}
		}

		#second column, first iteration; frontier i are from FIFO
		if ($$case[1] > 0 && $$case[2] == 0) {
			if (&contains(\@frontier_i, $dep_ref)) {
				return $fifo;
			} else {
				return $reuse;
			}
		}

		#steady state for t==0; frontier i-j are from FIFO
		if ($$case[1] > 0 && $$case[2] > 0) {
			if (&contains(\@frontier_ij, $dep_ref)) {
				return $fifo;
			} else {
				return $reuse;
			}
		}
	#t > 0	
	} else {
		#first iteration -- from FIFO if either i or j cross the tile boundary, other wise from PREV
		if ($$case[1] == 0 && $$case[2] == 0) {
			if ($$dep_ref[1] < 0 || $$dep_ref[2] < 0) {
				return $fifo;
			} else {
				return $prev;
			}
		}

		#first column; subset of frontier j crossing the tile boundary are from FIFO 
		if ($$case[1] == 0 && $$case[2] > 0) {
			if (&contains(\@frontier_j, $dep_ref)) {
				if ($$dep_ref[1] < 0) {
					return $fifo;
				} else {
					return $prev;
				}
			} else {
				return $reuse;
			}
		}

		#second column, first iteration; subset of frontier i crossing the tile boundary are from FIFO
		if ($$case[1] > 0 && $$case[2] == 0) {
			if (&contains(\@frontier_i, $dep_ref)) {
				if ($$dep_ref[2] < 0) {
					return $fifo;
				} else {
					return $prev;
				}
			} else {
				return $reuse;
			}
		}

		#steady state for t==0; frontier i-j are from PREV
		if ($$case[1] > 0 && $$case[2] > 0) {
			if (&contains(\@frontier_ij, $dep_ref)) {
				return $prev;
			} else {
				return $reuse;
			}
		}
	}
}



sub enumerateCases() {
	my @cases;
	foreach $t (0,1) {
		foreach $i (0,1) {
			foreach $j (0,1) {
				push @cases, [$t, $i, $j];
	}}}


	return @cases;
}

sub caseCondition() {
	my $case = $_[0];

	my $cond = "";
	foreach $d (0,1,2) {
		$cond .= ' && ' if (length ($cond) > 0);
		$cond .= $DIM_NAMES[$d];
		$cond .= ' == 0' if ($$case[$d] == 0);
		$cond .= ' >= 1' if ($$case[$d] == 1);
	}
	return $cond;
}



sub order() {
	my $order = 0;
	foreach $d (@_) {
		$order = max($order, abs $$d[1], abs $$d[2]);
	}
	return $order;
}

sub frontier_unrolled() {
	my $array_ref = $_[0];
	my $mode = $_[1];


	my @deps;
	my $checkDim;
	if ($mode eq 'j') {
		@deps = reverse &lexSort(@{$array_ref});
		$checkDim = 1;
	} elsif ($mode eq 'i') {
		@deps = reverse &lexSortShuffled($array_ref, [0,2,1]);
		$checkDim = 2;
	} elsif ($mode eq 'ij') {
		@deps = reverse &lexSortShuffled($array_ref, [0,2,1]);
		$checkDim = 2;
	} else {
printf("ERROR\n");
	}
	
	my @frontier;
	my $current = undef;
	my @init;

	my $count = 0;

	foreach $d (@deps) {
		my $sameIasCurrent = $$d[$checkDim] == $$current[$checkDim];
		if (!defined($current) || ($sameIasCurrent && $mode eq 'j' && $count<$UNROLL) || (!$sameIasCurrent)) {
			push @frontier, $d;
			$current = $d;
			$count+=1;
			$count = 1 if (!$sameIasCurrent);
		}
	}

	#post process for i and ij mode
	# removes frontiers convered by j
	if ($mode eq 'i' || $mode eq 'ij') {
		my @temp;
		foreach $dep (@frontier) {
			#those that are not at the i boundary, and is above the updated point in the j dim are covered by the j frontier
			next if ($$dep[1] < 0 && $$dep[2] >= $SKEW_SHIFT[2] && $mode eq 'i');
			#for i-j frontier, the it must be at the i boundary
			next if ($$dep[1] < 0 && $mode eq 'ij');
			push @temp, $dep;
		}
		@frontier = @temp;
	}

	@frontier = reverse @frontier;

	return @frontier;
}

sub frontier() {
	my $array_ref = $_[0];
	my $dim = $_[1];

	my @deps = reverse &lexSort(@{$array_ref});
	
	my @frontier;
	my $current = undef;

	foreach $d (@deps) {
		if (!defined($current) || $$d[1] != $$current[1]) {
			push @frontier, $d;
			$current = $d;
		}
	}

	return reverse @frontier;
}


sub lexSort() {

	return sort { 
		for (my $i =0; $i<$NDIM; $i++) { 
			return $$a[$i] <=> $$b[$i] if ($$a[$i] != $$b[$i]); 
		} return 0;
	} @_;
}

sub lexSortShuffled() {
	my @array = @{$_[0]};
	my @lexShuffle = @{$_[1]};

	return sort { 
		foreach $i (@lexShuffle) {
			return $$a[$i] <=> $$b[$i] if ($$a[$i] != $$b[$i]); 
		} return 0;
	} @array;
}


sub dump_deps() {

	foreach $a (@_) {
		my $s = &dump_dep($a);
		printf($s . "\n");
	}
}

sub dump_dep() {
	my $dep_ref = $_[0];
	my $delim = $_[1];

	$delim = ' ' if ($delim eq '');

	my $res="";
	foreach $e (@{$_[0]}) {
		$res .= $delim if (length($res) > 0);
		$res .= $e;
	}
	return $res;
}

sub readActor() {
	my $face = $_[0];

	my $fname = 'load_'.$VAR_NAME.'_'.$face;

	my ($d0, $l0, $lb0, $ub0, $s0) = &readActorLoop($face, 0);
	my ($d1, $l1, $lb1, $ub1, $s1) = &readActorLoop($face, 1);
	my ($d2, $l2, $lb2, $ub2, $s2) = &readActorLoop($face, 2);

	my $i0 = $DIM_NAMES[$d0];
	my $i1 = $DIM_NAMES[$d1];
	my $i2 = $DIM_NAMES[$d2];

	my @hskips = &readActorHoleSkips($face);

	my $skip;
	foreach $s (@hskips) {$skip .= "if ($s) continue;\n\t\t\t\t";}

	my $UNROLLm1 = $UNROLL-1;


	my $mainLoop;

	#it is simple when j is the innermost
	if ($d2 == 2) {
		$mainLoop = << "EOM";
	pack_t input;
	$l0 {
		$l1 {
			$l2 {
				//make sure numerator >0
				int index = ($i2 + $UNROLL - HALO_$d2) % $UNROLL;
				input.vals[index] = base[count];
				count++;
				if (index == $UNROLLm1) inputs.write(input); 
			}
		}
	}
EOM
	#otherwise need to keep track of multiple packs
	} elsif ($d1 == 2) {
		$mainLoop = << "EOM";
	pack_t input[$s2];
	$l0 {
		$l1 {
			$l2 {
				//make sure numerator >0
				int index = ($i1 + $UNROLL - HALO_$d1) % $UNROLL;
				input[$i2-$lb2].vals[index] = base[count];
				count++;
				if (index == $UNROLLm1) inputs.write(input[$i2]); 
			}
		}
	}
EOM
	#when the loop order doesn't follow lex order is tricky
	} elsif ($d0 == 2) {
		$mainLoop = << "EOM";
	int jSize = ($s0+$UNROLL-1) / $UNROLL; //ceil
	pack_t input[$s1][$s2][jSize];
	int jCount = 0;
	$l0 {
		//make sure numerator >0
		int index = ($i0 + $UNROLL - HALO_$d0) % $UNROLL;
		$l1 {
			$l2 {
				input[$i1-$lb1][$i2-$lb2][jCount].vals[index] = base[count];
				count++;
			}
		}
		if (index == $UNROLLm1) jCount++;
	}

	//must load the entire face to populate the fifo in lex order

	for ($DIM_NAMES[0]=HALO_0; $DIM_NAMES[0]<S0+HALO_0; $DIM_NAMES[0]++) {
		for ($DIM_NAMES[1]=0; $DIM_NAMES[1]<HALO_1; $DIM_NAMES[1]++) {
			for ($DIM_NAMES[2]=0; $DIM_NAMES[2]<jSize; $DIM_NAMES[2]++) {
				inputs.write(input[$DIM_NAMES[0]][$DIM_NAMES[1]][$DIM_NAMES[2]]); 
			}
		}
	}
EOM
	}


	return << "EOM";
void $fname($VAR_TYPE *base, fifo_t &inputs) {
	int $DIM_NAMES[0], $DIM_NAMES[1], $DIM_NAMES[2];
	int count = 0;
	

$mainLoop
}
EOM
}

sub readActorLoop() {
	my $face = $_[0];
	my $dim = $_[1];

	my %scan_region;
	$scan_region{'aux'} = [[0, '0', 'HALO_0'],    [1, '0', 'HALO_1'],         [2, '0', 'HALO_2']];
	$scan_region{'Pt'}  = [[1, '0', 'S1+HALO_1'], [2, 'HALO_2', 'S2+HALO_2'], [0, '0', 'HALO_0']];
	$scan_region{'Pi'}  = [[2, '0', 'S2+HALO_2'], [0, 'HALO_0', 'S0+HALO_0'], [1, '0', 'HALO_1']];
	$scan_region{'Pj'}  = [[0, '0', 'S0+HALO_0'], [1, 'HALO_1', 'S1+HALO_1'], [2, '0', 'HALO_2']];

	my $sr = ${$scan_region{$face}}[$dim];
	my $it = $DIM_NAMES[$$sr[0]];
	my $loop = "for ($it=$$sr[1]; $it<$$sr[2]; $it++)";

	#this is a symbolic simplification that works for the current SR only
	my $size = $$sr[2];
	$size =~ s/\+$$sr[1]//;

	return ($$sr[0], $loop, $$sr[1], $$sr[2], $size);
}

sub readActorHoleSkips() {
	my $face = $_[0];

	my %hole_bound;
	#aux is only at the lex negative corner
	$hole_bound{'aux'} = [[0, $DEP_MIN[1], $DEP_MIN[2]], [0, 0, 0]];
	#Pt covers both corners where j is positive
	$hole_bound{'Pt'}  = [[0, $DEP_MIN[1], 0], [0, $DEP_MAX[2], $DEP_MAX[2]]];
	#Pi covers both corners where i is negative 
	$hole_bound{'Pi'}  = [[0, $DEP_MIN[1], $DEP_MIN[2]], [0, 0, $DEP_MAX[2]]];
	#Pj only covers the i-positive j-negative corner, since it expands along t
	$hole_bound{'Pj'}  = [[0, 0, $DEP_MIN[2]], [0, $DEP_MAX[1], 0]];

	my @min = @{${$hole_bound{$face}}[0]};
	my @max = @{${$hole_bound{$face}}[1]};

	my @conds;

	foreach $hole (@HOLES) {
		if (($min[1] <= $$hole[1] && $$hole[1] <= $max[1]) && 
		    ($min[2] <= $$hole[2] && $$hole[2] <= $max[2])) {
			my $cond1 = "$DIM_NAMES[1] == ";
			my $cond2 = "$DIM_NAMES[2] == ";

			$cond1 .= (abs($DEP_MIN[1]) + $$hole[1]) if ($$hole[1] < 0);
			$cond1 .= 'S1+HALO_1-'.(1-$DEP_MAX[1]+$$hole[1]) if ($$hole[1] > 0);
			$cond2 .= (abs($DEP_MIN[2]) + $$hole[2]) if ($$hole[2] < 0);
			$cond2 .= 'S2+HALO_2-'.(1-$DEP_MAX[2]+$$hole[2]) if ($$hole[2] > 0);

			push @conds, "$cond1 && $cond2";
		}
		
	}

	return @conds;
}


sub inputShuffleActor() {

	my ($i0, $i1, $i2) = @DIM_NAMES;
	
	return << "EOM";
void input_shuffle(fifo_t &inputPt, fifo_t &inputPi, fifo_t &inputPj, fifo_t &inputAux, fifo_t &inputs) {

	int $i0, $i1, $i2;

	int offset = (HALO_2>$UNROLL)?(HALO_2 - $UNROLL):($UNROLL - HALO_2);

	for ($i0=0; $i0<S0+HALO_0; $i0++)
		for ($i1=0; $i1<S1+HALO_1; $i1++)
			for ($i2=offset; $i2<S2+HALO_2; $i2+=$UNROLL) {
				if ($i0 == 0) {
					if ($i1 < HALO_1 && $i2 < HALO_2) {
						inputs.write(inputAux.read());
					}
					if ($i1 >= HALO_1 && $i2 < HALO_2) {
						inputs.write(inputPj.read());
					}
					if ($i2 >= HALO_2) {
						inputs.write(inputPt.read());
					}
				} else {
					if ($i1 < HALO_1) {
						inputs.write(inputPi.read());
					}
					if ($i1 >= HALO_1 && $i2 < HALO_2) {
						inputs.write(inputPj.read());
					}
				}

			}
}

EOM
}


sub outputShuffleActor() {

	my ($i0, $i1, $i2) = @DIM_NAMES;
	
	return << "EOM";
void output_shuffle(fifo_t &outputs, fifo_t &outputPt, fifo_t &outputPi, fifo_t &outputPj, fifo_t &outputAux) {

	int $i0, $i1, $i2;

	int offset = (HALO_2>$UNROLL)?(HALO_2 - $UNROLL):($UNROLL - HALO_2);

	for ($i0=0; $i0<S0; $i0++)
		for ($i1=0; $i1<S1; $i1++)
			for ($i2=0; $i2<S2; $i2+=$UNROLL) {
				bool i1border = ($i1 >= S1-HALO_1);
				bool i2border = ($i2 + $UNROLL-1 >= S2-HALO_2);
				if ($i0 == S0-1) {
					pack_t output = outputs.read();

					outputPt.write(output);

					if (i1border && i2border) {
						outputAux.write(output);
					}
					if (i2border) {
						outputPj.write(output);
					}
				} else {
					if (i1border || i2border) {
						pack_t output = outputs.read();
						if (i1border) {
							outputPi.write(output);
						}
						if (i2border) {
							outputPj.write(output);
						}
					}
				}

			}
}

EOM
}

