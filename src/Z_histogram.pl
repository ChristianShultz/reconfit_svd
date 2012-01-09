#!/usr/bin/perl

#use warnings;

die "Z_histogram.pl <t0> <ops_list_file> <state1> [<state2> ...]" unless $#ARGV > 2;

#ops_list_file is a copy of ops_phases with the desired ops indicated by *
#you can use "all" to show all operators

$t0 = shift(@ARGV);
$ops_list_file = shift(@ARGV);
@states = @ARGV;
$Nstates = $#states + 1;
print "considering $Nstates states - but will normalise using the whole state set\n";

#### GET OP NAMES ####
%opnames = ();
&get_ops;

$dim  = `cat ops_phases | wc -l`;

print "OPERATORS :\n";
$count = 0;
foreach $key (sort {$a <=> $b} keys %opnames){
    print "$key => $opnames{$key}\n";
    $count++;
}
$num_ops = $count;
print "---> ${num_ops} operators in total\n";


## CAPTURE STATES ON THIS t0 ## 
$tmp = `ls t0${t0}/MassJackFiles/ | sed \'s/mass_t0_${t0}_reorder_state//\' | sed \'s/.jack//\' > /tmp/tmp`;
open(TMP, "< /tmp/tmp");
@list;
while($line = <TMP>){
    chomp $line;
    push(@list, $line);
}
system("rm /tmp/tmp");
@list = sort {$a <=> $b} @list;
$tot_state = $#list + 1;

print "ALL STATES : @list \n";


#### READ Z VALUES ####
%allZval = ();
%allZerr = ();
#keyed by operator number

foreach $key (sort {$a <=> $b} keys %opnames){
    @Zval = (); @Zerr = ();
    for($n = 0; $n < $dim; $n++ ){ push(@Zval, 0.0); push(@Zerr, 0.0); }
    &read_Z_op($key);

    # normalised for that op
    $allZval{ $key } = "@Zval";
    $allZerr{ $key } = "@Zerr";
}


#### MAKE THE HISTOGRAM DATA FILE ####
&write_data;

#### MAKE THE GNUPLOT SCRIPT ####
&make_gnu;

exit(0);

######################################################################
sub make_gnu {
    
    open(GNU, "> /tmp/histo_$random.gnu");
    print GNU "set yrange[-0.1:1.2]\n";
    $xmax = $Nstates +0.5;
    print GNU "set xrange[-1:$xmax]\n";
    print GNU "set key noenhanced\n";
    print GNU "set bars fullwidth\n";
    print GNU "set style data histograms\n";

    print GNU "set style fill solid border -1\n";
    print GNU "set style histogram errorbars gap 5 lw 1\n";

    $count = 0;
    foreach $key (sort {$a <=> $b} keys %opnames){
	$style = &ident_op( $opnames{$key} );
	$count++;
	print GNU "set style line $count lc rgb $style \n"; 
    }
    print GNU "set style increment user\n";

    print GNU "set xtics scale 0\n";
   
    print GNU "plot \'/tmp/histo_$random.data\'";
    $count = 0;
    $string = "";
    foreach $key (sort {$a <=> $b} keys %opnames){
	$aa = 2 + 2*$count; $bb = $aa + 1;
	if($count == 0){
	    $string = " using ${aa}:${bb}:xtic(1) ti col, \'\'";
	}
	else{
	    $string = $string . " using ${aa}:${bb} ti col, \'\'";
	}
	$count++;
    }
    chop $string; chop $string; chop $string; chop $string;

    print GNU $string;
    print GNU "\n pause -1\n";

# active this to keep a postscript file
#    print GNU "set term postscript enhanced solid color\n";
#    print GNU "set out \"histo.ps\"\n";
#    print GNU "replot\n";

    close(GNU);

    #make the window a sensible size
    $read = `xdpyinfo  | grep \'dimensions:\'`; chomp $read;
    ($a, $dims, $a, $a,$a) = split(' ', $read);
    ($width, $height) = split('x', $dims);
    $width = int( 0.9* $width);
    $height = int( 0.9* $height);

#    system("gnuplot -geometry 1500x800 -persist /tmp/histo_$random.gnu");
    system("gnuplot -geometry ${width}x${height} -persist /tmp/histo_$random.gnu");

# keep working files for debugging
#    system("cp /tmp/histo_${random}.data ./DATA");
#    system("cp /tmp/histo_${random}.gnu ./GNU");	

    system("rm /tmp/histo_${random}.*");
}


###############################################################################################################

sub ident_op {
    # attempts to color code operators 
    $op = $_[0];

    $hadron = "";
    $isospin = "";
    # add another dimension to these maps for flavour rep (eventually)

    # baryon colors
    %baryon_colors = ('J1o2' => {'conv' => '"black"', 'hyb' => '"#A9A9A9"', 'D2J1A'=> '"#DCDCDC"', 'local'=>'"#708090"', 'lapl'=>'"#F0FFFF"'},
	       'J3o2' => {'conv' => '"red"', 'hyb' => '"#FF1493"', 'D2J1A'=> '"#D2691E"', 'local'=>'"#A52A2A"'},
	       'J5o2' => {'conv' => '"#228B22"', 'hyb' => '"#7FFF00"', 'D2J1A' => '"#ADFF2F"' },
	       'J7o2' => {'conv' => '"blue"', 'hyb' => '"#1E90FF"', 'D2J1A'=>'"#20B2AA"'});

    # meson colors
    %meson_colors = ('J0' => {'conv' => '"black"', 'hyb' => '"#A9A9A9"'},
		     'J1' => {'conv' => '"red"', 'hyb' => '"#FF1493"'},
		     'J2' => {'conv' => '"#228B22"', 'hyb' => '"#7FFF00"'},
		     'J3' => {'conv' => '"blue"', 'hyb' => '"#1E90FF"'},
		     'J4' => {'conv' => '"gold"'},
		     'other' => {'other' => '"brown"'}    );
    
    %colors = (%baryon_colors, %meson_colors);
    

    @isoscalarnames = ('l', 's', 'octet', 'singlet', 'fl', 'fs', 'f1', 'f8', 'etal', 'etas', 'eta1', 'eta8', 'hl', 'hs', 'h1', 'h8', 'omegal', 'omegas', 'omega1', 'omega8');
    
    #baryon syntax   NucleonMG1g1MxD0J0S_J1o2_G1g1
    #meson syntax    b0xD1_J1__J1_T1 (old)   or   rho_rhoxD0_J0__J1_T1

    # need to add support for mesons and baryons in flight ...
    
    @break = split('_', $op);
    $irrep = $break[-1];
    $first = substr($irrep, 0, 1);
    if( ($first eq 'G') || ($first eq 'H') ){
	#its a baryon
	$hadron = "baryon";
	($opname, $spin, $irrep) = split('_', $op);
	($flavspin, $deriv) = split('x', $opname);
    }
    elsif( ($first eq 'A') || ($first eq 'B') || ($first eq 'T') || ($first eq 'E') ){
	#its a meson
	$hadron = "meson";
	($opstruct, $spinirrep) = split('__', $op);
	($spin, $irrep) = split('_', $spinirrep);
	($flavdirac, $deriv) = split('x', $opstruct);

	@flavdirac =  split('_', $flavdirac);

	$isospin = 1;

	if( scalar(@flavdirac) == 1 ){
	    # old op names
	    $dirac = $flavdirac[0];
	}
	else{
	    $flav = $flavdirac[0];
	    $dirac = $flavdirac[1];

	    if (grep {$_ eq $flav} @isoscalarnames) {
		$isospin = 0;
	    }
	}
	    
    }
    else{
	# will probably be called for two-particle opes
	print "don't know how to deal with $op yet! coloring it brown \n";
	$spin = $other;
    }

    if($hadron eq 'baryon'){

	# flavor reps not used yet
	if($flavspin =~ m/Sigma8/ ){$flavrep = 8;}
	elsif($flavspin =~ m/Sigma10/ ){$flavrep = 10;}
	#need extra coloration to show flavor
	#grow the map

	#baryon deriv types
	$type = "conv";
	if( ($deriv eq 'D2J1S') || ($deriv eq 'D2J1M') ){ $type = "hyb"; }
	elsif( $deriv eq 'D2J1A'){ $type = "D2J1A"; }
	elsif( $deriv eq 'D0J0S'){ $type = "local"; }
#    elsif ($deriv eq 'D2J0S'){ $type = "lapl"; }
    }


    elsif($hadron eq 'meson'){
	# flavour reps not used yet
	if($isospin == 1){
	    #don't need coloring by flavor
	}
	else{
	    #need extra coloration to show flavor
	    #grow the map
	}
   
        # meson deriv types
	$type = "conv";
	if( ( $deriv eq 'D2_J1') || ($deriv eq 'D3_J132_J2') ||  ($deriv eq 'D3_J131_J0') 
	    ||  ($deriv eq 'D3_J131_J1') ||  ($deriv eq 'D3_J131_J2') ) { $type = "hyb"; }   

    }
  

    else{
	print "DON'T KNOW WHAT TO DO WITH THIS OPERATOR : ${op} -> brown \n";
	$type = $other;
    }


    $color = $colors{$spin}{$type};
    return $color;
}


#################################################################################

sub get_masses {
    $state = $_[0];
    $line = `calc t0${t0}/MassJackFiles/mass_t0_${t0}_reorder_state${state}.jack`;
    my ($a, $m, $merr, $a) = split(' ',$line);

    $tmp = int(10000*sprintf("%5.4f", $merr));
    $string = sprintf("%5.4f(%d)", $m, $tmp);
    return $string;
}


sub get_ops {
    if($ops_list_file eq 'all'){
	open(FILE, "< ops_phases");
	while($line = <FILE>){
	    chomp $line;
	    ($num, $name, $r, $i) = split(' ',$line);
	    $opnames{$num} = $name;
	} 
    }
    else{
	open(FILE, "< ${ops_list_file}");
	while($line = <FILE>){
	    chomp $line;
	    @entries  = split(' ',$line);
	    if(($#entries == 4) && ($entries[4] eq "*")){$opnames{$entries[0]} = $entries[1];}
	}
    }
}

sub read_Z_op {
    $opnum = $_[0];

    for($n=0; $n < $tot_state; $n++){
	$state = $list[$n];
	$line = `calc t0${t0}/ZJackFiles/Z_t0_${t0}_reorder_state${state}_op${opnum}.jack`;
	my ($a, $z, $zerr, $a) = split(' ',$line); 
	$Zval[$state] = abs($z); 
	$Zerr[$state] = $zerr;
    }

    #use smart max finder that rejects noisy choices
    $max = &find_max_smart;

    for($n=0; $n < $tot_state; $n++){
	$state = $list[$n];	

	$Zval[$state] = sprintf("%5.4f", $Zval[$state] / $max);
	$Zerr[$state] = sprintf("%5.4f", $Zerr[$state] / $max);
    }    
}

sub find_max_smart {
    # rejects noisy max choices
    $locmax = 0;
    %tmp =();

    # will loop over all 'states' including those that aren't there - but they are zero 
    for($n = 0; $n < $dim; $n++){
	$tmp{$Zval[$n]} = $Zerr[$n];
    }
    foreach $key (sort {$b <=> $a} keys %tmp){
	if($tmp{$key} < 0.5*$key){ $locmax = $key; last; }
    }
    return $locmax;
}

sub write_data{
    $random = int(rand(1000));
    open(DATA, "> /tmp/histo_$random.data");

    print DATA "state ";
    foreach $key (sort {$a <=> $b} keys %opnames){
	print DATA "$opnames{$key} $opnames{$key} ";
    }
    print DATA "\n";

    foreach $state (@states){
	$mass = &get_masses($state);
	print DATA "[${state}]${mass}  ";
	foreach $key (sort {$a <=> $b} keys %opnames){
	    @tmp = split(' ', $allZval{ $key }); $z = $tmp[$state];
	    @tmp = split(' ', $allZerr{ $key } ); $zerr = $tmp[$state];
	    print DATA "$z $zerr   ";
	}
	print DATA "\n";
    }
    close(DATA);
}
