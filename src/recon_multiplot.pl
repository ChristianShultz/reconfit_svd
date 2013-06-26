#!/usr/bin/perl

use POSIX qw(ceil floor);

die "recon_multiplot.pl <t0> <ops> \n" unless $#ARGV >= 1;

$t0 = shift(@ARGV);
@states = @ARGV;

$random = int(rand(1000));

open(OUT, "> /tmp/multiplot_${random}.gnu");

print OUT "set lmargin 5\nset rmargin 0\n";
print OUT "set multiplot\n";
print OUT "unset key\n\n\n";

$num = $#states + 1;
$dim = ceil( sqrt($num) );

$s = 0.95/$dim;
$size = "$s, $s";

$count = 0;

open(OPS, "< ops_phases");
while($line = <OPS>){
    chomp $line;
    ($n, $op, $a, $b) = split(' ',$line); 
    push(@ops, $op);
}

foreach $state (@states){

    $filename = "t0${t0}/ReconPlots/recon_t0${t0}_${state}_${state}.ax";

    $x = $s*($count % $dim);
    $y = 1.0 -$s - $s * int($count / $dim);
    
    print OUT "set origin $x,$y\n";
    print OUT "set size $size\n";
 
    ($text, $ymin, $ymax) = &make_gnu_from_ax($filename);

    print OUT "set yrange[$ymin:$ymax]\n";
    
#    $xpos = 0.05*$xhigh ;
    $ypos = $ymin + 0.95 * abs($ymax - $ymin);
    $op = $ops[$state];
    print OUT "set label 1 \"C_${state}_${state} : $op\" at 5,$ypos \n";
#    $ypos = $ylow + 0.15 * abs($yhigh - $ylow);
#    print OUT "set label 2 \"m= $Zval +/- $Zerr\" at $xpos,$ypos \n";

    print OUT "$text\n\n";
    $count++;
}

print OUT "set nomultiplot\n";

print OUT "pause -1\n";

close(OUT);

#make the window a sensible size
##$read = `xdpyinfo  | grep \'dimensions:\'`; chomp $read;
##($a, $dims, $a, $a,$a) = split(' ', $read);
##($width, $height) = split('x', $dims);
$width = 900;
$height = 600;

system("gnuplot -geometry ${width}x${height} -persist /tmp/multiplot_${random}.gnu");

system("rm /tmp/multiplot_${random}.gnu");


exit(0);

sub make_gnu_from_ax{
    local($filename) = @_ ;

    open(RECON, "< $filename ");
    $cc = 0;
    while($line = <RECON>){
	chomp $line;
	$test = $line; chop $test; chop $test;
	
	if($test eq '#m'){
	    #print "new entry\n";
	    if($prev eq '#c0'){
		#print "  a fit contribution\n"; 
		$cc++;
	    }
	    else{
		#print "  the data\n";
		last;
	    }
	}
	
	$prev = $line;
    }
    close(RECON);
#    print "counted $count entries\n";
    $contributions = $cc - 6;

    my $text = "plot \'${filename}\' index 0 using 1:2 with lines lt 1 title \"\", \\\n";

    for($i = 2; $i <= $contributions; $i++){
	$index = $i - 1;
	$lt = $i;
	$text = $text . "\'${filename}\' index $index using 1:2 with lines lt $lt title \"\", \\\n";
}

    $line = `grep \"#y\" $filename`; chomp $line;
    ($yyy, $ymin, $ymax) = split(' ',$line);
    
#    $text = $text . "set yrange[${ymin}:${ymax}]\n";

    $c = $contributions;
    $text = $text . "\'${filename}\' index $c using 1:2 with lines linecolor rgb \"black\" lw 1.5 title \"\",\\\n"; $c++;
    $text = $text . "\'${filename}\' index $c using 1:2 with lines linecolor rgb \"black\" lw 1.0 title \"\",\\\n"; $c++;
    $text = $text . "\'${filename}\' index $c using 1:2 with lines linecolor rgb \"black\" lw 1.0 title \"\",\\\n"; $c++;
    
    $text = $text . "\'${filename}\' index $c using 1:2 with lines linecolor rgb \"gray\" title \"\",\\\n"; $c++;
    $text = $text . "\'${filename}\' index $c using 1:2 with lines linecolor rgb \"gray\" title \"\",\\\n"; $c++;
    $text = $text . "\'${filename}\' index $c using 1:2 with lines linecolor rgb \"gray\" title \"\",\\\n"; $c++; 
    
    $text = $text . "\'${filename}\' index $c using 1:2:3 with yerr linecolor rgb \"black\" title \"\"\\\n";
    

    @out = ($text, $ymin, $ymax);
    
    return @out;
}
