#!/usr/bin/perl

use POSIX qw(ceil floor);

die "multiplot.pl <states> \n" unless $#ARGV >= 0;

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

foreach $state (@states){

    $filename = "multi_t0_fits/PrinCorrPlots/multi_t0_prin_corr_fit_state${state}.ax";

    $x = $s*($count % $dim);
    $y = 1.0 -$s - $s * int($count / $dim);
    
    print OUT "set origin $x,$y\n";
    print OUT "set size $size\n";
 
    ($text, $xlow, $xhigh, $ylow, $yhigh, $chisq, $Zval, $Zerr) = &make_gnu_from_ax($filename);

    print OUT "set yrange[$ylow:$yhigh]\n";
    print OUT "set xrange[$xlow:$xhigh]\n";
    
    $xpos = 0.05*$xhigh ;
    $ypos = $ylow + 0.30 * abs($yhigh - $ylow);
    print OUT "set label 1 \"state$state, chisq=$chisq\" at $xpos,$ypos \n";
    $ypos = $ylow + 0.15 * abs($yhigh - $ylow);
    print OUT "set label 2 \"m= $Zval +/- $Zerr\" at $xpos,$ypos \n";

    print OUT "$text\n\n";

    $count++;
}

print OUT "set nomultiplot\n";

print OUT "pause -1\n";

close(OUT);

#make the window a sensible size
#$read = `xdpyinfo  | grep \'dimensions:\'`; chomp $read;
#($a, $dims, $a, $a,$a) = split(' ', $read);
#($width, $height) = split('x', $dims);
$width = 900;
$height = 600;

system("gnuplot -geometry ${width}x${height} -persist /tmp/multiplot_${random}.gnu");

system("rm /tmp/multiplot_${random}.gnu");

exit(0);


sub make_gnu_from_ax{
    local($filename) = @_ ;

    my $f = "\'${filename}\'";

    #below fit region
    my $text = "plot $f index 0 using 1:2 with lines ls 3, \\\n";
    $text = $text . "$f index 1 using 1:2 with lines ls 3, \\\n";
    $text = $text . "$f index 2 using 1:2 with lines ls 3, \\\n";

    #in fit region
    $text = $text . "$f index 3 using 1:2 with lines ls 1, \\\n";
    $text = $text . "$f index 4 using 1:2 with lines ls 1, \\\n";
    $text = $text . "$f index 5 using 1:2 with lines ls 1, \\\n";

    #above fit region
    $text = $text . "$f index 6 using 1:2 with lines ls 3, \\\n";
    $text = $text . "$f index 7 using 1:2 with lines ls 3, \\\n";
    $text = $text . "$f index 8 using 1:2 with lines ls 3, \\\n";

    #data
    #in fit region
    $text = $text . "$f index 9 using 1:2:3 with yerr ls 6, \\\n";
    #out of fit region
    $text = $text . "$f index 10 using 1:2:3 with yerr ls 3\n";

    #label text
    my $line = `cat $filename | grep 'gx'`;
    my ($junk, $chisq, $Z) = split('=', $line);
    chop $chisq; chop $chisq; chop $chisq;
    chop $Z;
    my($Zval, $Zerr) = split('\+-', $Z);
    chop $Zval; chop $Zerr;

    #get range
    $line =  `cat $filename | grep 'x '`;
    my ($a, $xlow, $xhigh) = split(' ', $line);

    $line = `cat $filename | grep 'y'`;
    my ($a, $ylow, $yhigh) = split(' ', $line);

    
    @out = ($text, $xlow, $xhigh, $ylow, $yhigh, $chisq, $Zval, $Zerr);

    return @out;
}
