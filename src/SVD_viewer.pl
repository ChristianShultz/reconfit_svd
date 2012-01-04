#!/usr/bin/perl

#use warnings;

die "SVD_viewer.pl <t0> <sing_val_num> " unless $#ARGV == 1;

$t0 = shift(@ARGV);
$num = shift(@ARGV);

$found = 0;

open(IN, "< SVDLogs/SVDULog_t0${t0}");
while($line = <IN>){
    chomp $line;
    ($a, $b) = split('=', $line);

    if($found){
	#print "$b\n";
	@valerrs = split(':', $b);
	@vals; @errs;
	foreach(@valerrs){
	    $this = $_;
	    ($val, $err) = split('\+/-', $this);
	    chop $val; chomp $err;
	    push(@vals, $val); push(@errs, $err); 
	#    print "val = $val, err = $err\n";
	}
    }

    if($a eq "sigma_${num} "){	
	$found = 1;
	($sv, $sverr) = split('\+/-', $b);
	$sv = "${sv}+/-${sverr}";
	$sv =~ s/\s+//g;
    }
    elsif($a eq "(reset) sigma_${num} "){
	$found = 1;
	($sv, $sverr) = split('\+/-', $b);
	$sv = "${sv}+/-${sverr}(reset)";
	$sv =~ s/\s+//g;
    }
    else{ $found = 0; }
}
close(IN);
    
#print "@vals\n";
#print "@errs\n";

# build a data file for gnuplot histogram
&write_data;
&make_gnu;

exit(0);


sub write_data{
    $random = int(rand(1000));
    open(DATA, "> /tmp/histo_$random.data");
    
    %opnames;
    open(FILE, "< ops_phases");
    while($line = <FILE>){
	chomp $line;
	($num, $name, $r, $i) = split(' ',$line);
	$opnames{$num} = $name;
    }

    print DATA "num ";
    $count = 0;
    foreach (@vals){
	# replace count with the operator names
	$name = $opnames{$count};
        print DATA "$name $name ";
	$count++;
    }
    print DATA "\n";

    $N = $#vals + 1;

    print DATA "$sv  ";
    for($n = 0; $n < $N; $n++){    
	$val = $vals[$n];
	$err = $errs[$n];
	print DATA "$val $err   ";
    }
    print DATA "\n";
    close(DATA);
}

sub make_gnu {
    
    open(GNU, "> /tmp/histo_$random.gnu");
    print GNU "set yrange[-1.2:1.2]\n";
    $xmax = 0.9;
    print GNU "set xrange[-0.5:$xmax]\n";
    print GNU "set key noenhanced\n";
    print GNU "set bars fullwidth\n";
    print GNU "set style data histograms\n";

    print GNU "set style fill solid border -1\n";
    print GNU "set style histogram errorbars gap 5 lw 1\n";

#    for($n = 0; $n < $N; $n++){
#	$nn = $n + 1;
#	print GNU "set style line $nn lc rgb \"black\"\n";
#    }
#    print GNU "set style increment user\n";

    print GNU "set xtics scale 0\n";
   
    print GNU "plot \'/tmp/histo_$random.data\'";

    $string = "";
    for($n = 0; $n < $N; $n++){    
        $aa = 2 + 2*$n; $bb = $aa + 1;
        if($n == 0){
            $string = " using ${aa}:${bb}:xtic(1) ti col, \'\'";
        }
        else{
            $string = $string . " using ${aa}:${bb} ti col, \'\'";
        }
    }
    chop $string; chop $string; chop $string; chop $string;

    print GNU $string;
    print GNU "\n pause -1\n";
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

    system("gnuplot -geometry ${width}x${height} -persist /tmp/histo_$random.gnu");

#    system("cp /tmp/histo_${random}.gnu ./GNU");
#    system("cp /tmp/histo_${random}.data ./DATA");

    system("rm /tmp/histo_${random}.*");
}
