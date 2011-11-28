#!/usr/bin/perl

# bit_test.pl -
#
# Sunday, November 27 2011
# 

#NB this is a comment
=for comment
 is a way of abusing POD to do long comments
=cut


use strict;
#name of the two directorys we are considering
my $pwd = qx(pwd);
my $path_old = "old_cho";
my $path_new = "new_cho";

#t0 considered
my $t = 8;
my $t0 = "t0${t}";

#assume no degeneracy and just sort by the value of the mass to form the hash
my $mname_old = "${path_old}/mass_t0${t}_state*_reordered.jack";
my $mname_new = "${path_new}/${t0}/MassJackFiles/mass_t0_${t}_reorder_state*.jack";

my @mold_names = ();
my @mnew_names = ();

@mold_names = qx (ls $mname_old);
@mnew_names = qx(ls $mname_new);

chomp @mold_names;
chomp @mnew_names;

my $num_states = $#mnew_names + 1;

my %mold = ();

foreach(@mold_names)
  {
    my @tmp = split(/\s+/, qx (calcbc $_));
    my ($trash,$state) = split(/state/,$_);
    my ($state,$trash) = split(/_/,$state);
    $mold{$tmp[1]} = $state;
  }

my %mnew = ();

foreach(@mnew_names)
  {
    my @tmp = split(/\s+/, qx (calcbc $_));
    my ($trash,$state) = split(/state/,$_);
    my ($state,$trash) = split(/\./,$state);
    $mnew{$tmp[1]} = $state;
  }

#check that we got what we wanted
=for comment #comment
my @keys = sort {$a <=> $b} keys %mnew;
foreach(@keys)
  {
    print "$_ $mnew{$_} $mold{$_} \n";
  }

my @keys = sort {$a <=> $b} keys %mold;
foreach(@keys)
  {
    print "$_ $mnew{$_} $mold{$_} \n";
  }
=cut #comment



my %hold = ();
my @keys = sort {$a <=> $b} keys %mold;
my $count = 0;
foreach(@keys)
  {
    $hold{$count} = $mold{$_};
    $count++;
  }

#keys should be the same
my @keys = sort {$a <=> $b} keys %mnew;
my %hnew = ();
my $count = 0;
foreach(@keys)
  {
    #print  "$mnew{$_} \n";
    $hnew{$count} = $mnew{$_};
    $count++;
  }

#abuse perl POD for comment, check the hash is working 
=for comment
my @keys = sort {$a <=> $b} keys %hnew;
foreach(@keys)
  {
    print "$_ $hnew{$_} \n";
  }


my @keys = sort {$a <=> $b} keys %hold;
foreach(@keys)
  {
    print "$_ $hold{$_} $hnew{$_}\n";
  }
=cut
#end comment


#begin the actual work
my $difflog;
my @keys = sort{$a <=> $b} keys %hold;

#check principal correlators
foreach(@keys)
  {
    my $string1 = "${path_old}/prin_corr_${t0}_state$hold{$_}_reordered.jack";
    my $string2 = " ${path_new}/${t0}/PrinCorrFiles/prin_corr_${t0}_state$hnew{$_}_reordered.jack";
    my $cmd = "diff $string1 $string2";
    my $tmp = qx($cmd);
    #print "$_ $cmd " . "\n";
    #print "$tmp \n";
    $difflog = $difflog . $cmd . "\n";
    $difflog = $difflog . $tmp;
  }

my $pcorr_diffs = $difflog;
$difflog = "";

open (DIFFLOG ,">pcorr_diffs");
print DIFFLOG $pcorr_diffs;
close (DIFFLOG);


#these were fit.. it doesn't make sense to check them
=for comment #comment
#check the mass jack files
foreach(@keys)
  {
    my $string1 = "${path_old}/mass_${t0}_state$hold{$_}_reordered.jack";
    my $string2 = " ${path_new}/${t0}/MassJackFiles/mass_t0_${t}_reorder_state$hnew{$_}.jack";
    my $cmd = "diff $string1 $string2";
    my $tmp = qx($cmd);
    #print "$_ $cmd " . "\n";
    #print "$tmp \n";
    $difflog = $difflog . $cmd . "\n";
    $difflog = $difflog . $tmp;
  }
my $mass_diffs = $difflog;
$difflog = "";

open (DIFFLOG ,">mass_diffs");
print DIFFLOG $mass_diffs;
close (DIFFLOG);

=cut
#end comment

#check generalized eigenvectors
foreach(@keys)
  {
    for(my $i = 0; $i < $num_states; $i++)
      {
	my $string1 = "${path_old}/v_${t0}_state$hold{$_}_op${i}_reordered.jack";
	my $string2 = " ${path_new}/${t0}/V_tJackFiles/V_t0_${t}_reordered_state$hnew{$_}_op${i}.jack";
	my $cmd = "diff $string1 $string2";
	my $tmp = qx($cmd);
	#print "$_ $cmd " . "\n";
	#print "$tmp \n";
	$difflog = $difflog . $cmd . "\n";
	$difflog = $difflog . $tmp;
      }
  }

my $vec_diffs = $difflog;
$difflog = "";
open (DIFFLOG ,">vec_diffs");
print DIFFLOG $vec_diffs;
close (DIFFLOG);

#these rely on the masses which were fit so we can't check them
=for comment

#check jackknife Z(t)
foreach(@keys)
  {
    for(my $i = 0; $i < $num_states; $i++)
      {
	my $string1 = "${path_old}/Z_${t0}_state$hold{$_}_op${i}_reordered.jack";
	my $string2 = " ${path_new}/${t0}/Z_tJackFiles/Zt_t0_${t}_reorder_state$hnew{$_}_op${i}.jack";
	my $cmd = "diff $string1 $string2";
	my $tmp = qx($cmd);
	#print "$_ $cmd " . "\n";
	#print "$tmp \n";
	$difflog = $difflog . $cmd . "\n";
	$difflog = $difflog . $tmp;
      }
  }

my $Z_diffs = $difflog;
$difflog = "";
open (DIFFLOG ,">Z_diffs");
print DIFFLOG $Z_diffs;
close (DIFFLOG);

=cut
#end comment

#these were fit so it doesnt really make sense to check them
=for comment #comment
#check jackknife Z fits
foreach(@keys)
  {
    for(my $i = 0; $i < $num_states; $i++)
      {
	my $string1 = "${path_old}/Z_fit_${t0}_state$hold{$_}_op${i}_reordered.jack";
	my $string2 = " ${path_new}/${t0}/ZJackFiles/Z_t0_${t}_reorder_state$hnew{$_}_op${i}.jack";
	my $cmd = "diff $string1 $string2";
	my $tmp = qx($cmd);
	#print "$_ $cmd " . "\n";
	#print "$tmp \n";
	$difflog = $difflog . $cmd . "\n";
	$difflog = $difflog . $tmp;
      }
  }

my $Zfit_diffs = $difflog;
$difflog = "";
open (DIFFLOG ,">Zfit_diffs");
print DIFFLOG $Zfit_diffs;
close (DIFFLOG);

=cut 
#end comment

#now go back and get rid of the stupid ones (gen eig vecs at t = t0, t=0)

my @files = ("pcorr_diffs","mass_diffs","vec_diffs","Z_diffs","Zfit_diffs");

my @stupid = (0,$t);

my @final = ();

foreach(@files)
  {
    open(FILE,"<$_");
    my @file = <FILE>;
    close(FILE);

    foreach(@stupid)
      {
	@file = reparset0($_,@file);
      }
    
    open (FILE,">${_}_reparse");
    print FILE @file;
    close(FILE); 

    push(@final,@file);
  }

open (DIFFLOG,">difflog");
print DIFFLOG @final;
close (DIFFLOG);






######
#subs#
######




sub reparset0 {
  my ($first,@lines) = @_;
  my @return;
  
  my $bound = $#lines +1;
  
  for(my $i = 0; $i < $bound; $i++)
    {
      if($lines[$i] =~ m/diff/)
	{
	  push (@return,$lines[$i]);
	}

      unless($lines[$i] =~ m/diff/)
	{
	  if($lines[$i] =~ m/c/)
	    {
	      my ($left,$right) = split(/c/,$lines[$i]);

	      if($left =~ m/\,/)
		{
		  push (@return, $lines[$i]);
		}
	      else
		{
		  if($left == $right)
		    {
		      my @tmp = ($first,$lines[$i],$lines[$i+1],$lines[$i+2],$lines[$i+3]);
		      if(parse4(@tmp))
			{
			  $i = $i + 3;
			}
		      else
			{
			  push (@return ,($lines[$i],$lines[$i+1],$lines[$i+2],$lines[$i+3]));
			  $i = $i + 3;
			}
		    }
		}
	    }
	  else
	    {
	      push(@return,$lines[$i]);
	    }
	}
    }

  return @return;

} #end reparse

sub parse4 {
  my ($num,@lines) = @_;
  my @tmp = split(/\s+/,$lines[1]);
  if ($num eq $tmp[1])
    {
      return 1;  #true
    }
  return; #false
}


