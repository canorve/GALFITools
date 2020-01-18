#!/usr/bin/perl
use Math::Trig;
use File::Copy;

## creates output files after DGCG run on galaxies.

$n = @ARGV[0]  or die "usage: $0 SersicIndex  \n";



$k=GetK($n);

print "k = $k for n = $n \n";



sub GetK
{
## solve the Sersic equation
# to get the dependence of K over 
# Sersic index

    my $n = @_[0];
   
    
    my ($resk);
    my ($fxa,$fxb);
    my ($lima,$limb);
    my ($fxres);

    my $count=1; 

    #limits
    $lima=0;
    $limb=100;
    

    $fxa = fx($n,$lima);
    $fxb = fx($n,$limb);

    $resk= ($lima + $limb)/2;
    
    $fxres=fx($n,$resk);


    if($fxa < 0 && $fxb > 0)
    {

	while($fxres > 0.00000001 || $fxres < -0.00000001 )
	{
	    
	    if($fxres > 0)
	    {
		$limb=$resk;
	    }
	    elsif($fxres < 0)
	    {
		$lima=$resk;
	    }
	    elsif($fxres==0)
	    {
		last;
	    }
	    $resk= ($lima + $limb)/2;
	    $fxres=fx($n,$resk);

	    $count++;

	    if ($count >= 10000) 
	    {
		last;
	    }

	}
    }
    elsif($fxa > 0 && $fxb < 0)
    {

	while($fxres > 0.00000001 || $fxres < -0.00000001 )
	{
	    if($fxres > 0)
	    {
		$lima=$resk;
	    }
	    elsif($fxres < 0)
	    {
		$limb=$resk;
	    }
	    elsif($fxres==0)
	    {
		last;
	    }
	    $resk= ($lima + $limb)/2;
	    $fxres=fx($n,$resk);
	    $count++;
	    
	    if ($count >= 10000)   # max iteration
	    {
		last;
	    }
	}
    }
    else
    {
	printf("no solution in the range: ($lima,$limb)\n");
    }
    
    $resk = substr($resk,0,10);
    return ($resk);
}

sub fx
{
# function to solve to get 
# relation between Sersic index and K
    my $n = @_[0];
    my $k = @_[1];

    my ($func);
    
    $func= exp(gammln(2*$n)) - 2 * exp(gammln(2*$n)) *gammp(2*$n,$k);

    return($func);
    
}

sub beta
{
# beta function
    my $z = @_[0];
    my $w = @_[1];
    
    return (exp(gammln($z) + gammln($w) - gammln($z + $w)));
}

sub gammln
{
# gamma function
# it returns the natural log of 
# gamma function
    my @cof=(76.18009172947146,-86.50532032941677, 24.01409824083091,-1.231739572450155,
	     0.1208650973866179e-2,-0.5395239384953e-5);
    my $xx = @_[0];
    
    my ($resul);
    my ($j);
    my ($x,$y,$tmp,$ser);

    $y=$x=$xx;
    $tmp=$x+5.5; $tmp -= ($x+0.5)*log($tmp);

    $ser=1.000000000190015;
    
    for ($j=0;$j<=5;$j++)
    {
	$ser += $cof[$j]/++$y;
    }
    $logarit=log(2.5066282746310005*$ser/$x);
    $resul= $logarit-$tmp;

    return ($resul);
}

sub  gammp
{
# incomplete gamma function
    my $a = @_[0];
    my $x = @_[1];
  
    my $ITMAX=100;
    my $EPS = 3.0e-7;
    my $FPMIN = 1.0e-30;
    my ($gamser,$gammcf,$gln);
    my ($n);
    my ($sum,$del,$ap);
    my ($i);
    my ($an,$b,$c,$d,$del2,$h);

    if ($x < 0.0 || $a <= 0.0)
    { 
	print "1\n" ;nrerror("Invalid arguments in routine gammp\n");
    }
    if ($x < ($a+1.0))
    {
	$gln=gammln($a);
	if ($x <= 0.0)
	{
	    if ($x < 0.0) 
	    {
		print "2\n" ; nrerror("x less than 0 in routine gser\n");
	    }
	    $gamser=0.0;
	    return ($gamser);
	}
	else
	{
	    $ap=$a;
	    $del=$sum=1.0/$a;
	    for ($n=1;$n<=$ITMAX;$n++)
	    {
		++$ap;
		$del *= $x/$ap;
		$sum += $del;
		if (abs($del) < abs($sum)*$EPS)
		{
		    $gamser=$sum*exp(-$x+$a*log($x)-($gln));
		    return ($gamser);
		}
	    }
	    print "3\n" ; nrerror("a too large, ITMAX too small in routine gser\n");
	    return;
	}
	return ($gamser);
    }
    else
    {
    
	$gln=gammln($a);         
	$b=$x+1.0-$a;
	$c=1.0/$FPMIN;
	$d=1.0/$b;
	$h=$d;
	for ($i=1;$i<=$ITMAX;$i++)
	{
	    $an = -$i*($i-$a);
	    $b += 2.0;
	    $d=$an*$d+$b;
	    if (abs($d) < $FPMIN) 
	    {
		$d=$FPMIN;
	    }
	    $c=$b+$an/$c;
	    if (abs($c) < $FPMIN)
	    { 
		$c=$FPMIN;
	    }
	    $d=1.0/$d;
	    $del2=$d*$c;
	    $h *= $del2;
	    if (abs($del2-1.0) < $EPS)
	    {
		last; 
	    }
	}
	if ($i > $ITMAX)
	{
	    print " 4\n" ;
	    nrerror("a too large, ITMAX too small in gcf\n");
	}
	$gammcf=exp(-$x+$a*log($x)-($gln))*$h;
	return (1.0-$gammcf);
    }
}

