
export LENGTHFILE=$2 &&
export MAPFILE=$3 &&

echo " -- LOCATION --" &&
samtools view $1 | perl -wlane '
BEGIN{
$LAST="";
if($ENV{MAPFILE} eq "NOMAP") {
$MAP=0;
} else {
print STDERR "      map file: $ENV{MAPFILE}";
$MAP=1;
open $in, "<", $ENV{MAPFILE} or die "ERROR & EXIT: Cannot load $ENV{MAPFILE}";
while ($line = <$in>) {
chomp $line;
@line = split "\t", $line;
$map{$line[0]}=$line[1];
}
}
};


$line = $_;
if($F[0] ne $LAST){
if (scalar keys %h2 > 1){$t="multi"} else {$t="unique"};
foreach $h (sort keys %h) {
print "$h\t$t\t$h{$h}{ID}\t$h{$h}{start}\t$h{$h}{length}";
}
%h=();
%h2=();
}

#####################################################
  unless ( $line =~ m/.*\s*NM:i:(\d+)\s.*/ )
  {
    die "ERROR & EXIT: Missing NM field";
  }
  $mm                = $1;
  $CIGAR             = $F[5];
  @CIGAR             = ( $CIGAR =~ /([0-9]+)([A-Z]+)/gi );
  $lengthBeforeMatch = 0;
  $matchLength       = 0;
  $match             = 0;
  $addToMatch        = 0;
  $addToMatchEnd     = 0;
  $atTheEnd = 0;
  for ( my $i = 0 ; $i < $#CIGAR ; $i += 2 )
  {
    if ( $CIGAR[ $i + 1 ] eq "M" || $CIGAR[ $i + 1 ] eq "=" )
    {
      $match++;
      $matchLength = $CIGAR[$i] + $addToMatch + $matchLength;
      $addToMatch  = 0;
    } else
    {
      if ( $match > 0 )
      {
       if ($atTheEnd && !($CIGAR[ $i + 1 ] eq "H" || $CIGAR[ $i + 1 ] eq "S")) {
        die "INTERNAL ERROR 3: No support for CIGAR letter $CIGAR[$i+1] ($i+1; $CIGAR[$i]$CIGAR[$i+1]). Please correct MOCAT2 accordingly. CIGAR=$CIGAR match=$match";
       }
        if ( ( $CIGAR[ $i + 1 ] eq "I" ) )
        {
          $addToMatch += $CIGAR[$i];
        } elsif ( $CIGAR[ $i + 1 ] eq "D" )
        {
          # no nothing
        } elsif ( $CIGAR[ $i + 1 ] eq "H" || $CIGAR[ $i + 1 ] eq "S") {
         $atTheEnd = 1;
         $addToMatchEnd += $CIGAR[$i]; 
        }else
        {
          die "INTERNAL ERROR 1: No support for CIGAR letter $CIGAR[$i+1] ($i+1; $CIGAR[$i]$CIGAR[$i+1]). Please correct MOCAT2 accordingly. CIGAR=$CIGAR match=$match";
        }
      }    # end match
      else
      {    # this is before matching begins
        if ( !( $CIGAR[ $i + 1 ] eq "H" || $CIGAR[ $i + 1 ] eq "S" || $CIGAR[ $i + 1 ] eq "I" ) )
        {
          die "INTERNAL ERROR 2: No support for CIGAR letter $CIGAR[$i+1] ($i+1; $CIGAR[$i]$CIGAR[$i+1]). Please correct MOCAT2 accordingly. CIGAR=$CIGAR match=$match";
        }
        $lengthBeforeMatch += $CIGAR[$i];
      }    # end not matched yet
    }    # end not matching M or =
  }    # end loop over cigar
  #$totalLength = $lengthBeforeMatch + $matchLength + $addToMatchEnd;
  $as = 100 - ( $mm / $matchLength ) * 100;
#####################################################

# SAVE TO HASH
if ($MAP == 1) {
die "ERROR & EXIT: Postprocess map file $ENV{MAPFILE} does not contain $F[2] $!" unless ($map{$F[2]});
$h2{$map{$F[2]}} = 1;
} else {
$h2{$F[2]} = 1;
}
$h{$F[2]}{start}=$F[3];
$h{$F[2]}{ID} = $as;
$h{$F[2]}{length}=$matchLength;
$LAST=$F[0];
# SAVE TO HASH

END{
if (scalar keys %h2 > 1){$t="multi"} else {$t="unique"};
foreach $h (keys %h) {
print "$h\t$t\t$h{$h}{ID}\t$h{$h}{start}\t$h{$h}{length}";
}
}
' | perl -F"\t" -wlane '
BEGIN{use Math::Round qw(:all);
open IN, "<$ENV{LENGTHFILE}" or die "MISSING \$ENV{LENGTHFILE}=$ENV{LENGTHFILE}!";
while (<IN>){
chomp;
@d=split "\t", $_;
$div{$d[0]} = $d[1];
}
close IN;
}
chomp(@F);
$x = round($F[3] / ($div{$F[0]}) * 200 );
$y = nearest(0.01, $F[3] / ($div{$F[0]}) * 200 );
if ($y =~ m/^\d+\.(\d+)$/) {
$y = $1 + 1;
} else {
$y = 1;
}
$h{$F[0]}{$F[1]}{nearest(.1, $F[2])}{$x}++;
if ($F[1] eq "unique" && nearest(.1, $F[2]) == 100) {
$h2{$F[0]}{$F[1]}{$x}{$y} = 1;
}

END{
@A=("unique", "multi");
foreach $a (sort keys %h){
foreach $b (@A) {
foreach $c (sort keys %{$h{$a}{$b}}) {
foreach $d (sort keys %{$h{$a}{$b}{$c}}) {
$t = scalar keys %{$h2{$a}{$b}{$d}};
print "$a\t$b\t$c\t$d\t$h{$a}{$b}{$c}{$d}\t$t";
}}}}}
' > $1.coverageLocation

echo "err $?"

if [ "$?" != "0" ]; then exit 666; fi

