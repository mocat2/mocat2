perl -lane '
BEGIN{$LAST=""};
$line=$_;
$a=scalar reverse sprintf "%b", $F[1];
$a=substr $a, 6,1;
if($a eq "0"){$a="2"};
 if ($F[0] eq $LAST || $LAST eq ""){
if ($F[2] ne "*") {push @{$H{$a}}, $line;};
} else {
for $r ((1,2)) {
foreach $str (@{$H{$r}}) {
$str =~ s/\t/\/$r\t/;
print $str;
}
}
%H=();
if ($F[2] ne "*") {push @{$H{$a}}, $line;};
}
$LAST=$F[0];
END{
for $r ((1,2)) {
foreach $str (@{$H{$r}}) {
$str =~ s/\t/\/$r\t/;
print $str;
}
}
}
'
