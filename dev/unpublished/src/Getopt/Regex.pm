# $Id: Regex.pm,v 1.2 2002/05/27 13:16:55 willijar Exp $

require Exporter;
package Getopt::Regex;
@ISA=qw(Exporter);
@EXPORT_OK=qw(GetOptions);
use strict; 

=head1 NAME

Getopt::Regex - handle command line options flexibly using regular expressions 

=head1 SYNOPSIS

    use Getopt::Regex;
    GetOptions(\@ARGV,[$regex,$ref,$takesarg],...);

B<\@ARGV> - reference to array of command line arguments
B<$regex> - regular expression for identifying the option
B<$ref> - reference to a variable to be set or a function which is passed the argument if it exists.
B<$takesarg> - if 1 subsequent command line argument is taken as its argument

=head1 DESCRIPTION

This package provides a flexible yet simple method for handling
command line options. It does not stamp over the callers namespace or,
currently, inforce any particular standard for the options - the user
can do this if they want. By using anonymous closures sophisticated
option specifications can be constructed.

The function B<GetOptions>, exported from the package takes a
reference to the argument list followed by a set of option
specifications which are references to arrays containing at least a
regular expression to match for the option and a reference to a
variable to be set or a function to be called. A third optional
argument for each option, if set to 1, pulles off the following
command line argument as an argument for that option.

The simplest use is to set a boolean variable if an argument is set

 GetOptions(\@ARGV,['-[v|V]',\$bool,0]);

If the option C<'-v'> or C<'-V'> occurs $bool will be set to 1, otherwise
it is left unchanged.

A subsequent command line argument may be used as as an option argument
as follows

 GetOptions(\@ARGV,['-f',\$fname,1]);
 
will set $fname to fname if '-f fname' is specified on the command line.

Processed arguments are removed from the argument list.
Only the first occurance of an argument is processed, if a
variable is being set as in the above examples

More complex argument specifications are possible using anonymous
functions as arguments. If the option takes an argument, the argument
is passed to the function.  Parts of the regular expression may also
be used in the anonymous closure being executed. e.g.

 GetOptions(\@ARGV,
  ['-D(.+)=(.+)',sub { diagnostic $1,$2; }    ,0],
  ['-D(.+)',     sub { diagnostic $1; $_[0]; },1]);

The first option specification will match options in the format
'-DDIAGNOSTIC=VALUE',the second will match occurances of the format
'-DDIAGNOSTIC VALUE'.  As can be seen use of the regex matching
variables can be made in handling the options.

When B<GetOptions> is called with a function reference the function is
called for all matching occurances of the regular expression, and the
proceesed options are removed from the argument list.

Another useful example is

 GetOptions(\@ARGV,['-(no-)*optimization',0,sub { $optimization=!$1; }]);

which identifies the option statements '-optimization' or
'-no-optimization' setting $optimization true or false approptiately.

The option '--' ends the search for matching options - further
arguments are not searched.

=head1 NOTE

Requires at least Perl 5.000 or Perl 5.001m if anonymous closures are
used.

=head1 HISTORY

$Log: Regex.pm,v $
Revision 1.2  2002/05/27 13:16:55  willijar
Added in bug fix given be Edmond Abrahamian <edmond@tripos.com> which
permits it to work if an argument is zero (but not "").

Revision 1.1  2002/05/27 13:14:53  willijar
Initial revision

Revision 1.3  1995/12/16 09:48:28  willijar
Rename package to be more consistant with guidlines

Revision 1.2  1995/12/12 22:01:35  willijar
Removed unnecessary default argument. Extended pod documentation.

Revision 1.1  1995/12/10 15:57:06  willijar
Initial revision

=head1 BUGS

Please let me know.

=head1 TODO

Possibly integrate default behaviour of other option functions in this
package.

=head1 AUTHOR

John A.R. Williams, <J.A.R.Williams@aston.ac.uk>

=cut

sub GetOptions {
  my ($argv,@options)=@_;
  foreach (@options) {
      my ($regex,$ref,$takesarg)=@{$_};
      my @args=@{$argv};
      @{$argv}=();
      my $arg;
    argloop:
      while (($arg=shift @args) ne "") {
	  if ($arg=~/$regex/) {
	      my $val=1;
	      if ($takesarg) { $val=shift @args; }
	      if (ref($ref) eq 'CODE') { &$ref($val); } 
	      else {  ${$ref}=$val; last argloop;}
	  }
	  else { 
	      push @{$argv},$arg;
	      if ($arg eq '--') { last argloop; }
	  }
      }	
      push @{$argv},@args; 
   }
}

1;
