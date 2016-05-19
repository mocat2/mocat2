#!/usr/bin/env perl
##############################
# File: Simple.pm
# Copyright (C) by Kai Wilker <kaiw@cpan.org>
# $Id: Simple.pm,v 1.7 2008/02/16 18:49:25 foo Exp foo $
##############################

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


package Config::File::Simple;

use 5.008;
use strict;
use warnings;
use Carp qw/ croak /;
use Tie::File;

our $VERSION = '1.00';

sub new {
    my ($class, $file) = @_;
    croak "Need a Configuration file! Usage: my \$object = new Config::File::Simple(\$config:file);" if !defined $file;
    my $self = bless { file => $file }, $class;
    return $self;
}

sub read {
    if (@_ > 2) { # If there are more than one variables given use multiple_read instead of read
        my $self = shift @_;
        my %values = $self->multiple_read(@_);
        return %values;
    }

    my ($self, $variable) = @_;
    my $value = 0; # This value will be returned

    # Tests: Is a variable given? Has the variable no special characters? Does the config file exist?
    croak "No variable is given! Usage: \$self->read(\$variable)" if !defined $variable;
    croak "The variable '$variable' has got special characters!" if $self->has_special_characters($variable);
    croak "The configuration file '$self->{'file'}' doesn't exist!" if ! -e $self->{'file'};

    open my $CONFIG, "<", $self->{'file'} or croak "Can't open file '$self->{'file'}': $!";
    while(my $line = <$CONFIG>) { # Now we parse the config file and search for the variable
        chomp $line;
        $line =~ s/ [^\\] \# .* //xms; # We don't need the comments
        $line =~ s/ ^ \s+ //xms; # Delete all space at the beginnig
        $line =~ s/ \s+ $//xms; # Delete all space at the end
        next if $line !~ / ^ $variable \s* = /xms; # Is this the right variable?
        $value = (split m/=/, $line)[1]; # We need the value
        $value =~ s/ ^ \s+ //xms; # Delete all space at the beginnig of the value
        $value =~ s/ \s+ $ //xms; # Delete all space at the end of the value
        $value =~ s/\\#/#/g; # Unescape the escaped hashs
    }
    
    close $CONFIG or croak "Can't close file '$self->{'file'}': $!";
    return $value;
}

sub multiple_read {
    my ($self, @variables) = @_;
    my %values; # This hash with the variables and values will be returned

    # Tests: Have the variables no special characters? Does the config file exist?
    foreach my $variable (@variables) { # Check all variables for special characters
        croak "The variable '$variable' has got special characters!" if $self->has_special_characters($variable);
    }
    croak "The configuration file '$self->{'file'}' doesn't exist!" if ! -e $self->{'file'};

    open my $CONFIG, "<", $self->{'file'} or croak "Can't open file '$self->{'file'}': $!";
    while(my $line = <$CONFIG>) {
        chomp $line;
        $line =~ s/ [^\\] \# .* //xms; # Delete all comments
        $line =~ s/ ^ \s+ //xms; # Delete all space at the beginnig
        $line =~ s/ \s+ $//xms; # Delete all space at the end

        foreach my $variable (@variables) {
            next if $line !~ / ^ $variable \s* = /xms; # Is this the right variable?
            my $value; # This value will be added to the hash: $values{$variable} = $value
            $value = (split m/=/, $line)[1]; # We need the value, not the variable
            $value =~ s/ ^ \s+ //xms; # Delete all space at the beginnig of the value
            $value =~ s/ \s+ $ //xms; # Delete all space at the end of the value
            $value =~ s/\\#/#/g; # Unescape the escaped hashs
            $values{$variable} = $value;
        }
    }

    close $CONFIG or croak "Can't close file '$self->{'file'}': $!";
    return %values;
}

sub variable_exists {
    my $self = shift @_;
    return 0 if ! -e $self->{'file'}; # If the config file doesn't exist, the variable doesn't exist, too.
    return $self->read(@_); # read() will return 0 if the variable doesn't exist, otherwise it'll return the value of the variable
}

sub has_special_characters {
    my ($self, $word) = @_;
    return ( $word !~ /^ \w+ $/xms );
}

sub set { # A wrapper for methods add() and change()
    my ($self, $variable, $value) = @_;

    croak "The variable '$variable' has got special charaters!"
        if $self->has_special_characters($variable);
    $value =~ s/#/\\#/g; # Escape the hashs in the value: # -> \#
    
    # If the variable exists change the value, otherwise add a new variable + value
    if($self->variable_exists($variable)) {
        $self->change($variable, $value);
    } else {
        $self->add($variable, $value);
    }
}

sub add {
    my ($self, $variable, $value) = @_;

    croak "The variable '$variable' has got special charaters!"
        if $self->has_special_characters($variable);

    # Add a new line at the end of the config file
    # variable = value
    open my $CONFIG, ">>", $self->{'file'} or croak "Can't open file $self->{'file'}: $!";
    print {$CONFIG} "$variable = $value\n"
        or croak "Can't write at file '$self->{'file'}': $!";
    close $CONFIG or croak "Can't close file $self->{'file'}: $!";
}

sub change { # Changes the value of a variable
    my ($self, $variable, $value) = @_;

    croak "The variable '$variable' has got special charaters!"
        if $self->has_special_characters($variable);

    my @config; # The content of the config file will be in this array
    tie @config, 'Tie::File', $self->{'file'}
        or croak "Can't tie file '$self->{'file'}': $!";
    # Search for the variable and change it's value
    foreach my $line (@config) {
        if($line =~ / ^ ( \s* $variable \s* = \s* ) /xms) {
            $line = "$1$value\n";
        }
    }
    untie @config or croak "Can't untie file '$self->{'file'}': $!";
}

sub add_comment { # Adds a comment at the end of the file
    my $self = shift @_;

    open my $CONFIG, ">>", $self->{'file'} or croak "Can't open file '$self->{'file'}': $!";
    foreach my $text (@_) {
        print {$CONFIG} "# $text\n" or croak "Can't write at file '$self->{'file'}': $!";
    }
    close $CONFIG or croak "Can't close file '$self->{'file'}': $!";
}

1;

__END__

=pod

=head1 NAME

Config::File::Simple - a module to parse and edit a configuration file

=head1 SYNOPSIS

    use Config::File::Simple;

    my $configuration_file = 'path/to/the/config/file';
    my $config = new Config::File::Simple($configuration_file);

=head1 DESCRIPTION

C<Config::File::Simple> is a OO module to parse and edit configuration files.
The configuration file's syntax is simple:

    VARIABLE1 = VALUE1
    variable2 = value \# 2! # This is a comment
    # This is another comment

Special characters mustn't contained in the variable's name. Hashs in value have to be escaped
with a backslash: \#, but Config::File::Simple will do this for you (see METHODS).

=head1 METHODS

=head2 Parse the config file

Get the value of a variable is simple, too. You can use C<read()> to do this.

Usage:

    my $value_of_foo = $config->read('foo');

If there's only one variable given, C<read()> will return it's value.

If you give more than one arguments, C<read()> will return a hash with the
variables and the values:

    my %hash = $config->read('foo', 'bar', 'quux');
    my $value_of_foo = $hash{'foo'};

Of course you can give arrays as arguments:

    my @array = qw/foo bar quux/;
    my %hash = $config->read(@array);
    my $value_of_foo = $hash{'foo'};

If there are escaped hash (\#) in the value, C<read()> will unescape them for you:

Example:

This is the config file:

    foo = bar
    quux = qu \#u ux

And this is the Perl script:

    print $config->read('quux');

The output will be:

    qu #u ux

=head2 Edit the config file

There are two methods, to create and edit the configuration file:
C<set()> and C<add_comment()>.

=head3 set()

C<set()> is a wrapper which uses the commands C<add()> and C<change()>. It changes the value of 
a variable. If the variable doesn't exist, it will add this variable.

Usage of C<set()>:

    my $variable = 'foo';
    my $value = 'bar';
    if ($config->set($variable, $value)) {
        print "Success!\n";
    }

Now the variable C<foo> has got the value C<bar>.

You can't give more than one arguments, and so you can't use it in list context. I hope this will
be fixed in a new version

=head3 add_comment()

C<add_comment()> adds a comment at the end of the file.

Usage:

    my $comment = 'This is a test!';
    if ($config->add_comment($comment)) {
        print "Success!\n";
    }

Now there is a new line at the end of the file:

    # This is a test!

In list context C<add_comment()> will add all comments:

    my @comments = ('Comment 1', 'Comment 2');
    if ($config->add_comment(@comments)) {
        print "Success!\n";
    }

Now there are these two lines at the end of the file:

    # Comment 1
    # Comment 2

=head1 AUTHOR

The module id developed by Kai Wilker. You can send questions to this module
to C<kaiw@cpan.org>.

=head1 BUGS

Feel free to send bug-reports to C<kaiw@cpan.org>.

=head1 SEE ALSO

See also L<Config::File>, L<Config::Simple>, L<Config::Any>.

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 by Kai Wilker. All rights reserved.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
