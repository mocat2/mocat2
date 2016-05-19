#!/usr/bin/perl

###########################################################
# Author:   KM Amada (kmamada@ifrec.osaka-u.ac.jp)
#
# Date      Changelog
###########################################################
# 07.26.13  Initial release
# 09.03.13  Added extensive warnings and error messages
#
###########################################################

use Getopt::Long;
use LWP::Simple;
use LWP::UserAgent;

my $BASEURL = "http://sysimm.ifrec.osaka-u.ac.jp/MAFFTash/REST/service.cgi";

my ( $WORKDIR, $PDBLIST, $OWNLIST, $HAT3FILE, $INSTRFILE );

GetOptions
(
    'd=s' => \$WORKDIR,
    'p=s' => \$PDBLIST,
    'o=s' => \$OWNLIST,
    'h=s' => \$HAT3FILE,
    'i=s' => \$INSTRFILE,
);


my $PDBLISTTMP = "/tmp/mafftash-rest-$$.pdb.inp";
unlink $PDBLISTTMP if -e $PDBLISTTMP;


######
# validation
&help("Required parameter : atleast one of either '-p' or '-o'") unless ( defined $PDBLIST || defined $OWNLIST);
&help("Required parameter : '-d'") if defined $OWNLIST && ! defined $WORKDIR;

$HAT3FILE = "hat3" unless defined $HAT3FILE;
$INSTRFILE = "instr" unless defined $INSTRFILE;
chop $WORKDIR if defined $WORKDIR && $WORKDIR =~ m/\/$/g;


######
# prepare inputs
my @files = ();
push(@files, "strweight" => "0.5");
push(@files, "premafft" => "1");


# pdb entries
if ( defined $PDBLIST )
{
    &bail("Error: Input file $PDBLIST does not exists!") unless -e $PDBLIST;


    if ( open(INPF,"<$PDBLIST") )
    {
        if ( open(OUTF,">$PDBLISTTMP") )
        {
            while(<INPF>)
            {
                chomp;
                if (/^(\w{5})$/)
                {
                    print OUTF ">PDBID\n$1\n";
                }
            }

            close OUTF;
        }
        else
        {
            close INPF;
            &bail("Error: Cannot open temporary file $PDBLISTTMP for writing!");
        }

        close INPF;
    }
    else
    {
        &bail("Error: Cannot open file $PDBLIST for reading!");
    }


    push(@files, "inputfile" => ["$PDBLISTTMP"]);

}



# upload own structures

my %ownids = ();

if ( defined $OWNLIST  )
{
    &bail("Error: Input file $OWNLIST does not exists!") unless -e $OWNLIST;

    if ( open(OWNINPF,"<$OWNLIST") )
    {
        while(<OWNINPF>)
        {
            chomp;

            if ( /^(\w{5})$/ )
            {
                my $fileref = "$WORKDIR/$1.pdb";

                unless (-e $fileref)
                {
                    close OWNINPF;
                    &bail("Error: File $fileref does not exists!");
                }


                push(@files, "inputownfile[]" => ["$fileref"]);
                $ownids{$1} = 1;
            }
        }

        close OWNINPF;
    }
    else
    {
        &bail("Error: Cannot open file $OWNLIST for reading!");
    }

}



######
# start rest service
my $browser = LWP::UserAgent->new;
$browser->timeout(0);


# post: running a mafftash job
my $postResponse = $browser->post ( $BASEURL, \@files, 'Content_Type' => 'form-data' );
#&bail(sprintf("[%d] %s\n", $postResponse->code, $postResponse->message)) unless($postResponse->is_success);
&bail(sprintf("[%d] %s\n", $postResponse->code, &parseError($postResponse->content))) unless($postResponse->is_success);


# get response from post request
my ($status, $mafftashid) = &parseResponse($postResponse->content);


# wait for results until it becomes available
while(1)
{
    sleep 5;

    # get: get results for mafftash job
    my $getResponse = $browser->get("$BASEURL/premafft/$mafftashid");

    if ( $getResponse->is_success )
    {
        # get response from get request
        ($status, $mafftashid) = &parseResponse($getResponse->content);


        # job is not yet done. wait
        if ( $status eq "done")
        {
            &bail("Error retrieving hat3 file!") unless ( getstore("$BASEURL/premafft/hat3/$mafftashid", $HAT3FILE) == 200 );
            &bail("Error retrieving instr file!") unless ( getstore("$BASEURL/premafft/instr/$mafftashid", $INSTRFILE) == 200 );
            last;
        }

        next;

    }
    else
    {
        #&bail(sprintf("[%d] %s\n", $getResponse->code, $getResponse->message));
        &bail(sprintf("[%d] %s\n", $getResponse->code, &parseError($getResponse->content)));
    }

}



# make sure outputs were generated
&bail("Error: Output file $HAT3FILE not found!") unless -e $HAT3FILE;
&bail("Error: Output file $INSTRFILE not found!") unless -e $INSTRFILE;


# warn if some ownids were ommitted
if ( scalar keys(%ownids) > 0 )
{
    my %instrids = ();

    if ( open(INSTRF,"<$INSTRFILE") )
    {
        while(<INSTRF>)
        {
            chomp;

            if ( /^>\d+_(\w{5})$/ )
            {
                $instrids{$1} = 1;
            }
        }

        close INSTRF;

        foreach my $id ( keys %ownids )
        {
            warn "Warning: Own structure $id was excluded from instr/hat3.\n" unless $instrids{$id};
        }
    }
    else
    {
        &bail("Error: Cannot open file $INSTRFILE for reading!");
    }

}



unlink $PDBLISTTMP if defined $PDBLISTTMP && -e $PDBLISTTMP;



####################
####################



sub parseResponse
{
    my $response = shift;


    #"status":"wait","mafftashid":"Ma8211432R"
    my $status =  ( $response =~ /\"status\":\"([^\s\"]+)\"/ ) ? $1 : "";
    my $mafftashid =  ( $response =~ /\"mafftashid\":\"([^\s\"]+)\"/ ) ? $1 : "";


    return ($status, $mafftashid);

}


sub parseError
{
    my $response = shift;

    #"error":"Invalid number of inputs found."
    my $errorstr = ( $response =~ /\"error\":\"([^\"]+)\"/ ) ? $1 : "";
    return $errorstr;
}



sub bail
{
    my $str = shift;
    print STDERR "$str\n" if defined $str;

    unlink $PDBLISTTMP if defined $PDBLISTTMP && -e $PDBLISTTMP;
    exit(1);
}



sub help
{
    my $str = shift;

    print <<'HELPME';

USAGE
  ./mafftash_premafft.pl -p [FILE]
  ./mafftash_premafft.pl -o [FILE] -d [DIRECTORY]
  ./mafftash_premafft.pl -p [FILE] -o [FILE] -d [DIRECTORY]


PARAMETERS
  -p [FILE]
     FILE contains a list of PDBIDs (one entry per line); make sure that the PDBIDs are in the standard 5-character pdbid+chain naming format

  -o [FILE] -d [DIRECTORY]
     FILE contains a list of IDs from your own structure/pdb files (one entry per line), IDs should be in the standard 5-character pdbid+chain naming format;
     for each ID in the list make sure that a corresponding structure file (same ID with .pdb extension) is stored in DIRECTORY

  -h [HATFILE]
     save the output hat3 file in HATFILE; if not set, the output is written to a file named 'hat3' in your current directory

  -i [INSTRFILE]
     save the output instr file in INSTRFILE; if not set, the output is written to a file named 'instr' in your current directory

HELPME

    &bail($str);
}



