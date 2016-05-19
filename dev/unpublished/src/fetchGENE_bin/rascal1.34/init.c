/* command line interface for Clustal W  */
/* DES was here MARCH. 1994 */
/* DES was here SEPT.  1994 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rascal.h"

/*
*	Prototypes
*/

static void init_alnout_options(ALNOUT_OPTPTR opt);
static void init_quickpw_options(QUICKPW_OPTPTR opt);
static void init_pw_options(PW_OPTPTR opt);
static void init_mult_options(MULT_OPTPTR opt);
static void init_gap_options(MULTGAP_OPTPTR opt);
static void init_prf_options(PRF_OPTPTR opt);
static void init_tree_options(TREE_OPTPTR opt);
/*
*	 Global variables
*/
static sint	debug;
ALNOUT_OPT alnout_opt;
QUICKPW_OPT quickpw_opt;
PW_OPT pw_opt;
MULT_OPT mult_opt;
MULTGAP_OPT gap_opt;
PRF_OPT prf_opt;
TREE_OPT tree_opt;


void init_options(OPTPTR opt)
{
	opt->explicit_type = 0;
        opt->alnout_opt=(&alnout_opt);
        opt->quickpw_opt=(&quickpw_opt);
        opt->pw_opt=(&pw_opt);
        opt->mult_opt=(&mult_opt);
        (*opt->mult_opt).gap_opt=(&gap_opt);
        opt->prf_opt=(&prf_opt);
        opt->tree_opt=(&tree_opt);

	init_alnout_options(opt->alnout_opt);
	init_quickpw_options(opt->quickpw_opt);
	init_pw_options(opt->pw_opt);
	init_mult_options(opt->mult_opt);
	init_gap_options((*opt->mult_opt).gap_opt);
	init_prf_options(opt->prf_opt);
	init_tree_options(opt->tree_opt);
}

static void init_tree_options(TREE_OPTPTR opt)
{
/*
   phylogenetic tree options
*/
	opt->treealgo = NJ;
	opt->use_ambiguities = FALSE;
	opt->tossgaps = FALSE;
	opt->kimura = FALSE;
	opt->output_tree_clustal   = FALSE;
	opt->output_tree_phylip    = TRUE;
	opt->output_tree_distances = FALSE;
	opt->output_tree_nexus     = FALSE;
}

static void init_gap_options(MULTGAP_OPTPTR opt)
{
/*
   multiple alignment gap options
*/
	strcpy(opt->hyd_residues,"GPSNDQEKR");
	opt->gap_dist = 4;
	opt->no_hyd_penalties = FALSE;
	opt->no_var_penalties = TRUE;
	opt->no_pref_penalties = FALSE;
	opt->use_endgaps = FALSE;
	opt->nendgappenalties = FALSE;
	opt->cendgappenalties = FALSE;

}

static void init_mult_options(MULT_OPTPTR opt)
{
/*
   multiple alignment parameters
*/
	opt->dna_gap_open = 15.0;
	opt->dna_gap_extend = 6.66;
	opt->transition_weight = 0.5;
	opt->dnamatnum = 1;
	strcpy(opt->dnamtrxname,"iub");
	strcpy(opt->dnausermtrxname,"");;
	opt->prot_gap_open = 10.0;
	opt->prot_gap_extend = 0.2;
	opt->matnum = 3;
	strcpy(opt->mtrxname,"gonnet");
	strcpy(opt->usermtrxname,"");
	opt->divergence_cutoff = 25;
	opt->no_weights = FALSE;
	opt->neg_matrix = FALSE;
	opt->reset_alignments_new  = FALSE;
	opt->reset_alignments_all  = FALSE;
	opt->propogate_motifs = FALSE;
	opt->use_motifs = FALSE;
	opt->use_ss_motifs = FALSE;
	strcpy(opt->motif_filename,"");
}

static void init_prf_options(PRF_OPTPTR opt)
{
	opt->profile_type = PROFILE;
}

static void init_pw_options(PW_OPTPTR opt)
{
/*
   slow pairwise alignment parameters
*/
	opt->dna_go_penalty = 15.0;
	opt->dna_ge_penalty = 6.66;
	opt->transition_weight = 0.5;
	opt->prot_go_penalty = 10.0;
	opt->prot_ge_penalty = 1.0;
	opt->matnum = 3;
	strcpy(opt->mtrxname,"gonnet");
	opt->dnamatnum = 1;
	strcpy(opt->dnamtrxname,"iub");
	strcpy(opt->usermtrxname,"");
	strcpy(opt->dnausermtrxname,"");
	opt->add_motif = FALSE;
	opt->quick_pairalign = FALSE;
	opt->treealgo = NJ;
}

static void init_quickpw_options(QUICKPW_OPTPTR opt)
{
/*
   quick pairwise alignment parameters
*/
	opt->dna_ktup      = 2;   /* default parameters for DNA */
	opt->dna_wind_gap  = 5;
	opt->dna_signif    = 4;
	opt->dna_window    = 4;

	opt->prot_ktup     = 1;   /* default parameters for proteins */
	opt->prot_wind_gap = 3;
	opt->prot_signif   = 5;
	opt->prot_window   = 5;
	opt->percent=TRUE;
}

static void init_alnout_options(ALNOUT_OPTPTR opt)
{
/*
   alignment output format parameters
*/
        opt->showaln        = TRUE; /* screen display for clustalw */
        opt->output_clustal = TRUE;
        opt->output_gcg     = FALSE;
        opt->output_phylip  = FALSE;
        opt->output_nbrf    = FALSE;
        opt->output_gde     = FALSE;
        opt->output_nexus   = FALSE;
        opt->output_gscope  = FALSE;
        opt->output_relacs  = FALSE;
        opt->output_tfa  = FALSE;
        opt->output_rsf  = FALSE;
        opt->lowercase       = TRUE; /* Flag for GDE output */
        opt->seq_numbers = FALSE;
        opt->output_order   = ALIGNED;
}

