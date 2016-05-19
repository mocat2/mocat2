#define MAXARGS 100

typedef struct {
	char *str;
	sint *flag;
	int type;
	char **arg;
} cmd_line_data;

/* 
   command line switches
*/
sint setoptions = -1;
sint sethelp = -1;
sint setinteractive = -1;
sint setbatch = -1;
sint setgapopen = -1;
sint setgapext = -1;
sint setpwgapopen = -1;
sint setpwgapext = -1;
sint setoutorder = -1;
sint setbootlabels = -1;
sint setmotifs = -1;
sint setpropogate = -1;
sint setssmotifs = -1;
sint setpwmatrix = -1;
sint setbootmatrix = -1;
sint setbootmatch = -1;
sint setmatrix = -1;
sint setpwdnamatrix = -1;
sint setdnamatrix = -1;
sint setnegative = -1;
sint setnoweights = -1;
sint setoutput = -1;
sint setoutputtree = -1;
sint setquicktree = -1;
sint settreealgo = -1;
sint setpwaddmotif = -1;
sint settype = -1;
sint setcase = -1;
sint setseqno = -1;
sint settransweight = -1;
sint setseed = -1;
sint setscore = -1;
sint setwindow = -1;
sint setktuple = -1;
sint setkimura = -1;
sint settopdiags = -1;
sint setpairgap = -1;
sint settossgaps = -1;
sint setnopgap = -1;
sint setnohgap = -1;
sint setnovgap = -1;
sint sethgapres = -1;
sint setvgapres = -1;
sint setuseendgaps = -1;
sint setmaxdiv = -1;
sint setgapdist = -1;
sint setdebug = -1;
sint setoutfile = -1;
sint setouttreefile = -1;
sint setinfile = -1;
sint setprofile1 = -1;
sint setprofile2 = -1;
sint setalign = -1;
sint setconvert = -1;
sint setnewtree = -1;
sint setusetree = -1;
sint setnewtree1 = -1;
sint setusetree1 = -1;
sint setnewtree2 = -1;
sint setusetree2 = -1;
sint setbootstrap = -1;
sint settree = -1;
sint setprofile = -1;
sint setsequences = -1;
sint setsecstr1 = -1;
sint setsecstr2 = -1;
sint setsecstroutput = -1;
sint sethelixgap = -1;
sint setstrandgap = -1;
sint setloopgap = -1;
sint setterminalgap = -1;
sint sethelixendin = -1;
sint sethelixendout = -1;
sint setstrandendin = -1;
sint setstrandendout = -1;

static char *type_arg[] = {
                "protein",
                "dna",
		""};

static char *bootlabels_arg[] = {
                "node",
                "branch",
		""};

static char *outorder_arg[] = {
                "input",
                "aligned",
                "intree",
		""};

static char *case_arg[] = {
                "lower",
                "upper",
		""};

static char *seqno_arg[] = {
                "off",
                "on",
		""};

static char *score_arg[] = {
                "percent",
                "absolute",
		""};

static char *treealgo_arg[] = {
                "nj",
                "bionj",
                "quicktree",
		""};

static char *output_arg[] = {
                "gcg",
                "gde",
                "pir",
                "phylip",
                "nexus",
                "gscope",
                "relacs",
                "rsf",
                "fasta",
		""};

static char *outputtree_arg[] = {
                "nj",
                "phylip",
                "dist",
                "nexus",
		""};

static char *outputsecstr_arg[] = {
                "structure",
                "mask",
                "both",
                "none",
		""};

/*
     command line initialisation

     type = 0    no argument
     type = 1    integer argument
     type = 2    float argument
     type = 3    string argument
     type = 4    filename
     type = 5    opts
*/
#define NOARG 0
#define INTARG 1
#define FLTARG 2
#define STRARG 3
#define FILARG 4
#define OPTARG 5


/* command line switches for DATA       **************************/
cmd_line_data cmd_line_file[] = {
     "infile",		&setinfile,		FILARG,	NULL,
     "profile1",	&setprofile1,		FILARG,	NULL,
     "profile2",	&setprofile2,		FILARG,	NULL,
     "",		NULL,			-1};
/* command line switches for VERBS      **************************/
cmd_line_data cmd_line_verb[] = {
     "help",		&sethelp,		NOARG,	NULL,
     "check",       &sethelp,    		NOARG,	NULL,
     "options",		&setoptions,		NOARG,	NULL,
     "align",		&setalign,		NOARG,	NULL,
     "newtree",		&setnewtree,		FILARG,	NULL,
     "usetree",		&setusetree,		FILARG,	NULL,
     "newtree1",	&setnewtree1,		FILARG,	NULL,
     "usetree1",	&setusetree1,		FILARG,	NULL,
     "newtree2",	&setnewtree2,		FILARG,	NULL,
     "usetree2",	&setusetree2,		FILARG,	NULL,
     "bootstrap",	&setbootstrap,		NOARG,	NULL,
     "tree",		&settree, 		NOARG,	NULL,
     "treealgo",	&settreealgo,		OPTARG,	treealgo_arg,
     "convert",		&setconvert,		NOARG,	NULL,
     "interactive",	&setinteractive,	NOARG,	NULL,
     "batch",		&setbatch,		NOARG,	NULL,
     "",		NULL,			-1};
/* command line switches for PARAMETERS **************************/
cmd_line_data cmd_line_para[] = {
     "type",		&settype,		OPTARG,	type_arg,
     "quicktree",	&setquicktree,		NOARG,	NULL,
     "profile",	&setprofile,	NOARG,	NULL,
     "sequences",	&setsequences,	NOARG,	NULL,
     "ssmotifs",	&setssmotifs,		NOARG,	NULL,
     "motifs",		&setmotifs,		FILARG,	NULL,
     "propogate",	&setpropogate,		NOARG,	NULL,
     "propagate",	&setpropogate,		NOARG,	NULL,
     "matrix",		&setmatrix,		FILARG,	NULL,
     "dnamatrix",	&setdnamatrix,		FILARG,	NULL,
     "negative",	&setnegative,		NOARG,	NULL,
     "noweights",	&setnoweights,		NOARG,	NULL,
     "gapopen", 	&setgapopen,		FLTARG,	NULL,
     "gapext",		&setgapext,		FLTARG,	NULL,
     "endgaps",		&setuseendgaps,		NOARG,	NULL,
     "nopgap",		&setnopgap,		NOARG,	NULL,
     "nohgap",		&setnohgap,		NOARG,	NULL,
     "novgap",		&setnovgap,		NOARG,	NULL,
     "hgapresidues",	&sethgapres,		STRARG,	NULL,
     "maxdiv",		&setmaxdiv,		INTARG,	NULL,
     "gapdist",		&setgapdist,		INTARG,	NULL,
     "pwmatrix",	&setpwmatrix,		FILARG,	NULL,
     "pwdnamatrix",	&setpwdnamatrix,	FILARG,	NULL,
     "pwgapopen",	&setpwgapopen,		FLTARG,	NULL,
     "pwgapext",	&setpwgapext,		FLTARG,	NULL,
     "pwmotif",		&setpwaddmotif,		NOARG,	NULL,
     "ktuple",		&setktuple,		INTARG,	NULL,
     "window",		&setwindow,		INTARG,	NULL,
     "pairgap",		&setpairgap,		INTARG,	NULL,
     "topdiags",	&settopdiags,		INTARG,	NULL,
     "score",		&setscore,		OPTARG,	score_arg,
     "transweight",	&settransweight,	FLTARG,	NULL,
     "seed",		&setseed,		INTARG,	NULL,
     "kimura",		&setkimura,		NOARG,	NULL,
     "tossgaps",	&settossgaps,		NOARG,	NULL,
     "bootlabels",	&setbootlabels,		OPTARG,	bootlabels_arg,
     "bootmatrix",      &setbootmatrix,         FILARG, NULL,
     "bootmatch",       &setbootmatch,          INTARG, NULL,
     "debug",		&setdebug,		INTARG,	NULL,
     "output",		&setoutput,		OPTARG,	output_arg,
     "outputtree",	&setoutputtree,		OPTARG,	outputtree_arg,
     "outfile",		&setoutfile,		FILARG,	NULL,
     "outtreefile",	&setouttreefile,	FILARG,	NULL,
     "outorder",	&setoutorder,		OPTARG,	outorder_arg,
     "case",		&setcase,		OPTARG,	case_arg,
     "seqnos",		&setseqno,		OPTARG,	seqno_arg,
     "nosecstr1",   &setsecstr1,		NOARG, NULL,
     "nosecstr2",   &setsecstr2,		NOARG, NULL,
     "secstrout",   &setsecstroutput,	OPTARG,  outputsecstr_arg,
     "helixgap",    &sethelixgap,		INTARG, NULL,
     "strandgap",   &setstrandgap,		INTARG, NULL,
     "loopgap",     &setloopgap,		INTARG, NULL,
     "terminalgap", &setterminalgap,	INTARG, NULL,
     "helixendin",  &sethelixendin,		INTARG, NULL,
     "helixendout", &sethelixendout,	INTARG, NULL,
     "strandendin", &setstrandendin,	INTARG, NULL,
     "strandendout",&setstrandendout,	INTARG, NULL,

     "",		NULL,			-1};


