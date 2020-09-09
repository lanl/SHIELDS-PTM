#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <argp.h>
#include <fcntl.h>
#include <Lgm_CTrans.h>
#include <Lgm_MagModelInfo.h>
#include <Lgm_QinDenton.h>


const  char *ProgramName = "ptm_lanlgeomag";
const  char *argp_program_version     = "0.99";
const  char *argp_program_bug_address = "<smorley@lanl.gov>";
static char doc[] = "\nGenerate field files for PTM.\n";

// Mandatory arguments
#define     nArgs   1
static char ArgsDoc[] = "OutDir";

/*
 *   Description of options accepted. The fields are as follows;
 *
 *   { NAME, KEY, ARG, FLAGS, DOC } where each of these have the following
 *   meaning;
 *
 *      NAME - Name of option's long argument (can be zero).
 *       KEY - Character used as key for the parser and it's short name.
 *       ARG - Name of the option's argument (zero if there isnt one).
 *     FLAGS - OPTION_ARG_OPTIONAL, OPTION_ALIAS, OPTION_HIDDEN, OPTION_DOC,
 *             or OPTION_NO_USAGE
 */

static struct argp_option Options[] = {
    {"IntModel",  'i', "internal_model", 0, "Internal Magnetic Field Model to use. Default is IGRF"},
    {"ExtModel",  'e', "external_model", 0, "External Magnetic Field Model to use. Default is T89."},
    {"StartDate", 'd', "yyyymmdd",       0, "Date "},
    {"Time",      't', "UTC",            0, "Time as HHMMSS"},
    {"Verbose",   'v', "Verbose",        0, "Verbose output? 1 - Yes; 0 - No"},
    { 0 }
    };


struct Arguments {
    char        *args[ nArgs ];

    char        IntModel[80];
    char        ExtModel[80];

    long int    StartDate;
    long int    Time;

    int         Verbose;
    };

/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state) {

    /* Get the input argument from argp_parse, which we
     *       know is a pointer to our arguments structure. */
        struct Arguments *arguments = state->input;
        switch( key ) {
            case 'd': // date
                arguments->StartDate = atol(arg);
                break;
            case 't':
                arguments->Time = atol(arg);
                break;
            case 'e': // external model
                strcpy( arguments->ExtModel, arg );
                break;
            case 'i': // internal model
                strcpy( arguments->IntModel, arg );
                break;
            case 'v':
                arguments->Verbose = atoi(arg);
                break;
            case ARGP_KEY_ARG:
                if (state->arg_num >= nArgs) {
                    /* Too many arguments. */
                    argp_usage (state);
                    }
                arguments->args[state->arg_num] = arg;
                break;
            case ARGP_KEY_END:
                if (state->arg_num < nArgs)
                /* Not enough arguments. */
                argp_usage (state);
                break;
            default:
                return ARGP_ERR_UNKNOWN;
            }
       return 0;
   }

/* Our argp parser. */
static struct argp argp = { Options, parse_opt, ArgsDoc, doc };


int main( int argc, char *argv[] ){

    struct              Arguments arguments;
    int                 verbose, UseTS07=0;
    char                IntModel[20], ExtModel[20];
    int                 i, j, k, sYear, sMonth, sDay, sDoy;
    int                 Year, Month, Day;
    long int            Date;

    double              sJD, JD, eJD, UTC;
    Lgm_Vector          u, B;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo(0);
    Lgm_QinDentonOne qd;
    char model_str[] = "T96";
    char buffer[256];
    char dirname[256];

    // grid extents to match ptm
    double xstart, xfinal;
    double ystart, yfinal;
    double zstart, zfinal;
    //number of points to match ptm
    double n = 121; //TODO: read on cmd line
    // Number of time steps
    int nt;
    int t_tot = 12;
    double deltat = 300;  // seconds between time snapshots
    double deltaJD = deltat/86400.0;

    /*
     *  Default option values.
     */
    arguments.StartDate = 20170907;
    strcpy( arguments.IntModel, "IGRF" );
    strcpy( arguments.ExtModel, "T96" );
    arguments.Time = 120000;
    arguments.Verbose = 1;

    /*
     *  Parse CmdLine arguments and options
     */
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
    verbose = arguments.Verbose;
    int hour = arguments.Time/10000;
    int minute = (arguments.Time - hour*10000)/100;
    int second = (arguments.Time - hour*10000) - minute*100;
    UTC = (float)hour + (float)minute/60.0 + (float)second/3600.0;
    Lgm_Set_Coord_Transforms(arguments.StartDate, UTC, mInfo->c);
    Lgm_Doy(arguments.StartDate, &sYear, &sMonth, &sDay, &sDoy);
    sJD = Lgm_JD(sYear, sMonth, sDay, UTC, LGM_TIME_SYS_UTC, mInfo->c);
    eJD = sJD + deltaJD*(t_tot-1);
    if (verbose > 1) {
        printf("Hour:Minute:Second = %02d:%02d:%02d (UTC=%f)\n", hour, minute, second, UTC);
        printf("JD = %f\n", JD);
        verbose = 1;
        }

    // Now set field options
    strcpy( IntModel,  arguments.IntModel );
    strcpy( ExtModel,  arguments.ExtModel );
    if ( !strcmp( ExtModel, "T87" ) ){
        mInfo->Bfield = Lgm_B_T87;
    } else if ( !strcmp( ExtModel, "CDIP" ) ){
        mInfo->Bfield = Lgm_B_cdip;
    } else if ( !strcmp( ExtModel, "EDIP" ) ){
        mInfo->Bfield = Lgm_B_edip;
    } else if ( !strcmp( ExtModel, "DUNGEY" ) ){
        mInfo->Bfield = Lgm_B_Dungey;
        strcpy(IntModel, "DUNGEY");
    } else if ( !strcmp( ExtModel, "IGRF" ) ){
        mInfo->Bfield = Lgm_B_igrf;
    } else if ( !strcmp( ExtModel, "TS04" ) ){
        mInfo->Bfield = Lgm_B_TS04;
    } else if ( !strcmp( ExtModel, "TS07" ) ){
        mInfo->Bfield = Lgm_B_TS07;
        UseTS07=1;
    } else if ( !strcmp( ExtModel, "T89c" ) ){
        mInfo->Bfield = Lgm_B_T89c;
    } else if ( !strcmp( ExtModel, "T96" ) ){
        mInfo->Bfield = Lgm_B_T96;
    } else if ( !strcmp( ExtModel, "OP77" ) ){
        mInfo->Bfield = Lgm_B_OP77;
    } else {
        // default
        mInfo->Bfield = Lgm_B_T89c;
    }
    if ( !strcmp( IntModel, "CDIP" ) ){
        mInfo->InternalModel = LGM_CDIP;
    } else if ( !strcmp( IntModel, "EDIP" ) ){
        mInfo->InternalModel = LGM_EDIP;
    } else if ( !strcmp( IntModel, "DUNGEY") ) {
        mInfo->InternalModel = LGM_DUNGEY;
    } else {
        // default
       mInfo->InternalModel = LGM_IGRF;
    }

    // Grid to match ptm_tec_interp.py
    // TODO: read on cmd line
    xstart = -15;
    xfinal = +15;
    ystart = -15;
    yfinal = +15;
    zstart = -15;
    zfinal = +15;

    int il = 1;  // Number of timestep. TODO: loop to generate multiple times if reqd.
    for (JD = sJD; JD <= eJD; JD += deltaJD ) {
        /*
         * Write to a single, new format, PTM output file
         * first line is:  Nx Ny Nz
         * second line is: X1, X2, ..., XN
         * third line is: Y1, Y2, ..., YN
         * fourth line is: Z1, Z2, ..., ZN
         * subsequent lines are: i j k Bx By Bz Ex Ey Ex
         * k increases fastest, then j, then i
         */

        /*
         *  Qin-Denton parameters can be obtained automatically by date/time ....
         */
        Lgm_JD_to_Date( JD, &Year, &Month, &Day, &UTC);
        Date = Year*10000 + Month*100 + Day;

        Lgm_Set_Coord_Transforms(Date, UTC, mInfo->c);
        Lgm_get_QinDenton_at_JD(JD, &qd, verbose , 0);  // Get (interpolate) the QinDenton vals
                                                        // from the values in the file at the
                                                        // given Julian Date.
        Lgm_set_QinDenton(&qd, mInfo);                  // Set params in mInfo structure.

        // open output files
        sprintf(dirname, "%s_%ld", model_str, Date);
        // First test for directory
        struct stat st = {0};
        if (stat(dirname, &st) == -1) {
            mkdir(dirname, 0700);
            }
        // E/B-Field output file
        sprintf(buffer, "%s_%ld/ptm_fields_%04d.txt", model_str, Date, il);
        FILE *fields = fopen(buffer, "w");

        // Write header
        fprintf(fields, "%4d %4d %4d\n", (int)n, (int)n, (int)n);
        for (i=0; i<n; i++) {
            fprintf(fields, "%12.5e ", xstart + (xfinal - xstart)/(n-1.0)*i);
            }
        fprintf(fields, "\n");
        for (i=0; i<n; i++) {
            fprintf(fields, "%12.5e ", ystart + (yfinal - ystart)/(n-1.0)*i);
            }
        fprintf(fields, "\n");
        for (i=0; i<n; i++) {
            fprintf(fields, "%12.5e ", zstart + (zfinal - zstart)/(n-1.0)*i);
            }
        fprintf(fields, "\n");

        // Loop over positions and write out to opened files
        for (i=0; i<n; i++){
            for (j=0; j<n; j++){
                for (k=0; k<n; k++){
                    //printf("%d, %d, %d\n", i, j, k);
                    u.x = xstart + (xfinal - xstart)/(n-1.0)*i;
                    u.y = ystart + (yfinal - ystart)/(n-1.0)*j;
                    u.z = zstart + (zfinal - zstart)/(n-1.0)*k;
                    // Call magnetic field function
                    mInfo->Bfield(&u, &B, mInfo);
                    // Write grid indices
                    fprintf(fields, "%4d %4d %4d ", i+1, j+1, k+1);
                    // Write magnetic field from semi-empirical model
                    // Write electric field as zero
                    fprintf(fields, " %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n", B.x, B.y, B.z, 0.0, 0.0, 0.0);
                }
            }
        }

        fclose(fields);
        il++;
    }

    // Now write tgrid.dat file
    sprintf(buffer, "%s_%ld/tgrid.dat", model_str, Date);
    FILE *tgrid = fopen(buffer, "w");
    for (nt=0; nt<t_tot; nt++) {
        fprintf(tgrid, "%.18e\n", 0.0 + deltat*nt);
        }
    fclose(tgrid);

    Lgm_FreeMagInfo( mInfo );

    exit(0);
}
