#include <stdio.h>
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
    {"Time",      't', "UTC",            0, "Time as fractional day"},
    { 0 }
    };


struct Arguments {
    char        *args[ nArgs ];       /* START_PA  END_PA  PA_INC */

    char        IntModel[80];
    char        ExtModel[80];

    long int    StartDate;
    double      Time;
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
                sscanf( arg, "%lf", &arguments->Time );
                break;
            case 'e': // external model
                strcpy( arguments->ExtModel, arg );
                break;
            case 'i': // internal model
                strcpy( arguments->IntModel, arg );
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
    int                 UseTS07=0;
    char                IntModel[20], ExtModel[20];
    int                 i, j, k;
    long int            JD, Date;

    double              UTC;
    Lgm_Vector          u, B;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo(0);
    Lgm_QinDentonOne qd;
    char model_str[] = "T96";

    // grid extents to match ptm
    double xstart, xfinal;
    double ystart, yfinal;
    double zstart, zfinal;
    //number of points to match ptm
    double n = 75; //TODO: read on cmd line

   /*
    *  Default option values.
    */
    arguments.StartDate = 20170907;
    strcpy( arguments.IntModel, "IGRF" );
    strcpy( arguments.ExtModel, "T96" );
    arguments.Time = 0;

    /*
     *  Parse CmdLine arguments and options
     */
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
    int hour = arguments.Time/10000;
    int minute = (arguments.Time - hour*10000)/100;
    int second = (arguments.Time - hour*10000) - minute*100;
    UTC = (float)hour + (float)minute/60.0 + (float)second/3600.0;
    Date = arguments.StartDate;
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    JD = Lgm_Date_to_JD( arguments.StartDate, UTC, mInfo->c );     // Compute JD.

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

    /*
     *  Qin-Denton parameters can be obtained automatically by date/time ....
     */
    Lgm_get_QinDenton_at_JD( JD, &qd, 1 , 0);        // Get (interpolate) the QinDenton vals
                                                    // from the values in the file at the
                                                    // given Julian Date.
    Lgm_set_QinDenton( &qd, mInfo );                 // Set params in mInfo structure.

    // Grid to match ptm_tec_interp.py
    // TODO: read on cmd line
    xstart = -15;
    xfinal = +15;
    ystart = -15;
    yfinal = +15;
    zstart = -15;
    zfinal = +15;

    // open output files
    int il = 1;
    char buffer[256];

    // B-Field
    sprintf(buffer, "./ptm_T89_data/bx3d_%s_%04d.txt", model_str, il);
    FILE *fBx = fopen(buffer, "w");
    sprintf(buffer, "./ptm_T89_data/by3d_%s_%04d.txt", model_str, il);
    FILE *fBy = fopen(buffer, "w");
    sprintf(buffer, "./ptm_T89_data/bz3d_%s_%04d.txt", model_str, il);
    FILE *fBz = fopen(buffer, "w");

    // E-Field
    sprintf(buffer, "./ptm_T89_data/ex3d_%s_%04d.txt", model_str, il);
    FILE *fEx = fopen(buffer, "w");
    sprintf(buffer, "./ptm_T89_data/ey3d_%s_%04d.txt", model_str, il);
    FILE *fEy = fopen(buffer, "w");
    sprintf(buffer, "./ptm_T89_data/ez3d_%s_%04d.txt", model_str, il);
    FILE *fEz = fopen(buffer, "w");

    // XYZ Grid
    sprintf(buffer, "./ptm_T89_data/xgrid.txt");
    FILE *fx = fopen(buffer, "w");
    sprintf(buffer, "./ptm_T89_data/ygrid.txt");
    FILE *fy = fopen(buffer, "w");
    sprintf(buffer, "./ptm_T89_data/zgrid.txt");
    FILE *fz = fopen(buffer, "w");

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
                // Write magnetic field from semi-empirical model
                fprintf(fBx, "%13g\n", B.x);
                fprintf(fBy, "%13g\n", B.y);
                fprintf(fBz, "%13g\n", B.z);
                // Write electric field as zero
                fprintf(fEx, "%13g\n", 0.0);
                fprintf(fEy, "%13g\n", 0.0);
                fprintf(fEz, "%13g\n", 0.0);
            }
        }
    }

    fclose(fBx);
    fclose(fBy);
    fclose(fBz);
    fclose(fEx);
    fclose(fEy);
    fclose(fEz);

    // write coordinates to data file
    for (i=0; i<n; i++) {
        u.x = xstart + (xfinal - xstart)/(n-1.0)*i;
        u.y = ystart + (yfinal - ystart)/(n-1.0)*i;
        u.z = zstart + (zfinal - zstart)/(n-1.0)*i;

        fprintf(fx, "%13g\n", u.x);
        fprintf(fy, "%13g\n", u.y);
        fprintf(fz, "%13g\n", u.z);
    }

    fclose(fx);
    fclose(fy);
    fclose(fz);

  Lgm_FreeMagInfo( mInfo );

  exit(0);
}
