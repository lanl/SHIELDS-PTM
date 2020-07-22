#include <Lgm_MagModelInfo.h>

/* BAL 04-Jan-2011 modified for more cases */

int main(){

    long int            Date;
    double              UTC;
    Lgm_Vector          u, B;
    Lgm_MagModelInfo    *mInfo;
    int                 i, j, k;
    char model_str[] = "T96";
    // grid extents to match ptm
    double xstart, xfinal;
    double ystart, yfinal;
    double zstart, zfinal;
    //number of points to match ptm
    double n = 75; //TODO: read on cmd line

    // TODO: read on cmd line
    Date = 20170907;
    UTC  = 0.0;

    mInfo = Lgm_InitMagInfo( );
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );

    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T96, mInfo );

    /*
     *  Qin-Denton parameters can be obtained automatically by date/time ....
     */
    JD = Lgm_Date_to_JD( Date, UTC, mInfo->c );     // Compute JD.
    Lgm_get_QinDenton_at_JD( JD, &p, 1 );           // Get (interpolate) the QinDenton vals
                                                    // from the values in the file at the
                                                    // given Julian Date.
    Lgm_set_QinDenton( &p, mInfo );                 // Set params in mInfo structure.

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
