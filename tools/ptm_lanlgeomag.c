#include <Lgm_MagModelInfo.h>

/* BAL 04-Jan-2011 modified for more cases */

int main(){

    long int            Date;
    double              UTC;
    Lgm_Vector          u, B;
    Lgm_MagModelInfo    *mInfo;
    int                 i, j, k;


    //Date = 20150831;
		Date = 20080311; 
    UTC  = 0.0;

	// for a given Date, need .par files in /packages2/.packages2/x86_64-pc-linux-gnu-rhel7/lanlgeomag/1.5.15/share/LanlGeoMag/Data/TS07D_FILES/Coeffs/ 
    mInfo = Lgm_InitMagInfo( );
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );


    /*
     * Can set model manually ...
     */
    mInfo->Bfield = Lgm_B_T89;
    //mInfo->Bfield = Lgm_B_TS04;
    //mInfo->Bfield = Lgm_B_TS04;
    //mInfo->Bfield = Lgm_B_T87;
    //mInfo->Bfield = Lgm_B_T96;
    //mInfo->Bfield = Lgm_B_TS07;


		char model_str[] = "T89";

    /*
     * Or there is a setter routine ...
     */
    //Lgm_MagModelInfo_Set_MagModel( LGM_EDIP, LGM_EXTMODEL_T89, mInfo );



    /*
     * For TS07, the coeffs need to be initialized for each new time...
     */
    //Lgm_SetCoeffs_TS07( Date, UTC, &mInfo->TS07_Info );

    /*
     *  Can set/over-ride Qin-Denton parameters manually ....
    mInfo->P      = 4.1011111111111118;
    mInfo->Dst    = 7.7777777777777777;
    mInfo->By     = 3.7244444444444444;
    mInfo->Bz     = -0.12666666666666665;
    mInfo->W[0]   = 0.12244444444444445;
    mInfo->W[1]   = 0.2514;
    mInfo->W[2]   = 0.089266666666666661;
    mInfo->W[3]   = 0.047866666666666668;
    mInfo->W[4]   = 0.22586666666666666;
    mInfo->W[5]   = 1.0461333333333334;
		*/


    /*
     *  Or Qin-Denton parameters can be obtained automatically by date/time ....
     */
//    JD = Lgm_Date_to_JD( Date, UTC, mInfo->c );     // Compute JD.
//    Lgm_get_QinDenton_at_JD( JD, &p, 1 );           // Get (interpolate) the QinDenton vals 
//                                                    // from the values in the file at the 
//                                                    // given Julian Date.
//    Lgm_set_QinDenton( &p, mInfo );                 // Set params in mInfo structure.


    printf("%13s", "Kp");
    printf("%13s", "Ugsmx (Re)");
    printf("%13s", "Ugsmy (Re)");
    printf("%13s", "Ugsmz (Re)");
    printf("%13s", "Bgsmx (nT)");
    printf("%13s", "Bgsmy (nT)");
    printf("%13s", "Bgsmz (nT)");
    printf("%13s", "Bmag (nT)\n");


	/*
    for (i=0; i<=5; i++) {
        mInfo->Kp = i;
        u.x = -6.6; u.y =  0.0;  u.z =  0.0;
        //Lgm_Convert_Coords( &u, &ugsm, GEO_TO_GSM, mInfo->c );
        Lgm_Convert_Coords( &u, &ugsm, SM_TO_GSM, mInfo->c );
        mInfo->Bfield( &ugsm, &B, mInfo );
        printf( "%13i", mInfo->Kp);
        printf( "%13g%13g%13g", ugsm.x, ugsm.y, ugsm.z );
        printf( "%13g%13g%13g", B.x, B.y, B.z );
        printf( "%13g\n", Lgm_Magnitude( &B ) );
    }

    for (j=0; j<100; j++){
        mInfo->Kp = 3;
        for (i=0; i<13; i++) {
            u.x = -1.0 - (double)i * 0.5;
            u.y =  0.0;  u.z =  0.0;
            mInfo->Bfield( &u, &B, mInfo );
            printf( "%13i", mInfo->Kp);
            printf( "%13g%13g%13g", u.x, u.y, u.z );
            printf( "%13g%13g%13g", B.x, B.y, B.z );
            printf( "%13g\n", Lgm_Magnitude( &B ) );
        }
    }
	*/

	// Grid to match ptm_tec_interp.py

	double n = 75;
	double xstart, xfinal;
	double ystart, yfinal;
	double zstart, zfinal;

	xstart = -15;
	xfinal = +15;

	ystart = -15;
	yfinal = +15;

	zstart = -15;
	zfinal = +15;

	// open output files
	int il = 1;
	char buffer[256];


	// B- Field 
	sprintf(buffer, "./ptm_T89_data/bx3d_%s_%04d.txt", model_str, il);
	FILE *fBx = fopen(buffer, "w");

	sprintf(buffer, "./ptm_T89_data/by3d_%s_%04d.txt", model_str, il);
	FILE *fBy = fopen(buffer, "w");

	sprintf(buffer, "./ptm_T89_data/bz3d_%s_%04d.txt", model_str, il);
	FILE *fBz = fopen(buffer, "w");

	// E- Field 
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

 	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			for (k=0; k<n; k++){
				//printf("%d, %d, %d\n", i, j, k);
				u.x = xstart + (xfinal - xstart)/(n-1.0)*i;
				u.y = ystart + (yfinal - ystart)/(n-1.0)*j;
				u.z = zstart + (zfinal - zstart)/(n-1.0)*k;

				mInfo->Bfield(&u, &B, mInfo);

				fprintf(fBx, "%13g\n", B.x);
				fprintf(fBy, "%13g\n", B.y);
				fprintf(fBz, "%13g\n", B.z);

				fprintf(fEx, "%13g\n", 0.0);
				fprintf(fEy, "%13g\n", 0.0);
				fprintf(fEz, "%13g\n", 0.0);

				//printf("%7g, %7g, %7g: %13g, %13g, %13g\n", u.x, u.y, u.z, B.x, B.y, B.z);
			}
		}
	}

	fclose(fBx);
	fclose(fBy);
	fclose(fBz);

	fclose(fEx);
	fclose(fEy);
	fclose(fEz);


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
