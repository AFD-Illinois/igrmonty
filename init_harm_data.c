
#include "decs.h"
#include "harm_model.h"
/*

get HARM simulation data from fname

checks for consistency of coordinates in data file with
values of coordinate parameters 

Uses standard HARM data file format

CFG 1 Sept 07

*/


void init_harm_data(char *fname)
{
	FILE *fp;
	double x[4];
	double rp, hp, V, dV, two_temp_gam;
	int i, j, k;

	/* header variables not used except locally */
	double t, tf, cour, DTd, DTl, DTi, dt;
	int nstep, DTr, dump_cnt, image_cnt, rdump_cnt, lim, failed;
	double r, h, divb, vmin, vmax, gdet;
	double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
	double J ;

	fp = fopen(fname, "r");

	if (fp == NULL) {
		fprintf(stderr, "can't open sim data file\n");
		exit(1);
	} else {
		fprintf(stderr, "successfully opened %s\n", fname);
	}

	/* get standard HARM header */
	fscanf(fp, "%lf ", &t);
	fscanf(fp, "%d ", &N1);
	fscanf(fp, "%d ", &N2);
	fscanf(fp, "%lf ", &startx[1]);
	fscanf(fp, "%lf ", &startx[2]);
	fscanf(fp, "%lf ", &dx[1]);
	fscanf(fp, "%lf ", &dx[2]);
	fscanf(fp, "%lf ", &tf);
	fscanf(fp, "%d ", &nstep);
	fscanf(fp, "%lf ", &a);
	fscanf(fp, "%lf ", &gam);
	fscanf(fp, "%lf ", &cour);
	fscanf(fp, "%lf ", &DTd);
	fscanf(fp, "%lf ", &DTl);
	fscanf(fp, "%lf ", &DTi);
	fscanf(fp, "%d ", &DTr);
	fscanf(fp, "%d ", &dump_cnt);
	fscanf(fp, "%d ", &image_cnt);
	fscanf(fp, "%d ", &rdump_cnt);
	fscanf(fp, "%lf ", &dt);
	fscanf(fp, "%d ", &lim);
	fscanf(fp, "%d ", &failed);
	fscanf(fp, "%lf ", &Rin);
	fscanf(fp, "%lf ", &Rout);
	fscanf(fp, "%lf ", &hslope);
	fscanf(fp, "%lf ", &R0);

	/* nominal non-zero values for axisymmetric simulations */
	startx[0] = 0.;
	startx[3] = 0.;

	stopx[0] = 1.;
	stopx[1] = startx[1] + N1 * dx[1];
	stopx[2] = startx[2] + N2 * dx[2];
	stopx[3] = 2. * M_PI;

	fprintf(stderr, "Sim range x1, x2:  %g %g, %g %g\n", startx[1],
		stopx[1], startx[2], stopx[2]);

	dx[0] = 1.;
	dx[3] = 2. * M_PI;

	/* Allocate storage for all model size dependent variables */
	init_storage();

	two_temp_gam =
	    0.5 * ((1. + 2. / 3. * (TP_OVER_TE + 1.) / (TP_OVER_TE + 2.)) +
		   gam);
	Thetae_unit = (two_temp_gam - 1.) * (MP / ME) / (1. + TP_OVER_TE);

	dMact = 0.;
	Ladv = 0.;
	bias_norm = 0.;
	V = 0.;
	dV = dx[1] * dx[2] * dx[3];
	for (k = 0; k < N1 * N2; k++) {
		j = k % N2;
		i = (k - j) / N2;
		fscanf(fp, "%lf %lf %lf %lf", &x[1], &x[2], &r, &h);

		/* check that we've got the coordinate parameters right */
		bl_coord(x, &rp, &hp);
		if (fabs(rp - r) > 1.e-5 * rp || fabs(hp - h) > 1.e-5) {
			fprintf(stderr, "grid setup error\n");
			fprintf(stderr, "rp,r,hp,h: %g %g %g %g\n",
				rp, r, hp, h);
			fprintf(stderr,
				"edit R0, hslope, compile, and continue\n");
			exit(1);
		}

		fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf",
		       &p[KRHO][i][j],
		       &p[UU][i][j],
		       &p[U1][i][j],
		       &p[U2][i][j],
		       &p[U3][i][j],
		       &p[B1][i][j], &p[B2][i][j], &p[B3][i][j]);


		fscanf(fp, "%lf", &divb);

		fscanf(fp, "%lf %lf %lf %lf",
		       &Ucon[0], &Ucon[1], &Ucon[2], &Ucon[3]);
		fscanf(fp, "%lf %lf %lf %lf", &Ucov[0],
		       &Ucov[1], &Ucov[2], &Ucov[3]);
		fscanf(fp, "%lf %lf %lf %lf", &Bcon[0],
		       &Bcon[1], &Bcon[2], &Bcon[3]);
		fscanf(fp, "%lf %lf %lf %lf", &Bcov[0],
		       &Bcov[1], &Bcov[2], &Bcov[3]);

		fscanf(fp, "%lf ", &vmin);
		fscanf(fp, "%lf ", &vmax);
		fscanf(fp, "%lf ", &vmin);
		fscanf(fp, "%lf ", &vmax);
		fscanf(fp, "%lf ", &gdet);

	        /* additional stuff: current */
		fscanf(fp, "%lf ", &J) ;
		fscanf(fp, "%lf ", &J) ;
		fscanf(fp, "%lf ", &J) ;
		fscanf(fp, "%lf ", &J) ;

		fscanf(fp, "%lf ", &J) ;
		fscanf(fp, "%lf ", &J) ;
		fscanf(fp, "%lf ", &J) ;
		fscanf(fp, "%lf\n", &J) ;

		bias_norm +=
		    dV * gdet * pow(p[UU][i][j] / p[KRHO][i][j] *
				    Thetae_unit, 2.);
		V += dV * gdet;

		/* check accretion rate */
		if (i <= 20)
			dMact += gdet * p[KRHO][i][j] * Ucon[1];
		if (i >= 20 && i < 40)
			Ladv += gdet * p[UU][i][j] * Ucon[1] * Ucov[0];

	}

	bias_norm /= V;
	dMact *= dx[3] * dx[2];
	dMact /= 21.;
	Ladv *= dx[3] * dx[2];
	Ladv /= 21.;
	fprintf(stderr, "dMact: %g, Ladv: %g\n", dMact, Ladv);


	/* done! */

}
