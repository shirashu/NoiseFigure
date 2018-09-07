/* ----------------------------------------------------
	ノイズ指数の計算
	2002/11/13  by Masahiro Nakayama 
	Ver. 1.0 (noise_figure1_0.cpp)
------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_0.h"
#include "current_spectrum1_0.h"
#include "common_function1_0.h" 
#include "noise_figure1_0.h"

/* ********************************** */
/* ----- ノイズ指数のMian関数 ----- */
void NoiseFigure(double sig[], double sid[], struct ComplexNumbers sigd[],
				 struct ComplexNumbers y11[], struct ComplexNumbers y21[],
				 double fr[], char *fp_name){

	/* --- 変数宣言 --- */
	int jf;
	double p[DATA], r[DATA], c[DATA], fmin[DATA];
	double y11_abs, y21_abs, sig_abs, sid_abs;
	double rn, gcor, a, b, gopt, x, y, z;


	/* --- 計算を行う --- */
	for(jf=0; ; jf++){
		if(fr[jf]<0){continue;}
		else if(fr[jf]>1.05e+11){break;}

		y11_abs = sqrt(y11[jf].re*y11[jf].re + y11[jf].im*y11[jf].im);
		y21_abs = sqrt(y21[jf].re*y21[jf].re + y21[jf].im*y21[jf].im);

		sig_abs = sig[jf];
		sid_abs = sid[jf];

		p[jf] = sid_abs / (4*BK*TA*y21_abs);
		r[jf] = sig_abs*y21_abs / (4*BK*TA*y11_abs*y11_abs);
		c[jf] = sigd[jf].im / sqrt(sig_abs*sid_abs);


		rn = sid_abs / (4*BK*TA*y21_abs*y21_abs);
		gcor = y11[jf].re - (y21[jf].re*sigd[jf].re + y21[jf].im*sigd[jf].im)/sid_abs;

		x =  (y11[jf].re*y21[jf].re + y11[jf].im*y21[jf].im)*sigd[jf].re;
		y = -(y11[jf].im*y21[jf].re - y11[jf].re*y21[jf].im)*sigd[jf].im;

		z = x+y;
		a = y21_abs*y21_abs*sig_abs/sid_abs + y11_abs*y11_abs - 2*z/sid_abs;

		b = y11[jf].im + (y21[jf].re*sigd[jf].im - y21[jf].im*sigd[jf].re)/sid_abs;
		b = b*b;

		gopt = sqrt(a-b);
		fmin[jf] = 10*log10(1+2*rn*(gcor+gopt));
	}

	/* --- ファイルに書き込む --- */
	WriteNoise(p, r, c, fmin, fr, fp_name);

}

/* ******************************** */
/* ----- ノイズ指数を書き込む ----- */
void WriteNoise(double p[], double r[], double c[], double fmin[],
				double fr[], char *fp_name){
    FILE *fp;
	int jf;

    fp = OpenFile(fp_name, 'w');

	fprintf(fp, "f[Hz]	P	R	C	fmin\n");
	
	for(jf=0; ; jf++){
		if(*(fr+jf)<0){continue;}
		else if(*(fr+jf)>1.05e+11){break;}

		fprintf(fp, "%1.6e	", fr[jf]);
		fprintf(fp, "%lf	%lf	%lf	",  p[jf], r[jf], c[jf]);
		fprintf(fp, "%lf\n",  fmin[jf]);
	}

	fclose(fp);
	printf("'%s' write\n", fp_name);

}
