/* ----------------------------------------------------
	ノイズ指数の計算
	2003/3/31 by Masahiro Nakayama 
	Ver. 1.1 (noise_figure1_1.cpp)
------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_1.h"
#include "current_spectrum1_1.h"
#include "common_function1_1.h" 
#include "noise_figure1_1.h"

/* ********************************** */
/* ----- ノイズ指数のMian関数 ----- */
void NoiseFigure(double sig[], double sid[], struct ComplexNumbers sigd[],
				 struct ComplexNumbers y11[], struct ComplexNumbers y12[],
				 struct ComplexNumbers y21[], struct ComplexNumbers y22[], 
				 double fr[], char *fp_name){

	/* --- 変数宣言 --- */
	double p[FIN], r[FIN], c[FIN], fmin[FIN], gn[FIN];
	struct ComplexNumbers zopt[FIN];


	/* --- 雑音パラメータの計算 --- */
	// Cappy提案の式
	if(NWAY==1){
		MinimumNF_Cappy(sig, sid, sigd,	y11, y12, y21, y22, 
						p, r, c, fmin, fr, gn, zopt);
	}
	// Gonzalez提案の式
  	else if(NWAY==0){
		MinimumNF_Gonzalez(sig, sid, sigd,	y11, y12, y21, y22, 
							p, r, c, fmin, fr, gn, zopt);
	}

	/* --- ファイルに書き込む --- */
	WriteNoise(p, r, c, fmin, fr, gn, zopt, fp_name);

}



/* ******************************************** */
/* --- ノイズ指数計算 Gonzalez提案の式による--- */
void MinimumNF_Gonzalez(double sig[], double sid[], struct ComplexNumbers sigd[],
				struct ComplexNumbers y11[], struct ComplexNumbers y12[],
				struct ComplexNumbers y21[], struct ComplexNumbers y22[], 
				double p[], double r[], double c[], double fmin[],
				double fr[], double gn[], struct ComplexNumbers zopt[]){

	int jf;
	double y11_abs, y21_abs, sig_abs, sid_abs;
	double rn, gcor, x, y, z, a, b, gopt;

	for(jf=0; ; jf++){
		if(fr[jf]<0){continue;} //周波数0Hzのインデックスを記録
		else if(fr[jf]==32168){break;} //最大周波数で計算をストップ

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
		fmin[jf] = 1+2*rn*(gcor+gopt);

		// 現状では求めることができない"zopt,gn"には"0"を代入
		zopt[jf].re = 0.0;
		zopt[jf].im = 0.0;
		gn[jf] = 0.0;
	}
}


/* ******************************************** */
/* --- ノイズ指数計算 Cappy提案の式による--- */
void MinimumNF_Cappy(double sig[], double sid[], struct ComplexNumbers sigd[],
				struct ComplexNumbers y11[], struct ComplexNumbers y12[],
				struct ComplexNumbers y21[], struct ComplexNumbers y22[], 
				double p[], double r[], double c[], double fmin[],
				double fr[], double gn[], struct ComplexNumbers zopt[]){
	int jf;
	double y11_abs, y21_abs, y12_abs, y22_abs;
	double rn, x, a, b;
	double gmat, bmat, d, e, f, g, h, i;
	double svg, svd;
	struct ComplexNumbers svgd, zc;

	for(jf=0; ; jf++){
		if(fr[jf]<0){continue;} //周波数0Hzのインデックスを記録
		else if(fr[jf]==32168){break;} //最大周波数で計算をストップ

		gmat = 	(y11[jf].re*y22[jf].re - y11[jf].im*y22[jf].im) - (y12[jf].re*y21[jf].re - y12[jf].im*y21[jf].im);
		bmat = 	(y11[jf].re*y22[jf].im + y11[jf].im*y22[jf].re) - (y12[jf].re*y21[jf].im + y12[jf].im*y21[jf].re);



		y11_abs = (y11[jf].re*y11[jf].re + y11[jf].im*y11[jf].im);
		y12_abs = (y12[jf].re*y12[jf].re + y12[jf].im*y12[jf].im);
		y21_abs = (y21[jf].re*y21[jf].re + y21[jf].im*y21[jf].im);
		y22_abs = (y22[jf].re*y22[jf].re + y22[jf].im*y22[jf].im);

		p[jf] = sid[jf] / (4*BK*TA*sqrt(y21_abs));
		r[jf] = sig[jf]*sqrt(y21_abs) / (4*BK*TA*y11_abs);
		c[jf] = sigd[jf].im / sqrt(sig[jf]*sid[jf]);

		a = gmat*gmat+bmat*bmat;
		b = y22_abs * sig[jf];
		x = y12_abs * sid[jf];
		d = y22[jf].re*y12[jf].re + y22[jf].im*y12[jf].im;
		e =-y22[jf].re*y12[jf].im + y22[jf].im*y12[jf].re;
		f = 2* (d*sigd[jf].re - e*sigd[jf].im);
		svg = 1/a*(b+x-f);

		a = gmat*gmat+bmat*bmat;
		b = y21_abs * sig[jf];
		x = y11_abs * sid[jf];
		d = y21[jf].re*y11[jf].re + y21[jf].im*y11[jf].im;
		e =-y21[jf].re*y11[jf].im + y21[jf].im*y11[jf].re;
		f = 2* (d*sigd[jf].re - e*sigd[jf].im);
		svd = 1/(a)*(b+x-f);

		a = gmat*gmat+bmat*bmat;
		b = ( y22[jf].re*y21[jf].re + y22[jf].im*y21[jf].im) * sig[jf];
		x = ( y12[jf].re*y11[jf].re + y12[jf].im*y11[jf].im) * sid[jf];
		d = y22[jf].re*y11[jf].re + y22[jf].im*y11[jf].im;
		e =-y22[jf].re*y11[jf].im + y22[jf].im*y11[jf].re;
		f = d*sigd[jf].re - e*sigd[jf].im;
		g = y12[jf].re*y21[jf].re + y12[jf].im*y21[jf].im;
		h =-y12[jf].re*y21[jf].im + y12[jf].im*y21[jf].re;
		i = g*sigd[jf].re - h*sigd[jf].im;
		svgd.re = 1/(a)*(-b - x + f + i);

		a = gmat*gmat+bmat*bmat;
		b = (-y22[jf].re*y21[jf].im + y22[jf].im*y21[jf].re) * sig[jf];
		x = (-y12[jf].re*y11[jf].im + y12[jf].im*y11[jf].re) * sid[jf];
		d = y22[jf].re*y11[jf].re + y22[jf].im*y11[jf].im;
		e =-y22[jf].re*y11[jf].im + y22[jf].im*y11[jf].re;
		f = d*sigd[jf].im + e*sigd[jf].re;
		g = y12[jf].re*y21[jf].re + y12[jf].im*y21[jf].im;
		h =-y12[jf].re*y21[jf].im + y12[jf].im*y21[jf].re;
		i = g*sigd[jf].im + h*sigd[jf].re;
		svgd.im = 1/(a)*(-b - x + f + i);

		a = y22[jf].re + y12[jf].re*svgd.re / svd  - y12[jf].im*svgd.im / svd;
		b = y22[jf].im + y12[jf].re*svgd.im / svd  + y12[jf].im*svgd.re / svd;
		x = (gmat)/(gmat*gmat+bmat*bmat);
		d =-(bmat)/(gmat*gmat+bmat*bmat);
		zc.re = a*x - b*d;
		zc.im = a*d + b*x;

		gn[jf] = svd / (4*BK*TA)*(gmat*gmat+bmat*bmat)/y21_abs;

		a = svg / (4*BK*TA);
		b = svd / (4*BK*TA);
		x = y12_abs / y21_abs;
		d = (svgd.re*svgd.re + svgd.im*svgd.im)/svd/svd;
		rn = a-b*x*d;

		fmin[jf] = 1+2*gn[jf]*(zc.re+sqrt(zc.re*zc.re+rn/gn[jf]));

		zopt[jf].re = sqrt(zc.re*zc.re+rn/gn[jf]);
		zopt[jf].im = -zc.im;

	}
}

/* ******************************** */
/* ----- ノイズ指数を書き込む ----- */
void WriteNoise(double p[], double r[], double c[], double fmin[],
				double fr[], double gn[], struct ComplexNumbers zopt[], 
				char *fp_name){
    FILE *fp;
	int jf;

    fp = OpenFile(fp_name, 'w');

	fprintf(fp, "f[Hz]	P	R	C	fmin	gn	ropt	xopt\n");
	
	for(jf=0; ; jf++){
		if(fr[jf]<=0){continue;}
		else if(fr[jf]>1.05e+11){break;}

		fprintf(fp, "%1.6e	", fr[jf]);
		fprintf(fp, "%lf	%lf	%lf	",  p[jf], r[jf], c[jf]);
		fprintf(fp, "%lf	",  10*log10(fmin[jf]));
		fprintf(fp, "%lf	%lf	%lf\n",  gn[jf], zopt[jf].re, zopt[jf].im);
	}

	fclose(fp);
	printf("'%s' write\n", fp_name);

}
