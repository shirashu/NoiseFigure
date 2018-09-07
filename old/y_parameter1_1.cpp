/* --------------------------------------------------------------------
	Yパラメータおよび各種物理定数の計算
	2003/3/31 by Masahiro Nakayama 
	Ver. 1.1 (y_parameter1_1.cpp)
-------------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_1.h"
#include "common_function1_1.h"
#include "y_parameter1_1.h"
#include "cooley_tukey_FFT1_1.h"

/* ************************************* 
     Yパラメータ計算のMian関数
************************************** */
void YParameter(double cur[][JT][IN0], double ss[][IN0], double rex[][IN0],
				struct ComplexNumbers y11[], struct ComplexNumbers y12[],
				struct ComplexNumbers y21[], struct ComplexNumbers y22[],
				double fr[], char *fn_yparam, char *fn_out){
	int jf;
	double dv[2];
	double g[IN0][IN0][FIN], b[IN0][IN0][FIN];
	double cgd[FIN], cgs[FIN];
	double ri[FIN], gm0[FIN], tau[FIN], cds[FIN], gds[FIN];

	// 変数の初期化
	dv[0] = DV0;	dv[1] = DV1;

	// ローパスフィルターをかける
	if(LOWPASS_Y==1){
		LowPassFilter_y(cur);
	}

	// Yパラメーター計算
	CalculateYPara_inte(cur, g, b,	rex, ss, dv, fr);
	// 三角形ウィンドによる平滑
	if(YAVE!=0){
		FlatYPara_inte(g, b, fr);
	}
	// 物理定数計算
	CalculateOtherPara(g, b, cgd, cgs, ri, gm0, tau, cds, gds, fr);

	// Yパラメーターを関数受け渡し用変数に代入
	for(jf=0; ; jf++){
		if(fr[jf]==32168){break;} //最大周波数で計算をストップ
		y11[jf].re = g[0][0][jf];
		y21[jf].re = g[1][0][jf];
		y11[jf].im = b[0][0][jf];
		y21[jf].im = b[1][0][jf];

		y12[jf].re = g[0][1][jf];
		y22[jf].re = g[1][1][jf];
		y12[jf].im = b[0][1][jf];
		y22[jf].im = b[1][1][jf];
	}

	// ファイルに書き込む
	WriteParameter(g, b, cgd, cgs, ri, gm0, tau, cds, gds, fr, fn_yparam, fn_out);

}


/* ********************************************************************** 
/* -------------------------------------------------------------------------
  以下各種関数
------------------------------------------------------------------------- */

/* ****************************************************** */
/* ----- 低域通過フィルタ ----- */
void LowPassFilter_y(double cur[][JT][IN0]){
	int jt, tau, ii;
	double temp[IN0][JTL0][IN0], rr0, rr1;
	double a;

	for(ii=0; ii<=1; ii++){
	
	for(jt=1; jt<=JTL0-1; jt++){
		temp[0][jt][ii]=0.0;	temp[1][jt][ii]=0.0;
		rr0=0.0;			rr1=0.0;
		for(tau=-FWIN; tau<=FWIN; tau++){
			if(tau==0){ a = 2*CUTFREQ*DT;}
			else{
				a = sin(2*PAI*CUTFREQ*DT*tau)/(PAI*tau)/2*(1+cos(PAI*tau/(FWIN+1)));
			}
			if(jt-tau<0){a=0;}
			else if(jt-tau>JTL0-1){a=0;}

			Sigma(cur[0][jt-tau][ii]*a, &temp[0][jt][ii], &rr0);
			Sigma(cur[1][jt-tau][ii]*a, &temp[1][jt][ii], &rr1);
		}
	}

	for(jt=0; jt<=JTL0-1; jt++){
		cur[0][jt][ii] = temp[0][jt][ii];
		cur[1][jt][ii] = temp[1][jt][ii];
	}
	}

}

/* ************************************** */
/* --- 数値積分法 ----------------------- */
/* --- Yパラメータを計算 ---------------- */
void CalculateYPara_inte(double cur[][JT][IN0], double g[][IN0][FIN],
						 double b[][IN0][FIN], double rex[][IN0],
						 double ss[IN0][IN0], double dv[], double fr[]){
	int i, ii, jt, jf;
	double w, s1, s2, r1, r2, t;

	for(jf=0; ; jf++){
		if(fr[jf]==32168){break;} //最大周波数で計算をストップ
		w = 2*PAI*fr[jf];
		for(ii=0; ii<=1; ii++){
		for(i=0 ; i<=IN0-1; i++){
			s1 = cur[i][JTL0-1][ii]*sin(w*double(JTL0-1)*DT)/2;
			r1 = 0.0;
			s2 = cur[i][JTL0-1][ii]*cos(w*double(JTL0-1)*DT)/2;
			r2 = 0.0;
			for(jt=1; jt<=JTL0-2; jt++){
		        t = DT*double(jt);
				Sigma(cur[i][jt][ii]*sin(w*t), &s1, &r1);
				Sigma(cur[i][jt][ii]*cos(w*t), &s2, &r2);
			}
			g[i][ii][jf] = (rex[i][ii] - ss[i][ii])/dv[ii] + w*s1*DT/dv[ii];
			b[i][ii][jf] = w*s2*DT/dv[ii];
		}
		}
	}
}


/* **************************************************** */
/* --- 数値積分法 ------------------------------------- */
/* --- Yパラメータを平滑化三角形ウィンドによる平滑　--- */
void FlatYPara_inte(double g[][IN0][FIN], double b[][IN0][FIN], double fr[]){
	int jf, k, i, ii, l[FIN], fmax;
	double t1[FIN], t2[FIN], rs1, rs2;

	for(jf=0; ; jf++){
		if(fr[jf]==32168){fmax=jf-1; break;} //最大周波数で計算をストップ
	}

	for(ii=0; ii<=1; ii++){
	for(i=0 ; i<=IN0-1; i++){

		for(jf=0; jf<=fmax; jf++){
			t1[jf] = 0; t2[jf] = 0;
			rs1 = 0;	rs2 = 0;

			l[jf] = YAVE;
			if(jf-l[jf]<0)    {do{l[jf]--;}while(jf-l[jf]<0);}
			else if(jf+l[jf]>fmax)    {do{l[jf]--;}while(jf+l[jf]>fmax);}

			for(k=-l[jf]; k<=l[jf]; k++){
				Sigma((l[jf]-abs(k))*g[i][ii][jf-k], &t1[jf], &rs1);
				Sigma((l[jf]-abs(k))*b[i][ii][jf-k], &t2[jf], &rs2);
			}
		}
		for(jf=0; jf<=fmax; jf++){
			if(l[jf] == 0){continue;}
			g[i][ii][jf] = t1[jf]/l[jf]/l[jf];
			b[i][ii][jf] = t2[jf]/l[jf]/l[jf];
		}

	}
	}
}

/* **************************** */
/* --- 各種物理定数を求める --- */
void CalculateOtherPara(double g[][IN0][FIN], double b[][IN0][FIN],
						double cgd[], double cgs[],	double ri[],
						double gm0[], double tau[], double cds[],
						double gds[], double fr[]){
	int jf;
	double g11, b11, g12, b12, g21, b21, g22, b22;
	double d, sd, sg, w;

	for(jf=0; ; jf++){
		if(fr[jf]<=0){continue;} //0GHz以下の周波数では計算しない
		if(fr[jf]==32168){break;} //最大周波数で計算をストップ
		w = 2*PAI*fr[jf];
		g11 = g[0][0][jf];		b11 = b[0][0][jf];
		g12 = g[0][1][jf];		b12 = b[0][1][jf];
		g21 = g[1][0][jf];		b21 = b[1][0][jf];
		g22 = g[1][1][jf];		b22 = b[1][1][jf];

		cgd[jf] = -b12/w;
		d   = b11+b12;
		sd  = d * d;
		sg  = g11*g11;
		cgs[jf] = d/w*(1+sg/sd);
		ri[jf]  = g11/(g11*g11 +b11*b11 +b12*b12 -2*b11*b12);
		gm0[jf] = sqrt((g21*g21+(b21-b12)*(b21-b12)) * (1+(w*cgs[jf]*ri[jf])*(w*cgs[jf]*ri[jf])));
		tau[jf] = 1/w * acos((g21 - w*b21*ri[jf]*cgs[jf] + w*b12*cgs[jf]*ri[jf]) / gm0[jf]);
		cds[jf] = (b22+b12)/w;
		gds[jf] = g22;
	}
}

/* ********************************** */
/* ----- パラメーターを書き込む ----- */
void WriteParameter(double g[][IN0][FIN], double b[][IN0][FIN],
					double cgd[], double cgs[], double ri[],
					double gm0[], double tau[], double cds[],
					double gds[], double fr[], char *fn_yparam, char *fn_out){
	int jf;
	FILE *fp_yparam, *fp_out;

	// --- ファイルオープン
    fp_yparam = OpenFile(fn_yparam, 'w');
    fp_out    = OpenFile(fn_out, 'w');

	fprintf(fp_yparam, "f[Hz]	Re[Y11]	Im[Y11]	Re[Y12]	Im[Y12]	Re[Y21]	Im[Y21]	Re[Y22]	Im[Y22]\n");
	fprintf(fp_out, "f[Hz]	Cgd[F/m]	Cgs[F/m]	Cds[F/m]	Ri[ohm/m]	gm0[S/m]	gds[S/m]	tau[s]\n");

	for(jf=0; ; jf++){
		if(fr[jf]<=0){continue;}
		else if(fr[jf]>1.05e+11){break;} //最大周波数で計算をストップ

		fprintf(fp_yparam, "%1.7e	" , fr[jf]);
		fprintf(fp_yparam, "%lf	%lf	" , g[0][0][jf], b[0][0][jf]);
		fprintf(fp_yparam, "%lf	%lf	" , g[0][1][jf], b[0][1][jf]);
		fprintf(fp_yparam, "%lf	%lf	" , g[1][0][jf], b[1][0][jf]);
		fprintf(fp_yparam, "%lf	%lf\n", g[1][1][jf], b[1][1][jf]);

		fprintf(fp_out, "%1.7e	%1.7e	%1.7e	%1.7e	", fr[jf], cgd[jf], cgs[jf], cds[jf]);
		fprintf(fp_out, "%1.7e	%lf	%lf	%1.7e\n", ri[jf], gm0[jf], gds[jf], tau[jf]);
	}

	// --- ファイルクローズ
	fclose(fp_yparam);
	fclose(fp_out);

	printf("'%s' write\n", fn_yparam);
	printf("'%s' write\n", fn_out);
}
