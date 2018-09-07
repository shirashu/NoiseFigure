/* --------------------------------------------------------------------
	Yパラメータおよび各種物理定数の計算
	2015 T.takahashi 
	Ver. 3.0 (y_parameter2_0.cpp)
-------------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "new.h"

#include "parameter2_0.h"
#include "common_function2_0.h"
#include "y_parameter2_0.h"

/* ************************************* 
     Yパラメータ計算のMian関数
************************************** */
void YParameter(struct YPara y[], double fr[],struct CircuitPara circ[], struct SetPara para){


	/* ================================================================= */
	/* ================================================================= */
	/*** --- 変数宣言 --- ***/
	double dv[2];	/*--Step印加電圧 */
	struct Current *cur, ss, rex;	/*--cur:電流値、ss,rex:Step前,後平均値 */
	struct Admittance *ad;			/*--アドミタンス */
//	struct CircuitPara *circ;		/*--等価回路パラメータ */

	/*** --- プロトタイプ宣言 --- ***/
	void ReadCurrent_Y(struct Current cur[], int ii, char *fn, struct SetPara para);
	void CurInit_Y(struct Current cur[], struct Current *ss, struct Current *rex, 
				   struct SetPara para);
	void LowPass_y(struct Current cur[], struct SetPara para);
	void CalculateYPara_inte(struct Current cur[], struct Admittance ad[],
							 struct Current *rex, struct Current *ss, double dv[], 
							 double fr[], int ii, struct SetPara para);
	void FlatYPara_inte(struct Admittance ad[], double fr[], struct SetPara para);
	void CalculateOtherPara(struct Admittance ad[],	struct CircuitPara circ[], 
							double fr[], struct SetPara para);
	void SetYPara(struct Admittance ad[], struct YPara y[], double fr[],
				  struct SetPara para);
	void WriteParameter(struct Admittance ad[], struct CircuitPara circ[], 
						double fr[], struct SetPara para);

	/*** ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	cur = new struct Current[para.data];
	ad = new struct Admittance[para.data];
//	circ = new struct CircuitPara[para.data];

	/* ================================================================= */
	/* ================================================================= */


	/* --- FILEから電流値を読み込む */
	ReadCurrent_Y(cur, 0, para.vgc, para);
	ReadCurrent_Y(cur, 1, para.vdc, para);

	/* --- 電流のふらつき、定常状態の電流値を計算 */
	CurInit_Y(cur, &ss, &rex, para);

	/* --- 変数の初期化 */
	dv[0] = para.dv0;	dv[1] = para.dv1;

	/* --- ローパスフィルターをかける */
	if(para.lowpass_y==1){
		LowPass_y(cur, para);
	}

	/* --- Yパラメーター計算 */
	CalculateYPara_inte(cur, ad, &rex, &ss, dv, fr, 0, para);

	/* --- 三角形ウィンドによる平滑 */
	if(para.yave!=0){
		FlatYPara_inte(ad, fr, para);
	}
	/* --- 物理定数計算 */
	CalculateOtherPara(ad, circ, fr, para);

	/* --- アドミタンス成分をYパラメーターに代入 */
	SetYPara(ad, y, fr, para);


	/* --- ファイルに書き込む */
	WriteParameter(ad, circ, fr, para);
	

	/* ----- メモリの開放 */
	delete [] cur;
	delete [] ad;
//	delete [] circ;

}


/* ********************************************************************** 
/* -------------------------------------------------------------------------
  以下各種関数
------------------------------------------------------------------------- */

/* ****************************************************** */
/* ----- 低域通過フィルタ(ラグウィンドウ) ----- */
void LowPass_y(struct Current cur[], struct SetPara para){

	/* ----- 変数宣言 */
	int jf;
	double	*x;

	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	x = new double [(para.jtl0 * 4)];


	/* ----- ローパスフィルターをかけるための処理 */
	for(jf=0; jf<=para.jtl0-1; jf++){
		x[jf] = cur[jf]._[0][0];
		x[para.jtl0*1 + jf] = cur[jf]._[0][1];
		x[para.jtl0*2 + jf] = cur[jf]._[1][0];
		x[para.jtl0*3 + jf] = cur[jf]._[1][1];
	}

	/* ----- ローパスフィルター処理 */
	LowPassFilter(x, para.jtl0, para.cutfreq_y, para.fwin_y, para);

	/* ----- ローパスフィルターをかけた後の処理 */
	for(jf=0; jf<=para.jtl0-1; jf++){
		cur[jf]._[0][0]= x[jf];
		cur[jf]._[0][1]= x[para.jtl0*1 + jf];
		cur[jf]._[1][0]= x[para.jtl0*2 + jf];
		cur[jf]._[1][1]= x[para.jtl0*3 + jf];
	}

	/* ----- メモリの開放 */
	delete [] x;

}

/* ************************************** */
/* --- Yパラメータを計算 ---------------- */
void CalculateYPara_inte(struct Current cur[], struct Admittance ad[],
						 struct Current *rex, struct Current *ss, 
						 double dv[], double fr[], int ii, struct SetPara para){
	int i, jt, jf;
	double w, t;
	struct Current s1, s2, r1, r2;

	for(jf=0; ; jf++){
		if(fr[jf]>para.fmax_y){break;} //最大周波数で計算をストップ
		w = 2*para.pai*fr[jf];
		for(ii=0;	ii<=1;	ii++){
		for(i=0;	i<=1;	i++){
			s1._[i][ii] = cur[para.jtl0-1]._[i][ii]*sin(w*double(para.jtl0-1)*para.dt)/2;
			r1._[i][ii] = 0.0;
			s2._[i][ii] = cur[para.jtl0-1]._[i][ii]*cos(w*double(para.jtl0-1)*para.dt)/2;
			r2._[i][ii] = 0.0;
			for(jt=1; jt<=para.jtl0-2; jt++){
		        t = para.dt*double(jt);
				Sigma(cur[jt]._[i][ii]*sin(w*t), &s1._[i][ii], &r1._[i][ii]);
				Sigma(cur[jt]._[i][ii]*cos(w*t), &s2._[i][ii], &r2._[i][ii]);
			}
			ad[jf].g[i][ii] = (rex->_[i][ii] - ss->_[i][ii])/dv[ii] + w*s1._[i][ii]*para.dt/dv[ii];
			ad[jf].b[i][ii] = w*s2._[i][ii]*para.dt/dv[ii];
		}
		}
	}

}


/* **************************************************** */
/* --- Yパラメータを平滑化三角形ウィンドによる平滑　--- */
void FlatYPara_inte(struct Admittance ad[], double fr[], struct SetPara para){
	int jf, i, ii, *l, fmax;
	double *t1, *t2, *x, *y;


	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	l = new int[para.data];
	t1 = new double[para.data];
	t2 = new double[para.data];
	x = new double [para.data];
	y = new double [para.data];

	for(jf=0; ; jf++){
		if(fr[jf]>para.fmax_y){fmax=jf-1;break;} //最大周波数で計算をストップ
	}

	for(ii=0; ii<=1; ii++){
	for(i=0 ; i<=1; i++){
		for(jf=0; jf<=fmax; jf++){
			x[jf] = ad[jf].g[i][ii];
			y[jf] = ad[jf].b[i][ii];
		}
		if(para.savegd!=0){
			TriangleWindow(x, fmax, para.yave, para.fin);
			TriangleWindow(y, fmax, para.yave, para.fin);
		}
		for(jf=0; jf<=fmax; jf++){
			ad[jf].g[i][ii]= x[jf];
			ad[jf].b[i][ii]= y[jf];
		}
	}
	}

	/* ----- メモリの開放 */
	delete [] l;
	delete [] t1;
	delete [] t2;

}

/* **************************** */
/* --- 各種物理定数を求める (すべて変更)--- */
void CalculateOtherPara(struct Admittance ad[],
						struct CircuitPara circ[], double fr[], struct SetPara para){
	int jf;
	double g11, b11, g12, b12, g21, b21, g22, b22;
	double a, b, e, f, g, gg,ggg, w ,ft;

	ft=0.0;

	for(jf=0; ; jf++){
		if(fr[jf]>para.fmax_y){break;} //最大周波数で計算をストップ
		w = 2*para.pai*fr[jf];
		g11 = ad[jf].g[0][0];		b11 = ad[jf].b[0][0];
		g12 = ad[jf].g[0][1];		b12 = ad[jf].b[0][1];
		g21 = ad[jf].g[1][0];		b21 = ad[jf].b[1][0];
		g22 = ad[jf].g[1][1];		b22 = ad[jf].b[1][1];

		circ[jf].cgd = fabs(-b12/w);

		a  = (b11-w*circ[jf].cgd)/w;
		b  = (b11-w*circ[jf].cgd)*(b11-w*circ[jf].cgd);
		e  = g11*g11;
		circ[jf].cgs = a*(1+e/b);

		f  = pow(b11-w*circ[jf].cgd,2.0);
		circ[jf].ri  = g11/(f + g11*g11);
		
		g  = pow(b21+w*circ[jf].cgd,2.0);
		gg  = pow(w*circ[jf].cgd*circ[jf].ri,2.0);
//		circ[jf].gm0 = sqrt(g21*g21 + g * (1 + e));
		circ[jf].gm0 = sqrt(g21*g21 + pow(b21+w*circ[jf].cgd, 2.0)) * sqrt(1 + pow(w*circ[jf].ri*circ[jf].cgs, 2.0));
		
//		Rgs = g11/(pow(b21+w*circ[jf].cgd ,2.0) + g11*g11);
		circ[jf].gm  =  sqrt(g21*g21 + pow(b21+w*circ[jf].cgd, 2.0)) * sqrt(1 + pow(w*circ[jf].ri*circ[jf].cgs, 2.0));
		ggg = w * circ[jf].cgs * circ[jf].ri *g21;

		circ[jf].tau = fabs(1/w * asin(-(b21+w*circ[jf].cgd+ggg)/circ[jf].gm));

		circ[jf].ft  = circ[jf].gm / (2*para.pai*(circ[jf].cgs+circ[jf].cgd));

		circ[jf].cds = fabs((b22+w*circ[jf].cgd)/w);
		circ[jf].gds = g22;
		
		
	}
//				ftの平均化
//	for(jf=1;jf<5 ; jf++){
//	ft += circ[jf].ft;
//	}
//	ft = ft/(jf-1);
//	for(jf=0; ; jf++){
//		if(fr[jf]>para.fmax_y){break;} //最大周波数で計算をストップ
//	circ[jf].ft = ft;		
//	}
}


/* ********************************** */
/* ----- パラメーターを書き込む ----- */
void SetYPara(struct Admittance ad[], struct YPara y[], double fr[],
			  struct SetPara para){
	int jf;

	for(jf=0; ; jf++){
		if(fr[jf]>para.fmax_y){break;} //最大周波数で計算をストップ
		
		y[jf]._11.re = ad[jf].g[0][0];
		y[jf]._11.im = ad[jf].b[0][0];
		y[jf]._21.re = ad[jf].g[1][0];
		y[jf]._21.im = ad[jf].b[1][0];

		y[jf]._12.re = ad[jf].g[0][1];
		y[jf]._12.im = ad[jf].b[0][1];
		y[jf]._22.re = ad[jf].g[1][1];
		y[jf]._22.im = ad[jf].b[1][1];

	}
}


/* ********************************** */
/* ----- パラメーターを書き込む ----- */
void WriteParameter(struct Admittance ad[],	struct CircuitPara circ[], 
					double fr[], struct SetPara para){
	int jf;
	FILE *fp_yparam, *fp_out;

	// --- ファイルオープン
    fp_yparam = OpenFile(para.yparam, 'w');
    fp_out    = OpenFile(para.output, 'w');

	fprintf(fp_yparam, "f[Hz]	Re[Y11]	Im[Y11]	Re[Y21]	Im[Y21]	Re[Y12]	Im[Y12]	Re[Y22]	Im[Y22]\n");
	fprintf(fp_out, "f[Hz]	Cgd[F/m]	Cgs[F/m]	Cds[F/m]	Ri[ohm/m]	gm0[S/m]	gds[S/m]	tau[s]\n");

	for(jf=0; ; jf++){
		if(fr[jf]>para.fmax_y){break;} //最大周波数で計算をストップ

		fprintf(fp_yparam, "%1.7e	" , fr[jf]);
		fprintf(fp_yparam, "%lf	%lf	" , ad[jf].g[0][0], ad[jf].b[0][0]);
		fprintf(fp_yparam, "%lf	%lf	" , ad[jf].g[1][0], ad[jf].b[1][0]);
		fprintf(fp_yparam, "%lf	%lf	" , ad[jf].g[0][1], ad[jf].b[0][1]);
		fprintf(fp_yparam, "%lf	%lf\n", ad[jf].g[1][1], ad[jf].b[1][1]);

		if(fr[jf] == 0.00){continue;}
		fprintf(fp_out, "%1.7e	%1.7e	%1.7e	%1.7e	", fr[jf], circ[jf].cgd, circ[jf].cgs, circ[jf].cds);
		fprintf(fp_out, "%1.7e	%lf	%lf	%1.7e\n", circ[jf].ri, circ[jf].gm0, circ[jf].gds, circ[jf].tau);
	}

	// --- ファイルクローズ
	fclose(fp_yparam);
	fclose(fp_out);

	printf("'%s' write\n", para.yparam);
	printf("'%s' write\n", para.output);
}
