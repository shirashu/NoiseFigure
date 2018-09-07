/* ----------------------------------------------------
	スペクトル計算
	2015/ T.takahashi 
	Ver. 3.0 (current_spectrum2_0.cpp)
------------------------------------------------------- */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "new.h"

#include "parameter2_0.h" 
#include "common_function2_0.h" 
#include "current_spectrum2_0.h" 

/* ***************************+************** */
/* ----- 電流のスペクトル計算のMian関数 ----- */
void CurrentSpectrum(struct Spectrum si[], double fr[], struct SetPara para){


	/* ================================================================= */
	/* ================================================================= */
	/*** --- 変数宣言 --- ***/
	struct Current *cur;		//-- ノイズ特性解析用データ
	struct Correlation *corr;	//-- 相関関数

	/*** --- プロトタイプ宣言 --- ***/
	void ReadCurrent_Spe(struct Current cur[], struct SetPara para);
	void CurInit_Spe(struct Current cur[], struct SetPara para);
	void LowPassFilter_S(struct Current cur[], struct SetPara para);
	void Spectrum_Cor(struct Current cur[], struct Correlation corr[],
					  struct Spectrum si[], double fr[], struct SetPara para);
	void Spectrum_FFT(struct Current cur[], struct Correlation corr[],
					  struct Spectrum si[], double fr[], struct SetPara para);

	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	cur = new struct Current[para.data];
	corr = new struct Correlation[para.data];

	/* ================================================================= */
	/* ================================================================= */


	/* --- FILEから電流値を読み込む */
	ReadCurrent_Spe(cur, para);
	/* --- 電流のふらつき、定常状態の電流値を計算 */
	CurInit_Spe(cur, para);

	/* --- 低域通過フィルタ */
	if(para.lowpass==1){
		LowPassFilter_S(cur, para);
	}

	/* --- FFTによりスペクトルを計算 */
	if(para.way==1){
		Spectrum_FFT(cur, corr, si, fr, para);
	}
	/* --- 相関関数からスペクトルを計算 */
	else if(para.way==0){
		Spectrum_Cor(cur, corr, si, fr, para);
	}

	/* ----- メモリの開放 */
	delete [] cur;
	delete [] corr;


}


/* ********************************************************************** 
/* -------------------------------------------------------------------------
  以下各種関数
------------------------------------------------------------------------- */

/* ****************************************************** */
/* ----- 低域通過フィルタ ----- */
void LowPassFilter_S(struct Current cur[], struct SetPara para){
	int jf;
	double *x;

	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	x = new double [(para.jtn0 * 4)];


	/* ----- ローパスフィルターをかけるための処理 */
	for(jf=0; jf<=para.jtl0-1; jf++){
		x[jf] = cur[jf]._[0];
		x[para.jtn0*1 + jf] = cur[jf]._[1];
		x[para.jtn0*2 + jf] = cur[jf]._[2];
		x[para.jtn0*3 + jf] = 0.0;
	}

	/* ----- ローパスフィルター処理 */
	LowPassFilter(x, para.jtn0, para.cutfreq_s, para.fwin_s, para);

	/* ----- ローパスフィルターをかけた後の処理 */
	for(jf=0; jf<=para.jtl0-1; jf++){
		cur[jf]._[0]= x[jf];
		cur[jf]._[1]= x[para.jtn0*1 + jf];
		cur[jf]._[2]= x[para.jtn0*2 + jf];
	}

	/* ----- メモリの開放 */
	delete [] x;




}
