/* ----------------------------------------------
	小信号パラメータ、雑音特性計算プログラム 
	2004/7/26 by Masahiro Nakayama 
	Ver.  2.0 (main_noise2_0.cpp)
---------------------------------------------- */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "new.h"

#include "parameter2_0.h"
#include "common_function2_0.h"
#include "y_parameter2_0.h"

/* ********************************************** */
/* ------------------ Main関数 ------------------ */
int main() {

	/* ================================================================= */
	/*** --- 変数宣言 --- ***/
	struct SetPara para;	//-- パラメータ一覧

	double *fr;				//-- 周波数インデックス
	struct YPara *y;		//-- Yパラメータ
	struct Spectrum *si;	//-- スペクトル密度
	struct CircuitPara *circ;		/*--等価回路パラメータ */

	/*** --- プロトタイプ宣言 --- ***/
	void YParameter(struct YPara y[], double fr[],struct CircuitPara circ[], struct SetPara para);
	void CurrentSpectrum(struct Spectrum si[], double fr[], struct SetPara para);
	void NoiseFigure(struct Spectrum si[], struct YPara y[],
					 double fr[],struct CircuitPara circ[], struct SetPara para);
	/* ================================================================= */

	/* ----- パラメータの読み込み */
	ReadParameter(&para);

	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	fr = new double[para.data];
	y = new struct YPara[para.data];
	si = new struct Spectrum[para.data];
	circ = new struct CircuitPara[para.data];

	/* ----- 求めたい周波数のインデックスを生成 */
	FrequencyIndex(fr, para);


	/* ********* Yパラメーター・等価回路定数計算 ********* */
	/* --- Yパラメーターの計算 */
	if(para.culc!=2){
		YParameter(y, fr,circ, para);
	}

	/* ********* 雑音スペクトル密度の計算 ********* */
	/* --- スペクトルの計算 */
	if(para.culc!=1){
		CurrentSpectrum(si, fr, para);
	}
	
	/* ********* 雑音特性の計算 ********* */
	/* --- ノイズパラメータ・雑音指数の計算 */
	if(para.culc==0){
		NoiseFigure(si, y, fr,circ, para);
	}


	/* ----- メモリの開放 */
	delete [] fr;
	delete [] y;
	delete [] si;

	/* ----- プログラムの停止処理 */
	if(para.p_stop == 1){
		StopPG();
	}

    return 0;
}
