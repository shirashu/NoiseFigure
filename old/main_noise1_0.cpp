/* ----------------------------------------------
	ノイズ指数計算プログラム 
	2002/11/13 by Masahiro Nakayama 
	Ver.  1.0 (main_noise1_0.cpp)
	パラメータに関しては，parameter.hを参照のこと．
---------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_0.h"
#include "current_spectrum1_0.h"
#include "cooley_tukey_FFT1_0.h"
#include "common_function1_0.h" 
#include "y_parameter1_0.h"
#include "noise_figure1_0.h"

/* ********************************************** */
/* ------------------ Main関数 ------------------ */
int main() {
	char fn_vdc[]=VDC, fn_vgc[]=VGC;
	char fn_yparam[]=YPARAM, fn_out[]=OUTPUT;
	char fn_spe[]=SPE ,fn_noise[]=NIOSE;
	double cur[IN0][JT][IN0], cur1[IN0][JT][IN0], ss[IN0][IN0], rex[IN0][IN0];
		/*--cur:電流値、cur1:ふらつき電流値*/
		/*--電流値 [0][*][*]:Gate,[1][*][*]:Drain,[2][*][*]:Source*/
		/*--       [*][*][0]:Gate Step電圧印加,[*][*][1]:Drain Step電圧印加*/
	struct ComplexNumbers y11[NUMS], y21[NUMS];//-- Yパラメータ
	double fr[NUMS];//周波数インデックス
	double corr[IN0][NUMS]; //-- 相関関数
	double sig[NUMS], sid[NUMS]; //-- ゲート，ドレインにおけるスペクトル
	struct ComplexNumbers sigd[NUMS];//-- ゲート・ドレイン相互におけるスペクトル


	/* ----- FILEから電流値を読み込み、変数に代入する */
	ReadCurrent(cur, 0, &fn_vgc[0]);
	ReadCurrent(cur, 1, &fn_vdc[0]);
	/* ----- 電流のふらつき、定常状態の電流値を計算 */
	CurInit(cur, cur1, ss, rex, &fn_spe[0]);
	/* ----- 求めたい周波数のインデックスを生成 */
	FrequencyIndex(fr);

	/* ----- Yパラメーター・等価回路定数の計算 */
	YParameter(cur1, ss, rex, y11, y21, fr, &fn_yparam[0], &fn_out[0]);

	/* ----- スペクトルからのノイズ計算 ----- */
	CurrentSpectrum(cur1, corr, sig, sid, sigd, fr, &fn_spe[0]);
	NoiseFigure(sig, sid, sigd, y11, y21, fr, &fn_noise[0]);

    return 0;
}