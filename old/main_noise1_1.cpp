/* ----------------------------------------------
	ノイズ指数計算プログラム 
	2003/3/31 by Masahiro Nakayama 
	Ver.  1.1 (main_noise1_1.cpp)
---------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_1.h"
#include "current_spectrum1_1.h"
#include "cooley_tukey_FFT1_1.h"
#include "common_function1_1.h" 
#include "y_parameter1_1.h"
#include "noise_figure1_1.h"

/* ********************************************** */
/* ------------------ Main関数 ------------------ */
int main() {
	char fn_vdc[]=VDC, fn_vgc[]=VGC, fn_vc[]=VC;
	char fn_yparam[]=YPARAM, fn_out[]=OUTPUT, fn_spe[]=SPE ,fn_noise[]=NIOSE;
	double fr[FIN];//周波数インデックス
	double cur[IN0][JT][IN0], cur1[IN0][JT][IN0], ss[IN0][IN0], rex[IN0][IN0];
		/*--cur:電流値、cur1:ふらつき電流値*/
		/*--電流値 [0][*][*]:Gate,[1][*][*]:Drain,[2][*][*]:Source*/
		/*--       [*][*][0]:Gate Step電圧印加,[*][*][1]:Drain Step電圧印加*/
	double cur2[IN0][JTN0]; //-- ノイズ特性解析用データ
	struct ComplexNumbers y11[FIN], y21[FIN];//-- Yパラメータ
	struct ComplexNumbers y12[FIN], y22[FIN];//-- Yパラメータ
	double corr[IN0+1][JTN0]; //-- 相関関数
	double sig[NUMN], sid[NUMN]; //-- ゲート，ドレインにおけるスペクトル
	struct ComplexNumbers sigd[NUMN];//-- ゲート・ドレイン相互におけるスペクトル


	// ----- 求めたい周波数のインデックスを生成
	FrequencyIndex(fr);

	// ********* Yパラメーター・等価回路定数計算 ********* 
	// --- FILEから電流値を読み込む 
	ReadCurrent(cur, 0, &fn_vgc[0]);
	ReadCurrent(cur, 1, &fn_vdc[0]);
	// --- 電流のふらつき、定常状態の電流値を計算 
	CurInit(cur, cur1, ss, rex);
	// --- Yパラメーターの計算
	YParameter(cur1, ss, rex, y11, y12, y21, y22, fr, &fn_yparam[0], &fn_out[0]);

	// ********* ノイズ特性解析計算 ********* 
	// --- FILEから電流値を読み込む 
	ReadCurrentVC(cur2, &fn_vc[0]);
	// --- 電流のふらつき、定常状態の電流値を計算 
	CurInitVC(cur2, &fn_spe[0]);
	// --- スペクトルの計算
	CurrentSpectrum(cur2, corr, sig, sid, sigd, fr, &fn_spe[0]);
	// --- ノイズパラメータ・指数の計算
	NoiseFigure(sig, sid, sigd, y11, y12, y21, y22, fr, &fn_noise[0]);

    return 0;
}