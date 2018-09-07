/* ----------------------------------------------
	パラメータ一覧
	2002/11/13  by Masahiro Nakayama 
	Ver. 1.0 (parameter1_0.h)
---------------------------------------------- */
#define JTL0 80000 /* ステップ電圧を加えた後の時間ステップ */
#define JTP0 40000 /* ステップ電圧を加える前の時間ステップ */
#define JT 120000   /* JTL0+JTP0 --- Simulationを行う全時間ステップ  */
#define NUMS 70000 /* 計算対象の時間ステップ数（電圧を加えた後のステップ数） */
#define NUMSBEF 20000 /* 計算対象の時間ステップ数（電圧を加える前のステップ数） */
#define EPP 1280000 // 超粒子1個あたりの電子数
#define DT 2.5e-15 // Δt：1ステップの時間
#define CUR0 0 //0:Ig=0、Is=Idと仮定，1:左記仮定なし、Ig=Is-Id

#define IN0 3 /* 電極 */
#define DV0 0.125	// ゲート電極に与えるステップ電圧
#define DV1 0.5		// ドレイン電極に与えるステップ電圧


/*************************************************/
/* --- Yパラメータおよびスペクトルの計算方法 --- */
#define WAY 0 //0:数値積分法，1:FFT法

/* --- 数値積分法のときに必要な設定 --- */
// --- Yパラメーター計算
#define YAVE 4 // Yパラメーターの平滑化を行うときの　重み付けのステップ数
//スペクトル計算
#define MSMALLS 1700 /* 相関関数の時間ステップm ---自己相関関数 */
#define MSMALLC 1700 /* 相関関数の時間ステップm ---相互相関関数 */
#define MLARGE 65000 /* M (MLARGE > MSMALL) */
#define STEP 1 /* 計算に用いる電流値、0：vgc(ゲートステップ)、1：vdc(ドレインステップ)*/

/* --- FFT法のときに必要な設定 --- */
// FFT計算
#define DATA  16384 /* FFT計算として参照するデータ数 */
#define JBIT  14 /* DATA = 2^JBIT　と設定すること */
//スペクトル計算
#define SWIN 0 // スペクトルウィンドウ、0:利用せず，1:利用する
#define SAVE 1 // 平滑化を行うときの　重み付けのステップ数


/*************************************************/
/*---------------------------------*/
#define PAI 3.1415926535897 /* 円周率π */
#define BK 1.38066e-23 /* ボルツマン定数 */
#define TA 300 /* 温度 */

#define VGC "vgc.txt" /*   ゲートステップ電圧の入力ファイル名 */
#define VDC "vdc.txt" /* ドレインステップ電圧の入力ファイル名 */
#define YPARAM "Yparam.txt" /* Yパラメーターのファイル名 */
#define OUTPUT "output.txt" /* アウトプットのファイル名 */
#define SPE "CurSpectrum.txt" /*   出力ファイル名 */
#define COR "correlation.txt" /*   出力ファイル名 */
#define NIOSE "noise.txt" /*   ノイズ指数出力ファイル名 */


/* ----- 複素数表現のための構造体 -----*/
struct ComplexNumbers{
  double re;
  double im;
};
