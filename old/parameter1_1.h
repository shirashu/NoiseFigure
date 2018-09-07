/* ----------------------------------------------
	パラメータ一覧
	2003/3/31 by Masahiro Nakayama 
	Ver. 1.1 (parameter1_1.h)
---------------------------------------------- */
#define DT 5.0e-15 // Δt：1ステップの時間
#define EPP 32000 // 超粒子1個あたりの電子数
#define FIN 10000 /* 周波数インデックスの配列数 */

/* ********************************************** */
// --- Yパラメータ計算
#define JTL0 12000 /* ステップ電圧を加えた後の時間ステップ */
#define JTP0 8000 /* ステップ電圧を加える前の時間ステップ */
#define JT 20000  /* JTL0+JTP0 --- Simulationを行う全時間ステップ  */
#define NUML 6000 /* 計算対象の時間ステップ数（電圧を加えた後のステップ数） */
#define NUMP 2000 /* 計算対象の時間ステップ数（電圧を加える前のステップ数） */
#define CUR0 1 //等価回路　　0:Ig=0、Is=Idと仮定，1:左記仮定なし、Ig=Is-Id
#define YAVE 30 // Yパラメーターの平滑化を行うときの　重み付けのステップ数
#define LOWPASS_Y 1 //0:行わない，1:行う

/* ********************************************** */
// --- 雑音解析 --- 
#define JTN0 300000 /* 全時間ステップ */
#define NUMN 285000 /* 計算対象の時間ステップ数（電圧を加えた後のステップ数） */
#define CUR1 1 //雑音解析　　0:Ig=0、Is=Idと仮定，1:左記仮定なし、Ig=Is-Id

// --- スペクトルの計算 ---
#define WAY 1 //0:数値積分法，1:FFT法
	// --- 数値積分法のときに必要な設定 --- 
	#define MSMALLS 20000 /* 相関関の時間ステップm ---自己相関関数 */
	#define MSMALLC 20000 /* 相関関数の時間ステップm ---相互相関関数 */
	#define MLARGE 270000 /* M (MLARGE > MSMALL) */

	// --- FFT法のときに必要な設定
	#define DATA 262144 /* FFT計算として参照するデータ数 */
	#define JBIT 18 /* DATA = 2^JBIT　と設定すること */
	#define TWIN 1 // スペクトルウィンドウ、0:利用せず，1:利用する
	#define FETAPP 1 // 、Sid一定,Sig f^2,Im[Sigd],f、0:利用せず，1:利用する
	#define SAVE 40 // 平滑化を行うときの　重み付けのステップ数(0のときは行わない)
	#define FMINAVE 3.0e+10 //2次曲線近似に使う周波数 - 最小値
	#define FMAXAVE 7.0e+10 //2次曲線近似に使う周波数 - 最大値

// --- 雑音指数解析 ---
#define NWAY 1 // ノイズ指数計算式 0:Gonzalez，1:Cappy

// --- 低域通過フィルタの設定 ---
#define LOWPASS 0 //0:行わない，1:行う
#define CUTFREQ 2.0e+11 /* カットオフ周波数 */
#define FWIN 1000 /* 周波数ウィンドウのフィルタ係数 */

/*************************************************/
/*---------------------------------*/
#define IN0 3 /* 電極 */
#define DV0 0.125	// ゲート電極に与えるステップ電圧
#define DV1 0.5		// ドレイン電極に与えるステップ電圧

#define VGC "vgc.txt" /*   ゲートステップ電圧入力の電流値ファイル名 */
#define VDC "vdc.txt" /* ドレインステップ電圧入力電流値ファイル名 */
#define VC "vc.txt" /* 雑音解析用電流データの入力ファイル名 */
#define YPARAM "Yparam.txt" /* Yパラメーターのファイル名 */
#define OUTPUT "output.txt" /* アウトプットのファイル名 */
#define SPE "CurSpectrum.txt" /*   出力ファイル名 */
#define COR "correlation.txt" /*   出力ファイル名 */
#define NIOSE "noise.txt" /*   ノイズ指数出力ファイル名 */

#define PAI 3.1415926535897 /* 円周率π */
#define BK 1.38066e-23 /* ボルツマン定数 */
#define TA 300 /* 温度 */

/* ----- 複素数表現のための構造体 -----*/
struct ComplexNumbers{
  double re;
  double im;
};
