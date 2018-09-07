/* ----------------------------------------------
	共通パラメータ一覧
	2004/7/26 by Masahiro Nakayama 
	Ver. 2.0 (parameter2_0.h)
---------------------------------------------- */

#define PARAMETER "__SetPara20.txt" /* パラメータ一覧ファイル名 */

void ReadParameter(struct SetPara *para);

/* **************************************************** */
/* ----- 複素数表現のための構造体 -----*/
struct ComplexNumbers{
	double re;
	double im;
};


/* **************************************************** */
/* ----- Yパラメータのための構造体 -----*/
struct YPara{
	struct ComplexNumbers _11;
	struct ComplexNumbers _12;
	struct ComplexNumbers _21;
	struct ComplexNumbers _22;
};


/* **************************************************** */
/* ----- 雑音スペクトル密度表現のための構造体 -----*/
struct Spectrum{
	double g;//-- ゲート，ドレインにおけるスペクトル
	double d;//-- ゲート，ドレインにおけるスペクトル
	struct ComplexNumbers gd;//-- ゲート・ドレイン相互におけるスペクトル
};



/* ******************************************* */
/* ----- 読み込みパラメータのための構造体 -----*/
struct SetPara{

	int culc;	// 計算対象、0:すべて、1:Yパラのみ2:スペクトルのみ

	double dt;	// Δt：1ステップの時間

	// --- Yパラメータ計算
	int jtp0;	// ステップ電圧を加える前の時間ステップ
	int jtl0;	// ステップ電圧を加えた後の時間ステップ
	int jt;		// JTL0+JTP0 --- Simulationを行う全時間ステップ
	int nump;	// 計算対象の時間ステップ数（電圧を加える前のステップ数）
	int numl;	// 計算対象の時間ステップ数（電圧を加えた後のステップ数）
	int cur0;	//等価回路　　0:Ig=0、Is=Idと仮定，1:左記仮定なし、Ig=Is-Id
	int yave;	// Yパラメーターの平滑化を行うときの　重み付けのステップ数

	double fmin_y;	// 出力ファイルの最小周波数
	double fmax_y;	// 出力ファイルの最小周波数の最大周波数
	double fstep_y;	// 出力ファイルの最小周波数の周波数幅（yパラ計算にのみ有効）

	int lowpass_y;	//0:行わない，1:行う
	double cutfreq_y;	// カットオフ周波数 */
	int fwin_y;		// 周波数ウィンドウのフィルタ係数 */

	double dv0;	// ゲート電極に与えるステップ電圧
	double dv1;	// ドレイン電極に与えるステップ電圧
	int fin;	// 周波数インデックスの配列数

	char vgc[128];		//   ゲートステップ電圧入力の電流値ファイル名 */
	char vdc[128];		// ドレインステップ電圧入力電流値ファイル名 */
	char yparam[128];	// Yパラメーターのファイル名 */
	char output[128];	// アウトプットのファイル名 */


	// --- スペクトルの計算 ---
	int jtn0;	// 全時間ステップ
	int numn;	// 計算対象の時間ステップ数（電圧を加えた後のステップ数）
	int fnum;
	int epp;	// 超粒子1個あたりの電子数

	int way;	//0:数値積分法，1:FFT法
		// --- FFT法のときに必要な設定
		int data;	// FFT計算として参照するデータ数
		int jbit;	// DATA = 2^JBIT　と設定すること
		int twin;	// スペクトルウィンドウ、0:利用せず，1:利用する

		int saved;	// ドレイン雑音に対しての平滑化の重み付け数(行わないときは1)
		int appr_sd;	// ドレイン雑音スペクトル密度[Sid = 一定]、0:利用せず，1:利用する
		double fmin_sd;	// 一定とする周波数の範囲 - 最小値
		double fmax_sd;	// 一定とする周波数の範囲 - 最大値
		double fmax_sd2;	// 一定とする周波数の範囲 - 最大値

		int saveg;	// ゲート雑音に対しての平滑化の重み付け数(行わないときは1)
		int appr_sg;	// ゲート雑音スペクトル密度[Sig ∝ f^2]、0:利用せず，1:利用する
		double fmin_sg;	// 一定とする周波数の範囲 - 最小値
		double fmax_sg;	// 一定とする周波数の範囲 - 最大値
		double fmax_sg2;	// 一定とする周波数の範囲 - 最大値

		int savegd;	// ドレイン・ゲート相互に対しての平滑化(行わないときは1)
		int appr_sgd;	// 相互スペクトル[Sigid ∝ f]、0:利用せず，1:利用する
		double fmin_sgd;	// 一定とする周波数の範囲 - 最小値
		double fmax_sgd;	// 一定とする周波数の範囲 - 最大値
		double fmax_sgd2;	// 一定とする周波数の範囲 - 最大値


		// --- 数値積分法のときに必要な設定 --- 
		int msmalls;	// 相関関数の時間ステップm ---自己相関関数
		int msmallc;	// 相関関数の時間ステップm ---相互相関関数
		int mlarge;		// M (MLARGE > MSMALL)


	// --- 低域通過フィルタの設定 ---
	int lowpass;	// 0:行わない，1:行う
	double cutfreq_s;	// カットオフ周波数 */
	int fwin_s;		// 周波数ウィンドウのフィルタ係数 */

	char vc[256];		// 雑音解析用電流データの入力ファイル名 */
	char spe[256];		// 雑音電流スペクトル密度出力ファイル名 */
	int cor_out;
	int current_out;
	int fouricomp_out;
	char cor[256];		// 相関関数出力ファイル名 */
	char current[256];		// 電流平均・雑音電流平均 */
	char fouricomp[256];		// ノイズ指数出力ファイル名 */

	/*---------------------------------*/

	// --- 雑音解析 --- 
	int nway;	// ノイズ指数計算式 0:Gonzalez，1:Cappy
	char noise[256];		// ノイズ指数出力ファイル名 */


	double pai;	// 円周率π
	double bk;	// ボルツマン定数
	int ta;		// 温度

	double wg;	// 雑音指数円の計算を行うゲート幅


	int p_stop;
};


