/* --------------------------------------------------------------------
	電流データ読み込み・計算処理（Yパラメータ計算用）
	2004/7/26 by Masahiro Nakayama 
	Ver. 2.0 (cur_proc_y2_0.cpp)
-------------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "new.h"

#include "parameter2_0.h"
#include "common_function2_0.h"
#include "y_parameter2_0.h"


/* **************************************************** */
/* ----- FILEから電流値を読み込み，等価回路特性用 ----- */
void ReadCurrent_Y(struct Current cur[], int ii, char *fn, struct SetPara para) {
	char buffer[2048];
	int jt;
    FILE *fp;

    fp = OpenFile(fn, 'r');

	/* ----- 1行ずつ変数に代入 ----- */
	for(jt=0; jt<=para.jt-1; jt++) {
		ReadLine(&buffer[0], fp);
		sscanf(&buffer[0], " %lf %lf %lf ", &cur[jt]._[1][ii],
			&cur[jt]._[0][ii], &cur[jt]._[2][ii]);
//							InAs用
//		sscanf(&buffer[0], " %lf %lf %lf ", &cur[jt]._[1][ii],
//			&cur[jt]._[2][ii], &cur[jt]._[0][ii]);
	}
	
	fclose(fp);
	printf("'%s' read\n", fn);

}

/* ************************************************************* */
/* 電流のふらつき、定常状態の電流値を計算、等価回路特性用------- */
void CurInit_Y(struct Current cur[], struct Current *ss, struct Current *rex, 
			   struct SetPara para){
	int jt, i, ii;
	struct Current *temp;
	double rs[3], rr[3];

	/* ----- 読み込んだパラメータに従い、配列のメモリを確保 */
	_set_new_handler(error);
	temp = new struct Current[para.jt];

	/* -----------------------------------------------------*/
	/* --- 定常状態における電流値の計算 */
	for(ii=0; ii<=1; ii++){
	for(i=0; i<=1; i++){
		/* --- ステップ入力前の平均電流値の計算 */
		for(jt=0; jt<=para.jt-1; jt++){
			temp[jt]._[i][ii] = cur[jt]._[i][ii];
		}
		ss->_[i][ii] = 0.0;	rs[i] = 0.0;
		for(jt=para.jtp0-para.nump; jt<=para.jtp0-1; jt++){
			Sigma(cur[jt]._[i][ii], &((*ss)._[i][ii]), &rs[i]);
		}
		ss->_[i][ii] =  ss->_[i][ii]/double(para.nump);

		/* --- ステップ入力後の平均電流値の計算 */
		rex->_[i][ii] = 0.0;	rr[i] = 0.0;
		for(jt=para.jt-para.numl; jt<=para.jt-1; jt++){
			Sigma(cur[jt]._[i][ii], &((*rex)._[i][ii]), &rr[i]);
		}
		rex->_[i][ii] = rex->_[i][ii]/double(para.numl);

		// --- 電流のふらつきを計算
		for(jt=0; jt<=para.jtl0-1; jt++){
			cur[jt]._[i][ii] = temp[jt+para.jtp0]._[i][ii] - rex->_[i][ii];
		}
	}
	}

	/* --------------------*/
	/* --- Ig=0に置き換える*/
	if(para.cur0==1){
		ss->_[0][0] = 0.0;	ss->_[0][1] = 0.0;
		rex->_[0][0] = 0.0;	rex->_[0][1] = 0.0;
	}

	/* ----- メモリの開放 */
	delete [] temp;

}