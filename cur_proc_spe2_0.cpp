/* --------------------------------------------------------------------
	電流データ読み込み・計算処理-スペクトル計算用-(分解能をあげるため三角ウィンドー→ハニング窓へ変更)
	2004/7/26 by Takuto Takkahashi
	Ver. 2.0 (cur_proc_spe2_0.cpp)
-------------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter2_0.h"
#include "common_function2_0.h"
#include "current_spectrum2_0.h" 

/* **************************************************** */
/* ----- FILEから電流値を読み込み、雑音特性解析用 ----- */
void ReadCurrent_Spe(struct Current cur[], struct SetPara para) {
	char buffer[2048];
	int jt;
    FILE *fp;

    fp = OpenFile(para.vc, 'r');

	for(jt=0; jt<=para.jtn0-1; jt++) {
		ReadLine(&buffer[0], fp);
		sscanf(&buffer[0], " %lf %lf %lf ", &cur[jt]._[1],
			&cur[jt]._[0], &cur[jt]._[2]);
//							InAs用
//			sscanf(&buffer[0], " %lf %lf %lf ", &cur[jt]._[1],
//						&cur[jt]._[2], &cur[jt]._[0]);

//		cur[jt]._[0]=cur[jt]._[0]*1.0;		//*gate width→[A^2/m]
//		cur[jt]._[1]=cur[jt]._[1]*1.0;
//		cur[jt]._[2]=cur[jt]._[2]*1.0;
	}

	fclose(fp);
	printf("'%s' read\n", para.vc);

}


/* ************************************************************* */
/* 電流のふらつき、定常状態の電流値を計算、雑音特性解析用------- */
void CurInit_Spe(struct Current cur[], struct SetPara para){
	int jt;
	double rex[3], rr1, rr2, rr3;
	double di_s,di_g,di_d;

	/* --- 定常状態における電流値の計算 */
	rex[0] = 0.0;	rex[1] = 0.0;	rex[2] = 0.0;
	rr1 = 0.0;	rr2 = 0.0;	rr3 = 0.0;
	di_s =0.0; di_g =0.0; di_d =0.0;
	for(jt=para.jtn0-para.numn; jt<=para.jtn0-1; jt++){
		Sigma(cur[jt]._[0], &rex[0], &rr1);
		Sigma(cur[jt]._[1], &rex[1], &rr2);
		Sigma(cur[jt]._[2], &rex[2], &rr3);
	}
	rex[0] = rex[0]/double(para.numn);
	rex[1] = rex[1]/double(para.numn);
	rex[2] = rex[2]/double(para.numn);

	/* --- 電流のふらつきを計算 */
	for(jt=0; jt<=para.jtn0-1; jt++){
		cur[jt]._[0] = cur[jt]._[0] - rex[0];
		cur[jt]._[1] = cur[jt]._[1] - rex[1];
		cur[jt]._[2] = cur[jt]._[2] - rex[2];
	}
	/* ---　電流揺らぎの計算 --- */
	for(jt=para.jtn0-para.numn; jt<=para.jtn0-1; jt++){
		di_g += pow(cur[jt]._[1],2.0);
		di_d += pow(cur[jt]._[0],2.0);
		di_s += pow(cur[jt]._[2],2.0);

//		di_g = di_g / para.numn/para.epp;
//		di_d = di_d / para.numn/para.epp;
//		di_s = di_s / para.numn/para.epp;
		di_g = di_g /(para.jtn0-1-para.numn);
		di_d = di_d /(para.jtn0-1-para.numn);
		di_s = di_s /(para.jtn0-1-para.numn);
	}
	/* --- 平均電流の書き込み */
	if(para.current_out ==1){
		FILE *fp;
		fp = OpenFile(para.current, 'w');
		fprintf(fp, "Id	Ig	Is\n");
		fprintf(fp, "%1.4e	%1.4e	%lf\n", rex[1], rex[0], rex[2]);
		fprintf(fp, "ΔId^2	ΔIg^2	ΔIs^2\n");
		fprintf(fp, "%1.4e	%1.4e	%1.4e\n", di_d, di_g, di_s);
		fclose(fp);
	}

}
