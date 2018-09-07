/* --------------------------------------------------------------------
	�d���f�[�^�ǂݍ��݁E�v�Z�����iY�p�����[�^�v�Z�p�j
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
/* ----- FILE����d���l��ǂݍ��݁C������H�����p ----- */
void ReadCurrent_Y(struct Current cur[], int ii, char *fn, struct SetPara para) {
	char buffer[2048];
	int jt;
    FILE *fp;

    fp = OpenFile(fn, 'r');

	/* ----- 1�s���ϐ��ɑ�� ----- */
	for(jt=0; jt<=para.jt-1; jt++) {
		ReadLine(&buffer[0], fp);
		sscanf(&buffer[0], " %lf %lf %lf ", &cur[jt]._[1][ii],
			&cur[jt]._[0][ii], &cur[jt]._[2][ii]);
//							InAs�p
//		sscanf(&buffer[0], " %lf %lf %lf ", &cur[jt]._[1][ii],
//			&cur[jt]._[2][ii], &cur[jt]._[0][ii]);
	}
	
	fclose(fp);
	printf("'%s' read\n", fn);

}

/* ************************************************************* */
/* �d���̂ӂ���A����Ԃ̓d���l���v�Z�A������H�����p------- */
void CurInit_Y(struct Current cur[], struct Current *ss, struct Current *rex, 
			   struct SetPara para){
	int jt, i, ii;
	struct Current *temp;
	double rs[3], rr[3];

	/* ----- �ǂݍ��񂾃p�����[�^�ɏ]���A�z��̃��������m�� */
	_set_new_handler(error);
	temp = new struct Current[para.jt];

	/* -----------------------------------------------------*/
	/* --- ����Ԃɂ�����d���l�̌v�Z */
	for(ii=0; ii<=1; ii++){
	for(i=0; i<=1; i++){
		/* --- �X�e�b�v���͑O�̕��ϓd���l�̌v�Z */
		for(jt=0; jt<=para.jt-1; jt++){
			temp[jt]._[i][ii] = cur[jt]._[i][ii];
		}
		ss->_[i][ii] = 0.0;	rs[i] = 0.0;
		for(jt=para.jtp0-para.nump; jt<=para.jtp0-1; jt++){
			Sigma(cur[jt]._[i][ii], &((*ss)._[i][ii]), &rs[i]);
		}
		ss->_[i][ii] =  ss->_[i][ii]/double(para.nump);

		/* --- �X�e�b�v���͌�̕��ϓd���l�̌v�Z */
		rex->_[i][ii] = 0.0;	rr[i] = 0.0;
		for(jt=para.jt-para.numl; jt<=para.jt-1; jt++){
			Sigma(cur[jt]._[i][ii], &((*rex)._[i][ii]), &rr[i]);
		}
		rex->_[i][ii] = rex->_[i][ii]/double(para.numl);

		// --- �d���̂ӂ�����v�Z
		for(jt=0; jt<=para.jtl0-1; jt++){
			cur[jt]._[i][ii] = temp[jt+para.jtp0]._[i][ii] - rex->_[i][ii];
		}
	}
	}

	/* --------------------*/
	/* --- Ig=0�ɒu��������*/
	if(para.cur0==1){
		ss->_[0][0] = 0.0;	ss->_[0][1] = 0.0;
		rex->_[0][0] = 0.0;	rex->_[0][1] = 0.0;
	}

	/* ----- �������̊J�� */
	delete [] temp;

}