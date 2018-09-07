/* ---------------------------------------------------------------
	�p�����[�^�t�@�C���ǂݍ��݊֐�
	2004/7/26 by Masahiro Nakayama 
	Ver. 2.0 (parameter2_0.cpp)
----------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "new.h"

#include "parameter2_0.h"
#include "common_function2_0.h" 

/* ******************************************* */
/* ----- �p�����[�^�̓ǂݍ��� ���C���֐� ----- */
void ReadParameter(struct SetPara *para){
	char buffer[1024], name[128], temp[128];
	FILE *fp;

	/*** --- �v���g�^�C�v�錾 --- ***/
	void SerchData_Int(char buffer[], struct SetPara *para);
	void SerchData_Double(char buffer[], struct SetPara *para);
	void SerchData_Char(char buffer[], struct SetPara *para);
	/* ================================================================= */

	/* ----- �p�����[�^�̊i�[���ꂽ�t�@�C���̃I�[�v������ */
    fp = OpenFile(PARAMETER, 'r');

	/* ********* �e�s�ɂ�����p�����[�^�̕ϐ��ǂݍ��ݏ��� *********	*/
	for(;;){
		ReadLine(&buffer[0], fp);		// ��s��ǂݍ���
		SerchData_Int(buffer, para);	// �����^�̕ϐ��Ƃ̔�r
		SerchData_Double(buffer, para);	// ���������^�̕ϐ��Ƃ̔�r
		SerchData_Char(buffer, para);	// ������^�̕ϐ��Ƃ̔�r

		/* �t�@�C���ŏI�s�̔��菈�� */
		sscanf(&buffer[0], " %s	", &name[0]);
		strcpy(temp, "end");
		if(strcmp(temp, name) == 0){break;}
	}

	fclose(fp);	// �t�@�C���N���[�Y
	printf("'%s' read\n", PARAMETER);


	/* ********* �ǂݍ��񂾃p�����[�^����v�Z�����e��p�����[�^�̏��� */
	/* ----- Y�p�����[�^�̑��X�e�b�v���̌v�Z */
	if(para->culc!=2){
		para->jt = para->jtl0 + para->jtp0;
	}

/* ----- �X�y�N�g�����x�v�Z(FFT���p��)�̑��X�e�b�v���̌v�Z */
	if(para->way==1){
		para->data = int(pow(2.0, para->jbit));	// data=2^JBIT ���v�Z
		printf(" 2^JBIT = %d\n", para->data );	// ��ʂɌv�Z���ʂ��o��

		/* �t�@�C���̍s����data�������Ă����ꍇ�̃G���[���� */
//		if(para->data > para->numn){
//			printf("Spectrum caluculation error : 2^JBIT > NUMN\n");
//			StopPG();
//			exit(1);
//		}
	}

	/* ----- ���g���C���f�b�N�X�ɕK�v�ȍs�񐔂̌v�Z */
	if(para->culc==1 || para->way == 0){	// Y�p�����[�^�̂݌v�Z�̂Ƃ�
		para->fin = int((para->fmax_y - para->fmin_y) / para->fstep_y) +
					10 + para->yave;	//
		para->fin = para->fin * 10;		//���f�s��Ȃǂ����邽�߁A���߂Ɋm��
	}
	else {	// ��L�ȊO
		para->fin = int((para->fmax_y - para->fmin_y) / (para->dt*para->data)) + 
					10 + para->yave;
		para->fin = para->fin * 10;
	}


}



/* *************************************************** */
/* ----- �p�����[�^�̓ǂݍ��� �����^�̕ϐ������� ----- */
void SerchData_Int(char buffer[], struct SetPara *para){
	char name[128];
	int value;

	/* ----- �ǂݍ��񂾈�s��ϐ��Ɋi�[ */
	sscanf(&buffer[0], " %s	%d	", &name[0], &value);

	/* ----- �ȉ��A�v���O�����Ŏg�p����ϐ��Ƃ̏ƍ� */

	/* --- ���� ---*/
	if(strcmp("CULC", name) == 0){para->culc = value; return ;}

	/* --- Y�p�����[�^�̌v�Z ---*/
	if(strcmp("JTL0", name) == 0){para->jtl0 = value; return ;}
	if(strcmp("JTP0", name) == 0){para->jtp0 = value; return ;}
	if(strcmp("NUML", name) == 0){para->numl = value; return ;}
	if(strcmp("NUMP", name) == 0){para->nump = value; return ;}
	if(strcmp("CUR0", name) == 0){para->cur0 = value; return ;}
	if(strcmp("YAVE", name) == 0){para->yave = value; return ;}
	if(strcmp("LOWPASS_Y", name) == 0){para->lowpass_y = value; return ;}
	if(strcmp("FWIN_Y", name) == 0){para->fwin_y = value; return ;}

	/* --- �X�y�N�g�����x�̌v�Z ---*/
	if(strcmp("EPP", name) == 0){para->epp = value; return ;}
	if(strcmp("FNUM", name) == 0){para->fnum = value; return ;}
	if(strcmp("JTN0", name) == 0){para->jtn0 = value; return ;}
	if(strcmp("NUMN", name) == 0){para->numn = value; return ;}
	if(strcmp("WAY", name) == 0){para->way = value; return ;}
	if(strcmp("MSMALLS", name) == 0){para->msmalls = value; return ;}
	if(strcmp("MSMALLC", name) == 0){para->msmallc = value; return ;}
	if(strcmp("MLARGE", name) == 0){para->mlarge = value; return ;}
	if(strcmp("DATA", name) == 0){para->data = value; return ;}
	if(strcmp("JBIT", name) == 0){para->jbit = value; return ;}
	if(strcmp("TWIN", name) == 0){para->twin = value; return ;}
	if(strcmp("SAVED", name) == 0){para->saved = value; return ;}
	if(strcmp("APPR_SD", name) == 0){para->appr_sd = value; return ;}
	if(strcmp("SAVEG", name) == 0){para->saveg = value; return ;}
	if(strcmp("APPR_SG", name) == 0){para->appr_sg = value; return ;}
	if(strcmp("SAVEGD", name) == 0){para->savegd = value; return ;}
	if(strcmp("APPR_SGD", name) == 0){para->appr_sgd = value; return ;}
	if(strcmp("LOWPASS", name) == 0){para->lowpass = value; return ;}
	if(strcmp("FWIN_S", name) == 0){para->fwin_s = value; return ;}
	if(strcmp("NWAY", name) == 0){para->nway = value; return ;}
	if(strcmp("COR_OUT", name) == 0){para->cor_out = value; return;}
	if(strcmp("CURRENT_OUT", name) == 0){para->current_out = value; return;}
	if(strcmp("FOURICOMP_OUT", name) == 0){para->fouricomp_out = value; return;}

	/* --- �ŏ��G���w���̌v�Z ---*/

	/* --- �G���w���~�̌v�Z ---*/

	// --- ���̑� ---
	if(strcmp("TA", name) == 0){para->ta = value; return;}
	if(strcmp("P_STOP", name) == 0){para->p_stop = value; return;}


}


/* ******************************************************* */
/* ----- �p�����[�^�̓ǂݍ��� ���������^�̕ϐ������� ----- */
void SerchData_Double(char buffer[], struct SetPara *para){
	char name[128];
	double value;

	/* ----- �ǂݍ��񂾈�s��ϐ��Ɋi�[ */
	sscanf(&buffer[0], " %s	%lf	", &name[0], &value);

	/* ----- �ȉ��A�v���O�����Ŏg�p����ϐ��Ƃ̏ƍ� */

	/* --- ���� ---*/
	if(strcmp("DT", name) == 0){para->dt = value; return;}

	/*--- Y�p�����[�^�̌v�Z ---*/
	if(strcmp("CUTFREQ_Y", name) == 0){para->cutfreq_y = value; return;}
	if(strcmp("FMIN_Y", name) == 0){para->fmin_y = value; return ;}
	if(strcmp("FMAX_Y", name) == 0){para->fmax_y = value; return ;}
	if(strcmp("FSTEP_Y", name) == 0){para->fstep_y = value; return ;}
	if(strcmp("DV0", name) == 0){para->dv0 = value; return;}
	if(strcmp("DV1", name) == 0){para->dv1 = value; return;}

	/* --- �X�y�N�g�����x�̌v�Z ---*/
	if(strcmp("FMIN_SD", name) == 0){para->fmin_sd = value; return;}
	if(strcmp("FMAX_SD", name) == 0){para->fmax_sd = value; return;}
	if(strcmp("FMAX_SD2", name) == 0){para->fmax_sd2 = value; return;}
	if(strcmp("FMIN_SG", name) == 0){para->fmin_sg = value; return;}
	if(strcmp("FMAX_SG", name) == 0){para->fmax_sg = value; return;}
	if(strcmp("FMAX_SG2", name) == 0){para->fmax_sg2 = value; return;}
	if(strcmp("FMIN_SGD", name) == 0){para->fmin_sgd = value; return;}
	if(strcmp("FMAX_SGD", name) == 0){para->fmax_sgd = value; return;}
	if(strcmp("FMAX_SGD2", name) == 0){para->fmax_sgd2 = value; return;}
	if(strcmp("CUTFREQ_S", name) == 0){para->cutfreq_s = value; return;}

	/* --- �ŏ��G���w���̌v�Z ---*/

	/* --- �G���w���~�̌v�Z ---*/
	if(strcmp("WG", name) == 0){para->wg = value; return;}

	/* --- ���̑� ---*/
	if(strcmp("PAI", name) == 0){para->pai = value; return;}
	if(strcmp("BK", name) == 0){para->bk = value;	return;}

}


/* ***************************************************** */
/* ----- �p�����[�^�̓ǂݍ��� ������^�̕ϐ������� ----- */
void SerchData_Char(char buffer[], struct SetPara *para){
	char name[128], value[128];

	/* ----- �ǂݍ��񂾈�s��ϐ��Ɋi�[ */
	sscanf(&buffer[0], " %s	%s	", &name[0], &value[0]);

	/* ----- �ȉ��A�v���O�����Ŏg�p����ϐ��Ƃ̏ƍ� */

	/* --- ���� ---*/


	/* --- Y�p�����[�^�̌v�Z ---*/
	if(strcmp("VGC", name) == 0){strcpy(para->vgc , value); return;}
	if(strcmp("VDC", name) == 0){strcpy(para->vdc , value); return;}
	if(strcmp("YPARAM", name) == 0){strcpy(para->yparam , value); return;}
	if(strcmp("OUTPUT", name) == 0){strcpy(para->output , value); return;}

	/* --- �X�y�N�g�����x�̌v�Z ---*/
	if(strcmp("VC", name) == 0){strcpy(para->vc , value); return;}
	if(strcmp("SPE", name) == 0){strcpy(para->spe , value); return;}
	if(strcmp("NIOSE", name) == 0){strcpy(para->noise , value);	return;}
	if(strcmp("COR", name) == 0){strcpy(para->cor , value);	return;}
	if(strcmp("CURRENT", name) == 0){strcpy(para->current , value);	return;}
	if(strcmp("FOURICOMP", name) == 0){strcpy(para->fouricomp , value);	return;}

	/* --- �ŏ��G���w���̌v�Z ---*/

	/* --- �G���w���~�̌v�Z ---*/

	/* --- ���̑� ---*/



}

