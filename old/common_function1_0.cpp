/* ---------------------------------------------------------------
	Main�֐��Ŏg�����ʓI�Ȋ֐�
	2002/11/13  by Masahiro Nakayama 
	Ver. 1.0 (common_function1_0.cpp)
----------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter1_0.h" 
#include "common_function1_0.h" 

/* ************************************************ */
/* ----- FILE����d���l��ǂݍ��݁A�ϐ��ɑ�� ----- */
void ReadCurrent(double cur[][JT][IN0], int ii, char *fp_name) {
	char buffer[128][64];
	int jt;
    FILE *fp;

    fp = OpenFile(fp_name, 'r');

	/* ----- 1�s���ϐ��ɑ�� ----- */
	for(jt=0; jt<=JT-1; jt++) {
		ReadLine(&buffer[0][0], fp);
		sscanf(&buffer[0][0], " %lf %lf %lf ", &cur[1][jt][ii],
			&cur[0][jt][ii], &cur[2][jt][ii]);
	}
	fclose(fp);
	printf("'%s' read\n", fp_name);

}

/* ********************************************** */
/* ----- ���߂������g���̃C���f�b�N�X�𐶐� ----- */
void FrequencyIndex(double fr[]){
	int jf;

	/* ------------------------------ */
	/* --- ���l�ϕ��@��p�����ꍇ --- */
	/* ------------------------------ */
	if(WAY==0){
		int fnum;
		double df, fmin, fmax;
		df   = 5.0e+9;		//���g���Ԋu
		fmin = -2.0e+10;	//�ŏ����g��
		fmax = 1.2e+11;		//�ő���g��
		fnum = (int)((fmax-fmin)/df); //���g���C���f�b�N�X�̌�

		for(jf=0; jf<=fnum; jf++){
			fr[jf] = df * double(jf) + fmin;
		}

		fr[fnum+1] = 32168; //�z��̏I��������
	}

	/* ------------------------- */
	/* --- FFT�@��p�����ꍇ --- */
	/* ------------------------- */
	if(WAY==1){
		for(jf=0; jf<=DATA/2-1; jf++){
			fr[jf] = jf/DT/DATA;
		}
		fr[DATA/2] = 32168; //�z��̏I��������
	}
}


/* *********************************************** */
/* �d���̂ӂ���A����Ԃ̓d���l���v�Z ------- */
void CurInit(double cur[][JT][IN0], double cur1[][JT][IN0], double ss[][IN0], 
			 double rex[][IN0], char *fp_name){
	int i, ii, jt;
	double rs1, rs2, rr1, rr2;

    FILE *fp;
    fp = OpenFile(fp_name, 'w');
	fprintf(fp, "Ig-vg-bef	Id-vg-bef	Ig-vg-aft	Id-vg-aft	");
	fprintf(fp, "Ig-vd-bef	Id-vd-bef	Ig-vd-aft	Id-vd-aft\n");

	// --- ����Ԃɂ�����d���l�̌v�Z
	for(ii=0; ii<=1; ii++){
		 ss[0][ii] = 0.0;	 ss[1][ii] = 0.0;
		rex[0][ii] = 0.0;	rex[1][ii] = 0.0;
		rs1 = 0.0;	rs2 = 0.0;
		rr1 = 0.0;	rr2 = 0.0;
		for(jt=JTP0-NUMSBEF; jt<=JTP0-1; jt++){
			Sigma(cur[0][jt][ii], &ss[0][ii], &rs1);
			Sigma(cur[1][jt][ii], &ss[1][ii], &rs2);
		}
		for(jt=JT-NUMS; jt<=JT-1; jt++){
			Sigma(cur[0][jt][ii], &rex[0][ii], &rr1);
			Sigma(cur[1][jt][ii], &rex[1][ii], &rr2);
		}
		 ss[0][ii] =  ss[0][ii]/double(NUMSBEF);
		 ss[1][ii] =  ss[1][ii]/double(NUMSBEF);
		rex[0][ii] = rex[0][ii]/double(NUMS);
		rex[1][ii] = rex[1][ii]/double(NUMS);
		fprintf(fp, "%lf	%lf	",  ss[0][ii],  ss[1][ii]);
		fprintf(fp, "%lf	%lf	", rex[0][ii], rex[1][ii]);
	}
	fprintf(fp, "\n");
	fclose(fp);


	/* -----------------------------------------------------*/
	/* --- �P����Ig=Is-Id�ɂ��������Čv�Z*/
	/* -----------------------------------------------------*/
	if(CUR0==1){
		// --- �d���̂ӂ�����v�Z
		for(i=0; i<=IN0-1; i++){
			for(ii=0; ii<=1; ii++){
				for(jt=0; jt<=JTL0; jt++){
					cur1[i][jt][ii] = cur[i][jt+JTP0][ii] - rex[i][ii];
				}
			}
		}
	}


	/* -----------------------------------------------------*/
	/* --- Ig=0�A�܂�Is=Id�Ɖ��肵�Čv�Z*/
	/* -----------------------------------------------------*/
	if(CUR0==0){
		// --- ����Ԃɂ�����d���l�̌v�Z
		for(ii=0; ii<=1; ii++){
			ss[1][ii] = 0.0;	rs1 = 0.0;
			for(jt=JTP0-NUMSBEF; jt<=JTP0-1; jt++){
				Sigma(cur[1][jt][ii]-cur[2][jt][ii], &ss[1][ii], &rs1);
			}
			ss[0][ii] = 0.0;
			ss[1][ii] =  ss[1][ii]/double(NUMSBEF)/2;
			ss[2][ii] = -ss[1][ii];
		}

		for(ii=0; ii<=1; ii++){
			rex[1][ii] = 0.0;	rex[0][ii] = 0.0;
			rr1 = 0.0;	rr2 = 0.0;
			for(jt=JT-NUMS; jt<=JT-1; jt++){
				Sigma(cur[0][jt][ii], &rex[0][ii], &rr1);
				Sigma(cur[1][jt][ii]-cur[2][jt][ii], &rex[1][ii], &rr2);
			}
			rex[0][ii] =  rex[0][ii]/double(NUMS);
			rex[1][ii] =  rex[1][ii]/double(NUMS)/2;
			rex[2][ii] = -rex[1][ii];
		}

		// --- �d���̂ӂ�����v�Z
		for(i=0; i<=IN0-1; i++){
			for(ii=0; ii<=1; ii++){
				for(jt=0; jt<=JTL0; jt++){
					cur1[i][jt][ii] = cur[i][jt+JTP0][ii] - rex[i][ii];
				}
			}
		}
		rex[0][0] = 0.0;
		rex[0][1] = 0.0;
	}


}

/* ****************************** */
/* ----- �t�@�C�����I�[�v�� ----- */
FILE* OpenFile(char *argv, char mode) {
   FILE *fp;
   if((fp = fopen(argv, &mode)) == NULL) {
           fprintf(stderr, "Error: File Not Found (%s)\n", argv);
           return 0;
   }
   return fp;
}

/* ************************ */
/* ----- ��s�ǂݍ��� ----- */
int ReadLine(char *buffer, FILE *fp) {
   int i=0;
   if(fgets(buffer, 256, fp) == NULL) {
           if(ferror(fp)) {
                   fprintf(stderr, "Error: Cannot Read File");
                   return 0;
           }
   if(feof(fp)) i = 1;
   }
   return i;
}

/* ********************************** */
/* ----- Siguma�v�Z(���l�v�Z�p) ----- */
void Sigma(double x,double *sum, double *reg)
{
	double temp;
	*(reg)  = *(reg) + x;
	temp = *(sum);
	*sum  = *(sum) + *(reg);
	temp = *(sum) - temp;
	*(reg)  = *(reg) - temp;
}
