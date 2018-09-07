/* ----------------------------------------------
	�p�����[�^�ꗗ
	2002/11/13  by Masahiro Nakayama 
	Ver. 1.0 (parameter1_0.h)
---------------------------------------------- */
#define JTL0 80000 /* �X�e�b�v�d������������̎��ԃX�e�b�v */
#define JTP0 40000 /* �X�e�b�v�d����������O�̎��ԃX�e�b�v */
#define JT 120000   /* JTL0+JTP0 --- Simulation���s���S���ԃX�e�b�v  */
#define NUMS 70000 /* �v�Z�Ώۂ̎��ԃX�e�b�v���i�d������������̃X�e�b�v���j */
#define NUMSBEF 20000 /* �v�Z�Ώۂ̎��ԃX�e�b�v���i�d����������O�̃X�e�b�v���j */
#define EPP 1280000 // �����q1������̓d�q��
#define DT 2.5e-15 // ��t�F1�X�e�b�v�̎���
#define CUR0 0 //0:Ig=0�AIs=Id�Ɖ���C1:���L����Ȃ��AIg=Is-Id

#define IN0 3 /* �d�� */
#define DV0 0.125	// �Q�[�g�d�ɂɗ^����X�e�b�v�d��
#define DV1 0.5		// �h���C���d�ɂɗ^����X�e�b�v�d��


/*************************************************/
/* --- Y�p�����[�^����уX�y�N�g���̌v�Z���@ --- */
#define WAY 0 //0:���l�ϕ��@�C1:FFT�@

/* --- ���l�ϕ��@�̂Ƃ��ɕK�v�Ȑݒ� --- */
// --- Y�p�����[�^�[�v�Z
#define YAVE 4 // Y�p�����[�^�[�̕��������s���Ƃ��́@�d�ݕt���̃X�e�b�v��
//�X�y�N�g���v�Z
#define MSMALLS 1700 /* ���֊֐��̎��ԃX�e�b�vm ---���ȑ��֊֐� */
#define MSMALLC 1700 /* ���֊֐��̎��ԃX�e�b�vm ---���ݑ��֊֐� */
#define MLARGE 65000 /* M (MLARGE > MSMALL) */
#define STEP 1 /* �v�Z�ɗp����d���l�A0�Fvgc(�Q�[�g�X�e�b�v)�A1�Fvdc(�h���C���X�e�b�v)*/

/* --- FFT�@�̂Ƃ��ɕK�v�Ȑݒ� --- */
// FFT�v�Z
#define DATA  16384 /* FFT�v�Z�Ƃ��ĎQ�Ƃ���f�[�^�� */
#define JBIT  14 /* DATA = 2^JBIT�@�Ɛݒ肷�邱�� */
//�X�y�N�g���v�Z
#define SWIN 0 // �X�y�N�g���E�B���h�E�A0:���p�����C1:���p����
#define SAVE 1 // ���������s���Ƃ��́@�d�ݕt���̃X�e�b�v��


/*************************************************/
/*---------------------------------*/
#define PAI 3.1415926535897 /* �~������ */
#define BK 1.38066e-23 /* �{���c�}���萔 */
#define TA 300 /* ���x */

#define VGC "vgc.txt" /*   �Q�[�g�X�e�b�v�d���̓��̓t�@�C���� */
#define VDC "vdc.txt" /* �h���C���X�e�b�v�d���̓��̓t�@�C���� */
#define YPARAM "Yparam.txt" /* Y�p�����[�^�[�̃t�@�C���� */
#define OUTPUT "output.txt" /* �A�E�g�v�b�g�̃t�@�C���� */
#define SPE "CurSpectrum.txt" /*   �o�̓t�@�C���� */
#define COR "correlation.txt" /*   �o�̓t�@�C���� */
#define NIOSE "noise.txt" /*   �m�C�Y�w���o�̓t�@�C���� */


/* ----- ���f���\���̂��߂̍\���� -----*/
struct ComplexNumbers{
  double re;
  double im;
};
