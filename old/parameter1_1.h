/* ----------------------------------------------
	�p�����[�^�ꗗ
	2003/3/31 by Masahiro Nakayama 
	Ver. 1.1 (parameter1_1.h)
---------------------------------------------- */
#define DT 5.0e-15 // ��t�F1�X�e�b�v�̎���
#define EPP 32000 // �����q1������̓d�q��
#define FIN 10000 /* ���g���C���f�b�N�X�̔z�� */

/* ********************************************** */
// --- Y�p�����[�^�v�Z
#define JTL0 12000 /* �X�e�b�v�d������������̎��ԃX�e�b�v */
#define JTP0 8000 /* �X�e�b�v�d����������O�̎��ԃX�e�b�v */
#define JT 20000  /* JTL0+JTP0 --- Simulation���s���S���ԃX�e�b�v  */
#define NUML 6000 /* �v�Z�Ώۂ̎��ԃX�e�b�v���i�d������������̃X�e�b�v���j */
#define NUMP 2000 /* �v�Z�Ώۂ̎��ԃX�e�b�v���i�d����������O�̃X�e�b�v���j */
#define CUR0 1 //������H�@�@0:Ig=0�AIs=Id�Ɖ���C1:���L����Ȃ��AIg=Is-Id
#define YAVE 30 // Y�p�����[�^�[�̕��������s���Ƃ��́@�d�ݕt���̃X�e�b�v��
#define LOWPASS_Y 1 //0:�s��Ȃ��C1:�s��

/* ********************************************** */
// --- �G����� --- 
#define JTN0 300000 /* �S���ԃX�e�b�v */
#define NUMN 285000 /* �v�Z�Ώۂ̎��ԃX�e�b�v���i�d������������̃X�e�b�v���j */
#define CUR1 1 //�G����́@�@0:Ig=0�AIs=Id�Ɖ���C1:���L����Ȃ��AIg=Is-Id

// --- �X�y�N�g���̌v�Z ---
#define WAY 1 //0:���l�ϕ��@�C1:FFT�@
	// --- ���l�ϕ��@�̂Ƃ��ɕK�v�Ȑݒ� --- 
	#define MSMALLS 20000 /* ���֊ւ̎��ԃX�e�b�vm ---���ȑ��֊֐� */
	#define MSMALLC 20000 /* ���֊֐��̎��ԃX�e�b�vm ---���ݑ��֊֐� */
	#define MLARGE 270000 /* M (MLARGE > MSMALL) */

	// --- FFT�@�̂Ƃ��ɕK�v�Ȑݒ�
	#define DATA 262144 /* FFT�v�Z�Ƃ��ĎQ�Ƃ���f�[�^�� */
	#define JBIT 18 /* DATA = 2^JBIT�@�Ɛݒ肷�邱�� */
	#define TWIN 1 // �X�y�N�g���E�B���h�E�A0:���p�����C1:���p����
	#define FETAPP 1 // �ASid���,Sig f^2,Im[Sigd],f�A0:���p�����C1:���p����
	#define SAVE 40 // ���������s���Ƃ��́@�d�ݕt���̃X�e�b�v��(0�̂Ƃ��͍s��Ȃ�)
	#define FMINAVE 3.0e+10 //2���Ȑ��ߎ��Ɏg�����g�� - �ŏ��l
	#define FMAXAVE 7.0e+10 //2���Ȑ��ߎ��Ɏg�����g�� - �ő�l

// --- �G���w����� ---
#define NWAY 1 // �m�C�Y�w���v�Z�� 0:Gonzalez�C1:Cappy

// --- ���ʉ߃t�B���^�̐ݒ� ---
#define LOWPASS 0 //0:�s��Ȃ��C1:�s��
#define CUTFREQ 2.0e+11 /* �J�b�g�I�t���g�� */
#define FWIN 1000 /* ���g���E�B���h�E�̃t�B���^�W�� */

/*************************************************/
/*---------------------------------*/
#define IN0 3 /* �d�� */
#define DV0 0.125	// �Q�[�g�d�ɂɗ^����X�e�b�v�d��
#define DV1 0.5		// �h���C���d�ɂɗ^����X�e�b�v�d��

#define VGC "vgc.txt" /*   �Q�[�g�X�e�b�v�d�����͂̓d���l�t�@�C���� */
#define VDC "vdc.txt" /* �h���C���X�e�b�v�d�����͓d���l�t�@�C���� */
#define VC "vc.txt" /* �G����͗p�d���f�[�^�̓��̓t�@�C���� */
#define YPARAM "Yparam.txt" /* Y�p�����[�^�[�̃t�@�C���� */
#define OUTPUT "output.txt" /* �A�E�g�v�b�g�̃t�@�C���� */
#define SPE "CurSpectrum.txt" /*   �o�̓t�@�C���� */
#define COR "correlation.txt" /*   �o�̓t�@�C���� */
#define NIOSE "noise.txt" /*   �m�C�Y�w���o�̓t�@�C���� */

#define PAI 3.1415926535897 /* �~������ */
#define BK 1.38066e-23 /* �{���c�}���萔 */
#define TA 300 /* ���x */

/* ----- ���f���\���̂��߂̍\���� -----*/
struct ComplexNumbers{
  double re;
  double im;
};
