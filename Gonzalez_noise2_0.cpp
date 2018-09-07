/* ----------------------------------------------------
	ノイズ指数の計算
	2004/7/26 by Masahiro Nakayama 
	Ver. 2.0 (Gonzalez_noise2_0.cpp)
------------------------------------------------------- */
#include <stdio.h>
#include <math.h>

#include "parameter2_0.h"
#include "common_function2_0.h" 
#include "y_parameter2_0.h"

/* ******************************************** */
/* --- ノイズ指数計算 Gonzalez提案の式による--- */
void MinimumNF_Gonzalez(struct Spectrum si[], 
				struct YPara y[], 
				double p[], double r[], double c[], double fmin[],
				double fr[], double gn[], struct ComplexNumbers zopt[], 
				struct ComplexNumbers zc[],double gass[],struct CircuitPara circ[], struct SetPara para){

	int jf;
	double y11_abs, y21_abs, sig_abs, sid_abs;
	double rn, gcor, x, w, z, a, b, gopt;
	double re_yopt,im_yopt,re_ya,im_ya,re_yb,im_yb,re_yc,re_yd,re_ye; 			/*Gass計算用パラメータ*/
//	double c2;

	for(jf=0; ; jf++){
		if(fr[jf]==32168){break;} /*最大周波数で計算をストップ */

		y11_abs = sqrt(y[jf]._11.re*y[jf]._11.re + y[jf]._11.im*y[jf]._11.im);
		y21_abs = sqrt(y[jf]._21.re*y[jf]._21.re + y[jf]._21.im*y[jf]._21.im);

		sig_abs = si[jf].g;
		sid_abs = si[jf].d;

		p[jf] = sid_abs / (4*para.bk*para.ta*y21_abs) * para.wg;
		r[jf] = sig_abs*y21_abs / (4*para.bk*para.ta*y11_abs*y11_abs) * para.wg;
		c[jf] = si[jf].gd.im / sqrt(sig_abs*sid_abs) * para.wg;


		rn = sid_abs / (4*para.bk*para.ta*y21_abs*y21_abs);   
		gcor = y[jf]._11.re - (y[jf]._21.re*si[jf].gd.re + y[jf]._21.im*si[jf].gd.im)/sid_abs;	/*符号変更*/

		x =  (y[jf]._11.re*y[jf]._21.re + y[jf]._11.im*y[jf]._21.im)*si[jf].gd.re;
		w = (-y[jf]._11.im*y[jf]._21.re + y[jf]._11.re*y[jf]._21.im)*si[jf].gd.im;

		z = x+w;
		a = y21_abs*y21_abs*sig_abs/sid_abs + y11_abs*y11_abs - 2*z/sid_abs;

		b = y[jf]._11.im + (y[jf]._21.re*si[jf].gd.im - y[jf]._21.im*si[jf].gd.re)/sid_abs;
		b = b*b;
	
		gopt = sqrt(a-b);
//Gonzaletzの方法でのNFminの求め方
		fmin[jf] = 1+2*rn*(gcor+gopt) * para.wg;

//Gonzaletzの方法と違うNFminの求め方
//	    c2    =pow(1-c[jf],2.0); 
//		fmin[jf] = 1 + 2*fr[jf]/circ[jf].ft*sqrt(p[jf]*r[jf]*c2);
//dBへ変換
		fmin[jf] = 10*log10(fmin[jf]);

	//Gassの計算(Yopt=Gopt+jBopt:Monte Carlo analysis of the dynamic behavior of III-V MOSFETs for low-noise RF applicatons)
//		re_yopt = 2*para.pai*fr[jf]*(circ[jf].cgs+circ[jf].cgd)*( (c[jf] * sqrt(r[jf]/p[jf]) )-1);		//参考にしている論文とGopt、Boptが逆！！
//		im_yopt = 2*para.pai*fr[jf]*(circ[jf].cgs+circ[jf].cgd)*sqrt(r[jf]/p[jf])*sqrt(1-c[jf]*c[jf]);
		re_yopt = 2*para.pai*fr[jf]*(circ[jf].cgs+circ[jf].cgd)*sqrt(r[jf]/p[jf])*sqrt(1-c[jf]*c[jf]);	//wCgg(root(R/P)*root(1-C^2))
		im_yopt = 2*para.pai*fr[jf]*(circ[jf].cgs+circ[jf].cgd)*( (c[jf] * sqrt(r[jf]/p[jf]) )-1);		//wCgg(C*root(R/P)-1)

		re_ya = y[jf]._21.re*y[jf]._12.re - y[jf]._21.im*y[jf]._12.im;		//Re{Y12Y21}
		im_ya = y[jf]._21.re*y[jf]._12.im + y[jf]._12.re*y[jf]._21.im;		//Im{Y12Y21}
		re_yb = re_yopt + y[jf]._11.re;							//Re{Yopt+Y11}
		im_yb = im_yopt + y[jf]._11.im;							//Im{Yopt+Y11}

		re_yc = (re_ya*re_yb + im_ya*im_yb) / (re_yb*re_yb +im_yb*im_yb);	// Y21Y12 / (Yopt+Y11)

	    re_yd = y[jf]._22.re - re_yc;			//Re[Y22-{Y21Y12/(Yopt+Y11)]

		re_ye= fabs( y21_abs*y21_abs / (re_yb*re_yb +im_yb*im_yb));	//|Y21/(Yopt+Y11)|**2

		gass[jf]= fabs((re_yopt / re_yd) * re_ye);
		gass[jf] = 10*log10(gass[jf]);

		/* 現状では求めることができない"zopt,gn"には"0"を代入 */
		zopt[jf].re = 0.0;
		zopt[jf].im = 0.0;
		zc[jf].re = 0.0;
		zc[jf].im = 0.0;
		gn[jf] = 0.0;
	}
}
