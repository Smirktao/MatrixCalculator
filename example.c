#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MatrixCalculator.h"
#include <time.h>

#define M_ROW 1 //矩阵行数
#define M_COLUMN 3 //矩阵列数



int main()
{
	//定义系统参数
	float	a1 = 0.5, a2 = 0.5, uk0 = 1, q = 0.024;
	float	y1 = 0.1, y0 = 0.2;
	int M = 80, N = 50; 

	//手动定义噪声 用于程序测试
	float d[80] = { -1.0149 , -0.4711, 0.1370,-0.2919, 0.3018, 0.3999,-0.9300,-0.1768,-2.1321, 1.1454,-0.6291,-1.2038,-0.2539,-1.4286,-0.0209,-0.5607, 2.1778, 1.1385,-2.4969, 0.4413,-1.3981,-0.2551, 0.1644, 0.7477,-0.2730, 1.5763,-0.4809, 0.3275, 0.6647, 0.0852, 0.8810, 0.3232,-0.7841,-1.8054, 1.8586,-0.6045, 0.1034, 0.5632, 0.1136,-0.9047,-0.4677,-0.1249, 1.4790,-0.8608, 0.7847, 0.3086,-0.2339,-1.0570,-0.2841,-0.0867,-1.4694, 0.1922,-0.8223,-0.0942, 0.3362,-0.9047,-0.2883, 0.3501,-1.8359, 1.0360, 2.4245, 0.9594,-0.3158, 0.4286,-1.0360, 1.8779, 0.9407, 0.7873,-0.8759, 0.3199,-0.5583,-0.3114,-0.5700,-1.0257,-0.9087,-0.2099,-1.6989, 0.6076,-0.1178, 0.6992 };
	
	//创建空矩阵 参数：（行数，列数）
	Matrix_t w0 = create_mat(M_ROW, M_COLUMN); 
	Matrix_t w1 = create_mat(M_ROW, M_COLUMN);
	Matrix_t x = create_mat(M_ROW, M_COLUMN);//注意 这里的x为了简化C代码，从matlab里的3行1列 改成了1行3列，所以下面运算的时候不再需要转置
											 
	//用数组定义矩阵初始值
	float	w0_value[M_ROW*M_COLUMN] = { 0.05,0.2,0.1 }; 
	float	w1_value[M_ROW*M_COLUMN] = { 0.2,-0.02,0.01 };
	float	x_value[M_ROW*M_COLUMN] = {	1,1,1};

	//把初始值赋值给矩阵
	set_mat_data(&w0, w0_value);
	set_mat_data(&w1, w1_value);
	set_mat_data(&x, x_value);

	//定义循环计算过程需要的中间变量
	int k = 0, i = 0;
	float d1 = 0, yp =	0;
	float uk = 0,bI = 0;
	float y2 = 0, y8 = 0,e1=0;
	Matrix_t w1_sub_w0 = create_mat(M_ROW, M_COLUMN);//w1-w0
	Matrix_t x_T = create_mat(M_ROW, M_COLUMN);
	Matrix_t w2 = create_mat(M_ROW, M_COLUMN);

	//定义结果变量
	Matrix_t y3 = create_mat(1, M);
	Matrix_t y4 = create_mat(1, M);
	Matrix_t yp1 = create_mat(1, M);
	Matrix_t w = create_mat(1, M);
	Matrix_t e = create_mat(1, M);

	for (k = 0; k < M; k++)
	{
		d1 = d[k];
		yp = 1- exp(-1 * (k+1) / 2.0); //因为C循环中K是从0到(M-1),而matlab是从1到M,所以要(K+1)
		for (i = 0; i < N; i++)
		{
			bI = dot_product(&w1, &x) - 0.1;
			uk = 1.0 / (1 + exp(-1 * bI));
			y2 = 0.33*y1 + 0.132*y0 + 0.5*uk + 0.038*uk0;
			y8 = 0.33*y1 + 0.132*y0 + 0.51*uk + 0.031*uk0 + q*d1;
			y0 = y1; y1 = y2; uk0 = uk;
			e1 = yp - y8;

			//下面几行分步计算 w2=w1+0.682*e1*x'+0.113*(w1-w0)
			sub_mat_2(&w1_sub_w0,&w1, &w0);//(w1-w0)
			scale_mat_2(&w1_sub_w0, 0.113);// 0.113*(w1-w0)
			copy_mat_data( &x_T, &x);//由于x改成了1行3列，所以不再需要转置，直接复制
			scale_mat_2(&x_T, 0.682*e1);// x' = 0.682*e1*x'
			add_mat_2(&w2, &w1, &x_T);// w2 = w1+ x'
		    add_mat_2(&w2,&w2, &w1_sub_w0); // w2 = w2 + w1_sub_w0 
			//至此 计算完毕 w2=w1+ 0.682*e1*x' + 0.113*(w1-w0) 

			copy_mat_data( &w0, &w1);//w0=w1;
			copy_mat_data( &w1, &w2);//w1=w2;
			
			if (e1 <= 0.001)
			{
				break;
			}
		}
		//因为C的数组下标从0开始，所以C语言里的数组第0行相当于matlb的数组第1行
		y3.data[0][k] = y2;
		y4.data[0][k] = y8;
		yp1.data[0][k] = yp;
		w.data[0][k] = q*d1;
		e.data[0][k] = e1;
	}
	show_mat("y3", &y3);
	show_mat("y4", &y4);
	//销毁所有矩阵
	free_mat(&w0);
	free_mat(&w1);
	free_mat(&x_T);
	free_mat(&w2);
	free_mat(&y3);
	free_mat(&y4);
	free_mat(&yp1);
	free_mat(&w);
	free_mat(&e);
	if (sizeof(d) == 4)
	{
		//如果d的占用内存是4，说明d是指针，需要释放
		free(d);
	}
	getchar();
	return 0;
}

