#include <iostream>
#include <ctime>
#include "MyGA.h"

double f(int x) 
{//Ŀ�꺯��
	return 15*x-x*x;
}

int main() {
	srand((unsigned)time(0));//�������������
	int x;
	double max;
	MyGA<int> A(500,500,0.7,0.001,0,15,f);
	A.GA(x, max);
	system("pause");
	return 0;
}