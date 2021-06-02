#include <iostream>
#include <ctime>
#include "MyGA.h"

double f(int x) 
{//目标函数
	return 15*x-x*x;
}

int main() {
	srand((unsigned)time(0));//产生随机数种子
	int x;
	double max;
	MyGA<int> A(500,500,0.7,0.001,0,15,f);
	A.GA(x, max);
	system("pause");
	return 0;
}