#pragma once
#include <iostream>
#include <ctime>

const double pi = 3.14159265;
const int len = 22;
typedef struct node {//染色体结构体
	bool chromo[len];
}node;

template <class T>
class MyGA {
public:
	MyGA(int n_generation, int Size, double pcross,double pmutate,double lower,double upper,double (*f)(T) )
	{
		_n_generation = n_generation;
		_Size = Size;
		//_bestval = bestval;
		_pcross = pcross;
		_pmutate = pmutate;
		_lower = lower;
		_upper = upper;
		bestchromo = NULL;
		group = new node[Size];
		temp = new node[Size];
		_f = f;
	}

	void GA(T& x, double& number);

private:
	void create(node& c);
	void decode(node& c, T& x);
	double fitness(node& c);
	void cross(node& c1, node& c2, int point);
	void mutate(node& c);
	double inline rand0() {//产生0到1的随机小数
		return rand() % 10000 / 10000.0;
	}
	void select();
	int getBest(T& x, double& number);
	
	double _upper;//解区间上界
	double _lower;//解区间下界
	int _n_generation;//进化代数
	int _Size;//种群规模
	double _bestval;//适应值最大值
	double _pcross;//交叉概率
	double _pmutate;//变异概率

	node* bestchromo;//记录最优个体
	node* group;//记录种群中的个体的数组
	node* temp;//记录种群中的个体的临时数组
	double (*_f)(T);
};



using namespace std;

template <class T>
void MyGA<T>::create(node& c) {//对单个染色体随机赋值
	srand((unsigned)time(0));//产生随机数种子
	for (int i = 0; i < len; i++) {
		c.chromo[i] = rand() % 2;
	}
}

template <class T>
void MyGA<T>::decode(node& c, T& x) {//二进制解码操作
	int num = pow(2,len);//即2的22次方
	double tem = 0;
	int n = 1;
	for (int i = 0; i < len; i++) {
		tem += c.chromo[i] * n;
		n *= 2;
	}
	x = (_upper-_lower)*(tem/ num)-_lower;
}




template <class T>
double MyGA<T>::fitness(node& c) {//适应度函数
	T x;
	decode(c, x);
	return _f(x);
}

template <class T>
void MyGA<T>::cross(node& c1, node& c2, int point) {//交叉操作
	node c3 = c1;
	for (int i = 0; i < len - point; i++) {
		c1.chromo[point + i] = c2.chromo[point + i];
	}
	for (int j = 0; j < len - point; j++) {
		c2.chromo[point + j] = c3.chromo[point + j];
	}
}

template <class T>
void MyGA<T>::mutate(node& c) {//变异操作
	int i = rand() % len;
	c.chromo[i] = !c.chromo[i];
}


template <class T>
void MyGA<T>::select() {//选择操作
	double* fitnessval=new double[_Size];
	double sum = 0;
	double* avgfitness= new double[_Size];
	int* id = new int[_Size];
	for (int i = 0; i < _Size; i++) {
		fitnessval[i] = fitness( group[i] );
		sum += fitnessval[i];//适应度总和
	}

	for (int i = 0; i < _Size; i++) {
		avgfitness[i] = fitnessval[i] / sum;
	}
	for (int i = 1; i < _Size; i++) {//适应度累加
		avgfitness[i] += avgfitness[i - 1];
	}
	for (int i = 0; i < _Size; i++) {//轮盘赌选择法
		double rannum = rand0();//产生0到1随机数
		int j;
		for (j = 0; j < _Size - 1; j++) {
			if (rannum < avgfitness[j]) {
				id[i] = j;
				break;
			}
		}
		if (j == _Size - 1) {
			id[i] = j;
		}
	}
	for (int i = 0; i < _Size; i++) {//将新个体替换旧个体
		temp[i] = group[i];
	}
	for (int i = 0; i < _Size; i++) {
		group[i] = temp[id[i]];
	}
}

template <class T>
int MyGA<T>::getBest(T& x, double& number) 
{//取得最优个体对应的位置
	double* fitnessval=new double[_Size];
	for (int i = 0; i < _Size; i++) {
		fitnessval[i] = fitness(group[i]);
	}
	int max_id = 0;
	for (int i = 1; i < _Size; i++) {
		if (fitnessval[i] > fitnessval[max_id]) {
			max_id = i;
		}
	}
	decode(group[max_id], x);
	number = _f(x);
	return max_id;
}

template <class T>
void MyGA<T>::GA(T& x, double& number) {//遗传算法流程
	for (int i = 0; i < _Size; i++) {
		create(group[i]);
	}


	bestchromo = &group[getBest(x, _bestval)];
	for (int i = 0; i < _n_generation; i++) {
		select();//选择操作
		int p = rand() % len;
		for (int j = 0, pre = -1; j < _Size; j++) {//根据概率交叉		
			if (rand0() < _pcross) {
				if (pre == -1)
					pre = j;
				else {
					cross(group[pre], group[j], p);
					pre = -1;
				}
			}
		}
		for (int k = 0, pre = -1; k < _Size; k++) {//根据概率进行变异
			if ((rand0() < _pmutate)) {
				mutate(group[k]);
			}
		}
		getBest(x, number);
		cout << "第" << i+1 << "代" << "最优x值为:" << x << "函数值为" << _f(x) << endl;//结果的输出
	}
}
