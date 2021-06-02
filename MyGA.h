#pragma once
#include <iostream>
#include <ctime>

const double pi = 3.14159265;
const int len = 22;
typedef struct node {//Ⱦɫ��ṹ��
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
	double inline rand0() {//����0��1�����С��
		return rand() % 10000 / 10000.0;
	}
	void select();
	int getBest(T& x, double& number);
	
	double _upper;//�������Ͻ�
	double _lower;//�������½�
	int _n_generation;//��������
	int _Size;//��Ⱥ��ģ
	double _bestval;//��Ӧֵ���ֵ
	double _pcross;//�������
	double _pmutate;//�������

	node* bestchromo;//��¼���Ÿ���
	node* group;//��¼��Ⱥ�еĸ��������
	node* temp;//��¼��Ⱥ�еĸ������ʱ����
	double (*_f)(T);
};



using namespace std;

template <class T>
void MyGA<T>::create(node& c) {//�Ե���Ⱦɫ�������ֵ
	srand((unsigned)time(0));//�������������
	for (int i = 0; i < len; i++) {
		c.chromo[i] = rand() % 2;
	}
}

template <class T>
void MyGA<T>::decode(node& c, T& x) {//�����ƽ������
	int num = pow(2,len);//��2��22�η�
	double tem = 0;
	int n = 1;
	for (int i = 0; i < len; i++) {
		tem += c.chromo[i] * n;
		n *= 2;
	}
	x = (_upper-_lower)*(tem/ num)-_lower;
}




template <class T>
double MyGA<T>::fitness(node& c) {//��Ӧ�Ⱥ���
	T x;
	decode(c, x);
	return _f(x);
}

template <class T>
void MyGA<T>::cross(node& c1, node& c2, int point) {//�������
	node c3 = c1;
	for (int i = 0; i < len - point; i++) {
		c1.chromo[point + i] = c2.chromo[point + i];
	}
	for (int j = 0; j < len - point; j++) {
		c2.chromo[point + j] = c3.chromo[point + j];
	}
}

template <class T>
void MyGA<T>::mutate(node& c) {//�������
	int i = rand() % len;
	c.chromo[i] = !c.chromo[i];
}


template <class T>
void MyGA<T>::select() {//ѡ�����
	double* fitnessval=new double[_Size];
	double sum = 0;
	double* avgfitness= new double[_Size];
	int* id = new int[_Size];
	for (int i = 0; i < _Size; i++) {
		fitnessval[i] = fitness( group[i] );
		sum += fitnessval[i];//��Ӧ���ܺ�
	}

	for (int i = 0; i < _Size; i++) {
		avgfitness[i] = fitnessval[i] / sum;
	}
	for (int i = 1; i < _Size; i++) {//��Ӧ���ۼ�
		avgfitness[i] += avgfitness[i - 1];
	}
	for (int i = 0; i < _Size; i++) {//���̶�ѡ��
		double rannum = rand0();//����0��1�����
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
	for (int i = 0; i < _Size; i++) {//���¸����滻�ɸ���
		temp[i] = group[i];
	}
	for (int i = 0; i < _Size; i++) {
		group[i] = temp[id[i]];
	}
}

template <class T>
int MyGA<T>::getBest(T& x, double& number) 
{//ȡ�����Ÿ����Ӧ��λ��
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
void MyGA<T>::GA(T& x, double& number) {//�Ŵ��㷨����
	for (int i = 0; i < _Size; i++) {
		create(group[i]);
	}


	bestchromo = &group[getBest(x, _bestval)];
	for (int i = 0; i < _n_generation; i++) {
		select();//ѡ�����
		int p = rand() % len;
		for (int j = 0, pre = -1; j < _Size; j++) {//���ݸ��ʽ���		
			if (rand0() < _pcross) {
				if (pre == -1)
					pre = j;
				else {
					cross(group[pre], group[j], p);
					pre = -1;
				}
			}
		}
		for (int k = 0, pre = -1; k < _Size; k++) {//���ݸ��ʽ��б���
			if ((rand0() < _pmutate)) {
				mutate(group[k]);
			}
		}
		getBest(x, number);
		cout << "��" << i+1 << "��" << "����xֵΪ:" << x << "����ֵΪ" << _f(x) << endl;//��������
	}
}
