#pragma comment(lib, "pthreadVC2.lib")
#include<stdio.h>
#include<cmath>
#include<algorithm>
#include<cstring>
#include<complex>
#include <windows.h>
#include <immintrin.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include<iostream>

using namespace std;

typedef struct {
	int	threadId;
} threadParm_t;
const int THREAD_NUM = 2;
pthread_mutex_t	mutex_task;
long long head, freq, tail;        
int next_task = 1;

const int cal_num= 64;
typedef complex<double> cd;
const int maxl = 1024*1024;
const double PI = 3.14159265358979;
char s1[maxl], s2[maxl];
int bit = 1, s = 2;
cd a[maxl], b[maxl];
int rev[2094153];
int output[maxl];


void getrev(int bit) {
	for (int i = 0; i<(1 << bit); i++) {
		rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit - 1));
	}
}
void out_Result() {
	
	for (int i = 0; i<s; i++) {
		output[i] += (int)(a[i].real() + 0.5);
		output[i + 1] += output[i] / 10;
		output[i] %= 10;
	}
	int i;
	for (i = 2*cal_num; !output[i] && i >= 0; i--);
	if (i == -1)printf("0");
	for (; i >= 0; i--) {
		printf("%d", output[i]);
	}
}

//optimization1 arrange W and reverse process
int dft;
void multi() {
	int i, A[cal_num], B[cal_num], C[2*cal_num];
	for (i = 0; i < cal_num; i++) {
		A[i] = (int)rand() % 10;
		B[i] = (int)rand() % 10;
	}
	for (int i = 0; i<cal_num; i++)
		for (int j = 0; j<cal_num; j++)
			C[i + j] += A[i] * B[j];
}
void fft_op(cd *a, int n) {
	for (int i = 0; i<n; i++) {
		if (i < rev[i])
			swap(a[i], a[rev[i]]);
	}
	for (int step = 1; step<n; step <<= 1) {
		cd wn = exp(cd(0, dft*PI / step));
		for (int j = 0; j<n; j += step << 1) {
			cd wnk(1, 0);
			for (int k = j; k<j + step; k++) {
				cd x = a[k];
				cd y = wnk*a[k + step];
				a[k] = x + y;
				a[k + step] = x - y;
				wnk *= wn;
			}
		}
	}
	if (dft == -1) {
		for (int i = 0; i<n; i++)
			a[i] /= n;
	}
}
void *arr_attitude(void *parm)
{
	threadParm_t *p = (threadParm_t *)parm;
	int r = p->threadId;
	int step = 0;

	while (1) {
		pthread_mutex_lock(&mutex_task);
		step = next_task;
		next_task <<= 1;
		pthread_mutex_unlock(&mutex_task);
		if (step >= s) break;
		cd wn = exp(cd(0, dft*PI / step));
		for (int j = 0; j<s; j += step << 1) {
			cd wnk(1, 0);
			for (int k = j; k<j + step; k++) {
				cd x = a[k];
				cd y = wnk*a[k + step];
				a[k] = x + y;
				a[k + step] = x - y;
				wnk *= wn;
			}
		}
	}

	pthread_exit(nullptr);
	return NULL;
}
void fft_pthread(cd *a, int n) {
	pthread_t thread[THREAD_NUM];
	threadParm_t threadParm[THREAD_NUM];
	for (int i = 0; i<n; i++) {
		if (i<rev[i])
			swap(a[i], a[rev[i]]);
	}
	for (int i = 0; i < THREAD_NUM; i++)
	{
		threadParm[i].threadId = i;
		pthread_create(&thread[i], nullptr, arr_attitude, (void *)&threadParm[i]);
	}

	for (int i = 0; i < THREAD_NUM; i++)
	{
		pthread_join(thread[i], nullptr);
	}
	if (dft == -1) {
		for (int i = 0; i<n; i++)
			a[i] /= n;
	}
}
void init(cd* s) {
	int i;
	for (i = 0; i < cal_num; i++) {
		s[i] = (int)rand() % 10;
		//s[i] = i;
	}
	for (; i < maxl; i++) {
		s[i] = 0;
	}
}
int main() {
	QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
	

	init(a); init(b);
	for (bit = 1; (1 << bit) < 2 * cal_num -1; bit++) {
		s <<= 1;
	}
	getrev(bit);
	dft = 1;
	QueryPerformanceCounter((LARGE_INTEGER *)&head);
	fft_op(a, s);
	fft_op(b, s);//dft
	for (int i = 0; i<s; i++)a[i] *= b[i];
	dft = -1;
	fft_op(a, s);//idft
	QueryPerformanceCounter((LARGE_INTEGER *)&tail);
	printf("fft: %lfms.\n", (tail - head) * 1000.0 / freq);
	//out_Result();

	init(a); init(b);
	mutex_task = PTHREAD_MUTEX_INITIALIZER;
	dft = 1;
	QueryPerformanceCounter((LARGE_INTEGER *)&head);
	fft_pthread(a, s);
	fft_pthread(b, s);//dft
	for (int i = 0; i<s; i++)a[i] *= b[i];
	dft = -1;
	fft_pthread(a, s);//idft
	QueryPerformanceCounter((LARGE_INTEGER *)&tail);
	printf("fft_pthread: %lfms.\n", (tail - head) * 1000.0 / freq);
	//out_Result();
	cout << endl;



	return 0;
}