# CSCI_596_FINAL_PROJECT

Intro
===

Fast Fourier Transform (FFT) is an efficient and fast algorithm for calculating polynomial multiplication. With the development of computer operation speed, the basic method of FFT, Discrete Fourier Algorithm (DFT), has attracted widespread attention. Indeed, when FFT calculation is used, the speed of polynomial multiplication can be greatly accelerated. However, with the emergence of some applications with greater computational demands, how to further optimize the FFT has become an urgent problem to be solved. We will study several optimization schemes for calculation based on the simplest polynomial multiplication, including the use of some regular conclusions to optimize the calculation process, the use of mathematical formula derivation to change the calculation thinking, and the use of parallel knowledge to execute the calculation process in multiple threads. We plan to add pthread multi-threaded calculation to it and refer to some more advanced paper research results, which have a theory of FFT Optimize planning and analysis.

# Naive polynomial multiplication

```c++
vector<double>ForceMul(vector<double>A,vector<double>B) {
    vector<double>ans;
    int aLen=A.size();
    int bLen=B.size();、
    int ansLen=aLen+bLen-1;
    for(int i=1;i<=ansLen;i++)
        ans.push_back(0);
    for(int i=0;i<aLen;i++)
        for(int j=0;j<bLen;j++)
            ans[i+j]+=A[i]*B[j]; 
    return ans;
}

```



# Basic FFT

A **fast Fourier transform** (**FFT**) includes 2 parts:  [discrete Fourier transform](https://en.wikipedia.org/wiki/Discrete_Fourier_transform) (DFT) and  its inverse (IDFT).  In fft_op, the first part is the binary flip, which is responsible for changing the original ordering based on the rules. The second part is to use different butterfly coefficients for different layers to calculate. The third part is for idft to divide each result by N.

```c++
void fft_op(cd *a, int n) {
    // part 1
	for (int i = 0; i<n; i++) {
		if (i < rev[i])
			swap(a[i], a[rev[i]]);
	}
    // part 2
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
    // part 3
	if (dft == -1) {
		for (int i = 0; i<n; i++)
			a[i] /= n;
	}
}
```
In main function, we use fft_op to finish DFT and IDFT:

```c++
	fft_op(a, s);
	fft_op(b, s);//dft
	for (int i = 0; i<s; i++)a[i] *= b[i];
	dft = -1;
	fft_op(a, s);//idft
```




Multithreaded FFT
===

After completing the optimized FFT algorithm, we need to consider how to implement it in multiple threads. Because when the amount of data is large, the multi-threaded implementation of FFT is very helpful to quickly complete the required calculations. 

In this experiment, we use Pthread to complete multithreading. First, we need to consider how to allocate tasks to each thread. For FFT, because of its particularity, when calculating one of **a[k] = x + y**, there is no need to consider the influence when other threads calculate the same **a[k]**. This provides great convenience for the distribution of our thread tasks. So we use step as a different parameter for each thread, which is equivalent to each thread executing a different divide and conquer layer of FFT, and add the calculated results to **a[k]**.

We implemented the Pthread-based FFT algorithm according to the above ideas. Here we show the processing part of each thread as follows:

```c++
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
```

Here we use `next_task` to pass parameters, so whenever a new thread is created, the value of `next_task` will be changed to continue running down.

```c++
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
```

In main function, we can use fft_pthread to finish DFT and IDFT:

```c++
	fft_pthread(a, s);
	fft_pthread(b, s);//dft
	for (int i = 0; i<s; i++)a[i] *= b[i];
	dft = -1;
	fft_pthread(a, s);//idft
```



Result
===

| | N=16 | N=256 | N=2048 | N=16384 | N=262144 |
| -------- | -------- | -------- | -------- | -------- | -------- |
| Basic Polynomial | 0.012718     | 0.260923 | 14.618655 | 1074.527939 | inf |
| FFT_OP | 1.172512 | 3.718561 | 38.532891 | 378.675382 | 8098.466381 |
| FFT_PTHREAD(2) | 1.186056 | 1.922460 | 12.603477 | 174.951241 | 1729.979298 |
| FFT_PTHREAD(4) | 3.309536 | 1.869127 | 7.257020 | 76.756450 | 1400.424697 |

*   N is the number of terms in a polynomial
