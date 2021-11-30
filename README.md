# CSCI_596_FINAL_PROJECT

Intro
===

Fast Fourier Transform (FFT) is an efficient and fast algorithm for calculating polynomial multiplication. With the development of computer operation speed, the basic method of FFT, Discrete Fourier Algorithm (DFT), has attracted widespread attention. Indeed, when FFT calculation is used, the speed of polynomial multiplication can be greatly accelerated. However, with the emergence of some applications with greater computational demands, how to further optimize the FFT has become an urgent problem to be solved. We will study several optimization schemes for calculation based on the simplest polynomial multiplication, including the use of some regular conclusions to optimize the calculation process, the use of mathematical formula derivation to change the calculation thinking, and the use of parallel knowledge to execute the calculation process in multiple threads. We plan to add pthread multi-threaded calculation to it and refer to some more advanced paper research results, which have a theory of FFT Optimize planning and analysis.


Multithreaded FFT
===

After completing the optimized FFT algorithm, we need to consider how to implement it in multiple threads. Because when the amount of data is large, the multi-threaded implementation of FFT is very helpful to quickly complete the required calculations. 

In this experiment, we use Pthread to complete multithreading. First, we need to consider how to allocate tasks to each thread. For FFT, because of its particularity, when calculating one of $a[k] = x + y$, there is no need to consider the influence when other threads calculate the same $a[k]$. This provides great convenience for the distribution of our thread tasks. So we use step as a different parameter for each thread, which is equivalent to each thread executing a different divide and conquer layer of FFT, and add the calculated results to $a[k]$.

We implemented the Pthread-based FFT algorithm according to the above ideas. Here we show the processing part of each thread as follows:

```
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

Result
===

| | N=16 | N=256 | N=2048 | N=16384 | N=262144 |
| -------- | -------- | -------- | -------- | -------- | -------- |
| Basic Polynomial | 0.012718     | 0.260923 | 14.618655 | 1074.527939 | inf |
| FET_OP | 1.172512 | 3.718561 | 38.532891 | 378.675382 | 8098.466381 |
| FET_PTHREAD(2) | 1.186056 | 1.922460 | 12.603477 | 174.951241 | 1729.979298 |
| FET_PTHREAD(4) | 3.309536 | 1.869127 | 7.257020 | 76.756450 | 1400.424697 |
