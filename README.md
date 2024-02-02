## Paper
Source code of the DATE '24 paper: "Efficient Spectral-Aware Power Supply Noise Analysis for Low-Power Design Verification" by Yinuo Bai, Xiaoyu Yang, Yicheng Lu, Dan Niu, Cheng Zhuo, Zhou Jin, Weifeng Liu

## Complier
>**g++ ./include/*.h ./parallel/*.hpp ./serial/*.hpp ./test/test.cpp -w -O3 -o main -fopenmp -D_GLIBCXX_PARALLEL**

## Test
>**./main /MTX_PATH/MTX.mtx SPARSE_RATIO THREADS_NUM** \
>**./main /MTX_PATH/MTX.mtx 0.02 8**

## Contact us
If you have any questions about running the code, please contact Xiaoyu Yang. \
E-mail: yangxiaoyu@student.cup.edu.cn