### ImportError: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.26‘ not found

```shell
# find whether this so.6 in your conda env
strings /home/app/anaconda3/lib/libstdc++.so.6 | grep GLIBCXX
# if so
export LD_LIBRARY_PATH=/home/haotian/software/miniconda3/lib
export LD_LIBRARY_PATH=/homeb/chenyu/miniconda/envs/ecloud/lib
```

### CUDA error: CUBLAS STATUS_INVALID VALUE when calling cublasGemmEx( handle ...

```shell
# Once I tried to use hugging-face transformer and pytorch-lightning to train something, this bug occurred again and again. Meanwhile, the model was deseperated from the cuda device once I straded training. Then, I found the problem is that my torch version is not compatible with the cuda version. Annoying, right? 
```

### After the power outage, the CUDA driver could not be found.

(1) Firstly, find the NVIDIA-Linux-x86 64-515.65.01.run at the ./Downloads

(2) `sudo bash NVIDIA-Linux-x86 64-515.65.01.run`

(3) Continue installation twice

(4) Q: Install NVIDIA's 32-bit compatibility libraries. A: No 

(5) Q: Would you like to run the nvidia-xconfig utility to automatically update your X configuration file xxx. A: No.
