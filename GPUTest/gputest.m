clc;
clear;
A = magic(5000);
f = ones(1,20),20;

tic;
B = filter(f, 1, A);
tCPU = toc;
disp(['Total time on CPU:    ' num2str(tCPU)])

tic;
AonGPU = gpuArray(A);
BonGPU = filter(f, 1, AonGPU);
BonGPU = gather(BonGPU);
wait(gpuDevice)
tGpu   = toc;
disp(['Total time on GPU:     ' num2str(tGpu)])
toc;
tic;
MATLABpool open gputest;
spmd
A = magic(5000);
f = ones(1,20),20;
AonGPU = gpuArray(A);
BonGPU = filter(f, 1, AonGPU);
BonGPU = gather(BonGPU);
wait(gpuDevice)
tGpu  = toc;
disp(['Parallel  time on GPU:     ' num2str(tGpu)])
end
matlabpool close
toc;