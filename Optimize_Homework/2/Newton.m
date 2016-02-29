function y=Newton()
clc;
clear;
syms x_1 x_2 x_3 x_4;
fun=(x_1-1)^2+(x_3-1)^2+100*(x_2-x_1^2)^2+100*(x_4-x_3^2)^2;
x=[-1.2;1;-1.2;1];
epsilon=0.000000001;
count=0;
g=gra(fun,x(1),x(2),x(3),x(4));
flag=norm(g);
while (flag>=epsilon) 
    hes=hess(fun,x(1),x(2),x(3),x(4));
    x=x-inv(hes)*g;
    x=double(x);
    g=gra(fun,x(1),x(2),x(3),x(4));
    flag=norm(g);
    count=count+1;
    fprintf('迭代次数：%d  ，结束条件：%.10f\n',count,flag)
    fprintf('x值：')
    fprintf('%.10f  ',x)
    fprintf('\n\n')
end
f=subs(fun,[x_1,x_2,x_3,x_4],[x(1),x(2),x(3),x(4)]);
f=double(f);
fprintf('最优值：%.10f  ',f)

function hes=hess(f,x1,x2,x3,x4)
syms x_1 x_2 x_3 x_4;
hes=hessian(f);
hes=subs(hes,{x_1,x_2,x_3,x_4},[x1,x2,x3,x4]);
hes=double(hes);

function g=gra(gf,x1,x2,x3,x4)  %计算梯度
syms x_1 x_2 x_3 x_4;
grad=gradient(gf);
g=subs(grad,{x_1,x_2,x_3,x_4},[x1,x2,x3,x4]);
g=double(g);