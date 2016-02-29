function y=BFGS()
clc;
clear;
syms x_1 x_2 x_3 x_4 lam;
fun=(x_1-1)^2+(x_3-1)^2+100*(x_2-x_1^2)^2+100*(x_4-x_3^2)^2;
x0=[-1.2;1;-1.2;1];
x1=[0;0;0;0];
epsilon=0.000000001;
count=0;
h=eye(4);

g0=gra(fun,x0(1),x0(2),x0(3),x0(4));
s=-h*g0;
a=x0+lam*s;
f=subs(fun,[x_1,x_2,x_3,x_4],[a(1),a(2),a(3),a(4)]);
Mini=nd(f,0.000001);
x1=x0+Mini*s;
g1=gra(fun,x1(1),x1(2),x1(3),x1(4));
flag=norm(g1);

while(flag>=epsilon)
    delta_x=x1-x0;
    delta_g=g1-g0;
    u=1+delta_g'*h*delta_g/(delta_x'*delta_g);
    h=h+(u*(delta_x)*delta_x'-h*delta_g*delta_x'-delta_x*delta_g'*h)/(delta_x'*delta_g);
    x0=x1;
    g0=g1;
    s=-h*g1;
    a=x1+lam*s;
    f=subs(fun,[x_1,x_2,x_3,x_4],[a(1),a(2),a(3),a(4)]);
    Mini=nd(f,0.0001);
    x1=x0+Mini*s;
    g1=gra(fun,x1(1),x1(2),x1(3),x1(4));
    flag=norm(g1);
    count=count+1;
    fprintf('迭代次数：%d  ，结束条件：%.10f\n',count,flag)
    fprintf('步长：%.10f\n',Mini)
    fprintf('x值：')
    fprintf('%.10f  ',x1)
    fprintf('\n\n')
end
f=subs(fun,[x_1,x_2,x_3,x_4],[x1(1),x1(2),x1(3),x1(4)]);
f=double(f);
fprintf('最优值：%.10f  ',f)

function g=gra(gf,x1,x2,x3,x4)  %计算梯度
syms x_1 x_2 x_3 x_4;
grad=gradient(gf);
g=subs(grad,[x_1,x_2,x_3,x_4],[x1,x2,x3,x4]);
g=double(g);

function Mi=nd(f,eps)
syms lam;
yijiedao=diff(f);
erjiedao=diff(yijiedao);
x=1;
flag=subs(yijiedao,lam,x);
flag=double(flag);
while(abs(flag)>=eps)
    a=subs(yijiedao,lam,x);
    a=double(a);
    b=subs(erjiedao,lam,x);
    b=double(b);
    x=x-a/b;
    flag=subs(yijiedao,lam,x);
    flag=double(flag);
end
Mi=x;