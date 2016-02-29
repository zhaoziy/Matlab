function y=zuisuxiajiang()
clc;
clear;
syms x_1 x_2 x_3 x_4;
fun=(x_1-1)^2+(x_3-1)^2+100*(x_2-x_1^2)^2+100*(x_4-x_3^2)^2;
x=[-1.2;1;-1.2;1];
epsilon=0.001;
count=0;
g=gra(fun,x(1),x(2),x(3),x(4));
flag=norm(g);
while (flag>=epsilon) 
    Mini=lamb(fun,x(1),x(2),x(3),x(4));
    a=x-Mini*g;
    f2=subs(fun,[x_1,x_2,x_3,x_4],[a(1),a(2),a(3),a(4)]);
    f2=double(f2);
    f1=subs(fun,[x_1,x_2,x_3,x_4],[x(1),x(2),x(3),x(4)]);
    f1=double(f1);
    if(f2>f1)
        x=x-0.01*g;
    else
        x=a;
    end
    g=gra(fun,x(1),x(2),x(3),x(4));
    flag=norm(g);
    count=count+1;
    fprintf('迭代次数：%d  ，结束条件：%.10f\n',count,flag)
    fprintf('步长：%.10f\n',Mini)
    fprintf('x值：')
    fprintf('%.10f  ',x)
    fprintf('\n\n')
end
f=subs(fun,[x_1,x_2,x_3,x_4],[x(1),x(2),x(3),x(4)]);
f=double(f);
fprintf('最优值：%.10f  ',f)

function Mini=lamb(gf,x1,x2,x3,x4)
syms x_1 x_2 x_3 x_4 lam;
grad=gradient(gf);
g=subs(grad,[x_1,x_2,x_3,x_4],[x1,x2,x3,x4]);
g=double(g);
a=[x1;x2;x3;x4]-lam*g;
f=subs(gf,[x_1,x_2,x_3,x_4],[a(1),a(2),a(3),a(4)]);
Mini=pingfen(f,0,1,0.000001);

function g=gra(gf,x1,x2,x3,x4)  %计算梯度
syms x_1 x_2 x_3 x_4;
grad=gradient(gf);
g=subs(grad,[x_1,x_2,x_3,x_4],[x1,x2,x3,x4]);
g=double(g);

function Min=pingfen(f,a1,b1,eps)
syms lam;
df=diff(f);
while((b1-a1)>=eps)
    x=(a1+b1)/2;
    sf=subs(df,lam,x);
    sf=double(sf);
    if(sf>0)
        b1=x;
    else
        a1=x;
    end
end
Min=0.5*(a1+b1);
