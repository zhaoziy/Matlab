function y=chengfa()
clc;
clear;
syms x_1 x_2;
count=1;
v=[x_1;x_2];
fun=4*x_1-x_2^2-12;
h=x_1^2+x_2^2-25;
g1=-10*x_1+x_1^2-10*x_2+x_2^2+34;
g2=-x_1;
g3=-x_2;

epsilon=0.001;
M=0.5;
c=2;
x=[1;2];
sg1=subs(g1,v,x);
sg1=double(sg1);
if(sg1>0)
   fg1=g1^2; 
else
   fg1=0; 
end
sg2=subs(g2,x_1,x(1));
sg2=double(sg2);
if(sg2>0)
    fg2=g2^2; 
else
    fg2=0;
end
sg3=subs(g3,x_2,x(2));
sg3=double(sg3);
if(sg3>0)
    fg3=g3^2; 
else
    fg3=0;
end
g_plus=h^2+fg1+fg2+fg3;
F=fun+M*g_plus;
x=BFGS(F,0.001);
t1=abs(subs(h,v,x));
t1=double(t1);
sg1=subs(g1,v,x);
sg2=subs(g2,x_1,x(1));
sg3=subs(g3,x_2,x(2));
sg=[sg1;sg2;sg3];
t2=max(sg);
t2=double(t2);
if(t1>t2)
    t=t1;
else
    t=t2;
end
while(t>=epsilon)
    M=c*M;
    sg1=subs(g1,v,x);
    sg1=double(sg1);
    if(sg1>0)
       fg1=g1^2; 
    else
       fg1=0; 
    end
    sg2=subs(g2,x_1,x(1));
    sg2=double(sg2);
    if(sg2>0)
        fg2=g2^2; 
    else
        fg2=0;
    end
    sg3=subs(g3,x_2,x(2));
    sg3=double(sg3);
    if(sg3>0)
        fg3=g3^2; 
    else
        fg3=0;
    end
    g_plus=h^2+fg1+fg2+fg3;
    F=fun+M*g_plus;
    x=BFGS(F,0.001);
    t1=abs(subs(h,v,x));
    t1=double(t1);
    sg1=subs(g1,v,x);
    sg2=subs(g2,x_1,x(1));
    sg3=subs(g3,x_2,x(2));
    sg=[sg1;sg2;sg3];
    t2=max(sg);
    t2=double(t2);
    if(t1>t2)
        t=t1;
    else
        t=t2;
    end
    count=count+1;
    fprintf('迭代次数：%d  ，结束条件：%.10f\n',count,t)
    fprintf('M：%d\n',M)
    fprintf('x值：')
    fprintf('%.10f  ',x)
    fprintf('\n\n')
end
f=subs(fun,v,x);
f=double(f);
fprintf('最优值：%.10f  ',f)


function y=BFGS(fun,epsilon)
syms x_1 x_2 lam;
x0=[10;10];
x1=[0;0];
count=0;
h=eye(2);

g0=gra(fun,x0(1),x0(2));
s=-h*g0;
a=x0+lam*s;
f=subs(fun,[x_1,x_2],[a(1),a(2)]);
Mini=nd(f,0.001);
x1=x0+Mini*s;
g1=gra(fun,x1(1),x1(2));
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
    f=subs(fun,[x_1,x_2],[a(1),a(2)]);
    Mini=nd(f,0.001);
    x1=x0+Mini*s;
    g1=gra(fun,x1(1),x1(2));
    flag=norm(g1);
    count=count+1;
end
y=x1;

function g=gra(gf,x1,x2)  %计算梯度
syms x_1 x_2;
grad=gradient(gf);
g=subs(grad,[x_1,x_2],[x1,x2]);
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
Mi=double(x);