function y=chengzi()
clc;
clear;
syms x_1 x_2 lambda;
count=1;
r=0.25;
alph=2;
epsilon=0.001;
c=2;
x=[2;2];
v=[x_1;x_2];
fun=4*x_1-x_2^2-12;
h=25-x_1^2-x_2^2;
g1=10*x_1-x_1^2+10*x_2-x_2^2-34;
g2=x_1;
g3=x_2;

lambda=0;
u1=0;
u2=0;
u3=0;

sym_11=0;
sym_22=0;
sym_33=0;

sg1=subs(g1,v,x);
sg2=subs(g2,v,x);
sg3=subs(g3,v,x);
sh0=subs(h,v,x);
sg1=double(sg1);
sg2=double(sg2);
sg3=double(sg3);
sh0=double(sh0);

sym_1=u1+c*g1;
sym_2=u2+c*g2;
sym_3=u3+c*g3;

panduan_a=subs(sym_1,v,x);
panduan_b=subs(sym_2,v,x);
panduan_c=subs(sym_3,v,x);
if(panduan_a<0)
    sym_11=sym_1;
end
if(panduan_b<0)
    sym_22=sym_2;
end
if(panduan_c<0)
    sym_33=sym_3;
end
fai=fun+lambda*h+(c/2)*(h^2)+1/(2*c)*((sym_11)^2-u1^2)+((sym_22)^2-u2^2)+((sym_33)^2-u3^2);
x=BFGS(fai,0.001);

sg1=subs(g1,v,x);
sg2=subs(g2,v,x);
sg3=subs(g3,v,x);
sh1=subs(h,v,x);
sg1=double(sg1);
sg2=double(sg2);
sg3=double(sg3);
sh1=double(sh1);
flag=sh1^2+(min([sg1;-(u1/c)]))^2+(min([sg2;-(u2/c)]))^2+(min([sg3;-(u3/c)]))^2;

while(flag>=epsilon^2)
    sym_11=0;
    sym_22=0;
    sym_33=0;
    u1=min([0;u1+c*sg1]);
    u2=min([0;u2+c*sg2]);
    u3=min([0;u3+c*sg3]);
    sym_1=u1+c*g1;
    sym_2=u2+c*g2;
    sym_3=u3+c*g3;
    lambda=lambda+c*sh1;
    panduan_a=double(subs(sym_1,v,x));
    panduan_b=double(subs(sym_2,v,x));
    panduan_c=double(subs(sym_3,v,x));
    if(panduan_a<0)
        sym_11=sym_1;
    end
    if(panduan_b<0)
        sym_22=sym_2;
    end
    if(panduan_c<0)
        sym_33=sym_3;
    end
    fai=fun+lambda*h+(c/2)*(h^2)+1/(2*c)*((sym_11)^2-u1^2)+((sym_22)^2-u2^2)+((sym_33)^2-u3^2);
    sh0=sh1;
    x=BFGS(fai,0.001);
    sg1=subs(g1,v,x);
    sg2=subs(g2,v,x);
    sg3=subs(g3,v,x);
    sh1=subs(h,v,x);
    sg1=double(sg1);
    sg2=double(sg2);
    sg3=double(sg3);
    sh1=double(sh1);
    flag=sh1^2+(min([sg1;-(u1/c)]))^2+(min([sg2;-(u2/c)]))^2+(min([sg3;-(u3/c)]))^2;
    if((abs(sh1)/abs(sh0))>r)
        c=alph*c;
    end
    count=count+1;
    fprintf('迭代次数：%d  ，结束条件：%.10f  ,c=%d\n\n',count,flag,c)
    fprintf('等式约束乘子：\nlambda=%f\n\n不等式约束乘子：\nu1=%f\nu2=%f\nu3=%f\n\n',lambda,u1,u2,u3)
    fprintf('x值：')
    fprintf('%.10f  ',x)
    fprintf('\n\n\n')
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