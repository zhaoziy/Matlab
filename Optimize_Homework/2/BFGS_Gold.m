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
Mini=Gold(f,0,1,0.000001);
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
    Mini=Gold(f,0,1,0.0001);
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

function Min=Gold(f,a1,b1,eps)
syms lam;
alph=0.618;
l=a1+(1-alph)*(b1-a1);
u=a1+alph*(b1-a1);
sl=subs(f,lam,l);
sl=double(sl);
su=subs(f,lam,u);
su=double(su);
while((b1-a1)>eps)
    if(sl>su)
        a1=l;
        b1=b1;
        l=u;
        u=a1+alph*(b1-a1);
    else
        a1=a1;
        b1=u;
        u=l;
        l=a1+(1-alph)*(b1-a1);
    end
    sl=subs(f,lam,l);
    sl=double(sl);
    su=subs(f,lam,u);
    su=double(su);
end
Min=0.5*(b1+a1);