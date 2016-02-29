function y=gonge()
clc;
clear;
syms lam;
epsilon=0.000000001;
count=0;
count_sum=0;
n=100;

for i=1:n
    x(i,1)=sym(['x',num2str(i)]);
end
sum=0;
for i=1:n
    sum=sum+x(i);
end
f=sum^2-sum;
a=zeros(n,1);

g0=gra(f,a,n);
s=-g0;
a=a+lam*s;
lamb=subs(f,x,a);
lam_s=nd(lamb ,0.0001);
a=subs(a,lam,lam_s);
a=double(a);
g1=gra(f,a,n);
flag=norm(g1);
while(flag>=epsilon)
    while(flag>=epsilon)
        if(count<n-1)
            u=((norm(g1))^2)/((norm(g0))^2);
            s=-g1+u*s;
        elseif(count==(n-1))
            count=0;
            break;
        end
        a=a+lam*s;
        lamb=subs(f,x,a);
        lam_s=nd(lamb ,0.0001);
        a=subs(a,lam,lam_s);
        a=double(a);
        g0=g1;
        g1=gra(f,a,n);
        flag=norm(g1);
        count=count+1;
        count_sum=count_sum+1;
        fprintf('迭代次数：%d  ，结束条件：%.10f\n',count_sum,flag)
        printf('步长：%.10f\n',lam_s)
        fprintf('x值：')
        fprintf('%.10f  ',a)
        fprintf('\n\n')
    end
    if(flag<epsilon)
        break;
    end
    g0=gra(f,a,n);
    s=-g0;
    a=a+lam*s;
    lamb=subs(f,x,a);
    lam_s=nd(lamb ,0.0001);
    a=subs(a,lam,lam_s);
    a=double(a);
    g1=gra(f,a,n)
    flag=norm(g1)
    count_sum=count_sum+1;
    fprintf('迭代次数：%d  ，结束条件：%.10f\n',count_sum,flag)
    fprintf('步长：%.10f\n',lam_s)
    fprintf('x值：')
    fprintf('%.10f  ',a)
    fprintf('\n\n')
end
f=subs(f,x,a);
f=double(f);
fprintf('最优值：%.10f  ',f)

function g=gra(gf,a,n)  %计算梯度
for i=1:n
x(i,1)=sym(['x',num2str(i)]);
end
grad=gradient(gf);
g=subs(grad,x,a);
g=double(g);

function Mi=nd(f,eps)
syms lam;
yijiedao=diff(f);
erjiedao=diff(yijiedao);
x=0.5;
flag=subs(yijiedao,lam,x);
flag=double(flag);
count1=0;
while(abs(flag)>=eps)
    a=subs(yijiedao,lam,x);
    b=subs(erjiedao,lam,x);
    x=x-a/b;
    flag=subs(yijiedao,lam,x);
    flag=double(flag)
    count1=count1+1
end
Mi=x;