function y=gonge()
clc;
clear;
syms lam;
epsilon=0.0000000000001;
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
lam_s=pingfen(lamb,0,1 ,0.000001);
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
        lam_s=pingfen(lamb,0,1,0.000001);
        a=subs(a,lam,lam_s);
        a=double(a);
        g0=g1;
        g1=gra(f,a,n);
        flag=norm(g1);
        count=count+1;
        count_sum=count_sum+1;
        fprintf('迭代次数：%d  ，结束条件：%.10f\n',count_sum,flag)
        fprintf('步长：%.10f\n',lam_s)
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
    lam_s=pingfen(lamb,0,1,0.000001);
    a=subs(a,lam,lam_s);
    a=double(a);
    g1=gra(f,a,n);
    flag=norm(g1);
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