function overallauc=positiontooverallauc1(beta,k1)
LOOCV(beta,k1);
load newglobalposition.mat;
load LD3;
newLD1=LD3;
load L_M.mat;
load M_D.mat;
load newinteraction.mat
[n,m]=size(newinteraction);
[n1,m1]=size(L_M);
[n2,m2]=size(M_D);
[pp,qq]=size(newLD1);


for i=1:pp
if newglobalposition(i)>n*m-pp+1
newglobalposition(i)=n*m-pp+1;
end
end
for k=1:n*m-pp+1
    tp=0;
    for t=1:pp
        if newglobalposition(1,t)<=k
            tp=tp+1;
        end
    end
    tpr(1,k)=tp/pp;
    fp=k*pp-tp;
    fpr(1,k)=fp/(pp*(n*m-pp));
end
if (k1==2)
    plot(fpr,tpr,'r');
else if(k1==3)
        plot(fpr,tpr,'b');
    else if(k1==4)
            plot(fpr,tpr,'g');
        end
    end
end

hold on
clear area;
area(1,1)=tpr(1,1)*fpr(1,1)/2;
for k=2:n*m-pp+1
    area(1,k)=[tpr(1,k-1)+tpr(1,k)]*[fpr(1,k)-fpr(1,k-1)]/2;
end
overallauc=sum(area);
end    
