function LOOCV(bata, k)
%bata=0.1;
load LD3
load LD4
newLD1=LD3;
L_M_D=LD4;
nlmd=length(L_M_D);
numld=length(newLD1);
L=unique(L_M_D(:,1));
M=unique(L_M_D(:,2));
D=unique(L_M_D(:,3));
nl=length(L);
nm=length(M);
nd=length(D);
L_M=zeros(nl,nm);
M_D=zeros(nm,nd);
newinteraction=zeros(nl,nd);
newglobalposition=[];
for i=1:nl
    indexL=find(strcmp(L_M_D(:,1),L(i)));
    for j=1:nd
        if(ismember(L_M_D(indexL,3),D(j)))
            newinteraction(i,j)=1;
        end
    end
end
save('newinteraction.mat','newinteraction');
for i=1:nl
    index1=find(strcmp(L_M_D(:,1),L(i)));
    miRNA=L_M_D(index1,2);
    for j=1:nm
        for t=1:length(index1)
            index2=find(strcmp(M,miRNA(t)));
            L_M(i,index2)=1;
        end
    end  
end
save('L_M.mat','L_M');
for i=1:nm
    index1=find(strcmp(L_M_D(:,2),M(i)));
    disease=L_M_D(index1,3);
    for j=1:nm
        for t=1:length(index1)
            index2=find(strcmp(D,disease(t)));
            M_D(i,index2)=1;
        end
    end
end
save('M_D.mat','M_D');

for cv=1:numld
    tic
    cv
    load L_M;
    load M_D;
    index_L_DcvInL=find(strcmp(L,newLD1(cv,1)));
    L_M(index_L_DcvInL,:)=0;
    index_L_DcvInD=find(strcmp(D,newLD1(cv,2)));
    M_D(:,index_L_DcvInD)=0;
        %Li 和 Lj的相似性SL(i,j)
    for i=1:nl
        indexLi_M=find(ismember(L_M(i,:),1));
        t1=length(indexLi_M);
        for j=1:nl
            if (i==j)
                SL(i,j)=1;
            else
                indexLj_M=find(ismember(L_M(j,:),1));
                t2=length(indexLj_M);
                %SL(i,j)=exp(-sum((sum(pdist2(M_D(indexLi_M,:),M_D(indexLj_M,:))./nd)))*(t1*t2)/(sum(sum(M_D))));
                  %SL(i,j)=exp(-sum(((sum(pdist2(M_D(indexLi_M,:),M_D(indexLj_M,:)).^2)./nd)))*(t1*t2)/(sum(sum(M_D))));
                  SL(i,j)=exp(-sum(((sum(pdist2(M_D(indexLi_M,:),M_D(indexLj_M,:)).^2))))./nd);
                  %SL(i,j)=exp(-sum((exp(-sum((pdist2(M_D(indexLi_M,:),M_D(indexLj_M,:)).^2)./nd)))));
                  %SL(i,j)=1/(1+exp(-15*SL(i,j)+log(9999)));
                %SL(i,j)=1/(1+exp(-pdist2(L_M(i,:),L_M(j,:)).^2));
            end
        end
    end
    %Di,Dj的相似性SD(i,j)
    for i=1:nd
        indexDi_M=find(ismember(M_D(:,i),1));
        t1=length(indexDi_M);
        for j=1:nd
            if(i==j)
                SD(i,j)=1;
            else
                indexDj_M=find(ismember(M_D(:,j),1));
                t2=length(indexDj_M);
                %SD(i,j)=exp(-sum((sum(pdist2(L_M(:,indexDi_M)',L_M(:,indexDj_M)')./nl)))*(t1*t2)/sum(sum(L_M)));
                %SD(i,j)=exp(-sum(((sum(pdist2(L_M(:,indexDi_M)',L_M(:,indexDj_M)').^2)./nl)))*(t1*t2)/sum(sum(L_M)));
                SD(i,j)=exp(-sum(((sum(pdist2(L_M(:,indexDi_M)',L_M(:,indexDj_M)').^2))))./nl);
                %SD(i,j)=exp(-sum(exp(-sum((pdist2(L_M(:,indexDi_M)',L_M(:,indexDj_M)').^2)./nl))));
               % SD(i,j)=1/(1+exp(-15*SD(i,j)+log(9999)));
                %SD(i,j)=1/(1+exp(-pdist2(M_D(:,i)',M_D(:,j)').^2));
            end
        end
    end
%     S=[SL,L_M*M_D;M_D'*L_M',SD];
%     F1=bata*S+bata.^2*S*S;
%     F=F1(1:nl,nl:nl+nd);
    A=L_M*M_D;
    if (k==2)
        F=bata*A+bata.^2*(SL*A+A*SD);
        else if (k==3)
            F=bata*A+bata.^2*(SL*A+A*SD)+bata.^3*(A*A'*A+SL*SL*A+SL*A*SD+A*SD*SD);
            else if(k==4)
                F=bata*A+bata.^2*(SL*A+A*SD)+bata.^3*(A*A'*A+SL*SL*A+SL*A*SD+A*SD*SD) +...
                    bata.^4*(SL^3*A +A*A'*SL*A+ SL*A*A'*A + A*SD*A'*A)+...
                    bata.^4*(A*A'*A*SD +SL^2*A*SD+ SL*A*SD^2+ A*SD^3);
                end
            end
    end
    
    index_cv_L=find(strcmp(L,L_M_D(cv,1)));
    index_cv_D=find(strcmp(D,L_M_D(cv,3)));
    finalscore=F(index_cv_L,index_cv_D);
    for i=1:nl
        for j=1:nd
            if newinteraction(i,j)==1
                F(i,j)=-10000;
            end
        end
    end
    % obtain the position of tested disease-microbe interaction as variable globalposition(1,cv),
[ll1,mm1]=size(find(F>=finalscore));
[ll2,mm2]=size(find(F>finalscore));
newglobalposition(1,cv)=ll2+1+(ll1-ll2-1)/2;
toc
end
save('newglobalposition.mat','newglobalposition');
end
