function [SL,SD]=similarity_calculate
   load AML;
   load AMD;
   L_M=AML';
   M_D=AMD;
  
    [nl,nm]=size(L_M);
    [nm,nd]=size(M_D);
    for i=1:nl
        index1=find(L_M(i,:));
        for j=1:nl
            if(i==j)
                SL(i,j)=1;
            else
                index2=find(L_M(j,:));
                 SL(i,j)=exp(-sum(sum(pdist2(M_D(index1,:),M_D(index2,:)).^2))./nd);
            end
        end
        
    end
    L_M=L_M';
    M_D=M_D';
    for i=1:nd
        index1=find(M_D(i,:));
        for j=1:nd
            if(i==j)
                SD(i,j)=1;
            else
                index2=find(M_D(j,:));
                SD(i,j)=exp(-sum((sum(pdist2(L_M(index1,:),L_M(index2,:)).^2))./nl));
            end
        end
        
    end
    save SD SD;
    save SL SL;
end