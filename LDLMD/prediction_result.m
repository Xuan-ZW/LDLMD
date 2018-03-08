function prediction_result
    load 'similarity_array.mat';
    load 'dataset1';
    prediction=[];
    [row,col]=size(S);
    for i=1:col
        prediction=[prediction;S(:,i)];
    end
    predictionlncRNA=[];
    for i=1:col
        predictionlncRNA=[predictionlncRNA,1:row];
    end
    prediction(:,2)=predictionlncRNA';
    predictiondisease=[];
    for i=1:col
        predictiondisease=[predictiondisease;i.*ones(row,1)];
    end
    prediction(:,3)=predictiondisease;
    result=sortrows(prediction,-1);
    save result;
    xlswrite('disease.xls',[1:length(disease)]','A1:A383');
    xlswrite('disease.xls',disease,'B1:B383');
    xlswrite('lncRNA.xls',[1:length(lncRNA)]','A1:A1127');
     xlswrite('lncRNA.xls',lncRNA,'B1:B1127');
    [a1,~,a]=xlsread('lncRNA');
    [b1,~,b]=xlsread('disease');
     m=length(result);
    
    l=zeros(m,1);
    d=zeros(m,1);
    ll={};
    dd={};
    for i=1:m
        if(~isempty(find(a1==result(i,2))))
            l(i,1)=find(a1==result(i,2));
            ll{i,1}=a{l(i,1),2};
        else
            ll{i,1}=result(i,2);
        end
    end
    
    for i=1:m
        if(~isempty(find(b1==result(i,3))))
            d(i,1)=find(b1==result(i,3));
            dd{i,1}=b{d(i,1),2};
        else
            dd{i,1}=result(i,3);
        end
    end
xlswrite('LMD final prediction result.xlsx',(1:m)','A1:A431641');
xlswrite('LMD final prediction result.xlsx',ll,'B1:B431641');
xlswrite('LMD final prediction result.xlsx',dd,'C1:C431641');
    
end