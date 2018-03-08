%% 导入lncRNA-disease 并进行相关处理
[~,~,lncRNA_disease1]=xlsread('L_D.xls','Sheet1');
[~, ~, lncRNA_disease2] = xlsread('newLD201704091.xls','A1:B2048');
lncRNA_disease=[lncRNA_disease1;lncRNA_disease2];
lncRNA_disease(2742,:)=[];

lncRNA_disease=lower(lncRNA_disease);
lncRNA_disease(:,2)=regexprep(lncRNA_disease(:,2),'carcinoma','cancer');
lncRNA_disease(:,2)=regexprep(lncRNA_disease(:,2),'neoplasms','cancer');
lncRNA_disease(:,2)=regexprep(lncRNA_disease(:,2),'tumor','cancer');
lncRNA_disease(:,2)=regexprep(lncRNA_disease(:,2),'cancers','cancer');
lncRNA_disease(:,2)=regexprep(lncRNA_disease(:,2),',','');
lncRNA_disease(:,2)=regexprep(lncRNA_disease(:,2),'\''s','');

% 构造lncRNA_disease 关联矩阵
N=size(lncRNA_disease,1);
lncRNA3=unique(lncRNA_disease(:,1));
disease3=unique(lncRNA_disease(:,2));

for i=1:N   
    lncRNA4=lncRNA_disease{i,1};
    disease4=lncRNA_disease{i,2};
    index1=find(strcmp(lncRNA3,lncRNA4));
    index2=find(strcmp(disease3,disease4));
    LD(i,1)=index1;
    LD(i,2)=index2;
    ALD(index1,index2)=1;
end

for i=1:size(unique(LD,'rows'),1)
    lncRNA_disease4{i,1}=lncRNA3{LD(i,1)};
    lncRNA_disease4{i,2}=disease3{LD(i,2)};
end
xlswrite('D:\我的文档\生物信息\LDLMD 重投\论文\SuppleTable3',lncRNA_disease4);



