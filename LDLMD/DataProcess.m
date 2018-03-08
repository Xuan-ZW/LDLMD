function [miRNA, lncRNA, disease,MD,ML,AMD,AML,LMD,miRNA_lncRNA,miRNA_disease,lncRNA_disease,LD3,LD4]=DataProcess
%% 处理数据，获得lncRNA-miRNA-disease,并且返回lncRNA, miRNA, disease 序列以及ID
%% 导入miRNA-lncRNA 并进行相关处理
%导入miRNA-lncRNA 
[~, ~, miRNA_lncRNA1] = xlsread('miRNA-lncRNA.xls','A2:A10213');
[~, ~, miRNA_lncRNA2] = xlsread('miRNA-lncRNA.xls','C2:C10213');
[~, ~, lncRNAmiRNAstarBasedatabase2015] = xlsread('D:\我的文档\生物信息\LDLMD 重投\LDLMD\lncRNA-miRNA_starBase_database_2015.xls','miRNA-lncRNA(no-hsa-mir)');
miRNA_lncRNA=[miRNA_lncRNA1,miRNA_lncRNA2];
miRNA_lncRNA=[miRNA_lncRNA;lncRNAmiRNAstarBasedatabase2015];
%全部变成小写
miRNA_lncRNA=lower(miRNA_lncRNA);
%miRNA 去掉-3p/-5p,
miRNA=miRNA_lncRNA(:,1);
miRNA_lncRNA(:,1)=regexprep(miRNA,'-[3,5]p','');
%获得miRNA,lncRNA 序列,并去除重复, 赋予ID
miRNA1=[sprintfc('%d',1:length(unique(miRNA_lncRNA(:,1))))',unique(miRNA_lncRNA(:,1))];

lncRNA=[sprintfc('%d',1:length(unique(miRNA_lncRNA(:,2))))',unique(miRNA_lncRNA(:,2))];


%% 导入miRNA_disease 并进行相关处理
%导入miRNA_disease
[~, ~, miRNA_disease] = xlsread('miRNA-disease_HMDD_2015.xls','A1:B5430');
miRNA_disease=lower(miRNA_disease);
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'carcinoma','cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'neoplasms','cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),',','');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer hepatocellular','hepatocellular cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer neuroendocrine','neuroendocrine cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer non-small-cell lung','non-small-cell lung cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer oral','oral cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer renal cell','renal cell cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer squamous cell','squamous cell cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer basal cell','basal cell cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer ductal breast','breast ductal cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer ehrlich tumor','ehrlich tumor cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer embryonal','embronal cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer endometrioid','endometrioid cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer germ cell and embryonal','germ cell and embryonal cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer glandular and epithelial','landular and epithlial cancer');
miRNA_disease(:,2)=regexprep(miRNA_disease(:,2),'cancer small cell','small cell cancer');

miRNA2=[sprintfc('%d',1:length(unique(miRNA_disease(:,1))))',unique(miRNA_disease(:,1))];

disease=[sprintfc('%d',1:length(unique(miRNA_disease(:,2))))',unique(miRNA_disease(:,2))];



%% 求miRNA1 和miRNA2 的交集，从而抽取出lncRNA-miRNA-disease,并且获得新的ID , 用ID表示的关联矩阵

[miRNA,r1,r2]=intersect(miRNA1(:,2),miRNA2(:,2));
r3=find(ismember(miRNA_lncRNA(:,1),miRNA));
r4=find(ismember(miRNA_disease(:,1),miRNA));
miRNA_lncRNA=miRNA_lncRNA(r3,:);
miRNA_disease=miRNA_disease(r4,:);
%lncRNA=unique(miRNA_lncRNA(:,2));
lncRNA=[sprintfc('%d',1:length(unique(miRNA_lncRNA(:,2))))',unique(miRNA_lncRNA(:,2))];

miRNA=[sprintfc('%d',1:length(unique([miRNA_lncRNA(:,1);miRNA_disease(:,1)])))',unique([miRNA_lncRNA(:,1);miRNA_disease(:,1)])];

disease=[sprintfc('%d',1:length(unique(miRNA_disease(:,2))))',unique(miRNA_disease(:,2))];
% 构造lncRNA_miRNA 关联矩阵
N=size(miRNA_lncRNA,1);
for i=1:N   
    miRNA1=miRNA_lncRNA{i,1};
    lncRNA2=miRNA_lncRNA{i,2};
    index1=find(strcmp(miRNA(:,2),miRNA1));
    index2=find(strcmp(lncRNA(:,2),lncRNA2));
    ML(i,1)=index1;
    ML(i,2)=index2;
    AML(index1,index2)=1;
end

%构造miRNA-disease 关联矩阵

N=size(miRNA_disease,1);
for i=1:N   
    miRNA1=miRNA_disease{i,1};
    disease2=miRNA_disease{i,2};
    index1=find(strcmp(miRNA(:,2),miRNA1));
    index2=find(strcmp(disease(:,2),disease2));
    MD(i,1)=index1;
    MD(i,2)=index2;
    AMD(index1,index2)=1;
end




ML=unique(ML,'rows');
MD=unique(MD,'rows');
LMD=AML'*AMD;

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
%% 寻找正样本，做留一交叉验证.即lncRNA-miRNA-disease.其中lncRNA 与disease是有关联的
co_lncRNA=intersect(miRNA_lncRNA(:,2),lncRNA3);
co_disease=intersect(miRNA_disease(:,2),disease3);
t1=find(ismember(lncRNA3,co_lncRNA));
t2=find(ismember(disease3,co_disease));
LD=unique(LD,'rows');
LD1=[];
for i=1:length(t1)
    LD1((i-1)*length(t2)+1:i*length(t2),1)=t1(i);
    LD1((i-1)*length(t2)+1:i*length(t2),2)=t2;
end

LD2=LD(find(ismember(LD,LD1,'rows')),:);

for i=1:size(LD2,1)
    LD3{i,1}=lncRNA3{LD2(i,1)};
    LD3{i,2}=disease3{LD2(i,2)};
end



t=1;

for i=1:size(LD3,1)
    L1=LD3{i,1};
    D1=LD3{i,2};
    miRNA1=miRNA_lncRNA(find(strcmp(L1,miRNA_lncRNA(:,2))),1);
    miRNA2=miRNA_disease(find(strcmp(D1,miRNA_disease(:,2))),1);
    miRNA3=intersect(miRNA1,miRNA2);
    if length(miRNA3)>0
        for j=1:length(miRNA3)
            LD4{t,1}=L1;
            LD4{t,2}=miRNA3{j};
            LD4{t,3}=D1;
            t=t+1;
        end
    end
end

save miRNA miRNA
save lncRNA lncRNA
save disease disease
save ML ML
save MD MD
save AML AML
save AMD AMD
save LMD LMD

end






