function [rankTopTenlncRNADisease] =ten_percent_of_predictive_result(lncRNA,disease,prediction)
%% 返回预测结果中， 每个lncRNA关联的前10个疾病。
ten_pairs_of_predictive_result=zeros(size(lncRNA,1)*10,3);
for i=1:size(lncRNA,1)
    ten_pairs_of_predictive_result(10*(i-1)+1:10*(i-1)+10,:)=prediction(376*(i-1)+1:376*(i-1)+10,:);
end
A= ten_pairs_of_predictive_result;
load disease;
load lncRNA;
for i=1:11020
    rankTopTenlncRNADisease{i,1}=A(i,1);
    rankTopTenlncRNADisease{i,2}=lncRNA(A(i,2),2);
    rankTopTenlncRNADisease{i,3}=disease(A(i,3),2);
end
save rankTopTenlncRNADisease rankTopTenlncRNADisease;
end











