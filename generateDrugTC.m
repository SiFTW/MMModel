function drugTC = generateDrugTC(maxVal,midPoint,slope,maxTime)
%this is inspired by https://peerj.com/articles/4251/

t=[1:maxTime];
drugTC=zeros(size(t));
for i=1:length(t)
    drugTC(i)=maxVal/(1+exp(-slope*(i-midPoint)));
end

end