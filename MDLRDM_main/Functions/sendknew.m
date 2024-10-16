function result = sendknew(data, time, kn)
[m,~] = size(data);
if nargin<3
    kn = 1;
elseif nargin<2   
    time = m;
end
datanew = data;
dataold = data;
for i = 1:time
    for j = 1:m
        listn = find(datanew(j,:)>0);%找到大于0的每一项
        newsum = sum(data(listn,:));%
        listsum = find(newsum>0);
        listk = setdiff(listsum,listn);%找到差别
        if(isempty(listk))
            continue;
        end
        datanew(j,listk) = kn^i;
    end
    if dataold == datanew
        break;
    else
        ffff = dataold - datanew;
    end
    dataold = datanew;
end
for i = 1:m
    dataold(i,i) = 0; 
end
result = dataold;