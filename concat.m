function [l0,count] = concat(bw,weight)
count=1;
for i=1:size(bw)
    s=size(cell2mat(bw(i)));
    %count+s(1)-1
    if(weight(i)==1)
   l0(count:count+s(1)-1,:)=cell2mat(bw(i)); 
    count=count+s(1);
    end
end

end