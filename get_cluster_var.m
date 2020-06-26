function [varst7,coefs,count]=get_cluster_var(init_w,end_w,count,mon_num,cnum,vercount)
varst7=[];
coefs=zeros(1,1);
tcount=1;
cluster_len=end_w-init_w+1;
for xx=init_w:end_w
    varst2=zeros(1,cnum+mon_num+2);
    varst2(mon_num+xx)=2;
    %varst2(mon_num+indj)=1; for weight
    varst7=[varst7;varst2]; %%%%second degree terms   
    coefs(tcount,1)=((cluster_len-1)/cluster_len)^2+1/cluster_len;
    count=count+1;
    tcount=tcount+1;
end

for xx=init_w:end_w-1
    for yy=xx+1:end_w
    varst2=zeros(1,cnum+mon_num+2);
    varst2(mon_num+xx)=1;
    varst2(mon_num+yy)=1; 
    varst7=[varst7;varst2]; %%%%second degree cross terms
    coefs(tcount,1)=-2/cluster_len-4/(cluster_len)^2;
    tcount=tcount+1;
    count=count+1;
    end
    
end
% disp('insside sup');
%     size(varst7)
%     disp('inside coef');
%  size(coefs)
%  count

end
