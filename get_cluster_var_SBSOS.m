%function [varst7,coefs,count]=get_cluster_var(init_w,end_w,count,mon_num,cnum,vercount)
varst7=[];
coefs=zeros(1,1);
tcount=1;
cluster_len=end_w-init_w+1;
%fvar=[];
for xx=init_w:end_w
    fi=zeros(1,vnum+1);
    fi(mon_num+xx)=2;
    fi(vnum+1)=(((cluster_len-1)/cluster_len)^2+1/cluster_len)/10;
    fvar=[fvar;fi]; %%%%second degree terms

    count=count+1;
    tcount=tcount+1;
end

for xx=init_w:end_w-1
    for yy=xx+1:end_w
    fi=zeros(1,vnum+1);
    fi(mon_num+xx)=1;
    fi(mon_num+yy)=1;
    fi(vnum+1)=(-2/cluster_len-4/(cluster_len)^2)/10;
    fvar=[fvar;fi]; %%%%second degree cross terms

    tcount=tcount+1;
    count=count+1;
    end
    
end
% fi=zeros(1,vnum);
% fi(vnum+1)=100;
% pop.G{ineq_count}=fvar;
% ineq_count=ineq_count+1; 



