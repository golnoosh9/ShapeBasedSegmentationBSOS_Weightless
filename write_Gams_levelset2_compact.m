function [coeffs,ws,exp_mat,R1,mv,l0,l0_new,mu,w,ee,c,indx1,indy1,cnum,cluster_num]= write_Gams_levelset2(im,BW,fixed,weight,csize,dsize,rnum,begin,programRelaxOrder)

%%%%%% the order of variables is as follows: t, lambdas, Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% points=gen_random_vars(cnum,100);%%%% order of variables in points is:
% weights, miu,e,c,weps,mcnum,cdiff(vase c5=c4)
% fpoly=[]

%load('four_random_inv.mat');
if begin==0
load('lip_random_inv.mat');
end
R1=0;ws=0;
% mset clear;
% mset clearmeas;

%im = imread(file);
sim=size(im);
sx=sim(1);
sy=sim(2);

max_dir=max(sx,sy);


maximum=max_dir/2;

% imr=im(:,:,1);
% img=im(:,:,2);
% im=(im2double(imr)./im2double(imr+img));
 mon_num=15;
%im=rgb2gray(im);


% BW2 = bwmorph(BW,'skel',Inf);
% BW3 = bwmorph(BW2,'branchpoints',1);
% while max(max(BW3))==1
% BW2=BW2-BW3
% BW3 = bwmorph(BW2,'branchpoints',1);
% end
% BW=BW2;

image_indexes=1:10:80;
s_points=size(image_indexes);
points_choose=floor(s_points(2)-1);
rand_indexes=randi([1,s_points(2)],points_choose,1)
bs = bwboundaries(BW,'noholes');
temp=cell2mat(bs(1));
sbs=size(temp)
max(image_indexes)

temp=temp(image_indexes(rand_indexes),:);
bs2{1}=temp;

%  temp=cell2mat(bs(2));
%  sbs=size(bs);
%  temp=temp([1,60],:);
%  bs2{2}=temp;
 bs=bs2';


l0=concat(bs,ones(size(bs)));
segs=size(bs);
for j=1:segs
temp=size(cell2mat(bs(j)));
sarray(j)=temp(1);
end
dsize=sarray*ones(segs);
cnum=dsize;

csize=ceil(dsize/cnum);
[xp,yp]=cluster_neighbors(bs,cnum,csize);


total_vars=[];
total_coefs=[];
size_bs=size(bs);
count=1;
init_w=1;
end_w=0;
cluster_size=3;


vx=45;
vercount=0;




for i=1:size_bs(1)
   
    tsize=size(cell2mat(bs(i)));
    cluster_num=tsize(1)/cluster_size;
    resid_cluster=mod(tsize(1),cluster_size);
    for j=1:cluster_num
      end_w=end_w+cluster_size;
      [temp_var,temp_coef,count]=get_cluster_var(init_w,end_w,count,mon_num,cnum,vercount);
      total_vars=[total_vars;temp_var];
      total_coefs=[total_coefs;temp_coef];
      init_w=init_w+cluster_size;
     
    end
    %%%%%%residue
    if resid_cluster>0
    end_w=end_w+resid_cluster;
      [temp_var,temp_coef,count]=get_cluster_var(init_w,end_w,count,mon_num,cnum,vercount);
         
      total_vars=[total_vars;temp_var];
      total_coefs=[total_coefs;temp_coef];
    init_w=init_w+resid_cluster;

    end
end

cluster_count=count;

im(:,:)=0;
s2=size(l0);
for(i=1:s2(1))
    im(l0(i,1),l0(i,2))=255;
end
% imshow(im);
% jjj

% BW=edge(im,'canny');
% bs2 = bwboundaries(BW);
% l02=concat(bs2,ones(size(bs2)));
% im(:,:)=0;
% s2=size(l02);
% for(i=1:s2(1))
%     im(l02(i,1),l02(i,2))=255;
% end

inds=find(im==0);
[indx,indy]=ind2sub(size(im),inds);
meanx=sx/2;
meany=sy/2;
nn=ones(size(l0));
h=[meanx;meany];
nn(:,1)=nn(:,1).*meanx;
nn(:,2)=nn(:,2).*meany;
l0_new=l0-nn;
tl0=l0;
tl0x=tl0(:,1);
tl0y=tl0(:,2);
m=size(l0);%total points
m=m(1)
%cluster_size=1;
%maximum=0;

indx1=(indx-ones(size(indx))*meanx)/maximum;
indy1=(indy-ones(size(indy))*meany)/maximum;
% for a=1:size(l0_new)
%      norm_factor=max(abs(l0_new(a,1)),abs(l0_new(a,2)));%normalize data
%       if(norm_factor>maximum)
%           maximum=norm_factor;
%       end
% end
 maximum=maximum;
    l0_new=l0_new/(maximum);
    tl0_new=l0_new;
    
    i_epsil=(1/(2*maximum));
    
%imtool(im);
% list(:,:)=find(BW==1);
% [x,y]=ind2sub(size(BW),list);
sl0=size(l0);
l0_s=size(l0_new);


 xs=l0_new(1:m,1);
 ys=l0_new(1:m,2);
 ixs=l0(1:m,1);
 iys=l0(1:m,2);
%   xs=[1;1;1;1];
%   ys=[1;1;1;1];

miu=0.8



%%%%%%% add the rectangle%%
%%%%%%%%%
num=1000;
l0=l0-nn;
minx=1;
miny=1;
maxx=sx;
maxy=sy;

sminx=cumsum(ones(maxy-miny,1));
added=ones(maxy-miny,1)*miny;
sminx=added+sminx;

sminx=(sminx-ones(size(sminx))*meany)/(maximum);
smaxx=sminx;

sminy=cumsum(ones(maxx-minx,1));
added=ones(maxx-minx,1)*minx;
sminy=added+sminy;
sminy=(sminy-ones(size(sminy))*meanx)/(maximum);
smaxy=sminy;

xdist=maxx-minx+1;
ydist=maxy-miny+1;
max_edge=max(xdist,ydist)
min_edge=min(xdist,ydist)

% maxy=maxy;
% maxx=maxx;
berr=0;

%[sminx,sminy,smaxx,smaxy]=build_rectangle_samples(minx,miny,maxx,maxy,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotcount=0;
miu=0.8;%percent
dots=1;
syms x y
%VEC = [1;x;y;x^2;x*y;y^2];
VEC=monomials([x;y],[0 1 2 3 4]);
segsize=size(bs);
min_edge/(2*maximum)
dsize=size(xs);
dsize=dsize(1)
csize=ceil(dsize/cnum)
pcount=1;

vars_vertical=zeros(1,cnum+mon_num+2);
vars_temp=zeros(1,cnum+mon_num+2);
vars_temp(mon_num+cnum+2)=1;
vars_vertical=[vars_vertical;vars_temp]


si=size(xs);
eval=0;obj=0;
sumres=0;count=1;scount=1;
err=0;perr=0;serr=0;varsum=0;pre_serr=0;

ineqs=1;

%%%%%5black_points
s=size(indx);
bnum=s(1);


vars_ineq=[]
coefs_ineq=[];

vars_cons=zeros(1,cnum+mon_num+2);
vars_all=[];
for w=1:cnum%%%%add constraints w^2-w=0
ineqPolySys{pcount}.typeCone = -1;
ineqPolySys{pcount}.dimVar =cnum+mon_num+2;
ineqPolySys{pcount}.degree =2;
ineqPolySys{pcount}.noTerms = 2;
vars1=zeros(1,cnum+mon_num+2);
vars1(mon_num+w)=2;

vars2=zeros(1,cnum+mon_num+2);
vars2(mon_num+w)=1;

varst7=[vars1;vars2];

ineqPolySys{pcount}.supports = [varst7];

coefs=zeros(2,1) ;
coefs(1)=1;
coefs(2)=-1;
%coefs(1)=1;


ineqPolySys{pcount}.coef     = coefs;


pcount=pcount+1;



end




for j=1:cnum
    aerr=0;
    xs=xp{j};

    ys=yp{j};

    csizes=size(xs);
    xs=(xs-(ones(csizes)*meanx))/maximum;
    ys=(ys-(ones(csizes)*meany))/maximum;
    csizes=csizes(1);
    csizes=1;
   
for i=1:csizes%cluster_size:size(xs)-cluster_size
    j
%     Derx=diff(VEC,x);
%     Derxy=diff(Derx,y)
%     Dereval=subs(Derxy,[x,y],[xs(i),ys(i)]);
%     Derevald=double(Dereval);
   eval=(subs(VEC,[x,y],[xs(i),ys(i)]));
   veval=(subs(VEC,[x,y],[xs(i),ys(i)]));
   
  
   %eps=-(eval*(power(w(i),2)))+e
 evald=double(eval);
%  
%    %%%%%%%ordering of variables: c1...c5, w1...wcnum,e,miu,eres 
   
% %%%%%k2(i)=((power(evald'*c,1)*w(j)>=-e);%%%%for now bised to outwards
 ineqPolySys{pcount}.typeCone = 1;
ineqPolySys{pcount}.dimVar =cnum+mon_num+2;
ineqPolySys{pcount}.degree =2;
ineqPolySys{pcount}.noTerms = 1+mon_num;
vars1=zeros(1,cnum+mon_num+2);
vars1(cnum+mon_num+1)=1;%%%for e


vars_deg1;
ineqPolySys{pcount}.supports = [vars1;varst7];

coefs=zeros(mon_num,1) ;
create_coef_4_1
%coefs(1)=1;


        ineqPolySys{pcount}.coef     = coefs;





pcount=pcount+1;

%%%%k2(i)=((power(evald'*c,1)*w(j)<=e);%%%%for now bised to outwards
 ineqPolySys{pcount}.typeCone = 1;
ineqPolySys{pcount}.dimVar =cnum+mon_num+2;
ineqPolySys{pcount}.degree =2;
ineqPolySys{pcount}.noTerms = 1+mon_num;
vars1=zeros(1,cnum+mon_num+2);
vars1(cnum+mon_num+1)=1;%%%for e

vars_deg1;
ineqPolySys{pcount}.supports = [vars1;varst7];

coefs=zeros(mon_num,1) ;
create_coef_4_1;
coefs(2:1+mon_num)=-coefs(2:1+mon_num);

        ineqPolySys{pcount}.coef     = coefs;





pcount=pcount+1;


% %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

end

%count=count+csize;
end




for i=1:maxy-miny
    eval=(subs(VEC,[x,y],[(minx-meanx)/(maximum),sminx(i)]));
    evald=double(eval);
     ineqPolySys{pcount}.typeCone = 1;
ineqPolySys{pcount}.dimVar =cnum+mon_num+2;
ineqPolySys{pcount}.degree =1;
ineqPolySys{pcount}.noTerms = mon_num+1;
vars1=zeros(1,cnum+mon_num+2);
vars_gen_weight;
ineqPolySys{pcount}.supports = [vars1;varst7];

coefs=zeros(6,1) ;
create_coef_4_1;
coefs(1)=-0.001 ;
        ineqPolySys{pcount}.coef     = coefs;
    % rminx(i)= (evald'*c-eres) >=0;
    % srminx(i)= evald'*c<=eres;

     eval=(subs(VEC,[x,y],[(maxx-meanx)/(maximum),smaxx(i)]));
  
    evald=double(eval);
    
         ineqPolySys{pcount+1}.typeCone = 1;
ineqPolySys{pcount+1}.dimVar =cnum+mon_num+2;
ineqPolySys{pcount+1}.degree =1;
ineqPolySys{pcount+1}.noTerms = mon_num+1;
vars1=zeros(1,cnum+mon_num+2);

vars_gen_weight;
ineqPolySys{pcount+1}.supports = [vars1;varst7];

coefs=zeros(6,1) ;
create_coef_4_1;
coefs(1)=-0.001 ;
        ineqPolySys{pcount+1}.coef     = coefs;
%      rmaxx(i)=( evald'*c) -eres>=0;
   %  srmaxx(i)= evald'*c<=eres;

    pcount=pcount+2;
end

for i=1:maxx-minx
 eval=(subs(VEC,[x,y],[sminy(i),(miny-meany)/(maximum)]));
    evald=double(eval);
         ineqPolySys{pcount}.typeCone = 1;
ineqPolySys{pcount}.dimVar =cnum+mon_num+2;
ineqPolySys{pcount}.degree =1;
ineqPolySys{pcount}.noTerms = mon_num+1;
vars1=zeros(1,cnum+mon_num+2);

vars_gen_weight;
ineqPolySys{pcount}.supports = [vars1;varst7];

coefs=zeros(6,1) ;
create_coef_4_1;
coefs(1)=-0.001 ;
        ineqPolySys{pcount}.coef     = coefs;
    % rminy(i)= ((evald'*c-eres)-eres) >=0;
  %   srminy(i)= evald'*c<=eres;

     eval=(subs(VEC,[x,y],[smaxy(i),(maxy-meany)/(maximum)]));
    evald=double(eval);
         ineqPolySys{pcount+1}.typeCone = 1;
ineqPolySys{pcount+1}.dimVar =cnum+mon_num+2;
ineqPolySys{pcount+1}.degree =1;
ineqPolySys{pcount+1}.noTerms = mon_num+1;
vars1=zeros(1,cnum+mon_num+2);
vars_gen_weight;
ineqPolySys{pcount+1}.supports = [vars1;varst7];

coefs=zeros(6,1) ;
create_coef_4_1;
coefs(1)=-0.001 ;
        ineqPolySys{pcount+1}.coef     = coefs;
   %  rmaxy(i)=(evald'*c-eres) >=0;
pcount=pcount+2;
end
                
dotsize=0;
weights=0;dots=m;cweights=0;one_formula=0;sec_formula=0;

%%%%%%%%%%%   
 %%%%%%%%%%%%%normalization of coeffs
 %%%%%%%%%%%%%%%%%%%%%%
              ineqPolySys{pcount}.typeCone = -1;
ineqPolySys{pcount}.dimVar =cnum+mon_num+2;
ineqPolySys{pcount}.degree =2;
ineqPolySys{pcount}.noTerms = 1+mon_num;
vars16=zeros(1,cnum+mon_num+2);
sec_degree

ineqPolySys{pcount}.supports = [varst7;vars16];

coefs=ones(1+mon_num,1);
%coefs(mon_num+1:(mon_num)*(mon_num-1)/2+mon_num)=2;
coefs(1+mon_num)=-1;
        ineqPolySys{pcount}.coef= coefs;
        
        pcount=pcount+1
             
  
        %%%%%%%%%%%%%%%%%
 
% %  %%%%%%invariant :exp
if begin==0
                ineqPolySys{pcount}.typeCone = 1;
                ineqPolySys{pcount}.dimVar =cnum+mon_num+2;
                ineqPolySys{pcount}.degree =2;
                ineqPolySys{pcount}.noTerms = 1+mon_num;
                
                vars_deg1
                %create_var2;
                varst2=zeros(1,cnum+mon_num+2);
                varst2(mon_num+cnum+2)=1;
                ineqPolySys{pcount}.supports = [varst2;varst7];
                %inv_coef=exp;
                coefs(2:1+mon_num)=inv_coeffs;%(1+mon_num+mon_num+mon_num*(mon_num-1)/2,1);
                coefs(1)=-1;
                %coefs(1+mon_num+mon_num+mon_num*(mon_num-1)/2+1)=1;
                ineqPolySys{pcount}.coef= coefs;
                
                pcount=pcount+1;
                
                ineqPolySys{pcount}.typeCone = 1;
                ineqPolySys{pcount}.dimVar =cnum+mon_num+2;
                ineqPolySys{pcount}.degree =2;
                ineqPolySys{pcount}.noTerms = 1+mon_num;
                
                vars_deg1
                %create_var2;
                varst2=zeros(1,cnum+mon_num+2);
                varst2(mon_num+cnum+2)=1;
                ineqPolySys{pcount}.supports = [varst2;varst7];
                %inv_coef=exp;
                coefs(2:1+mon_num)=inv_coeffs;%(1+mon_num+mon_num+mon_num*(mon_num-1)/2,1);
                coefs(1)=-1;
                %coefs(1+mon_num+mon_num+mon_num*(mon_num-1)/2+1)=1;
                ineqPolySys{pcount}.coef= -coefs;
                
                pcount=pcount+1;

        
        
end      
        

 
 
 
 
 


vars7=[];
for xx=1:cnum
    vars2=zeros(1,cnum+mon_num+2);
    vars2(mon_num+xx)=1;
vars7=[vars7;vars2];    
    
end



       % gen_weight;
        objPoly.typeCone = 1;
        objPoly.dimVar   = cnum+mon_num+2;
        objPoly.degree   = 2;
        objPoly.noTerms  = cluster_count-1+cnum+2;
        vars1=zeros(1,cnum+mon_num+2+vercount+cnum);
        vars1(cnum+mon_num+2)=2;%%inv
        
        vars2=zeros(1,cnum+mon_num+2);
        vars2(cnum+mon_num+1)=1;%%e
        
        vars1=zeros(1,cnum+mon_num+2);
        vars1(cnum+mon_num+2)=1;%%e_inv
        
       
 
         objPoly.supports =[total_vars;vars7;vars2;vars1];
  
      %  total_coeffs=zeros(size(total_coeffs));
    objPoly.coef     = [total_coefs/cnum;-ones(cnum,1)/cnum;10;10]
       
        lbd=-1*zeros(1,cnum+mon_num+2);
      
            lbd(mon_num+cnum+1)=-1;
         lbd(mon_num+1:cnum+mon_num)=0;
        lbd(mon_num+cnum+2)=-1;
         lbd(1:mon_num-1)=-2*5;
                
     
                
                
                
                
                
                
   
       
                
                

        ubd=ones(1,cnum+mon_num+2);
       
        ubd(1:mon_num)=2*5;
      
        ubd(mon_num+cnum+1)=1;
        ubd(mon_num+cnum+2)=1;

     
param.scalingSW=0;
%param.SDPsolverEpsilon=0.000000000001;
%param.SDPsolver='sdpa'
        param.relaxOrder =programRelaxOrder;
        % param.POPsolver='interior-point';
    %    save('before.mat');
      
 [param,SDPobjValue,POP,cpuTime,SDPsolverInfo,SDPinfo]=sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);


 
                allmonomials=full(SDPinfo.xIdxVec);
                allmonomials(1,:)=[];
                
                degrees=sum(allmonomials,2);
                IndexesDeg1=find(degrees==1);
                deg1monomials=allmonomials(IndexesDeg1,:);
                coeffs=SDPinfo.y(IndexesDeg1);
                coeffs=coeffs(1:mon_num);
                % coeffs=POP.xVect(1:mon_num);
                
                
                
                IndexesDeg2=find(degrees==2);
                deg2monomials=allmonomials(IndexesDeg2,:);
                mmSize=size(deg2monomials);
                varSize=mmSize(2);
                mon2Size=mmSize(1);
                deg2Values=SDPinfo.y(IndexesDeg2);
                
                IndexesDeg1=find(degrees==1);
                deg1monomials=allmonomials(IndexesDeg1,:);
                expectation=SDPinfo.y(IndexesDeg1);
                
                
                covarianceMat=zeros(varSize);
                for i=1:mon2Size
                indexes=find(deg2monomials(i,:)>0);
                sInd=size(indexes);
                if(sInd(2)==1)
                covarianceMat(indexes,indexes)=deg2Values(i);
                else
                x=indexes(1);
                y=indexes(2);
                covarianceMat(x,y)=deg2Values(i);
                covarianceMat(y,x)=deg2Values(i);
                
                end
                
                end
                
                variables_to_count=mon_num+cnum+2;
                coeffs_covar=covarianceMat(1:variables_to_count,1:variables_to_count);

                save('coeffs.mat','coeffs','coeffs_covar','mon_num','cnum','expectation');

e=1;
ws=1;w=1;
exp_mat=0;R1=0;mv=0;mu=0;ee=e;c=0;






