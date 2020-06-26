function [coeffs,ws,exp_mat,R1,mv,l0,l0_new,mu,w,ee,c,indx1,indy1,cnum,cluster_num]= write_Gams_levelset2_arplane(im,file,fixed,weight,csize,dsize,rnum,begin,programRelaxOrder)

%%%%%% the order of variables is as follows: t, lambdas, Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% points=gen_random_vars(cnum,100);%%%% order of variables in points is:
% weights, miu,e,c,weps,mcnum,cdiff(vase c5=c4)
% fpoly=[]
if begin==0
    
load('airplane_random_inv.mat');
end


R1=0;ws=0;

im = imread(file);
sim=size(im);
sx=sim(1);
sy=sim(2);

max_dir=max(sx,sy);


maximum=max_dir/2;

% imr=im(:,:,1);
% img=im(:,:,2);
% im=(im2double(imr)./im2double(imr+img));
 mon_num=28;
%im=rgb2gray(im);
BW=edge(im,'canny');
BW2 = bwmorph(BW,'skel',Inf);
BW3 = bwmorph(BW2,'branchpoints',1);
while max(max(BW3))==1
    BW2=BW2-BW3
    BW3 = bwmorph(BW2,'branchpoints',1);
end
BW=BW2
bs = bwboundaries(BW,'noholes');
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

bs2 = bwboundaries(BW);
l02=concat(bs2,ones(size(bs2)));
total_vars=[];
total_coefs=[];
size_bs=size(bs);
count=1;
init_w=1;
end_w=0;
cluster_size=randi([1,2])


vx=45;



vertical_points=(find(l0(:,2)==vx));
vertical_points=l0(vertical_points,:)
im(:,:)=0;
s2=size(vertical_points)
vercount=s2(1);
for(i=1:s2(1))
    vertical_points(i)
    im(vertical_points(i,1),vertical_points(i,2))=255;
end

%imshow(im);






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
s2=size(l02);
for(i=1:s2(1))
    im(l02(i,1),l02(i,2))=255;
end

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

%maximum=0;

indx1=(indx-ones(size(indx))*meanx)/maximum;
indy1=(indy-ones(size(indy))*meany)/maximum;

 maximum=maximum;
    l0_new=l0_new/(maximum);
    tl0_new=l0_new;
    
    i_epsil=(1/(2*maximum));
    

sl0=size(l0);
l0_s=size(l0_new);


 xs=l0_new(1:m,1);
 ys=l0_new(1:m,2);
 ixs=l0(1:m,1);
 iys=l0(1:m,2);

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
VEC=monomials([x;y],[0 1 2 3 4 5 6]);
segsize=size(bs);
min_edge/(2*maximum)
dsize=size(xs);
dsize=dsize(1)
csize=ceil(dsize/cnum)
pcount=1;

vars_vertical=zeros(1,cnum+mon_num+2+vercount+cnum);
vars_temp=zeros(1,cnum+mon_num+2+vercount+cnum);
vars_temp(mon_num+cnum+2)=1;
vars_vertical=[vars_vertical;vars_temp]


%e=100;
si=size(xs);
eval=0;obj=0;
sumres=0;count=1;scount=1;
err=0;perr=0;serr=0;varsum=0;pre_serr=0;

ineqs=1;

%%%%%5black_points
s=size(indx);
bnum=s(1);


vars_ineq=[];
coefs_ineq=[];

vars_cons=zeros(1,cnum+mon_num+2+vercount+cnum);
vars_all=[];
for w=1:cnum%%%%add constraints w^2-w=0
   ineqPolySys{pcount}.typeCone = -1;
    ineqPolySys{pcount}.dimVar =cnum+mon_num+2+vercount+cnum;
    ineqPolySys{pcount}.degree =2;
    ineqPolySys{pcount}.noTerms = 2;
    vars1=zeros(1,cnum+mon_num+2+vercount+cnum);
    vars1(mon_num+w)=2;

    vars2=zeros(1,cnum+mon_num+2+vercount+cnum);
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
ineqPolySys{pcount}.dimVar =cnum+mon_num+2+vercount+cnum;
ineqPolySys{pcount}.degree =2;
ineqPolySys{pcount}.noTerms = 1+mon_num;
vars1=zeros(1,cnum+mon_num+2+vercount+cnum);
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
ineqPolySys{pcount}.dimVar =cnum+mon_num+2+vercount+cnum;
ineqPolySys{pcount}.degree =2;
ineqPolySys{pcount}.noTerms = 1+mon_num;
vars1=zeros(1,cnum+mon_num+2+vercount+cnum);
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



% 
for i=1:maxy-miny
    eval=(subs(VEC,[x,y],[(minx-meanx)/(maximum),sminx(i)]));
    evald=double(eval);
     ineqPolySys{pcount}.typeCone = 1;
ineqPolySys{pcount}.dimVar =cnum+mon_num+2+vercount+cnum;
ineqPolySys{pcount}.degree =1;
ineqPolySys{pcount}.noTerms = mon_num+1;
vars1=zeros(1,cnum+mon_num+2+vercount+cnum);
vars_gen_weight;
ineqPolySys{pcount}.supports = [vars1;varst7];

coefs=zeros(6,1) ;
create_coef_4_1;
coefs(1)=-0.01 ;
        ineqPolySys{pcount}.coef     = coefs;
    % rminx(i)= (evald'*c-eres) >=0;
    % srminx(i)= evald'*c<=eres;

     eval=(subs(VEC,[x,y],[(maxx-meanx)/(maximum),smaxx(i)]));
  
    evald=double(eval);
    
         ineqPolySys{pcount+1}.typeCone = 1;
ineqPolySys{pcount+1}.dimVar =cnum+mon_num+2+vercount+cnum;
ineqPolySys{pcount+1}.degree =1;
ineqPolySys{pcount+1}.noTerms = mon_num+1;
vars1=zeros(1,cnum+mon_num+2+vercount+cnum);

vars_gen_weight;
ineqPolySys{pcount+1}.supports = [vars1;varst7];

coefs=zeros(6,1) ;
create_coef_4_1;
coefs(1)=-0.01 ;
        ineqPolySys{pcount+1}.coef     = coefs;
%      rmaxx(i)=( evald'*c) -eres>=0;
   %  srmaxx(i)= evald'*c<=eres;

    pcount=pcount+2;
end

for i=1:maxx-minx
 eval=(subs(VEC,[x,y],[sminy(i),(miny-meany)/(maximum)]));
    evald=double(eval);
         ineqPolySys{pcount}.typeCone = 1;
ineqPolySys{pcount}.dimVar =cnum+mon_num+2+vercount+cnum;
ineqPolySys{pcount}.degree =1;
ineqPolySys{pcount}.noTerms = mon_num+1;
vars1=zeros(1,cnum+mon_num+2+vercount+cnum);
vars_gen_weight;
ineqPolySys{pcount}.supports = [vars1;varst7];

coefs=zeros(6,1) ;
create_coef_4_1;
coefs(1)=-0.01 ;
        ineqPolySys{pcount}.coef     = coefs;
    % rminy(i)= ((evald'*c-eres)-eres) >=0;
  %   srminy(i)= evald'*c<=eres;

     eval=(subs(VEC,[x,y],[smaxy(i),(maxy-meany)/(maximum)]));
    evald=double(eval);
         ineqPolySys{pcount+1}.typeCone = 1;
ineqPolySys{pcount+1}.dimVar =cnum+mon_num+2+vercount+cnum;
ineqPolySys{pcount+1}.degree =1;
ineqPolySys{pcount+1}.noTerms = mon_num+1;
vars1=zeros(1,cnum+mon_num+2+vercount+cnum);
vars_gen_weight;
ineqPolySys{pcount+1}.supports = [vars1;varst7];

coefs=zeros(6,1) ;create_coef_4_1;
coefs(1)=-0.01 ;
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
ineqPolySys{pcount}.dimVar =cnum+mon_num+2+vercount+cnum;
ineqPolySys{pcount}.degree =2;
ineqPolySys{pcount}.noTerms = 1+mon_num;
vars16=zeros(1,cnum+mon_num+2+vercount+cnum);
sec_degree

ineqPolySys{pcount}.supports = [varst7;vars16];

coefs=-ones(1+mon_num,1);
%coefs(mon_num+1:(mon_num)*(mon_num-1)/2+mon_num)=2;
coefs(1+mon_num)=1;
        ineqPolySys{pcount}.coef= coefs;
        
                pcount=pcount+1;
        
 
% %  %%%%%%invariant :exp
if begin==0
  ineqPolySys{pcount}.typeCone = -1;
ineqPolySys{pcount}.dimVar =cnum+mon_num+2+vercount+cnum;
ineqPolySys{pcount}.degree =2;
ineqPolySys{pcount}.noTerms = 1+mon_num+mon_num+mon_num*(mon_num-1)/2;

full_deg2_mons;
%create_var2;
varst2=zeros(1,cnum+mon_num+2+vercount+cnum);
    varst2(mon_num+cnum+2)=1;
varst3=zeros(1,cnum+mon_num+2+vercount+cnum);
ineqPolySys{pcount}.supports = [varst2;varst7];
%inv_coef=exp;
coefs(2:1+mon_num+mon_num+mon_num*(mon_num-1)/2)=inv_coeffs;%(1+mon_num+mon_num+mon_num*(mon_num-1)/2,1);
%coefs(1+mon_num+mon_num+mon_num*(mon_num-1)/2+1)=-coefs(1+mon_num+mon_num+mon_num*(mon_num-1)/2+1);
%coefs(mon_num+1:(mon_num)*(mon_num-1)/2+mon_num)=2;
coefs(1)=-1;
%coefs(1+mon_num+mon_num+mon_num*(mon_num-1)/2+1)=1;
        ineqPolySys{pcount}.coef= coefs;
        
                pcount=pcount+1;
        
        
end
        

        
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
vars7=[];
for xx=1:cnum
    vars2=zeros(1,cnum+mon_num+2+vercount+cnum);
    vars2(mon_num+xx)=1;
vars7=[vars7;vars2];    
    
end

varsd=[];
for xx=1:cnum
    vars2=zeros(1,cnum+mon_num+2+vercount+cnum);
    vars2(mon_num+cnum+2+vercount+xx)=1;
    
varsd=[varsd;vars2];    
    
end

       % gen_weight;
        objPoly.typeCone = 1;
        objPoly.dimVar   = cnum+mon_num+2+vercount+cnum;
        objPoly.degree   = 2;
        objPoly.noTerms  = cluster_count-1+cnum+1+1+cnum;
        vars1=zeros(1,cnum+mon_num+2+vercount+cnum);
        vars1(cnum+mon_num+2)=2;%%inv
        
        vars2=zeros(1,cnum+mon_num+2+vercount+cnum);
        vars2(cnum+mon_num+1)=1;%%e
        
        vars3=zeros(1,cnum+mon_num+2+vercount+cnum);
        vars3(1)=1;
       % vars3(cnum+mon_num+2+vercount+cnum+1)=1;
        
        vars4=zeros(1,cnum+mon_num+2+vercount+cnum);
        vars4(2)=2;
        %vars4(cnum+mon_num+2+vercount+cnum+1)=1;
        
        vars5=zeros(1,cnum+mon_num+2+vercount+cnum);
        vars5(3)=2;
      %  vars5(cnum+mon_num+2+vercount+cnum+1)=1;
        
        vars6=zeros(1,cnum+mon_num+2+vercount+cnum);
        vars6(4)=1;
       % vars8=
%         vars6=zeros(1,cnum+mon_num+2+vercount+cnum+6);
%         vars6(4)=1;
%         vars6=zeros(1,cnum+mon_num+2+vercount+cnum);
%         vars6(cnum+mon_num+2+vercount+cnum+5)=1;

        
 
         objPoly.supports =[total_vars;vars7;vars2;vars1;varsd];
        size(total_coefs)
      %  total_coeffs=zeros(size(total_coeffs));
        objPoly.coef     = [total_coefs/(cnum);zeros(cnum,1)/cnum^2;1;500;zeros(cnum,1)*2]
       
        lbd=-1*zeros(1,cnum+mon_num+2+vercount+cnum);
 
            lbd(mon_num+cnum+1)=1/(24*maximum);
         lbd(mon_num+1:cnum+mon_num)=0;

         lbd(1:mon_num-1)=-2*5;

        start_point=1
        end_point=start_point+cluster_size
                iii=1;
                cluster_num=0;
        count=0;
        while iii==1
            count=count+1;
            if start_point>=cnum
                break
            end
            
            rchoice=randi([0,3]);
            rpoint=randi([start_point,end_point]);
            
            if rchoice==1%%%choose cluster
                lbd(mon_num+rpoint)=1;
                rchoice
            end
            start_point=start_point+cluster_size;
            end_point=end_point+cluster_size
            cluster_num=cluster_num+1
            if end_point>cnum
                end_point=cnum;
            end
            
        end
        

 


         lbd(mon_num+cnum+2)=-0.01;
 
         lbd(mon_num)=0;

        ubd=ones(1,cnum+mon_num+2+vercount+cnum);
       
                ubd(1:mon_num)=2*5;
             ubd(mon_num+cnum+1)=1/(2*maximum);
        ubd(mon_num+cnum+2)=0.01;
    param.scalingSW=0;
param.SDPsolver='sdpa';
        param.relaxOrder =programRelaxOrder;
        % param.POPsolver='interior-point';

  disp('formulated problem'); 
 [param,SDPobjValue,POP,cpuTime,SDPsolverInfo,SDPinfo]=sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);

 coeffs=SDPinfo.y(1:mon_num);
                save('coeffs_poly.mat','coeffs');
     
%weps1=(weps>=0);



%e=SDPinfo.y(cnum+mon_num+2+vercount+cnum-3);%%%%one w fixed
e=1;
ws=1;w=1;
exp_mat=0;R1=0;mv=0;mu=0;ee=e;c=0;






