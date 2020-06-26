%%var order: C, Wi, fixed epsilons for now
function[coeffs]= BSOS_ShapeSeg(begin,im,file)
if begin==0
load('air_random_inv.mat');
end

mon_num=28;
ReadImageDataGray;
sim=size(im);
sx=sim(1);
sy=sim(2);

max_dir=max(sx,sy);
maximum=max_dir/2;

image_indexes=[1:2,40:41,70:71,85:86,105:106,130:131];
s_points=size(image_indexes);
points_choose=floor(s_points(2)-1);
rand_indexes=randi([1,s_points(2)],points_choose,1)
bs = bwboundaries(BW,'noholes');
temp=cell2mat(bs(1));
sbs=size(temp)
max(image_indexes)

temp=temp(image_indexes,:);
bs2{1}=temp;

% temp=cell2mat(bs(2));
%  sbs=size(bs);
%  temp=temp([1,5],:);
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

vnum=mon_num+cnum;

size_bs=size(bs);
count=1;
init_w=1;
end_w=0;
cluster_size=3;
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

fvar=[];
ineq_count=1;
fderivative=[];

for i=1:size_bs(1)
tsize=size(cell2mat(bs(i)));
cluster_num=tsize(1)/cluster_size;
resid_cluster=mod(tsize(1),cluster_size);
for j=1:cluster_num
end_w=end_w+cluster_size;
get_cluster_var_SBSOS;
init_w=init_w+cluster_size;
end
%%%%%%residue
if resid_cluster>0
end_w=end_w+resid_cluster;
get_cluster_var_SBSOS;

init_w=init_w+resid_cluster;

end
end

im(:,:)=0;
s2=size(l0);
for(i=1:s2(1))
im(l0(i,1),l0(i,2))=255;
end
imshow(im);
%llll

cluster_count=count;

size_bs=size(bs);
count=1;
init_w=1;
end_w=0;
cluster_size=3;

syms x y
%VEC = [1;x;y;x^2;x*y;y^2];
VEC=monomials_gen([x;y],[0 1 2 3 4 5 6]);
sim=size(im);
sx=sim(1);
sy=sim(2);
VecDiff=diff(VEC);
VecDiff=VecDiff;

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

eval=(subs(VEC,[x,y],[xs(i),ys(i)]));
veval=(subs(VEC,[x,y],[xs(i),ys(i)]));


evalDiff=(subs(VecDiff,[x,y],[xs(i),ys(i)]));
evalDiffDouble=abs(double(evalDiff));
%eps=-(eval*(power(w(i),2)))+e
evald=double(eval);

gtotal=[];
for mon=1:mon_num%%-e<eval.w<e=>-1<eval.w.max<1=> 0<(eval.w.max+1)/2<1
gi=zeros(1,vnum+1);%%%vars plus multiplier
gi(mon_num+j)=1;%%%weight
gi(mon)=1;
gi(vnum+1)=evald(mon)*(maximum)/2;
gtotal=[gtotal;gi];

end
gi=zeros(1,vnum+1);
gi(vnum+1)=0.5;
pop.G{ineq_count}=[gtotal];

ineq_count=ineq_count+1;

fderivative=[];
for mon=1:mon_num
gi=zeros(1,vnum+1);%%%derivative elements
gi(mon)=1;
gi(vnum+1)=evalDiffDouble(mon)/10;
fderivative=[fderivative;gi];
end
gi=zeros(1,vnum+1);
gi(vnum+1)=-0.001;
pop.G{ineq_count}=[fderivative;gi];;

ineq_count=ineq_count+1;
% hhh

end


end

%BSOS_RectangleConstraints;
if(begin==0)%%%invariant constraint
gtotal=[];
for i=1:mon_num
gi=zeros(1,vnum+1);

gi(i)=1;
gi(vnum+1)=inv_coeffs(i)*100;
gtotal=[gtotal;gi];
end

pop.G{ineq_count}=gtotal;

ineq_count=ineq_count+1;

for i=1:mon_num
gi=zeros(1,vnum+1);

gi(i)=1;
gi(vnum+1)=-inv_coeffs(i)*100;
gtotal=[gtotal;gi];
end

pop.G{ineq_count}=gtotal;
ineq_count=ineq_count+1;
end

gnorm=[];
for i=1:mon_num
gi=zeros(1,vnum+1);
gi(i)=1;
gi(vnum+1)=1/(mon_num*2);
gnorm=[gnorm;gi];

end
gi=zeros(1,vnum+1);
gi(vnum+1)=-1;
pop.G{ineq_count}=[gnorm;gi];
ineq_count=ineq_count+1;
fsum=[];
for j=1:cnum
fi=zeros(1,vnum+1);
fi(mon_num+j)=1;
fi(vnum+1)=-1/cnum;%%%maximize sum=> minimize the negative
fsum=[fsum;fi];
gi2=zeros(1,vnum+1);
gi2(mon_num+j)=2;
gi2(vnum+1)=1;
gi1=zeros(1,vnum+1);
gi1(mon_num+j)=1;
gi1(vnum+1)=-1;
pop.G{ineq_count}=[gi2;gi1];

ineq_count=ineq_count+1;

gi1=zeros(1,vnum+1);
gi1(mon_num+j)=1;
gi1(vnum+1)=1;
pop.G{ineq_count}=[gi1];

ineq_count=ineq_count+1;

gi1=zeros(1,vnum+1);
gi1(mon_num+j)=1;
gi1(vnum+1)=-1;
gi0=zeros(1,vnum+1);

gi0(vnum+1)=1;
pop.G{ineq_count}=[gi1;gi0];

ineq_count=ineq_count+1;

gi2=zeros(1,vnum+1);
gi2(mon_num+j)=2;
gi2(vnum+1)=-1;
gi1=zeros(1,vnum+1);
gi1(mon_num+j)=1;
gi1(vnum+1)=1;
pop.G{ineq_count}=[gi2;gi1];
ineq_count=ineq_count+1;

end

pop.F = [fsum;fvar];

pop.n = vnum;
J={};
I={};




for i=1:ineq_count
J=[J,{i}];

end
pop.J=J

for i=1:cnum
I=[I,{[1:mon_num,mon_num+i]}];

end
cluster_point=1;
for i=1:cluster_count-1
I=[I,{mon_num+cluster_point:mon_num+cluster_point+cluster_size-1}];
cluster_point=cluster_point+cluster_size;
end

 for i=1:cnum+2
    pop.J{cnum+i}=1;
 end

%pop.J={1:ineq_count-1};
%pop.I={1:vnum};
pop.I=I;
pop.J=J;
pop.k=1;
pop.d=2;
save('pop_vars.mat','pop');

sdp = gendata2(pop,'SBSOS');

sol = csol(sdp,'sdpt3');

psol = postproc(pop,sdp,sol);
coeffs=psol.YY{1}(1+1:1+mon_num);
Ws=psol.YY{1}(1+mon_num+1:1+mon_num+cnum);
save('coeff_sol.mat','psol');

weights=psol.YY{1}(1+mon_num+1:1+mon_num+cnum);
ggg
end

%coeffs=psol.xsol(1:mon_num);
%DrawCoeffOnImage;
