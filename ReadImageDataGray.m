

%file='aircraft_1test.png';
 sample_num=0; 
 cnum=79;
im = imread(file);;
sim=size(im);
sx=sim(1);
sy=sim(2);
max_dir=max(sx,sy);
maximum=max_dir/2;2
%im=rgb2gray(im);
im= im2bw(im);


% imr=im(:,:,1);
% img=im(:,:,2);
% im=(im2double(imr)./im2double(imr+img));

BW=edge(im,'canny');




varsum=0;
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

im(:,:)=0;
s2=size(l0);
for(i=1:s2(1))
    im(l0(i,1),l0(i,2))=255;
end


meanx=sx/2;
meany=sy/2;
nn=ones(size(l0));
h=[meanx;meany];
nn(:,1)=nn(:,1).*meanx;
nn(:,2)=nn(:,2).*meany;

l0_new=l0-nn;
