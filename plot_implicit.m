function [cnt]=  plot_implicit(c)
% c2=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% c=c.*c2;


%[x,y]=meshgrid(-3:0.01:3);

syms x y

func=c(1)+c(2)*x+c(3)*y+c(4)*power(x,2)+c(5)*x*y+c(6)*power(y,2)+c(7)*power(x,3)+c(8)*power(x,2)*y+c(9)*x*power(y,2)+c(10)*power(y,3)+c(11)*power(x,4)+c(12)*power(x,3)*y+c(13)*power(x,2)*power(y,2)+c(14)*x*power(y,3)+c(15)*power(y,4);


m=ezplot(func,[-2,2]);

%cnt=get(m,'contourMatrix');

% c2=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% c=c.*c2;
% syms x y
% 
% func=c(1)+c(2)*x+c(3)*y+c(4)*power(x,2)+c(5)*x*y+c(6)*power(y,2)+c(7)*power(x,3)
% +c(8)*power(x,2)*y+c(9)*x*power(y,2)+c(10)*power(y,3)+c(11)*power(x,4)+c(12)*power(x,3)*y
% +c(13)*power(x,2)*power(y,2)+c(14)*x*power(y,3)+c(15)*power(y,4);
% 
% 
% 
% ezplot(func==0, [-1 1]);
% hold on
% xlabel('x axis')
% ylabel('y axis')
% [x,y]=meshgrid(-100:.01:100);
% func=c(1)+c(2)*x+c(3)*y+c(4)*power(x,2)+c(5)*x*y+c(6)*power(y,2)+c(7)*power(x,3)
% +c(8)*power(x,2)*y+c(9)*x*power(y,2)+c(10)*power(y,3)+c(11)*power(x,4)+c(12)*power(x,3)*y
% +c(13)*power(x,2)*power(y,2)+c(14)*x*power(y,3)+c(15)*power(y,4);
% 
% level=0.02;
% contour(x,y,func,1);
