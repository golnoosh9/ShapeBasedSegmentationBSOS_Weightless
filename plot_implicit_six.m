function  plot_implicit_six(c)
% c2=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% c=c.*c2;


%[x,y]=meshgrid(-3:0.01:3);

syms x y

func=c(1)+c(2)*x+c(3)*y+c(4)*power(x,2)+c(5)*x*y+c(6)*power(y,2)+...
    c(7)*power(x,3)+c(8)*power(x,2)*y+c(9)*x*power(y,2)+c(10)*power(y,3)...
    +c(11)*power(x,4)+c(12)*power(x,3)*y+c(13)*power(x,2)*power(y,2)+...
    c(14)*x*power(y,3)+c(15)*power(y,4)+c(16)*power(x,5)+c(17)*power(x,4)*y+...
    c(18)*power(x,3)*power(y,2)+c(19)*power(x,2)*power(y,3)+c(20)*x*power(y,4)+...
    c(21)*power(y,5)+c(22)*power(x,6)+c(23)*power(x,5)*y+c(24)*...
    power(x,4)*power(y,2)+c(25)*power(x,3)*power(y,3)+c(26)*power(x,2)*...
    power(y,4)+c(27)*x*power(y,5)+c(28)*power(y,6);


ezplot(func,[-2,2]);
