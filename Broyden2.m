format short
syms X Y;
f = X^2 + 2*X^2*Y^3+ Y^2;
x(1) = 1;
y(1) = 1;
x1=[ x(1) , y(1)] ;
e = 10^(-8); 
geradian = [subs(diff(f, X),[X,Y], [x(1),y(1)]) subs(diff(f, Y), [X,Y], [x(1),y(1)])];
% geradian
disp(geradian)
dx2 = diff(diff(f, X),X);
dy2 = diff(diff(f, Y),Y);
dxy = diff(diff(f, X),Y);
Hessian = [ subs(dx2, [X,Y], [x(1),y(1)]), subs(dxy, [X,Y], [x(1),y(1)]);...
    subs(dxy, [X,Y], [x(1),y(1)]), subs(dy2, [X,Y], [x(1),y(1)])]; % Hessian

inverse = inv(Hessian); 
disp(inverse)
i = 1;
k=1;

while k<4
    S = [x(i),y(i)]';
    p = S -inverse*geradian';
    x(i+1)=p(1,1);
    y(i+1) = p(2,1);
    g=geradian;
    a = [ x(i+1) , y(i+1) ]' - [x(i),y(i)]';
    i = i+1;
    geradian = [subs(diff(f, X),[X,Y], [x(i),y(i)]) subs(diff(f, Y), [X,Y], [x(i),y(i)])]; %geradian
    y0= (geradian - g)';
    m= y0- (Hessian*a);
    h=Hessian;
    Hessian= [subs(dx2, [X,Y], [x(i),y(i)]), subs(dxy, [X,Y], [x(i),y(i)]);...
        subs(dxy, [X,Y], [x(i),y(i)]), subs(dy2, [X,Y], [x(i),y(i)])];
    Hessian = h +  1/(norm(a,2) ^ 2 ) * m .*geradian';
    
    inverse = inv(Hessian);
    k=k+1;
end

disp(x(i))
disp(y(i))