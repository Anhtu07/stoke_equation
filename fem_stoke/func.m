function [ u1, u2, p, f1, f2, alpha, nu ] = func(x, y)
    alpha = 0;
    nu = 1;
    pi2 = 2*pi;
    u1=(y.*(y-1).^2 + y.^2.*(y-1)).*(cos(pi2*x)-1);
    u2=pi.*sin(pi2.*x).*y.^2.*(y-1).^2;
    %u1=-cos(pi2*x).*sin(pi2*y)+sin(pi2*y);
    %u2= sin(pi2*x).*cos(pi2*y)-sin(pi2*x);
    p=pi2*(cos(pi2*y)-cos(pi2*x));
    %f1=alpha*u1-nu*pi2*pi2*sin(pi2*y).*(2*cos(pi2*x)-1)+pi2*pi2*sin(pi2*x);
    %f2=alpha*u2+nu*pi2*pi2*sin(pi2*x).*(2*cos(pi2*y)-1)-pi2*pi2*sin(pi2*y);
    f1 = 4*pi^2*sin(2*pi*x) - (12*y - 6).*(cos(2*pi*x) - 1) + 4*pi^2*cos(2*pi*x).*(y.*(y - 1).^2 + y.^2.*(y - 1));
    f2 = 4*y.^2*pi^3.*sin(2*pi*x).*(y - 1).^2 - 2*y.^2*pi.*sin(2*pi*x) - 2*pi*sin(2*pi*x).*(y - 1).^2 - 4*pi^2*sin(2*pi*y) - 4.*y*pi.*sin(2*pi*x).*(2*y - 2);
end

