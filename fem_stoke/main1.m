tic;
Pen=10^10;
nx = 1024;
ny = 1024;
pex=0.317;pey=0.910;
[u1exy, u2exy, pexy, f1, f2, alpha, nu] = func(pex, pey);
[xy,th,pbx,pby,hx,hy]=kpde2dumsh(0,3,0,1,nx,ny);
ibc1=union(pbx(:,1),pbx(:,2)); ibc2= union(pby(:,1),pby(:,2));
ibcd=union(ibc1,ibc2);
np=size(xy,1); nt=size(th,1);
fprintf('So diem=%4d So tam giac=%4d \n',np,nt)
x=xy(:,1); y=xy(:,2); pi2=2*pi;
[u1e, u2e, pe, f1, f2, alpha, nu] = func(x, y);
[A, B1b, B2b, z, omega, tarea]=stok2dp1bmat(xy,th,alpha,nu);
[b, f1b, f2b]=stok2dp1brhs(xy,th,alpha,nu,f1,f2);
A(ibcd,ibcd) =A(ibcd,ibcd)+Pen*speye(length(ibcd));
A(np+ibcd,np+ibcd)=A(np+ibcd,np+ibcd)+Pen*speye(length(ibcd));
b(ibcd)=0; b(np+ibcd)=0;
A(2*np+1,2*np+1)=Pen; b(2*np+1)=0;
toc;
u=A\b;
u1=u(1:np); u2=u(np+1:2*np); p=u(2*np+1:3*np);
toc;
format long;


u1_xy = valinterp(xy, hx, hy, nx, u1, pex, pey);
u2_xy = valinterp(xy, hx, hy, nx, u2, pex, pey);
p_xy = valinterp(xy, hx, hy, nx, p, pex, pey);
fprintf('Sai so diem: u1 :');
peu1 = abs(u1exy - u1_xy)
fprintf('Sai so diem: u2');
peu2 = abs(u2exy - u2_xy)
fprintf('Sai so diem: p');
pep = abs(pexy - p_xy)



du1 = sparse(nt, 1);
du2 = sparse(nt, 1);
dp = sparse(nt, 1);

for i=1:3
    du1 = du1 + (u1(th(:, i)) - u1e(th(:, i))).^2/3;
    du2 = du2 + (u2(th(:, i)) - u2e(th(:, i))).^2/3;
    dp = dp + (p(th(:, i)) - pe(th(:, i))).^2/3;
end

l2e_u1 = sqrt(du1'*tarea);
l2e_u2 = sqrt(du2'*tarea);
l2e_p = sqrt(dp'*tarea);
fprintf('Sai so L2 cua u1: ');
l2e_u1
fprintf('Sai so L2 cua u2: ');
l2e_u2
fprintf('Sai so L2 cua p: ');
l2e_p

subplot(121);
trisurf(th,xy(:,1),xy(:,2),u1e,'facecolor','interp');

subplot(122);
trisurf(th,xy(:,1),xy(:,2),u1,'facecolor','interp');

