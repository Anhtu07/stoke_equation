function [b, f1b, f2b]=stok2dp1brhs(p,t,alpha,nu,f1,f2)
    %STOK2DP1BRHS P1-Bubble/P1 element:
    %Assembly of the right-hand side
    np=size(p,1);
    % 1. (f1,f2) at the center of triangles
    if (length(f1)==np), f1t=(f1(t(:,1))+f1(t(:,2))+f1(t(:,3)))/3;
    else f1t=f1; end
    if (length(f2)==np), f2t=(f2(t(:,1))+f2(t(:,2))+f2(t(:,3)))/3;
    else f2t=f2; end
    % 2. triangles area
    x21=p(t(:,2),1)-p(t(:,1),1); y12=p(t(:,1),2)-p(t(:,2),2);
    x32=p(t(:,3),1)-p(t(:,2),1); y23=p(t(:,2),2)-p(t(:,3),2);
    x13=p(t(:,1),1)-p(t(:,3),1); y31=p(t(:,3),2)-p(t(:,1),2);
    tarea=(x21.*y31-x13.*y12)/2;
    % 3. x^(T), y^(T), z and omega
    xt=[x32 x13 x21]; yt=[y23 y31 y12]; zt=(3/20)*alpha*tarea;
    omega=(81/40)*nu*(y23.^2+x32.^2+y31.^2+x13.^2+y23.*y31+x32.*x13)./tarea...
    +(81/280)*alpha*tarea;
    % 4. assembly of the right-hand side
    f1b=(9/20)*tarea.*f1t./omega; f2b=(9/20)*tarea.*f2t./omega;
    b1=sparse(np,1); b2=sparse(np,1); bb=sparse(np,1);
    for i=1:3
        b1=b1+sparse(t(:,i),1,(1/3)*f1t.*tarea-f1b.*zt,np,1);
        b2=b2+sparse(t(:,i),1,(1/3)*f2t.*tarea-f2b.*zt,np,1);
        bb=bb-sparse(t(:,i),1,(9/40)*f1b.*yt(:,i)+(9/40)*f2b.*xt(:,i),np,1);
    end
    b1=full(b1); b2=full(b2); bb=full(bb);
    % 5. right-hand side
    b=[b1; b2; bb];