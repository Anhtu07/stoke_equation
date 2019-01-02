function [A, B1b, B2b, z, omega, tarea]=stok2dp1bmat(p,t,alpha,nu)
    np=size(p,1);
    nt=size(t,1);
    % 1. triangles area
    x21=p(t(:,2),1)-p(t(:,1),1); y12=p(t(:,1),2)-p(t(:,2),2);
    x32=p(t(:,3),1)-p(t(:,2),1); y23=p(t(:,2),2)-p(t(:,3),2);
    x13=p(t(:,1),1)-p(t(:,3),1); y31=p(t(:,3),2)-p(t(:,1),2);
    tarea=(x21.*y31-x13.*y12)/2;
    % 2. x^(T), y^(T), z^(T) and omega
    xt=[x32 x13 x21]; yt=[y23 y31 y12]; zt=(3/20)*alpha*tarea;
    omega=(81/40)*nu*(y23.^2+x32.^2+y31.^2+x13.^2+y23.*y31+x32.*x13)./tarea...
    +(81/280)*alpha*tarea;
    % 3. assembly of stiffness and mass matrices
    Ah=sparse(np,np);
    for i=1:3
        for j=1:i
            Ah=Ah+sparse(t(:,i),t(:,j),alpha*tarea/12,np,np);
        end
    end
    Ah=Ah+Ah.';
    for i=1:3
        for j=1:3
            Ah=Ah+sparse(t(:,i),t(:,j),(nu/4)*xt(:,i).*xt(:,j)./tarea,np,np)...
            +sparse(t(:,i),t(:,j),(nu/4)*yt(:,i).*yt(:,j)./tarea,np,np);
        end
    end
    % 4. assembly of A-zz^t, divergence and pressure matrices
    B1=sparse(np,np); B2=sparse(np,np); E=sparse(np,np);
    B1b=sparse(np,nt); B2b=sparse(np,nt); z=sparse(np,nt);
    for i=1:3
        for j=1:3
            Ah=Ah-sparse(t(:,i),t(:,j),zt.*zt./omega,np,np);
            B1=B1-sparse(t(:,i),t(:,j),yt(:,j)/6+(9/40)*yt(:,i).*zt./omega,np,np);
            B2=B2-sparse(t(:,i),t(:,j),xt(:,j)/6+(9/40)*xt(:,i).*zt./omega,np,np);
            E=E-sparse(t(:,i),t(:,j),(81/1600)*yt(:,i).*yt(:,j)./omega,np,np)...
            -sparse(t(:,i),t(:,j),(81/1600)*xt(:,i).*xt(:,j)./omega,np,np);
        end
        B1b = B1b + sparse(t(:,i), 1:nt, (9/40)*yt(:,i), np, nt);
        B2b = B2b + sparse(t(:,i), 1:nt, (9/40)*xt(:,i), np, nt);
        z = z + sparse(t(:,i), 1:nt, zt, np, nt);
    end

    % 5. Final matrix
    Z=sparse(np,np);
    A=[Ah Z B1';
    Z Ah B2';
    B1 B2 E ];