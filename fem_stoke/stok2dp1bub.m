function [u1b,u2b]=stok2dp1bub(B1b, B2b, z, omega, f1b, f2b, p, u1, u2)
    u1b = (f1b - B1b'*p - z'*u1)./omega;
    u2b = (f2b - B2b'*p - z'*u2)./omega;
end

