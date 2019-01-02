function val = valinterp(xy, hx, hy, nx, u, x, y)
    idx_x = fix(x/hx) + 1;
    idx_y = fix(y/hy) + 1;
    idx_l = (idx_y-1)*nx + idx_x;
    idx_u = idx_y*nx + idx_x;
    if x + y - sum(xy(idx_u, :)) < 0
        idx_1=idx_l;idx_2=idx_l+1;idx_3=idx_u;
    else
        idx_1=idx_l+1;idx_2=idx_u+1;idx_3=idx_u;
    end
    idx=[idx_1, idx_2, idx_3];
    A = [[1,1,1]',xy(idx, :)];
    B = [u(idx_1), 0, 0; 0, u(idx_2), 0; 0, 0, u(idx_3)];
    c = A\B;
    val = sum([1, x, y]*c);
end

