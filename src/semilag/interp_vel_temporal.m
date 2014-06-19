function [ut, vt, wt] = interp_vel_temporal(u,v,w,t,tq,INTERP_TYPE)
%INTERP_VEL_TEMPORAL Temporal interpolation of velocity values.
%
% tq is the time that velocity values are queried.
% TODO: vectorize the 1D interpolation

xsize = size(u,1);
ysize = size(u,2);
zsize = size(u,3);
tsize = size(u,4);
ut = zeros(xsize, ysize, zsize, size(tq,2));
vt = ut;
wt = ut;
for i=1:xsize
    for j=1:ysize
        for k=1:zsize
            ulist = zeros(size(t));
            vlist = ulist;
            wlist = ulist;
            for cnt=1:tsize
                ulist(cnt) =u(i,j,k,cnt);
                vlist(cnt) =v(i,j,k,cnt);
                wlist(cnt) =w(i,j,k,cnt);
            end
            ut(i,j,k,:) = interp1(t,ulist,tq,INTERP_TYPE);
            vt(i,j,k,:) = interp1(t,vlist,tq,INTERP_TYPE);
            wt(i,j,k,:) = interp1(t,wlist,tq,INTERP_TYPE);
        end
    end
end
end