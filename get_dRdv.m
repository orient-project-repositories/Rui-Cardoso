function [dRdv] = get_dRdv(v)

if v~=[0 0 0]
    R=expm([0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0]);
else
    R=eye(3);
end

dRdv(:,:,1)=R*[0 0 0; 0 0 -1; 0 1 0];
dRdv(:,:,2)=R*[0 0 1; 0 0 0; -1 0 0];
dRdv(:,:,3)=R*[0 -1 0;1 0 0; 0 0 0];

end

