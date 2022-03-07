function [vn,in] = nonlinearSolver(p,K,vn, nleq)
b = 1;
iter = 0;
step = 1;
tol = 1e-5;

[in, jacobian] = nleq(vn);

F = p + K*in - vn;
while(norm(step) > tol && iter < 100)
    [in, jacobian] = nleq(vn);
    J = K*jacobian - eye(length(vn));

    step = J\F;
    vn_new = vn - b*step;
    [in_new, jacobian] = nleq(vn_new);
    Fnew = p + K*in_new - vn_new;
    if(norm(Fnew) < norm(F))
        F = Fnew;
        vn = vn_new;
        b = 1;
    else
        b = b/2;
    end
    
    iter = iter + 1;
end

% % Non-damped Newton-Raphson method
% tol = 1e-5;
% maxCount = 100;
% step = 1;
% counter = maxCount;
% 
% 
%     while norm(step) > tol && counter > 0
%         [in, jacobian] = nleq(vn);
%         J = K*jacobian - eye(length(vn));    
%     
%         ni = p + K*in - vn; 
%         step = J\ni;
%         vn = vn - step;
%         counter = counter - 1;
%     end
% end

