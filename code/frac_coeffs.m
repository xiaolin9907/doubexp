function T = frac_coeffs(N, mu, b, sgn , NN)

if nargin < 5
    NN = 100;
end
nn = N; 
K = kernel(mu, b, NN,sgn);
[ku, kd, kv] = svd(K);

kd = diag(kd);
dim_K = length(find(kd>1e-14));

if N < dim_K
    N = dim_K;
end

if N < length(kd)-dim_K 
    N = length(kd)-dim_K;
end


%%
% dim_K = 1; kv = zeros(100,1);kv(1)=1;
%%
% norm(ku(:,1:dim_K)*diag(kd(1:dim_K))*kv(:,1:dim_K)'-K)

kv = sparse([kv(:,1:dim_K);zeros(N+1+dim_K - length(kd),dim_K)]);
ku = [ku(:,1:dim_K)*diag(kd(1:dim_K));zeros(N+1+dim_K - length(kd),dim_K)];

ku = chebtech2.coeffs2vals(ku);


T = zeros(N+dim_K+1,N);
% 
II = Init_val(N+dim_K+1, mu);

Vold = zeros(N+dim_K+1, dim_K);
V = Vold;
Vnew = V;
if sgn == 1
    Vold(1,:) = II * kv;
    V(1,:) = 4*II * ultraS.multmat(size(kv,1),[1/2;1/2],0)*kv;
    V(2,:) = 4*Vold(1,:) - V(1,:);
    T(:,2) = cheb_mul_sum(Vold, ku);T(NN+2:end,2)=0;
    T(:,3) = cheb_mul_sum(V, ku);T(NN+3:end,3)=0;
    Vnew(3,:) = II * ultraS.multmat(size(kv,1),[9/4;-3;3/4],0)*kv;
    Vnew(2,:) = II * ultraS.multmat(size(kv,1),[3;0;-3],0)*kv;
    Vnew(1,:) = II * ultraS.multmat(size(kv,1),[15/4;3;9/4],0)*kv;
    T(:,4) = cheb_mul_sum(Vnew, ku); T(NN+3:end,4)=0;
else
    Vold(1,:) = II * kv;
    T(:,2) = cheb_mul_sum(Vold, ku);T(NN+2:end,2)=0;
    V(1,:) = -4*II * ultraS.multmat(size(kv,1),[1/2;1/2],0)*kv;
    V(2,:) = 4*Vold(1,:) + V(1,:);
    T(:,3) = cheb_mul_sum(V, ku);T(NN+3:end,3)=0;
end

% ODE初值 \int_0^1 U_{n-1}(2t-1) t^mu g_k(t) dt
Init = zeros(N, dim_K);
Init(1,:) = Vold(1,:);
multemp = ultraS.multmat(size(kv,1),[0;2],0);

Iold = kv;
I = multemp * Iold;
Init(2,:) = II * I;
for k = 3:N
    Inew = multemp*I - Iold;
    Init(k,:) = II * Inew;    
    Iold = I; I = Inew;
end
if sgn == -1
    Init = diag((-1).^(0:N-1).*(1:N))*Init;
else
     Init = diag((1:N))*Init;
end
% ODE
if sgn == -1

    Tra = spdiags([ones(N+3+dim_K,1),-ones(N+3+dim_K,1)],[-1;0],N+dim_K+2,N+dim_K+1);
    MD = ultraS.multmat(N+2+dim_K,[1;1],1)*ultraS.diffmat(N+2+dim_K,1);
    MDD = MD(1:end-1,:)*Tra;

    C = ultraS.convertmat(N+1+dim_K,0,0);
    CT = C*Tra(1:end-1,:);

    Tra = Tra(1:N+1+dim_K,:);
    MD = MD(1:end-1,1:end-1);


    for k = 2:N-1    
        right = 2*(k+1) .* (C*V) + (k+1)*k/(k-1) .* (C*Vold) +((k+1)/(k-1)).*(MD*Vold);
        right(1,:) = right(1,:) + k*Init(k+1,:);
        OP = MDD - k .* CT;
        Vnew = OP \ right;
        Vnew = Tra*Vnew;
        Vnew(1,:) = Vnew(1,:) + Init(k+1,:);
        T(:,k+2) = cheb_mul_sum(Vnew, ku);
        T(NN+k+2:end,k+2)=0;
        Vold = V;
        V = Vnew;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    Tra = spdiags([ones(N+3+dim_K,1),ones(N+3+dim_K,1)],[-1;0],N+dim_K+2,N+dim_K+1);
    MD = ultraS.multmat(N+2+dim_K,[1;-1],1)*ultraS.diffmat(N+2+dim_K,1);
    MDD = MD(1:end-1,:)*Tra;

    C = ultraS.convertmat(N+1+dim_K,0,0);
    CT = C*Tra(1:end-1,:);

    Tra = Tra(1:N+1+dim_K,:);
    MD = MD(1:end-1,1:end-1);

        Vold = V;
        V = Vnew;
    for k = 3:N-1    
        right = 2*(k+1) .* (C*V) - (k+1)*k/(k-1) .* (C*Vold) +((k+1)/(k-1)).*(MD*Vold);
        right(1,:) = right(1,:) - k*Init(k+1,:);
        OP = MDD + k .* CT;
        Vnew = OP \ right;
        Vnew = Tra*Vnew;
        Vnew(1,:) = Vnew(1,:) + Init(k+1,:);
        T(:,k+2) = cheb_mul_sum(Vnew, ku);
        T(NN+k+2:end,k+2)=0;
        Vold = V;
        V = Vnew;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if sgn == -1
    T = ultraS.multmat(N+1+dim_K,[1;1],0)*T;
    x = chebpts(size(K,1));
    y_mid = pi * sinh(b * x);
    ln_f = zeros(size(x));
    idx_pos = (y_mid >= 0);
    ln_f(idx_pos) = log(2) - log1p(exp(-y_mid(idx_pos)));
    idx_neg = ~idx_pos;
    ln_f(idx_neg) = log(2) - (-y_mid(idx_neg) + log1p(exp(y_mid(idx_neg))));
    y = exp(mu * ln_f);
    fc = chebtech2.vals2coeffs(y);
    for k = 1:size(T,2)
        T(1:size(K,1),k) = T(1:size(K,1),k) + (-1)^(k-1) * fc;
    end
else
    T = ultraS.multmat(N+1+dim_K,[-1;1],0)*T;
    x = chebpts(size(K,1));
    y_mid = pi * sinh(b * x);
    ln_f = zeros(size(x));
    idx_pos = (y_mid >= 0);
    ln_f(idx_pos) = log(2) - (y_mid(idx_pos) + log1p(exp(-y_mid(idx_pos))));
    idx_neg = ~idx_pos;
    ln_f(idx_neg) = log(2) - log1p(exp(y_mid(idx_neg)));
    y = exp(mu * ln_f);
    fc = chebtech2.vals2coeffs(y);
    for k = 1:size(T,2)
        T(1:size(K,1),k) = T(1:size(K,1),k) + fc;
    end
end
T = T ./ gamma(1+mu);
T = T(1:nn,1:nn);


end
% 


function II = Init_val(N, mu)
    II = zeros(1,N);

% gauss-jacobi 积分点和权重
    [GJx, GJw] = jacpts(N, 0, mu, [0,1]);

% II(k) = \int_0^1 t^mu T_{k-1}(2t-1) dt
    Iold = ones(N,1);
    d = 2*GJx-1;
    I = d;
    II(1) = dot(GJw, Iold);
    II(2) = dot(GJw, I);
    for k = 3:N
        Inew = 2*d.*I - Iold;
        II(k) = dot(GJw, Inew);
        Iold = I; I = Inew;
    end

end



function R = cheb_mul_sum(G, ku)
    G = chebtech2.coeffs2vals(G);
    M = G .* ku;
    R = sum(M,2);
    R = chebtech2.vals2coeffs(R);
end
% 
function Vc = kernel(mu,b,NN,sgn)
    n=NN;
    x = chebpts(n);
    t = chebpts(n,[0,1]);


    V = compute_G_matrix(x, t, b, mu, sgn);

    Vc1 = chebtech2.vals2coeffs(V');
    Vc = chebtech2.vals2coeffs(Vc1');
    
end


function G_matrix = compute_G_matrix(x, t, b, mu, sgn)
    x = x(:);
    t = t(:)';

    if sgn == -1
        s = x - (1 + x) .* t;
        w_factor = 1 + x;          % 对应 (x-s)/t
        idx_boundary = (x == -1);  % 边界点在 -1
    elseif sgn == 1
        s = x + (1 - x) .* t;
        w_factor = 1 - x;          % 对应 (s-x)/t
        idx_boundary = (x == 1);   % 边界点在 1
    else
        error('sgn 必须为 1 或 -1');
    end
    
    % 提取无奇异性的中间变量
    w = (b / 2) .* w_factor .* t;
    cosh_sx = cosh( (b / 2) .* (x + s) );
    z = pi .* cosh_sx .* sinh(w);
    
    % 1. 分子提取出来的稳定项
    term1 = log( (pi * b / 2) .* w_factor ); % N x 1，自动广播
    term2 = logcosh_stable( (b / 2) .* (x + s) );
    term3 = log_sinhc_stable(z);
    term4 = log_sinhc_stable(w);
    
    % 2. 原公式分母的两个 cosh 项
    den1 = logcosh_stable( (pi / 2) .* sinh(b .* x) ); % N x 1
    den2 = logcosh_stable( (pi / 2) .* sinh(b .* s) ); % N x M
    
    % 3. 组合对数并乘上指数 mu
    ln_G = term1 + term2 + term3 + term4 - den1 - den2;
    
    % 4. 还原为实际数值
    G_matrix = exp(mu .* ln_G);
    
    % 边界物理意义修正：
    % 强制对应的端点导数为 0 (从而整个 G_matrix 在该行为 0)
    if any(idx_boundary)
        G_matrix(idx_boundary, :) = 0;
    end
end

% ================= 辅助稳定函数 (防 Inf/NaN 核心) ================= %

function y = log_sinhc_stable(v)
    % 直接计算 ln(sinh(v)/v)
    y = zeros(size(v));
    idx1 = abs(v) > 1e-8;
    
    % 大范围计算
    y(idx1) = logabs_sinh_stable(v(idx1)) - log(abs(v(idx1)));
    
    % 极小范围泰勒展开
    idx2 = ~idx1 & (v ~= 0);
    v2 = v(idx2).^2;
    y(idx2) = v2 / 6 - v2.^2 / 180;
end

function y = logabs_sinh_stable(x)
    % 稳定的 ln|sinh(x)|，防止 x > 710 时 sinh 溢出
    abs_x = abs(x);
    y = zeros(size(x));
    
    idx_large = abs_x > 20;
    y(idx_large) = abs_x(idx_large) - log(2) + log1p(-exp(-2 .* abs_x(idx_large)));
    y(~idx_large) = log(abs(sinh(x(~idx_large))));
end

function y = logcosh_stable(x)
    % 稳定的 ln(cosh(x))
    abs_x = abs(x);
    y = abs_x - log(2) + log1p(exp(-2 .* abs_x));
end









