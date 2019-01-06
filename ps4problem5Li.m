clear;
clc;
infnorms1 = zeros(10, 1); % preallocate array of infinite norms (unstable)
infnorms2 = zeros(10, 1); % preallocate array of infinite norms (stable)
 
for n=10:10:100 
%{
a. Create the Hilbert matrix Hn of size n (using hilb(n)) and consider its
columns h1, . . . , hn as a basis of Rn. The matrix Hn is non-singular, and
thus its columns indeed form a basis, but it is very close to singular 
(i.e. its columns are close to being linearly dependent), and this leads 
to numerical problems.
%}
    W = hilb(n);

%{
b. Implement the basic Gram–Schmidt algorithm to construct an orthogonal
basis v1, . . . , vn of Rn from h1, . . . , hn. Please don’t use any 
advanced built-in function for orthogonalization (such as orthog), just 
basic matrix operations. At the end of the process normalize your vectors 
so that the basis is orthonormal.
%}

    V = zeros(n, n); % preallocated orthonormal basis 
    V(:,1) = W(:,1) / norm(W(:,1)); % set first vector
    for k=2:n
        sum = 0;
        for i=1:(k-1)
            sum = sum + dot(W(:,k), V(:,i)) / norm(V(:,i))^2 * V(:,i);
        end
        V(:,k) = W(:,k) - sum; % find remaining vectors (orthogonal basis)
        V(:,k) = V(:,k) / norm(V(:,k)); % normalize vectors 
    end
  
%{
c. If vectors v1, . . . , vn obtained in (b) are orthonormal, then V = 
[v1, . . . , vn] must be orthogonal. As a measure of orthogonality, 
compute the infinite norm  deltaV(n) = ||In  VT V ||inf, which is a matrix 
norm (use function norm(A,Inf)). The closer deltaV(n) to zero, the closer 
the columns of V are to being orthogonal.
%}

    V = sym(V); 
    infnorms1(n/10) = norm(V, inf); % compute infinite norm

%{
d. Repeat (b) and (c), to construct an orthonormal basis {u1, . . . , un} 
using the modified Gram–Schmidt algorithm (which is numerically more 
stable) described in lecture 10 (page 46), and compute ?U (n).
%} 
    
    U = zeros(n, n); % preallocated orthonormal basis 
    U(:,1) = W(:,1) / norm(W(:,1)); % set and normalize first vector 
    for k=2:n
        for i=k:n
            % modify W
            W(:,i) = W(:,i) - dot(W(:,i), U(:,k-1)) * U(:,k-1);
        end
        U(:,k) = W(:,k) / norm(W(:,k)); % set and normalize next vector
    end
    infnorms2(n/10) = norm(U, inf); % compute infinite norm
end

%{
e. To compare the basic and modified Gram–Schmidt algorithms, plot deltaV(n)
and deltaU(n) versus n.
%}

x = 10:10:100;
figure(1);
plot(x, infnorms1, 'o', x, infnorms2, 'o');
legend('unstable implementation', 'stable implementation', 'Location', ...
    'northwest');
xlabel('n');
ylabel('Infinite norm');