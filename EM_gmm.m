function [ label, model ] = EM_gmm( X, k )
% EM.
% Input--------------------------------------------------------------------
%   X: d x n data matrix
%   k: Number of Gaussian Components
% Output-------------------------------------------------------------------
%   label: 1 x n cluster label
%   model: trained model structure(mu,Sigma,w)

%Variable Initialisations--------------------------------------------------
maximum_iter = 500; % maximum Number of iterations to go for
%The likelihood function for each iteration initialised with -inf
likelihood = -inf(1,maximum_iter);
d = size(X,1);%Dimension of each vector
n = size(X,2);%Number of vectors
label = randi([1 k],1,n);%Initialising the labels
%Wni is n x k matrix. Each row corresponds to one vector.
Z = full(sparse(1:n,label,1,n,k,n));
tol = 1e-6;%tolerance

%Iterations----------------------------------------------------------------
for i = 2:maximum_iter
    [~,label(1,:)] = max(Z,[],2);%Extracting labels from Z
    Z = Z(:,unique(label));%Removing the empty clusters
    k = size(Z,2);%Finding the updated number of clusters
    %Finding Previous mu,Sigma and w
    nk = sum(Z,1);%Number of vectors in each cluster
    w = nk/n;%Weights for each gaussian component (1 x k)
    mu = bsxfun(@times, X*Z, 1./nk);%Finding mean of clusters (d x k)
    Sigma = zeros(d,d,k);%Sigma Matrices for each Cluster (d x d x k)
    z = sqrt(Z);%Finding sqrt of Z
    for j = 1:k
        Xo = bsxfun(@minus,X,mu(:,j));%Subtracting mean
        Xo = bsxfun(@times,Xo,z(:,j)');%Multiplying each column with r
        Sigma(:,:,j) = Xo*Xo'/nk(j)+eye(d)*(1e-6);%Variance Limiting
    end
    model.mu = mu;%Assigning mean to the model
    model.Sigma = Sigma;%Assigning Sigma to the model
    model.w = w;%Assigning weights to the model
    %Finding the new Z
    for j = 1:k
        %Saving the log values of pdf to manage small value product issues
        Z(:,j) = loggausspdf(X,mu(:,j),Sigma(:,:,j));
    end
    Z = bsxfun(@plus,Z,log(w));%Adding log(w) to multiply
    %Finding sum of loglikelihood over all components for each vector
    T = logsumexp(Z,2);
    likelihood(i) = sum(T)/n; %loglikelihood calculation
    Z = exp(bsxfun(@minus,Z,T));%getting the Z back
    %Checking if tolerance reached
    if abs(likelihood(i)-likelihood(i-1)) < tol*abs(likelihood(i))
        break;
    end
    %END
 end

%function for calculating the log of gauss pdf-----------------------------
function y = loggausspdf(X, mu, Sigma)
d = size(X,1);
X = bsxfun(@minus,X,mu);
[U,p]= chol(Sigma);
if p ~= 0
    error('ERROR: Sigma is not PD.');
end
Q = U'\X;
q = dot(Q,Q,1);% quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(U)));% normalization constant
y = -(c+q)/2;


%function to compute log(sum(exp)) to avoid numerical underflow------------
function s = logsumexp(X, dim)
% subtract the largest in each dim
y = max(X,[],dim);
s = y+log(sum(exp(bsxfun(@minus,X,y)),dim));   
i = isinf(y);
if any(i(:))
    s(i) = y(i);
end



