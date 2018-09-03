function [ Y ] = loglikelihood_cal( X , parameter, N ,M, n )
% Input---------------------------------------------------------------
% X : coeff matrix
% parameters : gmm model
% Output--------------------------------------------------------------
% Y : row vector having log likelihood for al the M samples
Y = zeros(1,M);%Initialising the vector
% size_mu = size(parameter.mu)
% k = size_mu(2);
% mu = parameter.mu;
% sigma = parameter.Sigma;
% w = parameter.w;
% size_X = size(X)
% for i = 1:k
%     pdf_matrix(:,i) = mvnpdf(X',(mu(:,i))',sigma(:,:,i));
% end
 obj = gmdistribution((parameter.mu)',parameter.Sigma,parameter.w);
 pdf_val = pdf(obj,X');
 log_val = log(pdf_val);
 for i = 1:M
    Y(1,i) = sum(log_val(i:i+n-1));
 end
