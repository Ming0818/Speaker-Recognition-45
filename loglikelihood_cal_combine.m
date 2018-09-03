function [ Y ] = loglikelihood_cal_combine(X,parameter_solo,parameter_whsp,parameter_fast,M,n)
% Input---------------------------------------------------------------
% X : coeff matrix
% parameter_solo : gmm model for solo speech
% parameter_whsp : gmm model for whisper speech
% parameter_fast : gmm model for fast speech
% Output--------------------------------------------------------------
% Y : row vector having log likelihood for all the M samples
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
 obj = gmdistribution((parameter_solo.mu)',parameter_solo.Sigma,parameter_solo.w);
 pdf_val = pdf(obj,X');
 log_val(1,:) = log(pdf_val);
 obj = gmdistribution((parameter_whsp.mu)',parameter_whsp.Sigma,parameter_whsp.w);
 pdf_val = pdf(obj,X');
 log_val(2,:) = log(pdf_val);
 obj = gmdistribution((parameter_fast.mu)',parameter_fast.Sigma,parameter_fast.w);
 pdf_val = pdf(obj,X');
 log_val(3,:) = log(pdf_val);
 log_val = max(log_val,[],1);
 for i = 1:M
    Y(1,i) = sum(log_val(i:i+n-1));
 end



end

