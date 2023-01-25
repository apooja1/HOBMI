function [H,S,D]=sobi_2(X,n,p)

[m,N,ntrials]=size(X);
if nargin<1 || nargin > 3
  help sobi
elseif nargin==1
 n=m; % Source detection (hum...)
 p=min(100,ceil(N/3)); % Number of time delayed correlation matrices to be diagonalized 
 % Authors note: For noisy data, use at least p=100 the time-delayed covariance matrices.
elseif nargin==2
 p=min(100,ceil(N/3)); % Default number of correlation matrices to be diagonalized
                       % Use < 100 delays if necessary for short data epochs
end

%
% Make the data zero mean
%
X(:,:)=X(:,:)-kron(mean(X(:,:)')',ones(1,N*ntrials)); 

%
% Pre-whiten the data based directly on SVD
%
[UU,S,VV]=svd(X(:,:)',0);
Q= pinv(S)*VV';
X(:,:)=Q*X(:,:);

% Alternate whitening code
% Rx=(X*X')/N;
% if m<n % assumes white noise
%   [U,D]=eig(Rx); 
%   [puiss,k]=sort(diag(D));
%   ibl= sqrt(puiss(n-m+1:n)-mean(puiss(1:n-m)));
%    bl = ones(m,1) ./ ibl ;
%   BL=diag(bl)*U(1:n,k(n-m+1:n))';
%   IBL=U(1:n,k(n-m+1:n))*diag(ibl);
% else    % assumes no noise
%    IBL=sqrtm(Rx);
%    Q=pinv(IBL);
% end
% X=Q*X;

%
% Estimate the correlation matrices
%
 k=1;
 pm=p*m; % for convenience   p is the number of covariance matrix to be optimized 
 for u=1:m:pm
   k=k+1; 
   for t = 1:ntrials 
       if t == 1
           Rxp=X(:,k:N,t)*X(:,1:N-k+1,t)'/(N-k+1)/ntrials;  %pXp
       else
           Rxp=Rxp+X(:,k:N,t)*X(:,1:N-k+1,t)'/(N-k+1)/ntrials;
       end
   end
   M(:,u:u+m-1)=norm(Rxp,'fro')*Rxp;  % Frobenius norm =
 end                                  % sqrt(sum(diag(Rxp'*Rxp)))

%
% Perform joint diagonalization
%
epsil=1/sqrt(N)/100; 
encore=1; 
V=eye(m);
while encore 
 encore=0;
 for p=1:m-1
  for q=p+1:m
   % Perform Givens rotation
   g=[   M(p,p:m:pm)-M(q,q:m:pm)  ;
         M(p,q:m:pm)+M(q,p:m:pm)  ;
      i*(M(q,p:m:pm)-M(p,q:m:pm)) ];
	  [vcp,D] = eig(real(g*g')); 
          [la,K]=sort(diag(D));
   angles=vcp(:,K(3));
   angles=sign(angles(1))*angles;
   c=sqrt(0.5+angles(1)/2);
   sr=0.5*(angles(2)-j*angles(3))/c; 
   sc=conj(sr);
   oui = abs(sr)>epsil ;
   encore=encore | oui ;
   if oui  % Update the M and V matrices 
    colp=M(:,p:m:pm);
    colq=M(:,q:m:pm);
    M(:,p:m:pm)=c*colp+sr*colq;
    M(:,q:m:pm)=c*colq-sc*colp;
    rowp=M(p,:);
    rowq=M(q,:);
    M(p,:)=c*rowp+sc*rowq;
    M(q,:)=c*rowq-sr*rowp;
    temp=V(:,p);
    V(:,p)=c*V(:,p)+sr*V(:,q);
    V(:,q)=c*V(:,q)-sc*temp;
   end%% if
  end%% q loop
 end%% p loop
end%% while

%
% Estimate the mixing matrix 
%
H = pinv(Q)*V; 

%
% Estimate the source activities
%
if nargout>1
  S=V'*X(:,:); % estimated source activities
end