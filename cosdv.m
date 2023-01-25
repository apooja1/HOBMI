function [metric,numDomains] = cosdv(x, y)
%COSF - a copula based test of dependence called the copula statistic. See
%       <TODO: insert paper reference here>
% Inputs:
%  x - realizations of the first random variable
%  y - realizations of the second random variable
% Outputs
%  metric - a real-valued number between 0 and 1 which represents the
%           statistical dependence between x and y
%
% TODO:
%   [ ] - Make this multivariate.
%**************************************************************************
%*                                                                        *
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>                *
%*                                                                        *
%* This program is free software: you can redistribute it and/or modify   *
%* it under the terms of the GNU General Public License as published by   *
%* the Free Software Foundation, either version 3 of the License, or      *
%* (at your option) any later version.                                    *
%*                                                                        *
%* This program is distributed in the hope that it will be useful,        *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of         *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
%* GNU General Public License for more details.                           *
%*                                                                        *
%* You should have received a copy of the GNU General Public License      *
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
%*                                                                        *
%**************************************************************************

general_debug_print = 0;
lambda_debug_print = 0;

u = pobs(x);
v = pobs(y);

[u,I] = sort(u);
v = v(I);

n = length(x);
E_copula = zeros(1,n);

for ii=1:n
    E_copula(ii) = sum(u(ii)>=u & v(ii)>=v)/(n+1);
end
M_copula = min(u,v);
W_copula = max(u+v-1,0);
Pi_copula = u.*v;

score = 0; wanc = 1;

if(E_copula(1)==Pi_copula(1))
    wanc = 0;
end
if( (E_copula(1)==M_copula(1)) || (E_copula(1)==W_copula(1)) || (E_copula(1)==E_copula(2)) )
    wanc = 1;
end
indep = (mean(abs(diff(E_copula)))>0.12);
ss_minIdx = 1;
numDomains = 0;
for ii=2:n
    ss = E_copula(ss_minIdx:ii);
    if(~issorted(ss) && ~issorted(fliplr(ss)))      % just a way of testing if they are tied?
        
        if(general_debug_print)
            fprintf('[cosdv] -- Domain -> %d:%d\n',ss_minIdx,ii-1);
        end
        
        if(E_copula(ii-1)==E_copula(ii-2))
            jj = ii-2;
        else
            jj = ii-1;
        end
        
        % compute the relative distance function
        if (E_copula(jj)>Pi_copula(jj))
            w = (E_copula(jj)-Pi_copula(jj))/(M_copula(jj)-Pi_copula(jj));
            if(lambda_debug_print)
                fprintf('[cosdv] -- >> a\n');
            end
        else
            if (E_copula(jj)<Pi_copula(jj))
                w=(Pi_copula(jj)-E_copula(jj))/(Pi_copula(jj)-W_copula(jj));
                if(lambda_debug_print)
                    fprintf('[cosdv] -- >> b\n');
                end
            else
                w=0;
                if(lambda_debug_print)
                    fprintf('[cosdv] -- >> c\n');
                end
            end
        end
      
        if (( (round(E_copula(ii-1),2)==round(M_copula(ii-1),2)) || ...
              (round(E_copula(ii-1),2)==round(W_copula(ii-1),2)) || ...
              (E_copula(ii-1)==(1/(n+1)))) && (length(ss)>4) )
            if(lambda_debug_print)
                fprintf('[cosdv] -- >> d\n');
            end
            w=1;
        end
      
        if (( (E_copula(ii)==E_copula(ii-2)) || ...
              (E_copula(ii-1)==E_copula(ii-2))) && (length(ss)>4) ) 
            if(lambda_debug_print)
                fprintf('[cosdv] -- >> e\n');
            end
            w=1;
        end
      
        if (( (round(Pi_copula(ii-1),2)==round(M_copula(ii-1),2)) || ...
              (round(Pi_copula(ii-1),2)==round(W_copula(ii-1),2))) && (length(ss)>4))
            if(lambda_debug_print)
                fprintf('[cosdv] -- >> f\n');
            end
            w=1;
        end
        condd = (round(E_copula(ii-1),2)==round(Pi_copula(ii-1),2));
        if(condd && indep)
            if(lambda_debug_print)
                fprintf('[cosdv] -- >> g\n');
            end
            w = 0;
        end
        
        if(general_debug_print)
            fprintf('[cosdv] -- lambda_min=%0.02f lambda_max=%0.02f\n',wanc,w);
        end
        
        score = score+(length(ss)-1)*(w+wanc)/2;
        wanc = w;
        ss_minIdx = ii;
        numDomains = numDomains + 1;
    else
        if(ii==n)
            if(general_debug_print)
                fprintf('[cosdv] -- Domain -> %d:%d\n',ss_minIdx,ii);
            end
            % compute the relative distance function
            if (E_copula(ii)>Pi_copula(ii))
                if(lambda_debug_print)
                    fprintf('[cosdv] -- >> a\n');
                end
                w=(E_copula(ii)-Pi_copula(ii))/(M_copula(ii)-Pi_copula(ii));
            else
                if (E_copula(ii)<Pi_copula(ii))
                    if(lambda_debug_print)
                        fprintf('[cosdv] -- >> b\n');
                    end
                    w=(Pi_copula(ii)-E_copula(ii))/(Pi_copula(ii)-W_copula(ii));
                else
                    if(lambda_debug_print)
                        fprintf('[cosdv] -- >> c\n');
                    end
                    w=0;
                end
            end
            
            if ( ((E_copula(ii-1)==M_copula(ii-1)) || (E_copula(ii-1)==W_copula(ii-1))) && ...
                  (length(ss)>=4 ) )
                if(lambda_debug_print)
                    fprintf('[cosdv] -- >> d\n');
                end
                w=1;
            end
            if(general_debug_print)
                fprintf('[cosdv] -- lambda_min=%0.02f lambda_max=%0.02f\n',wanc,w);
            end
            score=score+(length(ss)-1)*(w+wanc)/2;
            numDomains = numDomains + 1;
        end
    end
end

metric = (score+1)/n;
Rr=corrcoef(u,v);
if Rr(1,2)<0
metric=-metric;
end

end