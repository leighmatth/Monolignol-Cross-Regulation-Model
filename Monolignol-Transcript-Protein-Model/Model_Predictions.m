% Copyright (c) 2019, North Carolina State University. All rights reserved.
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation under version 2 of the License.
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

function Ypred=Model_Predictions(Xact,Xtarg,Xfulltrans)
if nargin<3
    prevhalf_flag=0;
else
    prevhalf_flag=1;
    if size(Xfulltrans)~=[size(Xtarg,1)/2 size(Xtarg,2)]
        error('Wrong dimensions in Xfulltrans (full transcript profile for previous model).')
    end
end

if size(Xact)~=size(Xtarg)
    error('Xact and Xtarg must have same dimensions.')
end

load ../ModelDetails.mat
load ../TranscriptProteinAbundances.mat Ywt_avg

Ywt_avg=Ywt_avg{:,:}';

[M,Nact]=size(Xtarg);

% For Old Model:
% 'full' uses experimental transcript values for just the targeted and avg.
% WT for the un-targeted transcripts
muprevfull=[Ywt_avg(1:M/2,1); zeros(M/2,1)];
Xadj_prevfull=Xtarg;

if prevhalf_flag
    % 'half' uses experimental transcript values for all transcripts, not just targeted
    muprevhalf=zeros(M,1);
    Xadj_prevhalf=[Xfulltrans;zeros(M/2,Nact)];
end

% For New Model:
% Set proteins of targeted genes to percentage of their wild-type abundance
Ytranswtavg=Ywt_avg(1:M/2,1);
Yprotwtavg=Ywt_avg(1+M/2:M,1);
Yprotnewtarg=Yprotwtavg.*(Xtarg(1:M/2,:)./Ytranswtavg);

Xadj_full=[Xtarg(1:M/2,:); Yprotnewtarg]; % adjust targeted values to include the proteins

% Allocate space for variables
Ypred.prevhalf=zeros(M,Nact);
Ypred.prevfull=zeros(M,Nact);
Ypred.full=zeros(M,Nact);

for j=1:Nact
    
    Kappa=diag(~([Xact(1:M/2,j); Xact(1:M/2,j)])); %for experiment j, diagonal matrix with 0 for targeted transcripts/proteins, 1 for others
    Kappa_prev=diag([~Xact(1:M/2,j);ones(M/2,1)]);
    
    % predictions from previous model
    Ypred.prevhalf(:,j)=(eye(M)-Kappa_prev*Boldmodel)\(Kappa_prev*muprevhalf+Xadj_prevhalf(:,j));
    Ypred.prevfull(:,j)=(eye(M)-Kappa_prev*Boldmodel)\(Kappa_prev*muprevfull+Xadj_prevfull(:,j));
    
    % predictions from our proposed model
    Ypred.full(:,j)=(eye(M)-Kappa*B)\(Kappa*mue+Xadj_full(:,j));
end
