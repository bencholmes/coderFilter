function y = coderFilter(u, fs, cutoff, Q,  type)
% CODERFILTER a MATLAB Coder compatible second order low/hi/bandpass filter
% directly modelled from an analogue prototype of LCR.

% Author:   Ben Holmes
% Date:     2018/05/28
% License:  BSD

% Input
% where N is the number of cascaded sections...
%   - u:        Input signal to be filtered.
%   - fs:       The sampling frequency of the signal.
%   - cutoff:   Vector of length N containing cutoff frequencies of each 
%   cascaded filter.
%   - Q:        Vector of length N containing Q factor of each cascaded
%   filter.
%   - type:     Cell array of length N containing character arrays
%   indicating 'high', 'band', or 'low' pass filters.

%% Input checking
% If only a string is provided for the type wrap it into a cell array.
if ischar(type)
    type = {type};
end

if length(type) ~= length(cutoff) || length(type) ~= length(Q)
    error('Type, cutoff, and Q must all be the same length');
end

nCascades = length(type);
%% Calculate SSM

% c notates cell here, designed to facilitate Coder compatibility.
Ac = cell(1,nCascades);
Bc = cell(1,nCascades);
Cc = cell(1,nCascades);
Dc = cell(1,nCascades);

% Calculate state-space matrices
for nn=1:nCascades
    [Ac{nn}, Bc{nn}, Cc{nn}, Dc{nn}] = analogue2nd_calculateSSM(2*pi*cutoff(nn),...
                                    Q(nn), fs, type{nn});
end

% Find numeric matrices from cells.
if nCascades > 1
    [A,B,C,D] = cascadeSSM(Ac,Bc,Cc,Dc);
else
    A = Ac{1};
    B = Bc{1};
    C = Cc{1};
    D = Dc{1};
end

%% Run model
Ns = size(u,2);
x  = zeros(size(A,1),Ns+1);
y  = zeros(size(D,1),Ns);

for nn=1:Ns
    x(:,nn+1) = A*x(:,nn) + B*u(:,nn);
    y(:,nn) = C*x(:,nn) + D*u(:,nn);
end

end

function [A,B,C,D] = analogue2nd_calculateSSM(omega_c, Q, fs, type)

if strcmp(type,'band')
    Q = 1/Q;
end

R1 = 1e3;
C1 = Q/(R1*omega_c);
L1 = R1/(omega_c*Q);
T = 1/fs;

% Reactive components
if strcmp(type,'band')
    Nr = [  0  0  1];
    
    Nx = [  1 -1  0;... L
            0  1 -1];
        
    Nv = [  1  0  0];
    
    No = [  0  0  1];
else
    Nr = [  0  1];
    
    Nx = [  1 -1;... L
            0  1];
        
    Nv = [  1  0];
    
    No = [  0  1];
end

Gr = 1/R1;
GxVec = [T/(2*L1) (2*C1)/T];
ZVec = [-1 1];

if strcmp(type,'high') 
    GxVec = fliplr(GxVec);
    ZVec = fliplr(ZVec);
end

Gx = diag(GxVec);
Z = diag(ZVec);

% Acquire sizes
nStates     = size(Nx,1);   % number of states
nInputs     = size(Nv,1);   % number of voltage sources
nNodes      = size(No,2);   % number of nodes (sum)
nOutputs    = size(No,1);   % number of outputs

% Find constant system matrix
S = [Nr.'*Gr*Nr+...
    Nx.'*Gx*Nx,...
    Nv.';...
    Nv,...
    zeros(nInputs,nInputs)];

A = 2*Z*Gx...
    *[Nx zeros(nStates,nInputs)]...
    *(S\[Nx zeros(nStates,nInputs)].')...
     - Z;

B = 2*Z*Gx...
    *[Nx zeros(nStates,nInputs)]...
    *(S\[zeros(nInputs,nNodes) eye(nInputs)].');


C = [No zeros(nOutputs,nInputs)]...
    *(S\[Nx zeros(nStates,nInputs)].');

D = [No zeros(nOutputs,nInputs)]...
    *(S\[zeros(nInputs,nNodes) eye(nInputs)].');
        
end