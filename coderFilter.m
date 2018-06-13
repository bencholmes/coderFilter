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

function [Ac,Bc,Cc,Dc] = cascadeSSM(A,B,C,D)
% CASCADESSM    A function to cascade linear state-space matrices for the
% purpose of cascading filters etc.
%
% This function is derived from the answer
% https://math.stackexchange.com/questions/2201067/cascade-of-state-space-models-for-linear-systems
% From user
% https://math.stackexchange.com/users/506682/arash
%
% Inputs:
%  - A,B,C,D are cell matrices where each cell is a matrix to cascade
%
% Outputs:
%  - Ac, Bc, Cc, Dc are the resultant cascaded matrices of default numeric
%  type.
%
% Notes:
%  The function is designed to be useable with MATLAB Coder, so has
%  coder.varsize defining the maximum size of the cascaded matrices. This
%  is set by default at [100 100] but can be changed by you if you like.

% Make sure that each cell matrix is of the same length, i.e. there are the
% same number of matrices to cascade of A, B, C, and D.
if length(A) ~= length(B) || length(A) ~= length(C) || length(A) ~= length(D)
    error('There must be the same number of each matrix');
end

nCasc = length(A);

% If there are more than one elements to cascade,
if nCasc > 1
    coder.varsize('Ac',[100 100]);
    coder.varsize('Bc',[100 100]);
    coder.varsize('Cc',[100 100]);
    coder.varsize('Dc',[100 100]);

    % Cascade the first two,
    [Ac, Bc, Cc, Dc] = cascadeOnce({A{1} A{2}}, {B{1} B{2}},...
                                   {C{1} C{2}}, {D{1} D{2}});

    % Then loop through the remaining matrices.
    for nn=3:nCasc
        % Append the result with the next matrix to cascade.
        [Ac, Bc, Cc, Dc] = cascadeOnce({Ac, A{nn}}, {Bc, B{nn}}, {Cc, C{nn}}, {Dc, D{nn}});
    end
else
    % Don't cascade, just output.
    % Useful for when you want to sometimes cascade and sometimes not.
    Ac = A{1};
    Bc = B{1};
    Cc = C{1};
    Dc = D{1};
end

end

function [Ac, Bc, Cc, Dc] = cascadeOnce(A, B, C, D)

Ac = [A{1} zeros(size(A{1},1), size(A{2},2)); B{2}*C{1} A{2}];
Bc = [B{1}; B{2}*D{1}];
Cc = [D{2}*C{1} C{2}];
Dc = D{2}*D{1};

end