function [Ac,Bc,Cc,Dc] = cascadeSSM(A,B,C,D)
% CASCADESSM    A function to cascade linear state-space matrices for the
% purpose of cascading filters etc.
%
% Author:  Ben Holmes
% Date:    2018/05/20
% License: GPL V3
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