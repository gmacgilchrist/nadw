function Dop=ariane_Dop(opchar,D,operation)
% function Dop = ariane_Dop(opchars,D,operation)
%
% Perform an operation (sum or average) on the characteristics in each bin
% of a distribution, defined by the distribution structure array D.
%
% INPUT
%   opchar      ariane output of characteristic on which to operate
%   D           structure array of distribution


%
% OUTPUT
%   S           contitional array of trajectories to be included
%
% EXAMPLE
% S = ariane_subchars({final_section,final_age},{7,[0 31536000]})
% selects trajectories that 'exit' across section number 7, and that are
% less than a year old. S can then be used to isolate other characteristics
% of these trajectories e.g. init_lon(S) gives the initial position
% particles in this subset.
%

% G.A. MacGilchrist (27/03/18) gmacgilchrist@gmail.com
n = length(fieldnames(D));

switch operation
    case{'sum'}
        disp('ariane_Dop: Summing values in each bin')
        if n==4;
            Dop=ariane_Dsum(opchar,D);
        elseif n==6;
            Dop=ariane_D2sum(opchar,D);
        end
    case{'mean'}
        disp('ariane_Dop: Averaging (mean) characteristics in each bin')
        if n==4;
            Dop=ariane_Dmean(opchar,D);
        elseif n==6;
            Dop=ariane_D2mean(opchar,D);
        end
    case{'median'}
        disp('ariane_Dop: Averaging (median) characteristics in each bin')
        if n==4;
            Dop=ariane_Dmedian(opchar,D);
        elseif n==6;
            Dop=ariane_D2median(opchar,D);
        end
end
    