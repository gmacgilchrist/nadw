function S=ariane_S(subchars,values)
% function S = ariane_S(subchars,values)
%
% Isolate a subset of Ariane trajectores based on their values.
%
% INPUT
%   subchars    cell array of unique ariane outputs
%   values      corresponding values at which to extract subset (can be
%               single values, or 2-element array to define range
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

disp('ariane_S: subsetting ariane output based on characteristic values')
S = ones(size(subchars{1}));

for c=1:length(subchars)
    switch length(values{c})
        case{1}
            S = S & subchars{c}==values{c};
        case{2}
            S = S & subchars{c}>values{c}(1) & subchars{c}<values{c}(2);
    end
end
