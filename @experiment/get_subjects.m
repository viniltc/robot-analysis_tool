function [subjs, subjnames] = get_subjects(exp)
% Gets SUBJECT objects (and their names) associated with a given experiment
% (C) V. Sanguineti (2008)

% keyboard
subjs = exp.subj;
subjnames = exp.subjnames;