function [eofs,baseval,depth,weights,xi,f,num_eofs,num_depth] = eb_read_eeof(filename,obj_bool)
%eb_read_eeof loads all arguments of EEOF files as doubles
% 
% input 1: filename
% input 2 (optional): boolean to load separately (default) or under one structure

if nargin == 1
    obj_bool = false;
end

% loads from file as doubles
eofs        = double(ncread(filename,'eofs')); 
baseval     = double(ncread(filename,'baseval'));
num_eofs    = double(ncread(filename,'num_eofs'));
num_depth   = double(ncread(filename,'num_depth'));
depth       = double(ncread(filename,'depth'));
weights     = double(ncread(filename,'weights'));
xi          = double(ncread(filename,'pdf_val'));
f           = double(ncread(filename,'pdf_freq'));

if obj_bool
    % assign to object
    eeof_obj.eofs       = eofs;
    eeof_obj.baseval    = baseval;
    eeof_obj.num_eofs   = num_eofs;
    eeof_obj.num_depth  = num_depth;
    eeof_obj.depth      = depth;
    eeof_obj.weights    = weights;
    eeof_obj.xi         = xi;
    eeof_obj.f          = f;
    eofs = eeof_obj;
end
    
end

