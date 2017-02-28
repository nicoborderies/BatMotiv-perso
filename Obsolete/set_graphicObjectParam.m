 function [f] = set_graphicObjectParam(f,param)
 
    % default parameters
        if nominal(class(f))=='matlab.ui.Figure'
            default.Position = [1 0 0.9 0.9];
            default.Units = 'normalized';
        elseif  nominal(class(f))=='matlab.graphics.axis.Axes'
            default.fontSize = 20;
        elseif nominal(class(f))=='matlab.graphics.chart.primitive.Line'
            default.Colors = [0 0 0];
            default.linesize = 2;
            default.markerSize = 10;
        end
        
        if nargin<2
            param = default;
        else
            arg = setdiff(fieldnames(default),fieldnames(param));
            nA = numel(arg);
            for iA = 1:nA
                param.(arg{iA}) =  default.(arg{iA});
            end
        end
 
    % set object properties
        arg = fieldnames(param);
        np = numel(arg);
        for ip = 1:np
            try
            f.(arg{ip}) =  param.(arg{ip});
            end
        end
 
 end