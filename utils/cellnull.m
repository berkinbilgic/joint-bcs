function [ C ] = cellnull( C )

for s = 1:length(C)  
   for w = 1:length(C{s})
       C{s}{w} = zeros(size(C{s}{w}));
   end
end

end

