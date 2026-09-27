function [ C ] = cellmult( C,T )

for s = 1:length(C)  
   for w = 1:length(C{s})
       C{s}{w} = C{s}{w}.*T{s}{w};
   end
end

end

