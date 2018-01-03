function order = genSliceOrder( n, type)
%
% Matthias Moosmann -Moosoft-
%
switch type
case 1, %verschachtelt
   %n gerade  
   if n/2==round(n/2),     
      for i=1:n,              
         if i/2==round(i/2),
            order(i)=(n-i+2)/2 ;            
         else 
            order(i)= (2*n-i+1)/2;
         end                 
      end
   %n ungerade      
   else
      for i=1:n,
         if i/2==round(i/2), 
            order(i)=(2*n-i+2)/2;            
         else 
            order(i)=(n-i+2)/2;           
         end
      end      
   end
   order = fliplr( order);
      
case 0, % ascending( MRT: first slice top)
   order = 1 : n;
end

