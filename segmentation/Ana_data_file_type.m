function [the_other, bit_per_pixel,something_else]=Ana_data_file_type(datatype,filetype)
% Given the datatype ot the fyletype gives the other, that is:
% Given the filetype gives the datatype (in Save_Ana_file.m ) or
% Given the datatype gives the (real) filetype (in Read_Ana_File).
% One of the inputs is emtpty=[].
% Matlab does not contains all the types of C !!!!.

% rgpm 4-04-2003 Felicidades Rodolfo y Rolando.

typedata=[   0,        1,       2,     8,      16,    32,        4];          %according to analyze users
typefile={'Unknown','Binary','uint8','int32','float', 'double', 'short'};     %according to analyze users
realfile={'Unknown','uint8', 'uint8','int32','single','double', 'int16' };    %according to Matlab
bitxpixel=[  0,      8,        8,      32,     32,     64,       16];         %according to both

if isempty(datatype), % DataType to save
    kk=strmatch('double',typefile);
    the_other=typedata(kk),
    bit_per_pixel=bitxpixel(kk),
    something_else=realfile(kk),
    
elseif isempty(filetype), % FileType to read
    kk=find(typedata==datatype);
    the_other=char(realfile(kk)),
    bit_per_pixel=bitxpixel(kk),
    something_else=typefile(kk),
    
else
    myerror(['Option not available /',mfilename])
end

