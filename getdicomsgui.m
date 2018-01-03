%extension=char(input('What extension do you use in your image files?\nPlease type between single quotes (case sensitive)\n'))
extension=evalin('base','extension');
tofind=['.*\.' extension];

tensorvol=[];
changeto=(uigetdir);
eval(['cd ' changeto ''])
dicomdir=dir(changeto);
dicomdir=struct2cell(dicomdir);
dicomdir=dicomdir(1,:);
dicomdir=dicomdir';
%dicomdir=char(dicomdir);
counter=0;


for num=1:length(dicomdir)
    filename=char(dicomdir(num,:));
    f=regexp(filename,tofind);    
    if (f==1) 
        counter=counter+1;
        %st=char(st);
        ima=dicomread(filename);
        %eval (['slice' num2str(counter) '=ima']);
        disp(filename)
        tensorvol=cat(3,tensorvol,ima);
    end
end
disp('found')
disp (counter)
disp ('image files in this directory')
%tensorvol=tensorvol';
tensorvol=reshape(tensorvol,256,256,(counter/7),7);
assignin('base','counter',counter)
assignin('base','tensorvol',tensorvol)
%disp ('they are now in a volume called tensorvol')


