

Anatomy = '';
Comments = '';
DateOfAcquisition = '';
DateOfModification = '';
Name = '';
Tesselation = '';
TriConn = '';
VertConn = '';

SUB( 1).id = 'c01'; SUB( 1).name = '';  SUB( 1).date = '26-Jan-1997';
SUB( 2).id = 'c02'; SUB( 2).name = '';  SUB( 2).date = '';
SUB( 3).id = 'c03'; SUB( 3).name = '';  SUB( 3).date = '';
SUB( 4).id = 'c04'; SUB( 4).name = '';  SUB( 4).date = '';
SUB( 5).id = 'c05'; SUB( 5).name = '';  SUB( 5).date = '';
SUB( 6).id = 'c06'; SUB( 6).name = '';  SUB( 6).date = '';
SUB( 7).id = 'c07'; SUB( 7).name = '';  SUB( 7).date = '';
SUB( 8).id = 'c08'; SUB( 8).name = '';  SUB( 8).date = '';
SUB( 9).id = 'c09'; SUB( 9).name = '';  SUB( 9).date = '';
SUB(10).id = 'c10'; SUB(10).name = '';  SUB(10).date = '';

subjpath  = 'd:\matlab\brainstorm_v1\subjects';
studypath = 'd:\matlab\brainstorm_v1\studies';

for s = 1:10,
    
    fprintf('...processing %s\n',SUB(s).id);
    
    subjfile  = sprintf('%s\\%s%02d\\%s%02d_brainstormsubject.mat',subjpath, char(g),s,char(g),s);
    studyfile = sprintf('%s\\%s%02d\\%s%02d_brainstormstudy.mat',  studypath,char(g),s,char(g),s);
    
    %load(subjfile)
    
    Comments = 'PTSDPET PROJECT';
    
    DateOfAcquisition = SUB(s).date;
    
    DateOfModification = date;
    
    Name = SUB(s).name;
    
    Anatomy = sprintf('%s%02d\\%s%02d_subjectimage.mat',char(g),s,char(g),s);
    
    Tesselation = sprintf('%s%02d\\%s%02d_subjecttess.mat',char(g),s,char(g),s);
    
    TriConn = '';
    VertConn = '';
    
    save(subjfile,'Comments','DateOfAcquisition','DateOfModification',...
        'Name','Anatomy','Tesselation','TriConn','VertConn');
    
    
    % Setup and save the brainstormstudy.mat file in studies directory
    
    BrainStormSubject = sprintf('%s%02d\\%s%02d_brainstormsubject.mat',char(g),s,char(g),s);
    
    DateOfStudy = DateOfAcquisition;
    
    Session = 'ERP';
    
    save(studyfile,'BrainStormSubject','DateOfModification','DateOfStudy',...
        'Name','Session');
    
end
