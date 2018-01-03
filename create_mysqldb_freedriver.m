
% CREATE MySQL DATABASE for 'PTSDPET' study
% Assume that we have the free mysql function
% (see http://www.math.toronto.edu/almgren/mysql/)
% and a mysql server running on localhost

mysql('open')
error = mysql('status');
if error,
  error('cannot connect to mysql database on localhost');
end
mysql('use ptsdpet')


for i=1:20,
    if i < 11,
        ID{i} = sprintf('c%02d',i);
        GP(i) = 0;
    else
        ID{i} = sprintf('p%02d',i-10);
        GP(i) = 1;
    end
end

%              CONTROLS...                                                 PTSD...
%                F     M     F     M     M     M     M     M     M     F     M     M     F     M     M     F     M     M     M     F
DEMO.SEX   = [   1.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0   1.0   0.0   0.0   1.0   0.0   0.0   0.0   1.0 ];
DEMO.AGE   = [  50.0  52.0  51.0  49.0  47.0  46.0  42.0  49.0  44.0  44.0  52.0  44.0  54.0  56.0  43.0  51.0  48.0  51.0  52.0  57.0 ];
DEMO.EDUC  = [  12.0  11.0  11.0  11.0  11.0  14.0  12.0  16.0   8.0  10.0  11.0  11.0  12.0  11.0  12.0  13.0  13.0   8.0  10.0  10.0 ];
DEMO.OCC   = [   4.6   5.6   5.3   4.9   4.9   3.9   4.3   4.5   5.5   5.5   4.2   4.1   3.9   2.0   4.3   4.1   4.0   4.6   5.0   5.3 ];
DEMO.IQ    = [ 105.0 112.0 109.0 114.0 108.0 118.0 111.0 110.0  96.0 111.0 105.0 110.0 111.0 110.0 104.0 109.0 107.0 113.0 111.0 106.0 ];

% create a new demographics table
sqlcom = ['create table demo ( ',...
          'id char(3) not null primary key, ',...
          'gp smallint, ',...
          'sex smallint, ',...
          'age smallint, ',...
          'educ smallint, ',...
          'occ double, ',...
          'IQ smallint);'];

str = mysql('show tables');
if strmatch('demo',str),
  mysql('drop table demo;');
end
mysql(sqlcom);

for i = 1:length(GP),
	values = [GP(i),DEMO.SEX(i),DEMO.AGE(i),DEMO.EDUC(i),DEMO.OCC(i),DEMO.IQ(i)];
	valueString = num2str(values,', %5.3f');
	valueString = ['''',ID{i},'''', valueString];
	mysql(['insert into demo values(',valueString,')']);
end



CLIN.BDI   = [  10     9    12     2     5     0     3     0     2     9    23    37    29     9    32    29    26    16    34    13   ];
CLIN.STANX = [  43    26    36    22    29    27    30    22    39    20    47    68    51    51    43    66    49    49    59    45   ];
CLIN.TRANX = [  41    29    42    36    35    26    22    29    41    32    59    69    64    27    67    61    60    50    73    58   ];
CLIN.GHQ   = [  12     6     1     0     0     0     0     2     0     1    10    30    26     7    30    26    29     5    27    16   ];

IES.AVOID  = [   0     0     0     0     0     0     0     0     0     0    13    29    26     1    27    29    26     6    16    32   ];
IES.INTR   = [   0     0     0     0     0     0     0     0     0     0    31    33    30     6    31    35    26    12    33    21   ];
IES.TOT    = [   0     0     0     0     0     0     0     0     0     0    44    62    56     7    58    64    52    18    49    53   ];

CAPS.B     = [   0     0     0     0     0     0     0     0     0     0    22    15    25     2     6    19    15    12    31     9   ];
CAPS.C     = [   0     0     0     0     0     0     0     0     0     0    32    45    27    10    33    33    37    46    34    20   ];
CAPS.D     = [   0     0     0     0     0     0     0     0     0     0    32    36    18    20    23    29    35    34    26    23   ];
CAPS.TOTF  = [   0     0     0     0     0     0     0     0     0     0    44    47    41    19    28    43    40    46    48    26   ];
CAPS.TOTI  = [   0     0     0     0     0     0     0     0     0     0    42    49    29    13    34    38    47    46    43    26   ];
CAPS.TOT   = [   0     0     0     0     0     0     0     0     0     0    86    96    70    32    62    81    87    92    91    52   ];

% create clin create table SQL command
sqlcom = ['create table clin ( ',...
          'id char(3) not null primary key, ',...
          'gp smallint, ',...
          'bdi smallint, ',...
          'stanx smallint, ',...
          'tranx smallint, ',...
          'ghq smallint, ',...
          'ies_a smallint, ',...
          'ies_i smallint, ',...
          'ies_t smallint, ',...
          'caps_b smallint, ',...
          'caps_c smallint, ',...
          'caps_d smallint, ',...
          'caps_tf smallint, ',...
          'caps_ti smallint, ',...
          'caps_t  smallint ',...
          ');'];

			
str = mysql('show tables');
if strmatch('clin',str),
  mysql('drop table clin;');
end
mysql(sqlcom);

for i = 1:length(GP),
	values = [GP(i) CLIN.BDI(i) CLIN.STANX(i) CLIN.TRANX(i) CLIN.GHQ(i) IES.AVOID(i) IES.INTR(i) IES.TOT(i) CAPS.B(i) CAPS.C(i) CAPS.D(i) CAPS.TOTF(i) CAPS.TOTI(i) CAPS.TOT(i) ];
	valueString = num2str(values,', %5.3f');
	valueString = ['''',ID{i},'''', valueString];
	mysql(['insert into clin values(',valueString,')']);
end



RT.OAC     = [ 417   NaN   566   NaN   481   NaN   371   NaN   537   NaN   NaN   NaN   NaN   NaN   NaN   492   556   NaN   476   NaN   ];
RT.OUC     = [ NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   100   NaN   NaN   NaN   NaN   NaN   413   NaN   NaN   ];
RT.OAT     = [ 454   470   572   477   451   429   403   464   536   521   644   736   564   NaN   574   468   641   539   520   620   ];
RT.OUT     = [ NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   538   NaN   NaN   NaN   NaN   NaN   NaN   522   NaN   NaN   NaN   NaN   ];
RT.OFN     = [   0     0     0     1     2     3     0     0     1     0     6     3     1   NaN     9     0     6     1     2     2   ];
RT.OFP     = [   2     0     1     0     4     0     2     0     8     0     0     1     0   NaN     0     2     2     1     1     0   ];
RT.OFPRT   = [ 417   NaN   566   NaN   481   NaN   371   NaN   537   NaN   NaN   100   NaN   NaN   NaN   507   556   413   476   NaN   ];

RT.TAC     = [ NaN   NaN   NaN   637   NaN   NaN   394   405   548   694   573   NaN   NaN   NaN   NaN   NaN   697   466   561   789   ];
RT.TUC     = [ NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   ];
RT.TAT     = [ 486   408   576   496   453   429   403   454   575   584   672   751   572   NaN   565   464   704   515   491   675   ];
RT.TUT     = [ NaN   341   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   ];
RT.TAD     = [ NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   448   NaN   NaN   NaN   NaN   NaN   739   906   NaN   796   616   ];
RT.TUD     = [ 725   NaN   NaN   646   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   846   NaN   NaN   NaN   NaN   NaN   ];
RT.TFN     = [   0     2     2     5     1     1     0     0     4     0    11     7     0   NaN    16    10    19     7     4    16   ];
RT.TFP     = [   1     1     0     6     0     0     2     1     2     3     1     0     0   NaN     1     1     5     2     3    19   ];
RT.TFPRT   = [ 725   341   NaN   639   NaN   NaN   394   405   548   612   573   NaN   NaN   NaN   846   739   801   466   640   734   ];


% create rt create table SQL command
sqlcom = ['create table behav ( ',...
          'id char(3) not null primary key, ',...
          'gp smallint, ',...
          'rt_oac smallint, ',...
          'rt_ouc smallint, ',...
          'rt_oat smallint, ',...
          'rt_out smallint, ',...
          'rt_ofn smallint, ',...
          'rt_ofp smallint, ',...
          'rt_ofprt smallint, ',...
          'rt_tac smallint, ',...
          'rt_tuc smallint, ',...
          'rt_tat smallint, ',...
          'rt_tut smallint, ',...
          'rt_tad smallint, ',...
          'rt_tud smallint, ',...
          'rt_tfn smallint, ',...
          'rt_tfp smallint, ',...
          'rt_tfprt smallint ',...
          ');'];

str = mysql('show tables');
if strmatch('behav',str),
  mysql('drop table behav;');
end
mysql(sqlcom);

for i = 1:length(GP),
	values = [GP(i) RT.OAC(i) RT.OUC(i) RT.OAT(i) RT.OUT(i) RT.OFN(i) RT.OFP(i) RT.OFPRT(i) RT.TAC(i) RT.TUC(i) RT.TAT(i) RT.TUT(i) RT.TAD(i) RT.TUD(i) RT.TFN(i) RT.TFP(i) RT.TFPRT(i) ];
	valueString = num2str(values,', %5.3f');
	valueString = ['''',ID{i},'''', valueString];
	valueString = strrep(valueString,'NaN','NULL');
	mysql(['insert into behav values(',valueString,')']);
end


mysql('close')
