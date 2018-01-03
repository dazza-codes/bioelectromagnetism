function [T] = coregister_fiducials(P,Q)

% coregister_fiducials - Determine fiducial points coregistration
% 
% [T] = coregister_fiducials(P,Q)
% 
% P and Q are both a (3x3) matrix, three points in Cartesian XYZ.
% So P and Q can be column vectors for each fiducial location,
% where P is defined as:
%
% P(:,1) is nasion XYZ (3,1)
% P(:,2) is left preauricular XYZ (3,1)
% P(:,3) is right preauricular XYZ (3,1)
%
% Similarly for the other set of fiducials Q. The order of these
% points (nasion, left, right) is essential!
% 
% T is a struct, with fields T.P2Q and T.Q2P, each of which contains the
% rotations and translations in a 4x4 matrix
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% History:  02/2004, Darren.Weber_at_radiology.ucsf.edu
%                    adapted from notes by Tom Ferree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


error('under construction');


DB = 1; % debug, 0 off, 1 on

if DB,
  % These 2 point sets can be coregistered with a rotation of
  % 90 degrees anticlockwise around Z.  Each point set contains 3
  % points along orthogonal XYZ axes, they already have the same
  % origin (0,0,0).
  
  
  P = [ 1 0 0;  0 1 0; 0 -1 0]'; % nasion +X, left +Y, right -Y
  Q = [ 0 1 0; -1 0 0; 1  0 0]'; % nasion +Y, left -X, right +X
  
  % add a small translation in X, so P has a different origin from Q 
  P(1,:) = P(1,:) + 0.5;
  
else
  if ~exist('P','var') | ~exist('Q','var'),
    msg = sprintf('coregister_fiducials: no input P or Q\n');
    error(msg);
  elseif isempty(P) | isempty(Q),
    msg = sprintf('coregister_fiducials: empty input P or Q\n');
    error(msg);
  elseif ~(isequal(size(P),size(ones(3,3)))),
    msg = sprintf('coregister_fiducials: P must be 3x3 matrix\n');
    error(msg);
  elseif ~(isequal(size(Q),size(ones(3,3)))),
    msg = sprintf('coregister_fiducials: Q must be 3x3 matrix\n');
    error(msg);
  end
end



%-------------------------------------------------------------
% Plot the inputs

fig = findobj('tag','fiducials start');
if fig,
  % replace it!
  close(fig)
end

% when plotting, the X values are all in row1, y in row2, z in row3
figure('name','fiducials start','tag','fiducials start')
Color = eye(3); % red, green, blue
scatter3(P(1,:),P(2,:),P(3,:),80,Color,'filled');
text(P(1,1),P(2,1),P(3,1),'P nasion');
text(P(1,2),P(2,2),P(3,2),'P left');
text(P(1,3),P(2,3),P(3,3),'P right');
view(3); hold on
Color = [ 0 1 1; 1 0 1; 1 1 0 ]; % cyan, magenta, yellow
scatter3(Q(1,:),Q(2,:),Q(3,:),80,Color,'filled');
text(Q(1,1),Q(2,1),Q(3,1),'Q nasion');
text(Q(1,2),Q(2,2),Q(3,2),'Q left');
text(Q(1,3),Q(2,3),Q(3,3),'Q right');
legend('P Nasion','P Left','P Right','Q Nasion','Q Left','Q Right',-1)


%---------------------------------------------------------
% We assume P & Q are column vectors for each fiducial location,
% nasion, left, right, in that order, so
% P(:,1) = nasion
% P(:,2) = left
% P(:,3) = right

PN = P(:,1);
PL = P(:,2);
PR = P(:,3);

fprintf('...P Nasion = [ %8.4f %8.4f %8.4f ]\n',PN);
fprintf('...P Left   = [ %8.4f %8.4f %8.4f ]\n',PL);
fprintf('...P Right  = [ %8.4f %8.4f %8.4f ]\n',PR);

QN = Q(:,1);
QL = Q(:,2);
QR = Q(:,3);

fprintf('...Q Nasion = [ %8.4f %8.4f %8.4f ]\n',QN);
fprintf('...Q Left   = [ %8.4f %8.4f %8.4f ]\n',QL);
fprintf('...Q Right  = [ %8.4f %8.4f %8.4f ]\n',QR);


%-------------------------------------------------------------
% First determine the line passing through the nasion, which
% intersects the line between the left and right at a right angle.
% The point of intersection will be defined as the origin, Po & Qo

PL2PR = PR - PL;
PL2PN = PN - PL;
%PR2PL = -1 * PL2PR;
%PN2PL = -1 * PL2PN;

QL2QR = QR - QL;
QL2QN = QN - QL;
%QR2QL = -1 * QL2QR;
%QN2QL = -1 * QL2QN;

strang = 1;
if strang,
  
  Po = PL + ( ( dot(PL2PR,PL2PN) / ( norm(PL2PR) ^2 ) ) * PL2PR );
  Qo = QL + ( ( dot(QL2QR,QL2QN) / ( norm(QL2QR) ^2 ) ) * QL2QR );
  
else
  
  % Imagine a right angle triange formed by the left ear, nasion and some
  % unknown point along the line from the left ear to the right ear (PO).
  
  % Find the angle at the left ear between the line from the left to right
  % ear and the line from the left ear to the nasion
  Pcos_theta = dot(PL2PR,PL2PN) / ( norm(PL2PR) * norm(PL2PN) );
  Qcos_theta = dot(QL2QR,QL2QN) / ( norm(QL2QR) * norm(QL2QN) );
  
  % The length of the line from the left ear to the nasion is the hypotenuse.
  % Since cos_theta = adjacent / hypotenuse, we now find the length of the
  % adjacent line, which lies along the vector from left to right ear
  Padjacent_length = norm( PL2PN ) * Pcos_theta;
  Qadjacent_length = norm( QL2QN ) * Qcos_theta;
  
  % We now calculate the unit vector from left to right ear
  PL2PR_unit_vector = PL2PR / norm( PL2PR );
  QL2QR_unit_vector = QL2QR / norm( QL2QR );
  
  % So, multiplication of this unit vector by the desired length will give
  % the required origin, after we add this vector to the left ear location
  Po = PL + (PL2PR_unit_vector * Padjacent_length);
  Qo = QL + (QL2QR_unit_vector * Qadjacent_length);
  
end


% Add these offsets to the output T struct
T.Po = Po;
T.Qo = Qo;

if DB,
  scatter3(Po(1),Po(2),Po(3),60,'k','filled');
  scatter3(Qo(1),Qo(2),Qo(3),60,'k','filled');
end

%-------------------------------------------------------------
% define the +X vector in the direction of the nasion
% define the +Y vector in the direction of the left ear
% define the +Z vector as their cross product, using the right hand rule
PX = PN - Po;
PY = PL - Po;
PZ = cross( PX, PY );

QX = QN - Qo;
QY = QL - Qo;
QZ = cross( QX, QY );

% Now define the unit vectors
PXunit = PX / norm(PX);
PYunit = PY / norm(PY);
PZunit = cross ( PXunit , PYunit );

QXunit = QX / norm(QX);
QYunit = QY / norm(QY);
QZunit = cross ( QXunit , QYunit );

angle = acos( dot(PZunit , QZunit) ) * (180 / pi)

fprintf('...P Nasion unit vector = [ %8.4f %8.4f %8.4f ]\n',PXunit);
fprintf('...P Left unit vector   = [ %8.4f %8.4f %8.4f ]\n',PYunit);
fprintf('...P Right unit vector  = [ %8.4f %8.4f %8.4f ]\n',PZunit);

fprintf('...Q Nasion unit vector = [ %8.4f %8.4f %8.4f ]\n',QXunit);
fprintf('...Q Left unit vector   = [ %8.4f %8.4f %8.4f ]\n',QYunit);
fprintf('...Q Right unit vector  = [ %8.4f %8.4f %8.4f ]\n',QZunit);

% Given the calculations above, these errors should never occur
if dot(PXunit,PYunit), error('PX unit vector is not orthogonal to PY unit vector'); end
if dot(PXunit,PZunit), error('PX unit vector is not orthogonal to PZ unit vector'); end
if dot(PYunit,PZunit), error('PY unit vector is not orthogonal to PZ unit vector'); end
if dot(QXunit,QYunit), error('QX unit vector is not orthogonal to QY unit vector'); end
if dot(QXunit,QZunit), error('QX unit vector is not orthogonal to QZ unit vector'); end
if dot(QYunit,QZunit), error('QY unit vector is not orthogonal to QZ unit vector'); end

%-------------------------------------------------------------
% Now obtain the translation vectors



% Scaling, compare vector magnitude of Nasion in P & Q


%-------------------------------------------------------------
% Now obtain the rotation matrices

% The rotation matrices are simply column vectors of the unit vectors
Pmatrix = [ PXunit PYunit PZunit ];
Qmatrix = [ QXunit QYunit QZunit ];

T.rot.P2Q = Qmatrix;
T.rot.Q2P = Qmatrix';




%-------------------------------------------------------------
% Plot the final results

fig = findobj('tag','fiducials end');
if fig,
  % replace it!
  close(fig)
end

P2Q = (T.rot.P2Q * (P - repmat(T.Po,1,size(P,2)) ) ) + repmat(T.Qo,1,size(P,2))
Q2P = (T.rot.Q2P * (Q - repmat(T.Qo,1,size(P,2)) ) ) + repmat(T.Po,1,size(P,2))

% when plotting, the X values are all in row1, y in row2, z in row3
figure('name','fiducials end','tag','fiducials end');
Color = eye(3); % red, green, blue
scatter3(P2Q(1,:),P2Q(2,:),P2Q(3,:),80,Color,'filled');
text(P2Q(1,1),P2Q(2,1),P2Q(3,1),'P2Q nasion');
text(P2Q(1,2),P2Q(2,2),P2Q(3,2),'P2Q left');
text(P2Q(1,3),P2Q(2,3),P2Q(3,3),'P2Q right');
view(3); hold on
Color = [ 0 1 1; 1 0 1; 1 1 0 ]; % cyan, magenta, yellow
scatter3(Q2P(1,:),Q2P(2,:),Q2P(3,:),80,Color,'filled');
text(Q2P(1,1),Q2P(2,1),Q2P(3,1),'Q2P nasion');
text(Q2P(1,2),Q2P(2,2),Q2P(3,2),'Q2P left');
text(Q2P(1,3),Q2P(2,3),Q2P(3,3),'Q2P right');
legend('P2Q Nasion','P2Q Left','P2Q Right','Q2P Nasion','Q2P Left','Q2P Right',-1)



return
