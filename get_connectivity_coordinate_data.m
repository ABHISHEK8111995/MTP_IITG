% % Number of Elements
nel = nx*ny ;
% 
% % Number of Nodes
nno = (nx+1)*(ny+1) ;
% 
% % Number of dofs per element 
ndoel = 8 ;
% 
% % Number of dofs in the FE mesh, Number of dofs per element
 ndof = 2*nno ; ndoelo = 8 ;
% 
% filename = 'Input/fe_data.txt' ;
% fid = fopen(filename,'w') ;
% fprintf(fid,'%g \t %g \t %g \t %g \t %g',ny,nx,nel,ndoel,ndof);
% fclose(fid);

% Element Connectivity
%
% Note that the below code works only for rectangular bodies. You can
% either write your own code or take input from comemrcial FE packages like
% ANSYS. In that case change the below code accordingly.

CON = zeros(nel,4) ; % since 4 nodes per element
iel = 0 ;
for j = 1:ny
    for i = 1:nx 
        iel = iel + 1 ; 
        CON(iel,:) = [ (j-1)*(nx+1)+i  (j-1)*(nx+1)+i+1  j*(nx+1)+i+1  j*(nx+1)+i ] ;    
    end
end

  
% filename = 'Input/connectivity.txt' ;
% fid = fopen(filename,'w') ;
% for iel = 1:nel
%     if iel < nel
%         fprintf(fid,'%g \t %g \t %g \t %g \n',CON(iel,:));   
%     else
%         fprintf(fid,'%g \t %g \t %g \t %g',CON(iel,:));
%     end
% end
% fclose(fid);


% Initial Node Coordinates
 Xn = zeros(nno,2) ; % 2 dof per node i.e. x coordinate and y coordinate.
% ino = 0 ;
% for j = 0:ny
%     for i = 0:nx 
%         ino = ino + 1 ;
%         Xn(ino,:) =  [ Xo + i*Lx/nx  Yo + j*Ly/ny ] ;
%     end
% end
% 

% condition for taking different data from different txt files containing
% coordinate data
if iteration == 1
       Xn = load('coordnew_1.txt');
else if iteration == 2
       Xn = load('coordnew_2.txt');
else if iteration == 3
       Xn = load('coordnew_3.txt');
end 
end
end


% filename = 'Input/coordinate.txt' ;
% fid = fopen(filename,'w') ;
% for ino = 1:nno
%     if ino < nno
%         fprintf(fid,'%20.15f \t %20.15f \n',Xn(ino,:));   
%     else
%         fprintf(fid,'%20.15f \t %20.15f',Xn(ino,:));   
%     end
% end
% fclose(fid);


% Current node coordinates initialized as initial coordinates to start
% with.
xn = Xn ; % Xn contains the coordinates of the reference configuration all the time
          % xn contains the coordinates of the current configuration all
          % the time        