% SHAHROOZ IRANMANESH
% STUDENT NUMBER =97150504
% FINAL PROJECT FINITE ELEMENT
% ----------------------------

%----------------------------
%INPUT DATA
%----------------------------

nel=8;                        %number of element
nnel=3;                       %number of nodes per element
ndof=2;                       %number of dofs per node
nnode=10;                     % total number of nodes
sdof=nnode*ndof;              %tatal dofs
edof=nnel*ndof;               %degree of freedom per element
emodule=1000000;              %modulus
poisson=0.3;                  %ratio

%----------------------------
%DATA FOR COORDINATE OF NODAL
%----------------------------

gcoord=[0 0;0 1;1 0;1 1;2 0;
        2 1;3 0;3 1;4 0;4 1];
    
%---------------------------
% DATA FOR NODAL CONNECTIVITY FOR EACH ELEMENT
%---------------------------

nodes=[1 3 4;1 4 2;3 5 6;3 6 4;5 7 8;
       5 8 6;7 9 10;7 10 8;];
   
 %-------------------------
 %DATA FOR BOUNDRY CONDITION
 %-------------------------
 
 bcdof=[1 2 3];              %vector containing dofs with boundry condition
 bcval=[0 0 0];              %vector containing boundry condition with dofs in bcdof
 
 %-----------------------
 %INITIALIZATION
 %-----------------------
 
 ff=zeros(sdof,1);            %force vector
 kk=zeros(sdof,sdof);         %sys matrix
 disp=zeros(sdof,1);          %displacement vector
 eldisp=zeros(edof,1);        %element displacement 
 stress=zeros(nel,3);         %matrix contain stress
 strain=zeros(nel,3);         %matrix contain strain
 index=zeros(edof,1);         %index vector
 kinmtx=zeros(3,edof);        %kinematic matrix
 matmtx=zeros(3,3);           %constitutive matrix
 
 %------------------------
 %FORCE VECTOR
 %------------------------
 
 ff(17)=500;                  %force applied at node 9
 ff(19)=500;                  %force applied at node 10
 
 %-----------------------
 %COMPUTE ELEMENT AND ASEMBLE
 %-----------------------
 
 matmtx=fematiso(1,emodule,poisson);    %constitutive matrix
 
 for iel=1:nel                          %for total number of element
     
 nd(1)=nodes(iel,1);                    %1st conected node
 nd(2)=nodes(iel,2);                    %2nd conected node
 nd(3)=nodes(iel,3);                    %3rd conected node
 
 x1=gcoord(nd(1),1);y1=gcoord(nd(1),2); %coordinate 1st node
 x2=gcoord(nd(2),1);y2=gcoord(nd(2),2); %coordinate 2nd node
 x3=gcoord(nd(3),1);y3=gcoord(nd(3),2); %coordinate 3rd node
 
 index=feeldof(nd,nnel,ndof);           %extract dofs for element
 
 %-------------------------
 %DERIVATIVE OF SHAPE
 %-------------------------
 
 area=0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2);  %area traiangle
 area2=area*2;
 dhdx=(1/area2)*[(y2-y3)*(y3-y1)*(y1-y2)];          %derivative x
 dhdy=(1/area2)*[(x3-x2)*(x1-x3)*(x2-x1)];          %derivative y
 
 kinmtx2=fekine2d(nnel,dhdx,dhdy);                %kinematic matrix
 k=kinmtx2*matmtx*kinmtx2*area;                   %stifness matrix
 kk=feasmbl1(kk,k,index);                         %assemble element
 
 end                                              %end loop
 
 %---------------------
 %APPLY BOUNDRY CONDITION
 %--------------------
 
 [kk,ff]=feaolyc2(kk,ff,bcdof,bcval);
 
 %----------------------
 %SOLVE MATRIX
 %----------------------
 
 disp=kk\ff;
 
 %----------------------
 %ELEMENT STRESS COMPUTE
 %----------------------
 
 for ielp=1:nel                      %loop for total number of element
 
 nd(1)=nodes(ielp,1);                    %1st conected node
 nd(2)=nodes(ielp,2);                    %2nd conected node
 nd(3)=nodes(ielp,3);                    %3rd conected node   
 
 x1=gcoord(nd(1),1);y1=gcoord(nd(1),2);  %coordinate 1st node
 x2=gcoord(nd(2),1);y2=gcoord(nd(2),2);  %coordinate 2nd node
 x3=gcoord(nd(3),1);y3=gcoord(nd(3),2);  %coordinate 3rd node
 
index=feeldof(nd,nnel,ndof);             %extract sys dofs for element

%-----------------------
%EXRACT ELEMENT DISPLACEMENT
%-----------------------

for i=1:edof

eldisp(i)=disp(index(i));
end

erea=0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2);   %area triangle
area2=area*2;
 dhdx=(1/area2)*((y2-y3)(y3-y1)(y1-y2));          %derivative x
 dhdy=(1/area2)*((x3-x2)(x1-x3)(x2-x1));          %derivative y
 
 kinmtx2=fekine2d(nnel,dhdx,dhdy);                %kinematic matrix
 estrain=kinmtx2*eldisp;                          %compute strain
 estress=matmtx*estrain;                          %compute stress
 
 for i=1:3
 strain(ielp,i)=estrain(i);                       %store for each element
 stress(ielp,i)=estress(i);                       %store for each element
 
 end
 
 end
 
 %------------------
 %PRINT FEM SOLUTIONS
 %------------------
 
 num=1:1:sdof;
 displace=[num disp]                             %print total displacement
 
 for i=1:nel
 stresses=[i stresses(i,:)]                      %print stresses
 
 end
 
 
