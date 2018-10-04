%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    2003 Program for Women in Mathematics,  Mathematical Biology     %%
%% Institute for Advanced Study, Princeton University, May 12-22, 2003 %%
%%                  L.J. Fauci, Tulane University,                     %%
%%        K.A. Rejniak - Mathematical Biosciences Institute, OSU       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CellGrowth


  %--------------------------------------------------------------------%
  % model parameters - fixed                                           %
  %--------------------------------------------------------------------%
  example =1 
  Ng=64;                               % fluid grid size
  Ngg=32;                              % grid size for pictures only
  xmin=-1; xmax=1; cen=(xmax+xmin)/2;  % fluid domain (square)
  hg=(xmax-xmin)/Ng;                   % fluid mesh width

  Src=1;                               % source strength
  Spr=100;                             % spring stiffness
  rho=1;                               % fluid density
  mu=1;                                % fluid viscosity

  dt=0.05;                             % time step

  Nb=64;                               % number of boundary points
  len=0.2;                             % body radius/half length
  shape=1;                             % body shape: 1-circle, 2-square
  connect=0;                           % spring connections:
                                       % 0-adjacent only, 1-adjacent and secondary
                                       % 2-adjacent and center, 3-adjcent and opposite

  scaleF=1.0;                          % scalling parameter for forces
  scaleV=1.0;                          % scalling parameter for velocities

  NumLoop=40;                          % number of steps
  mod_num=5;                           % frequency
  %Mov=moviein(NumLoop);

  %-- define parameters for cell shape and cytoskeleton --%
  if (example==1)
    mu=1; Src=5; scaleF=0.01; scaleV=1;
    shape=1; connect=3;  % 1-circle 3-adjacent and opposite forces
  else
    mu=0.01; Src=5; scaleF=0.01; scaleV=0.5;
    shape=2; connect=2;  % 2-square 2-adjacent and center forces
  end

  %-- define fluid grid --%
  for ii=1:Ng+1
   for jj=1:Ng+1
      xg(ii,jj,1)=xmin+(ii-1)*hg;
      xg(ii,jj,2)=xmin+(jj-1)*hg;
   end
  end



  %----------------------%
  %-- development loop --%
  for repeat_num=1:1   % 5 runs: cells shape; forces; velocities

    %-- cell shape --%
    xb=DefineCellShape(shape,Nb,len);
    hb=sqrt((xb(1,1)-xb(1,2))^2+(xb(2,1)-xb(2,2))^2);

    %-- distributed sources and sinks --%
    sb(1,1)=cen;  sb(2,1)=cen;  sbb(1,1)= Src;     % a source at the center
    sb(1,2)=xmin; sb(2,2)=xmin; sbb(1,2)=-Src;     % a sink in the corner
    Nbs=2;

    %-- initial velocity --%
    ug=zeros(Ng+1,Ng+1,2);
    ub=zeros(2,Nb);


    for loop_num=1:NumLoop

      %-- new position of boundary points --%
      xb=xb+ub*dt;

      %-- boundary forces --%
      fb=zeros(size(xb));
      fadj=AdjacentForces(xb,Nb,hb,Spr);            % adjacent
      fsec=SecondaryForces(xb,Nb,hb,Spr,connect);   % secondary
      fcen=CenterForces(xb,Nb,cen,len,Spr,connect); % center      
      fopp=OppositeForces(xb,Nb,len,Spr,connect);   % opposite
      fb=fadj+fsec+fcen+fopp;    % add all forces

      %-- grid sources --%
      sg=BoundToGrid1(sb,sbb,Nbs,Ng,hg,hg,0.5*hg,xmin,xmax);
      
      %-- grid forces --%
      fg=zeros(size(ug));
      fg=BoundToGrid2(xb,fb,Nb,Ng,hg,hg,0.5*hg,xmin,xmax);
      
      %-- compute grid velocity from NavierStokes --%
      vg=NavierStokes(ug,fg,sg,Ng,rho,mu,dt,hg);
      ug=vg;

      %-- boundary velocities --%
      ub=GridToBound(xb,Nb,vg,Ng,hg,hg,xmin,xmax);


      %-- draw the results --%
      %clf
      %axis([xmin,xmax,xmin,xmax])
      %axis equal
      %hold on

      %-- draw cell shape --%
      %plot(xb(1,:),xb(2,:))                        % cell membrane
      %plot([xb(1,Nb),xb(1,1)],[xb(2,Nb),xb(2,1)])
      %plot(xb(1,:),xb(2,:),'*r')                   % boundary points


      %if (repeat_num==1)
      %  DrawCellCytoskeleton(xb,Nb,cen,connect)
      %  title('(1) Growing Cell: cell boundary-red stars, cell cytoskeleton-grey lines','FontSize',15)
      %elseif (repeat_num==2)
      %  scaleFG=scaleF*2;
      %  quiver(xb(1,:),xb(2,:),scaleFG*fb(1,:),scaleFG*fb(2,:),0,'k')
      %  title('(2) Boundary forces-black lines,  cell boundary-red stars','FontSize',15)
      %elseif (repeat_num==3)
      %  scaleFG=scaleF/3;
      %  quiver(xg(:,:,1),xg(:,:,2),scaleFG*fg(:,:,1),scaleFG*fg(:,:,2),0,'b')
      %  title('(3) Fluid forces-blue lines,  cell boundary-red stars','FontSize',15)
      %elseif (repeat_num==4)
      %  scaleFV=scaleV/1.5;
      %  quiver(xg(1:2:end,1:2:end,1),xg(1:2:end,1:2:end,2),...
      %  scaleFV*vg(1:2:end,1:2:end,1),scaleFV*vg(1:2:end,1:2:end,2),0,'b')
      %  title('(4) Fluid velocities-blue lines,  cell boundary-red stars','FontSize',15)
      %elseif (repeat_num==5)
      %  scaleFV=scaleV/1.5;
      %  quiver(xb(1,:),xb(2,:),scaleFV*ub(1,:),scaleFV*ub(2,:),0,'k')
      %  title('(5) Boundary velocities-black lines,  cell boundary-red stars','FontSize',15)
      %end
      %axis([xmin,xmax,xmin,xmax])
      %axis equal
      %axis([xmin,xmax,xmin,xmax]*1.1)
      %pause(0.2)

    end   % for loop_num
  end  % for repeat_num


end % end function
%--------------------------------------------------------------------%


%--------------------------------------------------------------------%
% draw cell pseudo-skeleton                                          %
%--------------------------------------------------------------------%
function DrawCellCytoskeleton(xb,Nb,cen,connect)
  if (connect==0)
    for ii=1:Nb-1
      plot([xb(1,ii),xb(1,ii+1)],[xb(2,ii),xb(2,ii+1)],'Color',[0.5,0.5,0.5])
    end
    plot([xb(1,Nb),xb(1,1)],[xb(2,Nb),xb(2,1)],'Color',[0.5,0.5,0.5])
  elseif (connect==1)
    for ii=1:Nb-2
      plot([xb(1,ii),xb(1,ii+2)],[xb(2,ii),xb(2,ii+2)],'Color',[0.5,0.5,0.5])
    end
    plot([xb(1,Nb-1),xb(1,1)],[xb(2,Nb-1),xb(2,1)],'Color',[0.5,0.5,0.5])
    plot([xb(1,Nb)  ,xb(1,2)],[xb(2,Nb)  ,xb(2,2)],'Color',[0.5,0.5,0.5])
  elseif (connect==2)
    for ii=1:Nb/2
      plot([xb(1,ii),cen],[xb(2,ii),cen],'Color',[0.5,0.5,0.5])
    end
  elseif (connect==3)
    Nbb=Nb/4;
    for ii=1:Nbb+1
      plot([xb(1,ii),xb(1,3*Nbb+2-ii)],[xb(2,ii),xb(2,3*Nbb+2-ii)],'Color',[0.5,0.5,0.5])
    end
  end
end % function DrawCellCytoskeleton
%--------------------------------------------------------------------%

%--------------------------------------------------------------------%
% define adjacent forces between adjacent boundary points:           %
% (n-th to n-1 and n+1)                                              %
%--------------------------------------------------------------------%
function fbb=AdjacentForces(xb,Nb,hb,Spr)
  fbb=zeros(size(xb));
  for ii=1:Nb
    Lrest=hb;
    if (ii==1)
      dl1=xb(1,Nb)  -xb(1,ii);
      dl2=xb(2,Nb)  -xb(2,ii);
      dr1=xb(1,ii+1)-xb(1,ii);
      dr2=xb(2,ii+1)-xb(2,ii);
    elseif (ii==Nb)
      dl1=xb(1,ii-1)-xb(1,ii);
      dl2=xb(2,ii-1)-xb(2,ii);
      dr1=xb(1,1)   -xb(1,ii);
      dr2=xb(2,1)   -xb(2,ii);
    else
      dl1=xb(1,ii-1)-xb(1,ii);
      dl2=xb(2,ii-1)-xb(2,ii);
      dr1=xb(1,ii+1)-xb(1,ii);
      dr2=xb(2,ii+1)-xb(2,ii);
    end
    ndl=sqrt(dl1^2+dl2^2);
    ndr=sqrt(dr1^2+dr2^2);

    fbb(1,ii)=fbb(1,ii)+(Spr*(ndl-Lrest)*dl1/ndl)+(Spr*(ndr-Lrest)*dr1/ndr);
    fbb(2,ii)=fbb(2,ii)+(Spr*(ndl-Lrest)*dl2/ndl)+(Spr*(ndr-Lrest)*dr2/ndr);
  end
end % function AdjacentForces
%----------------------------------------------------------------------%

%----------------------------------------------------------------------%
% define secondary adjacent forces between second to adjacent boundary %
% points: (n-th to n-2 and n+2)                                        %
%----------------------------------------------------------------------%
function fbb=SecondaryForces(xb,Nb,Lrest,Spr,connect)
  fbb=zeros(size(xb));
  for ii=1:Nb
    if (connect==1)
      Lrest=2*hb;
      if (ii==1)
        dl1=xb(1,Nb-1)-xb(1,ii);
        dl2=xb(2,Nb-1)-xb(2,ii);
        dr1=xb(1,ii+2)-xb(1,ii);
        dr2=xb(2,ii+2)-xb(2,ii);
      elseif (ii==2)
        dl1=xb(1,Nb)-xb(1,ii);
        dl2=xb(2,Nb)-xb(2,ii);
        dr1=xb(1,ii+2)-xb(1,ii);
        dr2=xb(2,ii+2)-xb(2,ii);
      elseif (ii==Nb-1)
        dl1=xb(1,ii-2)-xb(1,ii);
        dl2=xb(2,ii-2)-xb(2,ii);
        dr1=xb(1,1)-xb(1,ii);
        dr2=xb(2,1)-xb(2,ii);
      elseif (ii==Nb)
        dl1=xb(1,ii-2)-xb(1,ii);
        dl2=xb(2,ii-2)-xb(2,ii);
        dr1=xb(1,2)-xb(1,ii);
        dr2=xb(2,2)-xb(2,ii);
      else
        dl1=xb(1,ii-2)-xb(1,ii);
        dl2=xb(2,ii-2)-xb(2,ii);
        dr1=xb(1,ii+2)-xb(1,ii);
        dr2=xb(2,ii+2)-xb(2,ii);
      end
      ndl=sqrt(dl1^2+dl2^2);
      ndr=sqrt(dr1^2+dr2^2);

      fbb(1,ii)=fbb(1,ii)+(Spr*(ndl-Lrest)*dl1/ndl)+(Spr*(ndr-Lrest)*dr1/ndr);
      fbb(2,ii)=fbb(2,ii)+(Spr*(ndl-Lrest)*dl2/ndl)+(Spr*(ndr-Lrest)*dr2/ndr);
    end
  end
end % SecondaryForces
%--------------------------------------------------------------------%

%----------------------------------------------------------------------%
% define central forces between boundary points and the cell nucleus   %
% only holf of all boundary points are connected to cell nucleus to    %
% allow the whole cell to grow                                         %
%----------------------------------------------------------------------%
function fbb=CenterForces(xb,Nb,cen,len,Spr,connect)
  fbb=zeros(size(xb));
  for ii=1:Nb/2
    if (connect==2)
      Lrest=len;

      dl1=cen-xb(1,ii);
      dl2=cen-xb(2,ii);
      ndl=sqrt(dl1^2+dl2^2);

      fbb(1,ii)=fbb(1,ii)+(Spr*(ndl-Lrest)*dl1/ndl);
      fbb(2,ii)=fbb(2,ii)+(Spr*(ndl-Lrest)*dl2/ndl);
    end
  end
end % function CenterForces
%--------------------------------------------------------------------%

%--------------------------------------------------------------------%
% define opposite forces between opposite boundary points            %
%--------------------------------------------------------------------%
function fbb=OppositeForces(xb,Nb,len,Spr,connect)
  fbb=zeros(size(xb));
  for ii=1:Nb
    if (connect==3)
      Lrest=2*len;
      Nb2=Nb/4;

      if (ii <= (Nb2+1))
        dl1=xb(1,3*Nb2+2-ii)-xb(1,ii);
        dl2=xb(2,3*Nb2+2-ii)-xb(2,ii);
        dr1=-dl1;
        dr2=-dl2;

        ndl=sqrt(dl1^2+dl2^2);
        ndr=sqrt(dr1^2+dr2^2);

        fbb(1,ii)=fbb(1,ii)+(Spr*(ndl-Lrest)*dl1/ndl);
        fbb(2,ii)=fbb(2,ii)+(Spr*(ndl-Lrest)*dl2/ndl);

        fbb(1,3*Nb2+2-ii)=fbb(1,3*Nb2+2-ii)+(Spr*(ndr-Lrest)*dr1/ndr);
        fbb(2,3*Nb2+2-ii)=fbb(2,3*Nb2+2-ii)+(Spr*(ndr-Lrest)*dr2/ndr);
      end
    end
  end
end % function OppositeForces
%--------------------------------------------------------------------%

%--------------------------------------------------------------------%
% defne cell shape depending on the value of parameters shape:       %
% shape=1  cell is circular; shape==2  cell is a square;  Nb defines %
% the number of boundary points; len defines a radius for a circular %
% cell; or a length of a side of the square                          %
%--------------------------------------------------------------------%
function xb=DefineCellShape(shape,Nb,len)
  if (shape==1)          %% circle
    hb=2*pi/Nb;
    for ii=1:Nb
      xb(1,ii)=len*cos((ii-1)*hb);
      xb(2,ii)=len*sin((ii-1)*hb);
    end
  elseif (shape==2)       %% square
    Nbb=Nb/4; hb=2*len/Nbb;
    for ii=1:Nbb
      xb(1,ii)     =-len+(ii-1)*hb; xb(2,ii)     =-len;
      xb(1,Nbb+ii)  = len;           xb(2,Nbb+ii)  =-len+(ii-1)*hb;
      xb(1,2*Nbb+ii)= len-(ii-1)*hb; xb(2,2*Nbb+ii)= len;
      xb(1,3*Nbb+ii)=-len;           xb(2,3*Nbb+ii)= len-(ii-1)*hb;
    end
  end   % cell shape
end  % function DefineCellShape
%-----------------------------------------------------------------%

%-----------------------------------------------------------------%
% spreads the material values sb(Nb) (forces, sources) defined at %
% material points xb(2,Nb) to the fluid grid sg(Ng+1,Ng+1) in the %
% square domain [xmin,xmax]^2 with mesh width hg, material points %
% separation hb and a radius of the discrete delta function hdl.  %
%-----------------------------------------------------------------%
function sg=BoundToGrid1(xb,sb,Nb,Ng,hdl,hg,hb,xmin,xmax)
  sg=zeros(Ng+1,Ng+1);

  % passive value - do nothing
  pas=-100;

  for n3=1:Nb
    % move points into the domain
    xbb1=IntoDom(xb(1,n3),xmin,xmax);
    xbb2=IntoDom(xb(2,n3),xmin,xmax);

    % determine indeces of the nearest lower-down grid point
    Nx=1+floor((xbb1-xmin)/hg);
    Ny=1+floor((xbb2-xmin)/hg);

    % tests all 16 possible grid points
    for ii=-1:2
      for jj=-1:2
        % compute the interpolation Delta function
        llx=xmin+(Nx-1)*hg+ii*hg;
        rr=abs(xbb1-llx);
        dx=DeltaFun(rr,hdl);
        lly=xmin+(Ny-1)*hg+jj*hg;
        rr=abs(xbb2-lly);
        dy=DeltaFun(rr,hdl);

        % determine indices of the grid points to update
        [x1,x2]=IndDel(llx,ii,Nx,Ng,xmin,xmax);
        [y1,y2]=IndDel(lly,jj,Ny,Ng,xmin,xmax);

        % update the values if poits are not pasive
        if (dx*dy > 0)
          sg(x1,y1)  =sg(x1,y1)  +sb(1,n3)*dx*dy*hb;
          if (x2 ~= pas)
            sg(x2,y1) =sg(x2,y1) +sb(1,n3)*dx*dy*hb;
          end
          if (y2 ~= pas)
            sg(x1,y2) =sg(x1,y2) +sb(1,n3)*dx*dy*hb;
          end
          if ((x2 ~= pas) & (y2 ~= pas))
            sg(x2,y2) =sg(x2,y2) +sb(1,n3)*dx*dy*hb;
          end
        end

      end  % for jj
    end % for ii
  end % for n3
end % function BoundToGrid1
%-------------------------------------------------------------------%

%-------------------------------------------------------------------%
% spreads the material values sb(2,Nb) (forces, sources) defined at %
% material points xb(2,Nb) to the fluid grid sg(Ng+1,Ng+1,2) in the %
% square domain [xmin,xmax]^2 with mesh width hg, material points   %
% separation hb and a radius of the discrete delta function hdl.    %
%-------------------------------------------------------------------%
function sg=BoundToGrid2(xb,sb,Nb,Ng,hdl,hg,hb,xmin,xmax)
  sg=zeros(Ng+1,Ng+1,2);

  % passive value - do nothing
  pas=-100;

  for n3=1:Nb
    % move points into the domain
    xbb1=IntoDom(xb(1,n3),xmin,xmax);
    xbb2=IntoDom(xb(2,n3),xmin,xmax);
 
    % determine indeces of the nearest lower-down grid point
    Nx=1+floor((xbb1-xmin)/hg);
    Ny=1+floor((xbb2-xmin)/hg);

    % tests all 16 possible grid points
    for ii=-1:2
      for jj=-1:2
        % compute the interpolation Delta function
        llx=xmin+(Nx-1)*hg+ii*hg;
        rr=abs(xbb1-llx);
        dx=DeltaFun(rr,hdl);
        lly=xmin+(Ny-1)*hg+jj*hg;
        rr=abs(xbb2-lly);
        dy=DeltaFun(rr,hdl);

        % determine indices of the grid points to update
        [x1,x2]=IndDel(llx,ii,Nx,Ng,xmin,xmax);
        [y1,y2]=IndDel(lly,jj,Ny,Ng,xmin,xmax);

        % update the values if poits are not pasive
        if (dx*dy > 0)
          sg(x1,y1,1)=sg(x1,y1,1)+sb(1,n3)*dx*dy*hb;
          sg(x1,y1,2)=sg(x1,y1,2)+sb(2,n3)*dx*dy*hb;

          if (x2 ~= pas)
            sg(x2,y1,1)=sg(x2,y1,1)+sb(1,n3)*dx*dy*hb;
            sg(x2,y1,2)=sg(x2,y1,2)+sb(2,n3)*dx*dy*hb;
          end
          if (y2 ~= pas)
            sg(x1,y2,1)=sg(x1,y2,1)+sb(1,n3)*dx*dy*hb;
            sg(x1,y2,2)=sg(x1,y2,2)+sb(2,n3)*dx*dy*hb;
          end
          if ((x2 ~= pas) & (y2 ~= pas))
            sg(x2,y2,1)=sg(x2,y2,1)+sb(1,n3)*dx*dy*hb;
            sg(x2,y2,2)=sg(x2,y2,2)+sb(2,n3)*dx*dy*hb;
          end
        end

      end % jj
    end % ii
  end % for n3
end % function BoundToGride
%---------------------------------------------------------------------%

%---------------------------------------------------------------------%
% Uses the Fast Fourier Method to solve the Navier-Stokes equations   %
% for updating grid velocity ug due to grid forces fg and grid source %
% distribution sg, rho and mu are fluid constants, dt is a time step, %
% hg is a mesh width.                                                 %
%---------------------------------------------------------------------%
function  vg=NavierStokes(ug,fg,sg,Ng,rho,mu,dt,hg)
  vg=zeros(size(fg));

  % stage n terms: force density fg, source distribution sg and current
  % velocity ug
  for n1=1:Ng+1
    for n2=1:Ng+1
      for ik=1:2
        % upwind scheme for the advection term
        if (ug(n1,n2,1) < 0)
          in1=PeriodInd(n1,Ng,1);
          pom=ug(in1,n2,ik)-ug(n1,n2,ik);
        else
          in1=PeriodInd(n1,Ng,-1);
          pom=ug(n1,n2,ik)-ug(in1,n2,ik);
        end
        vg(n1,n2,ik)=ug(n1,n2,1)*pom;

        if (ug(n1,n2,2) < 0)
          in2=PeriodInd(n2,Ng,1);
          pom=ug(n1,in2,ik)-ug(n1,n2,ik);
        else
          in2=PeriodInd(n2,Ng,-1);
          pom=ug(n1,n2,ik)-ug(n1,in2,ik);
        end
        vg(n1,n2,ik)=vg(n1,n2,ik)+ug(n1,n2,2)*pom;
        vg(n1,n2,ik)=-dt*vg(n1,n2,ik)/hg;

        % central difference for the grad of source term
        if (ik == 1)
          in1=PeriodInd(n1,Ng,1);
          in2=PeriodInd(n1,Ng,-1);
          pom=sg(in1,n2)-sg(in2,n2);
        elseif (ik == 2)
          in1=PeriodInd(n2,Ng,1);
          in2=PeriodInd(n2,Ng,-1);
          pom=sg(n1,in1)-sg(n1,in2);
        end
        vg(n1,n2,ik)=vg(n1,n2,ik)+dt*mu*pom/(6*hg*rho*rho);

        % current vlocity and force terms
        vg(n1,n2,ik)=vg(n1,n2,ik)+ug(n1,n2,ik)+dt*fg(n1,n2,ik)/rho;
      end % for ik
    end % for n2
  end % for n1

  % the Fast Fourier transforms of source distribution sg and stage n term vg
  disp(vg(:,:,1))
  fsg =fft2(sg(1:Ng,1:Ng));
  fug1=fft2(vg(1:Ng,1:Ng,1));
  fug2=fft2(vg(1:Ng,1:Ng,2));
  disp(fug1(5,5))
  % determines fug - the Fourier Transform of the velocity field at the stage n+1
  Eps=0.0000001;
  for n1=0:Ng-1
    for n2=0:Ng-1
      B1=sin(2*pi*n1/Ng);
      B2=sin(2*pi*n2/Ng);
      Bb=B1^2+B2^2;
      Aa=1+4*mu*dt*(sin(pi*n1/Ng)^2+sin(pi*n2/Ng)^2)/(rho*hg*hg);

      if (Bb < Eps)
        fvg1(n1+1,n2+1,1)=real(fug1(n1+1,n2+1))/Aa;
        fvg1(n1+1,n2+1,2)=imag(fug1(n1+1,n2+1))/Aa;
        fvg2(n1+1,n2+1,1)=real(fug2(n1+1,n2+1))/Aa;
        fvg2(n1+1,n2+1,2)=imag(fug2(n1+1,n2+1))/Aa;
      else
        Bv=B1*real(fug1(n1+1,n2+1))+B2*real(fug2(n1+1,n2+1));
        fvg1(n1+1,n2+1,1)=(Bb*real(fug1(n1+1,n2+1))-B1*Bv)/(Aa*Bb);
        fvg2(n1+1,n2+1,1)=(Bb*real(fug2(n1+1,n2+1))-B2*Bv)/(Aa*Bb);

        Bv=B1*imag(fug1(n1+1,n2+1))+B2*imag(fug2(n1+1,n2+1));
        fvg1(n1+1,n2+1,2)=(Bb*imag(fug1(n1+1,n2+1))-B1*Bv)/(Aa*Bb);
        fvg2(n1+1,n2+1,2)=(Bb*imag(fug2(n1+1,n2+1))-B2*Bv)/(Aa*Bb);

        fvg1(n1+1,n2+1,1)=fvg1(n1+1,n2+1,1)+hg*B1*imag(fsg(n1+1,n2+1))/(Bb*rho);
        fvg1(n1+1,n2+1,2)=fvg1(n1+1,n2+1,2)-hg*B1*real(fsg(n1+1,n2+1))/(Bb*rho);
        fvg2(n1+1,n2+1,1)=fvg2(n1+1,n2+1,1)+hg*B2*imag(fsg(n1+1,n2+1))/(Bb*rho);
        fvg2(n1+1,n2+1,2)=fvg2(n1+1,n2+1,2)-hg*B2*real(fsg(n1+1,n2+1))/(Bb*rho);
      end
    end % for n2
  end % for n1
  %disp(fvg1(5,5,1))
  %disp(fvg1(5,5,2))
  % the inverse Fast Fourier Method of fvg
  fvgg=fvg1(1:Ng,1:Ng,1)+i*fvg1(1:Ng,1:Ng,2);  
  vg1=ifft2(fvgg);  
  fvgg=fvg2(1:Ng,1:Ng,1)+i*fvg2(1:Ng,1:Ng,2);
  vg2=ifft2(fvgg);
  for ii=1:Ng
    for jj=1:Ng
      vg(ii,jj,1)=real(vg1(ii,jj));
      vg(ii,jj,2)=real(vg2(ii,jj));
    end
  end
  vg(Ng+1,:,:)=vg(1,:,:);
  vg(:,Ng+1,:)=vg(:,1,:);

end % function NavierStokes
%-----------------------------------------------------------------------%

%-----------------------------------------------------------------------%
% interpolates the grid values fg(Ng,Ng,2) (velocities) to the material %
% values fb(2,Nb) defined at the material points xb(2,Nb) in the square %
% domain [xmn,xmx]^2 with mesh width hg and a radius of the discrete    %
% Dirac delta hdl.                                                      %
%-----------------------------------------------------------------------%
function fb=GridToBound(xb,Nb,fg,Ng,hdl,hg,xmn,xmx)
  fb=zeros(size(xb));

  for n3=1:Nb
    % moves points into the domain
    xbb1=IntoDom(xb(1,n3),xmn,xmx);
    xbb2=IntoDom(xb(2,n3),xmn,xmx);

    %computes the indices of the nearest down-left grid point
    Nx=1+floor((xbb1-xmn)/hg);
    Ny=1+floor((xbb2-xmn)/hg);

    % test all 16 neighboring points
    for ii=-1:2
      for jj=-1:2
        % determine the value of the interpolation Delta-function
        llx=xmn+(Nx-1)*hg+ii*hg;        
        rr=abs(xbb1-llx);
        dx=DeltaFun(rr,hdl);
        lly=xmn+(Ny-1)*hg+jj*hg;
        rr=abs(xbb2-lly);
        dy=DeltaFun(rr,hdl);

        % determine the indices of grid points gaining positive impact
        [x1,x2]=IndDel(llx,ii,Nx,Ng,xmn,xmx);
        [y1,y2]=IndDel(lly,jj,Ny,Ng,xmn,xmx);

        % update the values if inside the impact domain
        if (dx*dy > 0)
          fb(1,n3)=fb(1,n3)+fg(x1,y1,1)*dx*dy*hg*hg;
          fb(2,n3)=fb(2,n3)+fg(x1,y1,2)*dx*dy*hg*hg;
        end
      end % for jj
    end % for ii
  end % for n3
end % function GridToBound
%-------------------------------------------------------------------%

%-------------------------------------------------------------------%
% determines value of the unitary bellshaped discrete approximation %
% to the Dirac delta function of radius h for a point located at    %
% distance r from the center.                                       %
%-------------------------------------------------------------------%
function dist=DeltaFun(r,h)
  if (abs(r) < (2*h))
    dist=0.25*(1+cos(0.5*pi*r/h))/h;
  else
    dist=0;
  end
end % function DeltaFun
%-------------------------------------------------------------------%

%-------------------------------------------------------------------%
% determines indeces in1 and in2 of grid elements which need  to be %
% updated in the interpolation of point ll (shifted by ij according %
% to the reference point Nij). Depending on the location of point   %
% ll, there may be at most two active indeces for each element, due %
% to periodicity of fluid domain.                                   %
%   0 elements to update, if ll inside the domain,                  %
%   1 element  to update, if ll is close to the boundary,           %
%   2 elements to update, if ll is close to the corner.             %
%-------------------------------------------------------------------%
function [in1,in2]=IndDel(ll,ij,Nij,Nmx,xmn,xmx)
  % passive value - do nothing
  pas=-100;

  in2=pas;
  if (ll < xmn)
    in1=Nmx+1+ij;
  elseif ((ll == xmn) | (ll == xmx))
    in1=Nmx+1;
    in2=1;
  elseif (ll > xmx)
    in1=ij;
  else
    in1=Nij+ij;
  end
end % function IndDel
%--------------------------------------------------------------------%

%--------------------------------------------------------------------%
% transforms the real coordinate xy of the body into the correspon-  %
% ding coordinate pom inside the periodic domain [xmn,xmx]x[xmn,xmx] %                             %
%--------------------------------------------------------------------%
function pom=IntoDom(xy,xmin,xmax)
  len=xmax-xmin;
  pom=xy;
  while (pom>xmax)
    pom=pom-len;
  end
  while (pom<xmin)
    pom=pom+len;
  end
end % function IntoDom
%-------------------------------------------------------------------%

%-------------------------------------------------------------------%
% computes shifting of the index nn in the square periodic domain   %
% of size Ng+1 where the 1st and the (Ng+1)st columns and rows are  %
% identical according to the value of which:                        %
%    which = 1 -- nn shifted one element to the right               %
%    which =-1 -- nn shifted one element to the left                %
%-------------------------------------------------------------------%
function ind=PeriodInd(nn,Ng,which)
  ind=0;
  if (which == (-1))
    if (nn == 1)
      ind=Ng;
    else
      ind=nn-1;
    end
  end
  if (which == 1)
    if (nn == Ng+1)
      ind=2;
    else
      ind=nn+1;
    end
  end
end % function PeriodInd
%------------------------------------------------------------------%
