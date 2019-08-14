
%==========================================================================
% Plot Stress (contours)
%==========================================================================


%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

h = gcf;
set(h,'Color','w');

if ~exist('stress_scaleFrac')
    stress_scaleFrac = 1/2;
end

if ~exist('stress_component')
    stress_component = 1; % {1,2,3}
end

if ~exist('stress_axisAngle')
    stress_axisAngle = 0;
end

if ~exist('stress_valueSign')
    stress_valueSign = 0; % -1 (negative), [] (all), 1 (positive)
end

if ~exist('stress_plotBox')
    stress_plotBox = [];
end

fprintf('\n\tstress_scaleFrac = %i\n',stress_scaleFrac)
fprintf('\tstress_component = %i\n',stress_component)
fprintf('\tstress_axisAngle = %i\n',stress_axisAngle)
fprintf('\tstress_valueSign = %i\n',stress_valueSign)

dotsz = 8; % dot/marker size

% what to plot: compression or tension ??

%--------------------------------------------------------------------------
% Evaluate von Mises at Gauss points (if not done already)
%--------------------------------------------------------------------------

if ~exist('mStres_std','var') || isempty(mStres_std) ...
        || ~exist('mGsCrd_std','var') || isempty(mGsCrd_std)
    
    [mStres_std,mGsCrd_std] = Stress_std(mNDspS,mNdCrd,vElStd, ...
     mLNodS,vElPhz,nPhase,cDMatx,mPrLod,mGsStd_omgShp,mGsStd_omgDrv);
    
    mStres_std = single(mStres_std);
    mGsCrd_std = single(mGsCrd_std);
    
    [mStres_enr,mGsCrd_enr] = Stress_enr(mNDspS,mNDspE,mNdCrd,vElEnr,mLNodS, ...
     cLNodE,vElPhz,nPhase,cDMatx,mPrLod,cGsEnr_omgShS,cGsEnr_omgDvS,cGsEnr_omgDvE);
    
    mStres_enr = single(mStres_enr);
    mGsCrd_enr = single(mGsCrd_enr);
    
end

c = cos(stress_axisAngle); s = sin(stress_axisAngle); % dir. cos. of local axis
T = [c^2,s^2,2*s*c;s^2,c^2,-2*s*c;-s*c,s*c,c^2-s^2]; % global to local mapping

vStres = T(stress_component,:)*[mStres_std,mStres_enr];
mGsCrd = [mGsCrd_std;mGsCrd_enr];

%--------------------------------------------------------------------------
% Only plot what's within a bounding box (if at all specified)
%--------------------------------------------------------------------------

if exist('plotBox_stress','var') ...
        && ~isempty(plotBox_stress) ...
        && norm(diff(plotBox_stress,1))>0
    
    q = mGsCrd(:,1) > stress_plotBox(1,1) & mGsCrd(:,1) < stress_plotBox(2,1) & ...
        mGsCrd(:,2) > stress_plotBox(1,2) & mGsCrd(:,2) < stress_plotBox(2,2) ;
    
    vStres = vStres(q);
    mGsCrd = mGsCrd(q,:);

end

%--------------------------------------------------------------------------
% Adjust contrast
%--------------------------------------------------------------------------

if isempty(stress_valueSign) || stress_valueSign==0
    q = true(length(vStres),1);
else % get values of specified sign only
    q = sign(vStres)==sign(stress_valueSign); 
end

vStres = vStres(q);
mGsCrd = mGsCrd(q,:);

%--------------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------------

clf(h); % better to clear the figure just before plotting
hold on;
axis equal;

% plot von Mises stress
scatter(mGsCrd(:,1),mGsCrd(:,2),dotsz, ...
    sign(vStres).*(abs(vStres).^stress_scaleFrac),'s','fill')

% plot cracks
for ii = 1:nCrack
    plot(cCkCrd{ii}(:,1),cCkCrd{ii}(:,2),'w','linewidth',1)
end

%--------------------------------------------------------------------------
% Format axis
%--------------------------------------------------------------------------

xmin = min(mGsCrd(:,1));
xmax = max(mGsCrd(:,1));
ymin = min(mGsCrd(:,2));
ymax = max(mGsCrd(:,2));

m = max(xmax-xmin,ymax-ymin)*0.05;
axis([xmin-m,xmax+m,ymin-m,ymax+m])
set(gca,'layer','top','box','on')

if exist('lengthUnits','var')
    xlabel(['x (',lengthUnits,')'])
    ylabel(['y (',lengthUnits,')'])
    title('Stress (scaled)')
end

% set color sinc.
colormap(jet) % bone, hot, cool
colorbar('location','EastOutside');

% clear var.
clear vStres mGsCrd c s T
clear xmin ymax ymin ymax m

figure(h);

%--------------------------------------------------------------------------
