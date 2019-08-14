
%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

h = gcf; hold on;
set(h,'Color','w');

if ~exist('deformed_scale','var')
    deformed_scale = 10; end
if ~exist('deformed_dotsz','var')
    deformed_dotsz = 3; end

deformed_scale
deformed_dotsz

color_enr = [1.0,0.0,0.0];
color_std = [0.9,0.9,0.9];

%--------------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------------

% % non-enriched elements
% if ~isempty(vElStd) 
%     
%     [mGsDsp,mGsCrd] = Disps_std(vElStd);
%     
%     if exist('deformed_plotBox','var')
%         q = mGsCrd(:,1) > deformed_plotBox(1,1) & mGsCrd(:,1) < deformed_plotBox(2,1) & ...
%             mGsCrd(:,2) > deformed_plotBox(1,2) & mGsCrd(:,2) < deformed_plotBox(2,2) ;
%         
%         mGsCrd = mGsCrd(q,:);
%         mGsDsp = mGsDsp(q,:);
%         
%     end
%     
%     mGsCrd = mGsCrd+mGsDsp*deformed_scale;
%     
%     scatter(mGsCrd(:,1),mGsCrd(:,2),deformed_dotsz,color_std,'o','fill')
%     
% end

if 0
    
    % enriched elements
    if ~isempty(vElEnr)
        
        [mGsDsp,mGsCrd] = Disps_enr(vElEnr);
        
        if exist('deformed_plotBox','var')
            q = mGsCrd(:,1) > deformed_plotBox(1,1) & mGsCrd(:,1) < deformed_plotBox(2,1) & ...
                mGsCrd(:,2) > deformed_plotBox(1,2) & mGsCrd(:,2) < deformed_plotBox(2,2) ;
            
            mGsCrd = mGsCrd(q,:);
            mGsDsp = mGsDsp(q,:);
            
        end
        
        mGsCrd = mGsCrd+mGsDsp*deformed_scale;
        scatter(mGsCrd(:,1),mGsCrd(:,2),deformed_dotsz,color_enr,'o','fill')
        
    end

else
    
    elems_cut = unique([cell2mat(cBrStd_eTp(:)); ...
        cell2mat(cHvBln_elm(:));cell2mat(cHvStd_elm(:))]);
    
    if ~isempty(elems_cut)
        
        [mGsDsp,mGsCrd] = Disps_enr(elems_cut);
        
        if exist('deformed_plotBox','var')
            q = mGsCrd(:,1) > deformed_plotBox(1,1) & mGsCrd(:,1) < deformed_plotBox(2,1) & ...
                mGsCrd(:,2) > deformed_plotBox(1,2) & mGsCrd(:,2) < deformed_plotBox(2,2) ;
            
            mGsCrd = mGsCrd(q,:);
            mGsDsp = mGsDsp(q,:);
            
        end
        
        mGsCrd = mGsCrd+mGsDsp*deformed_scale;
        scatter(mGsCrd(:,1),mGsCrd(:,2),deformed_dotsz,color_enr,'o','fill')
        
    end

end

%--------------------------------------------------------------------------
% Format axis
%--------------------------------------------------------------------------

axis equal

if exist('lengthUnits','var')
    xlabel(['x (',lengthUnits,')'])
    ylabel(['y (',lengthUnits,')'])
else
    xlabel('x')
    ylabel('y')
end

clear elems_cut
clear mGsCrd mGsDsp

%--------------------------------------------------------------------------
