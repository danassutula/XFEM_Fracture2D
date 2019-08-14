
%--------------------------------------------------------------------------
% Figure sizes (for LaTeX) 
%   n.b. see 'Auxiliary_Scope' for sizes for plotting on the go
%--------------------------------------------------------------------------

% N.B. to resize figure, use: FigResize(szfig_std,szfnt_std,h)
% to save figure in pdf, use: SavePDF(figure_name,h)

% full screen
scrsz = get(0,'ScreenSize');
scrsz = scrsz([3,4]);

% figure positions (centered)

szfig_std   = [scrsz(1)/2-280,  scrsz(2)/2-210,     560,       420];
szfig_big   = [scrsz(1)/2-560,  scrsz(2)/2-420,    1020,       840];
szfig_Big   = [scrsz(1)/2-840,  scrsz(2)/2-630,    1680,      1260];
szfig_wde   = [scrsz(1)/2-560,  scrsz(2)/2-210,    1020,       420];
szfig_Wde   = [scrsz(1)/2-840,  scrsz(2)/2-210,    1680,       420];
szfig_tll   = [scrsz(1)/2-210,  scrsz(2)/2-280,     420,       560];

szfig_480p  = [scrsz(1)/2-360,  scrsz(2)/2-240,     720,       480];
szfig_720p  = [scrsz(1)/2-640,  scrsz(2)/2-360,    1280,       720];
szfig_1080p = [scrsz(1)/2-960,  scrsz(2)/2-540,    1920,      1080];

szfig_full  = [             1,              1, scrsz(1), scrsz(2)];

szfig_ppt   = [0,  0,    1020,       840]*1.2;

% font size (for LaTeX)
szfnt_sml =  9;
szfnt_std = 11;
szfnt_big = 22;
szfnt_Big = 44;

%--------------------------------------------------------------------------