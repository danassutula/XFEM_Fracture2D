function [vGsPnt,vGsWgt] = Gauss_DomLin(nGauss)

switch nGauss
    case 1
        vGsPnt = 0;
        vGsWgt = 2;
    case 2
        vGsPnt = [-1;1]./sqrt(3);
        vGsWgt = [1;1];
    case 3
        vGsPnt = [-1;0;1].*sqrt(3/5);
        vGsWgt = [5;8;5]./9;
    case 4
        vGsPnt = [-sqrt((3-2*sqrt(6/5))/7);-sqrt((3+2*sqrt(6/5))/7);sqrt((3+2*sqrt(6/5))/7);sqrt((3-2*sqrt(6/5))/7)];
        vGsWgt = [(18+sqrt(30))/36;(18-sqrt(30))/36;(18-sqrt(30))/36;(18+sqrt(30))/36];
    case 5
        vGsPnt = [-sqrt(5-2*sqrt(10/7))/3;-sqrt(5+2*sqrt(10/7))/3;0;sqrt(5+2*sqrt(10/7))/3;sqrt(5-2*sqrt(10/7))/3];
        vGsWgt = [(322+13*sqrt(70))/900;(322-13*sqrt(70))/900;128/225;(322-13*sqrt(70))/900;(322+13*sqrt(70))/900];
    otherwise
        if nGauss > 5
           [vGsPnt,vGsWgt] = LgWeight(nGauss,-1,1);
        end
end