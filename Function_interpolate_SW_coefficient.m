%% The actual SW coefficient is obtained through three-dimensional interpolation by using SSP and SVF
% LUT_3D
% SSP
% SVF
% Band_name: The name of the band, Band_name = 1,2,3,4,5

% author HZW

function [A] = Function_interpolate_SW_coefficient(LUT_3D, SSP, SVF, Band_name)

    LUT = LUT_3D;
    row_num = 121; % Look up the row_num row in the table, which represents a channel
    SSP;
    SVF;

    % Channel 10
    if Band_name == 1 
        x=LUT(1:row_num,1)./10;
        y=LUT(1:row_num,2)./10;
        v0=LUT(1:row_num,5);
        v1=LUT(1:row_num,6);
        v2=LUT(1:row_num,7);
        v3=LUT(1:row_num,8);
        v4=LUT(1:row_num,9);
        
        F0 = scatteredInterpolant(x,y,v0);
        A0 = F0(SSP, SVF);
        F1 = scatteredInterpolant(x,y,v1);
        A1 = F1(SSP, SVF);
        F2 = scatteredInterpolant(x,y,v2);
        A2 = F2(SSP, SVF);
        F3 = scatteredInterpolant(x,y,v3);
        A3 = F3(SSP, SVF);
        F4 = scatteredInterpolant(x,y,v4);
        A4 = F4(SSP, SVF);
    
        A = [A0, A1, A2, A3, A4];
    end

        % Channel 11
    if Band_name == 2 
        x=LUT(row_num+1:row_num*2,1)./10;
        y=LUT(row_num+1:row_num*2,2)./10;
        v0=LUT(row_num+1:row_num*2,5);
        v1=LUT(row_num+1:row_num*2,6);
        v2=LUT(row_num+1:row_num*2,7);
        v3=LUT(row_num+1:row_num*2,8);
        v4=LUT(row_num+1:row_num*2,9);
        
        F0 = scatteredInterpolant(x,y,v0);
        A0 = F0(SSP, SVF);
        F1 = scatteredInterpolant(x,y,v1);
        A1 = F1(SSP, SVF);
        F2 = scatteredInterpolant(x,y,v2);
        A2 = F2(SSP, SVF);
        F3 = scatteredInterpolant(x,y,v3);
        A3 = F3(SSP, SVF);
        F4 = scatteredInterpolant(x,y,v4);
        A4 = F4(SSP, SVF);

        A = [A0, A1, A2, A3, A4];
    end

    % Channel 12
    if Band_name == 3 
        x=LUT(2*row_num+1:row_num*3,1)./10;
        y=LUT(2*row_num+1:row_num*3,2)./10;
        v0=LUT(2*row_num+1:row_num*3,5);
        v1=LUT(2*row_num+1:row_num*3,6);
        v2=LUT(2*row_num+1:row_num*3,7);
        v3=LUT(2*row_num+1:row_num*3,8);
        v4=LUT(2*row_num+1:row_num*3,9);

        F0 = scatteredInterpolant(x,y,v0);
        A0 = F0(SSP, SVF);
        F1 = scatteredInterpolant(x,y,v1);
        A1 = F1(SSP, SVF);
        F2 = scatteredInterpolant(x,y,v2);
        A2 = F2(SSP, SVF);
        F3 = scatteredInterpolant(x,y,v3);
        A3 = F3(SSP, SVF);
        F4 = scatteredInterpolant(x,y,v4);
        A4 = F4(SSP, SVF);

        A = [A0, A1, A2, A3, A4]; 
    end

    % Channel 13
    if Band_name == 4 
        x=LUT(3*row_num+1:row_num*4,1)./10;
        y=LUT(3*row_num+1:row_num*4,2)./10;
        v0=LUT(3*row_num+1:row_num*4,5);
        v1=LUT(3*row_num+1:row_num*4,6);
        v2=LUT(3*row_num+1:row_num*4,7);
        v3=LUT(3*row_num+1:row_num*4,8);
        v4=LUT(3*row_num+1:row_num*4,9);
        
        F0 = scatteredInterpolant(x,y,v0);
        A0 = F0(SSP, SVF);
        F1 = scatteredInterpolant(x,y,v1);
        A1 = F1(SSP, SVF);
        F2 = scatteredInterpolant(x,y,v2);
        A2 = F2(SSP, SVF);
        F3 = scatteredInterpolant(x,y,v3);
        A3 = F3(SSP, SVF);
        F4 = scatteredInterpolant(x,y,v4);
        A4 = F4(SSP, SVF);
    
        A = [A0, A1, A2, A3, A4];
    end

    % Channel 14
    if Band_name == 5 
        x=LUT(4*row_num+1:row_num*5,1)./10;
        y=LUT(4*row_num+1:row_num*5,2)./10;
        v0=LUT(4*row_num+1:row_num*5,5);
        v1=LUT(4*row_num+1:row_num*5,6);
        v2=LUT(4*row_num+1:row_num*5,7);
        v3=LUT(4*row_num+1:row_num*5,8);
        v4=LUT(4*row_num+1:row_num*5,9);

        F0 = scatteredInterpolant(x,y,v0);
        A0 = F0(SSP, SVF);
        F1 = scatteredInterpolant(x,y,v1);
        A1 = F1(SSP, SVF);
        F2 = scatteredInterpolant(x,y,v2);
        A2 = F2(SSP, SVF);
        F3 = scatteredInterpolant(x,y,v3);
        A3 = F3(SSP, SVF);
        F4 = scatteredInterpolant(x,y,v4);
        A4 = F4(SSP, SVF);

        A = [A0, A1, A2, A3, A4]; 
    end
    
end