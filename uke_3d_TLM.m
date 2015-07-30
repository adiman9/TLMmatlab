clear all;
gridSize = 49;
iterations = 20000;

fcalc1 = 232.1901;
fcalc2 = 751.034;
fcalc3 = 1482.003;
fcalc4 = 2138.4419;
fcalc5 = 2422.2298;
fcalc6 = 2809.2132;

delta_L = 0.01;
c = 331;
root_3 = sqrt(3);
delta_t = delta_L/(root_3*c);

helmholtzCutoffFreq = c / (10 * delta_L);

reflectedLeft = zeros(gridSize,gridSize,gridSize);
reflectedRight = zeros(gridSize,gridSize,gridSize);
reflectedFront = zeros(gridSize,gridSize,gridSize);
reflectedBack = zeros(gridSize,gridSize,gridSize);
reflectedTop = zeros(gridSize,gridSize,gridSize);
reflectedBottom = zeros(gridSize,gridSize,gridSize);

incidentLeft = zeros(gridSize,gridSize,gridSize);
incidentRight = zeros(gridSize,gridSize,gridSize);
incidentFront = zeros(gridSize,gridSize,gridSize);
incidentBack = zeros(gridSize,gridSize,gridSize);
incidentTop = zeros(gridSize,gridSize,gridSize);
incidentBottom = zeros(gridSize,gridSize,gridSize);

currentPressure = zeros(gridSize,gridSize,gridSize);

XPressureDFT1 = zeros(gridSize,gridSize);
YPressureDFT1 = zeros(gridSize,gridSize);
ZPressureDFT1 = zeros(gridSize,gridSize);

XPressureDFT2 = zeros(gridSize,gridSize);
YPressureDFT2 = zeros(gridSize,gridSize);
ZPressureDFT2 = zeros(gridSize,gridSize);

XPressureDFT3 = zeros(gridSize,gridSize);
YPressureDFT3 = zeros(gridSize,gridSize);
ZPressureDFT3 = zeros(gridSize,gridSize);

XPressureDFT4 = zeros(gridSize,gridSize);
YPressureDFT4 = zeros(gridSize,gridSize);
ZPressureDFT4 = zeros(gridSize,gridSize);

XPressureDFT5 = zeros(gridSize,gridSize);
YPressureDFT5 = zeros(gridSize,gridSize);
ZPressureDFT5 = zeros(gridSize,gridSize);

XPressureDFT6 = zeros(gridSize,gridSize);
YPressureDFT6 = zeros(gridSize,gridSize);
ZPressureDFT6 = zeros(gridSize,gridSize);

Xpressuregraph = zeros(gridSize,gridSize, iterations);
Ypressuregraph = zeros(gridSize,gridSize, iterations);
Zpressuregraph = zeros(gridSize,gridSize, iterations);

measuredPressure = (1:iterations);
alpha1 = (1:iterations);
alpha2 = (1:iterations);
alpha3 = (1:iterations);
alpha4 = (1:iterations);
alpha5 = (1:iterations);
alpha6 = (1:iterations);

ukeback = imread('uke back 1 49.bmp');
ukefront = imread('uke front 1 49.bmp');
ukeside = imread('uke side 1 49.bmp');

reflection = (1-sqrt(3))/(1+sqrt(3));

% figure('units','normalized','outerposition',[0 0 1 1])

iterstring = num2str(iterations);


for x = 1:iterations
    alpha1(x) = exp(-2*i*pi*fcalc1*delta_t*x);
    alpha2(x) = exp(-2*i*pi*fcalc2*delta_t*x);
    alpha3(x) = exp(-2*i*pi*fcalc3*delta_t*x);
    alpha4(x) = exp(-2*i*pi*fcalc4*delta_t*x);
    alpha5(x) = exp(-2*i*pi*fcalc5*delta_t*x);
    alpha6(x) = exp(-2*i*pi*fcalc6*delta_t*x);
end


for x = 1:iterations
    
%     f=0.05;
%     
%     level=sin(2*pi*f*x);
    
    if x == 1
        level=1;
    else
        level=0;
    end;
    
    incidentBack(10, 25, 15) = incidentBack(10, 25, 15) + level;
    incidentFront(10, 25, 15) = incidentFront(10, 25, 15)+ level;
    incidentLeft(10, 25, 15) = incidentLeft(10, 25, 15)+ level;
    incidentRight(10, 25, 15) = incidentRight(10, 25, 15) + level;
    incidentTop(10, 25, 15) = incidentTop(10, 25, 15) + level;
    incidentBottom(10, 25, 15) = incidentBottom(10, 25, 15) + level;
    
    
    for i = 1:gridSize
       
       for j = 1:gridSize
              
           for z = 1:gridSize
               
               if (z == 11 && ukefront(j, i) == 0) || (z == 5 && ukeback(j, i) == 0) || (z > 5 && z < 11 && ukeside(j, i) == 0)
                
                   reflectedFront(z, j, i) = incidentFront(z, j, i);
                   reflectedBack(z, j, i) = incidentBack(z, j, i);
                   reflectedLeft(z, j, i) = incidentLeft(z, j, i);
                   reflectedRight(z, j, i) = incidentRight(z, j, i);
                   reflectedTop(z, j, i) = incidentTop(z, j, i);
                   reflectedBottom(z, j, i) = incidentBottom(z, j, i);
               
               else
                   
                   reflectedFront(z, j, i) = (1/3) * (incidentRight(z, j, i) - 2*incidentFront(z, j, i) + incidentLeft(z, j, i) + incidentBack(z, j, i) + incidentTop(z, j, i) + incidentBottom(z, j, i));
                   reflectedBack(z, j, i) = (1/3) * (incidentRight(z, j, i) + incidentFront(z, j, i) + incidentLeft(z, j, i) - 2*incidentBack(z, j, i) + incidentTop(z, j, i) + incidentBottom(z, j, i));
                   reflectedLeft(z, j, i) = (1/3) * (incidentRight(z, j, i) + incidentFront(z, j, i) - 2*incidentLeft(z, j, i) + incidentBack(z, j, i) + incidentTop(z, j, i) + incidentBottom(z, j, i));
                   reflectedRight(z, j, i) = (1/3) * (-2*incidentRight(z, j, i) + incidentFront(z, j, i) + incidentLeft(z, j, i) + incidentBack(z, j, i) + incidentTop(z, j, i) + incidentBottom(z, j, i));
                   reflectedTop(z, j, i) = (1/3) * (incidentRight(z, j, i) + incidentFront(z, j, i) + incidentLeft(z, j, i) + incidentBack(z, j, i) -2* incidentTop(z, j, i) + incidentBottom(z, j, i));
                   reflectedBottom(z, j, i) = (1/3) * (incidentRight(z, j, i) + incidentFront(z, j, i) + incidentLeft(z, j, i) + incidentBack(z, j, i) + incidentTop(z, j, i) - 2*incidentBottom(z, j, i));
                   
               end
           end
           
       end
      
    end 
   
    
    for i = 1:gridSize
       
       for j = 1:gridSize
              
           for z = 1:gridSize
               
                   %bottom corners, edges and face
                   if j == 1
                       
                       if i == 1
                           %bottom front left corner
                           if z == 1
                               incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                               incidentBack(z, j, i) = reflection * reflectedBack(z, j, i);
                               incidentLeft(z, j, i) = reflection * reflectedLeft(z, j, i);
                               incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                               incidentTop(z, j, i) = reflection * reflectedTop(z, j, i);
                               incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                            
                           %bottom back left corner
                           elseif z == gridSize
                               incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                               incidentBack(z, j, i) = reflection * reflectedBack(z, j, i);
                               incidentLeft(z, j, i) = reflection * reflectedLeft(z, j, i);
                               incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                               incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                               incidentBottom(z, j, i) = reflection * reflectedBottom(z, j, i);
                           %bottom left edge
                           else
                               incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                               incidentBack(z, j, i) = reflection * reflectedBack(z, j, i);
                               incidentLeft(z, j, i) = reflection * reflectedLeft(z, j, i);
                               incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                               incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                               incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                           end
                           
                       elseif i == gridSize
                           
                           %bottom front right corner
                           if z == 1
                               incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                               incidentBack(z, j, i) = reflection * reflectedBack(z, j, i);
                               incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                               incidentRight(z, j, i) = reflection * reflectedRight(z, j, i);
                               incidentTop(z, j, i) = reflection * reflectedTop(z, j, i);
                               incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                           %bottom back right corner
                           elseif z == gridSize
                               incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                               incidentBack(z, j, i) = reflection * reflectedBack(z, j, i);
                               incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                               incidentRight(z, j, i) = reflection * reflectedRight(z, j, i);
                               incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                               incidentBottom(z, j, i) = reflection * reflectedBottom(z, j, i);
                           %bottom right edge
                           else
                               incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                               incidentBack(z, j, i) = reflection * reflectedBack(z, j, i);
                               incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                               incidentRight(z, j, i) = reflection * reflectedRight(z, j, i);
                               incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                               incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                           end
                       
                       %bottom front edge
                       elseif z == 1
                           incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                           incidentBack(z, j, i) = reflection * reflectedBack(z, j, i);
                           incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                           incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                           incidentTop(z, j, i) = reflection * reflectedTop(z, j, i);
                           incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                       %bottom back edge
                       elseif z == gridSize
                           incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                           incidentBack(z, j, i) = reflection * reflectedBack(z, j, i);
                           incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                           incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                           incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                           incidentBottom(z, j, i) = reflection * reflectedBottom(z, j, i);
                       %bottom face
                       else
                           incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                           incidentBack(z, j, i) = reflection * reflectedBack(z, j, i);
                           incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                           incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                           incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                           incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                       end
                   
                   %top corners, edges and face
                   elseif j == gridSize
                       
                       
                       if i == 1
                           %top front left corner
                           if z == 1
                               incidentFront(z, j, i) = reflection * reflectedFront(z, j, i);
                               incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                               incidentLeft(z, j, i) = reflection * reflectedLeft(z, j, i);
                               incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                               incidentTop(z, j, i) = reflection * reflectedTop(z, j, i);
                               incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                            
                           %top back left corner
                           elseif z == gridSize
                               incidentFront(z, j, i) = reflection * reflectedFront(z, j, i);
                               incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                               incidentLeft(z, j, i) = reflection * reflectedLeft(z, j, i);
                               incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                               incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                               incidentBottom(z, j, i) = reflection * reflectedBottom(z, j, i);
                           %top left edge
                           else
                               incidentFront(z, j, i) = reflection * reflectedFront(z, j, i);
                               incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                               incidentLeft(z, j, i) = reflection * reflectedLeft(z, j, i);
                               incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                               incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                               incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                           end
                           
                       elseif i == gridSize
                           
                           %top front right corner
                           if z == 1
                               incidentFront(z, j, i) = reflection * reflectedFront(z, j, i);
                               incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                               incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                               incidentRight(z, j, i) = reflection * reflectedRight(z, j, i);
                               incidentTop(z, j, i) = reflection * reflectedTop(z, j, i);
                               incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                           %top back right corner
                           elseif z == gridSize
                               incidentFront(z, j, i) = reflection * reflectedFront(z, j, i);
                               incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                               incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                               incidentRight(z, j, i) = reflection * reflectedRight(z, j, i);
                               incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                               incidentBottom(z, j, i) = reflection * reflectedBottom(z, j, i);
                           %top right edge
                           else
                               incidentFront(z, j, i) = reflection * reflectedFront(z, j, i);
                               incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                               incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                               incidentRight(z, j, i) = reflection * reflectedRight(z, j, i);
                               incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                               incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                           end
                       
                       %top front edge
                       elseif z == 1
                           incidentFront(z, j, i) = reflection * reflectedFront(z, j, i);
                           incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                           incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                           incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                           incidentTop(z, j, i) = reflection * reflectedTop(z, j, i);
                           incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                       %top back edge
                       elseif z == gridSize
                           incidentFront(z, j, i) = reflection * reflectedFront(z, j, i);
                           incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                           incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                           incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                           incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                           incidentBottom(z, j, i) = reflection * reflectedBottom(z, j, i);
                       %top face
                       else
                           incidentFront(z, j, i) = reflection * reflectedFront(z, j, i);
                           incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                           incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                           incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                           incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                           incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                       end
                       
                   elseif i == 1
                       
                       %front left edge
                       if z == 1
                           incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                           incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                           incidentLeft(z, j, i) = reflection * reflectedLeft(z, j, i);
                           incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                           incidentTop(z, j, i) = reflection * reflectedTop(z, j, i);
                           incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                          
                       %back left edge
                       elseif z == gridSize
                           incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                           incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                           incidentLeft(z, j, i) = reflection * reflectedLeft(z, j, i);
                           incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                           incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                           incidentBottom(z, j, i) = reflection * reflectedBottom(z, j, i);
                       %left face
                       else
                           incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                           incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                           incidentLeft(z, j, i) = reflection * reflectedLeft(z, j, i);
                           incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                           incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                           incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                       end
                       
                   elseif i == gridSize
                       
                       %front right edge
                       if z == 1
                           incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                           incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                           incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                           incidentRight(z, j, i) = reflection * reflectedRight(z, j, i);
                           incidentTop(z, j, i) = reflection * reflectedTop(z, j, i);
                           incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                       %back right edge
                       elseif z == gridSize
                           incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                           incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                           incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                           incidentRight(z, j, i) = reflection * reflectedRight(z, j, i);
                           incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                           incidentBottom(z, j, i) = reflection * reflectedBottom(z, j, i);
                       %right face
                       else
                           incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                           incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                           incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                           incidentRight(z, j, i) = reflection * reflectedRight(z, j, i);
                           incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                           incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                       end
                     
                   %front face
                   elseif z == 1
                       incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                       incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                       incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                       incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                       incidentTop(z, j, i) = reflection * reflectedTop(z, j, i);
                       incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                     
                   %back face
                   elseif z == gridSize
                       incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                       incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                       incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                       incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                       incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                       incidentBottom(z, j, i) = reflection * reflectedBottom(z, j, i);
                       
                   %rest of the freespace
                   else
                       incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                       incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                       incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                       incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                       incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                       incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
                   end
               
               currentPressure(z, j, i) = reflectedRight(z, j, i) + reflectedLeft(z, j, i) + reflectedFront(z, j, i) + reflectedBack(z, j, i) + reflectedTop(z, j, i) + reflectedBottom(z, j, i);
    
                            
               if j == 25
                   
                   Ypressuregraph(z, i, x) = currentPressure(z, j, i);
                   
                   YPressureDFT1(z, i) = Ypressuregraph(z, i, x) * alpha1(x) + YPressureDFT1(z, i);
                   YPressureDFT2(z, i) = Ypressuregraph(z, i, x) * alpha2(x) + YPressureDFT2(z, i);
                   YPressureDFT3(z, i) = Ypressuregraph(z, i, x) * alpha3(x) + YPressureDFT3(z, i);
                   YPressureDFT4(z, i) = Ypressuregraph(z, i, x) * alpha4(x) + YPressureDFT4(z, i);
                   YPressureDFT5(z, i) = Ypressuregraph(z, i, x) * alpha5(x) + YPressureDFT5(z, i);
                   YPressureDFT6(z, i) = Ypressuregraph(z, i, x) * alpha6(x) + YPressureDFT6(z, i);
                   
                   if (z == 5 && i > 13 && i < 37) || ((i == 14 || i == 36) && z > 4 && z < 12) || (z == 11 && ((i > 31 && i < 37) || (i > 13 && i < 28))) 
                       Ypressuregraph(z, i, x) = 0.09;
                   end

               end
               
               if i == 25
                   
                   Zpressuregraph(z, j, x) = currentPressure(z, j, i);
                   
                   ZPressureDFT1(z, j) = Zpressuregraph(z, j, x) * alpha1(x) + ZPressureDFT1(z, j);
                   ZPressureDFT2(z, j) = Zpressuregraph(z, j, x) * alpha2(x) + ZPressureDFT2(z, j);
                   ZPressureDFT3(z, j) = Zpressuregraph(z, j, x) * alpha3(x) + ZPressureDFT3(z, j);
                   ZPressureDFT4(z, j) = Zpressuregraph(z, j, x) * alpha4(x) + ZPressureDFT4(z, j);
                   ZPressureDFT5(z, j) = Zpressuregraph(z, j, x) * alpha5(x) + ZPressureDFT5(z, j);
                   ZPressureDFT6(z, j) = Zpressuregraph(z, j, x) * alpha6(x) + ZPressureDFT6(z, j);
                   
                    if ((z == 5 || z == 11) && j > 17 && j < 32) || ((j == 18 || j == 31) && z > 5 && z < 11)
                       Zpressuregraph(z, j, x) = 0.09;
                    end
                   
               end
               
               if z == 7
                  
                   Xpressuregraph(j, i, x) = currentPressure(z, j, i);
                   
                   XPressureDFT1(j, i) = Xpressuregraph(j, i, x) * alpha1(x) + XPressureDFT1(j, i);
                   XPressureDFT2(j, i) = Xpressuregraph(j, i, x) * alpha2(x) + XPressureDFT2(j, i);
                   XPressureDFT3(j, i) = Xpressuregraph(j, i, x) * alpha3(x) + XPressureDFT3(j, i);
                   XPressureDFT4(j, i) = Xpressuregraph(j, i, x) * alpha4(x) + XPressureDFT4(j, i);
                   XPressureDFT5(j, i) = Xpressuregraph(j, i, x) * alpha5(x) + XPressureDFT5(j, i);
                   XPressureDFT6(j, i) = Xpressuregraph(j, i, x) * alpha6(x) + XPressureDFT6(j, i);

                   if ukeside(j, i) == 0
                       Xpressuregraph(j, i, x) = 0.09;
                   end
               end    
               
           end    
       end  
    end 
    
   measuredPressure(x) = currentPressure(20, 25, 30);
  
%    Ypressuregraph(20, 30, x) = 0.12;
%    
%    
%    
%    
%    subplot(131)
%    pcolor(Xpressuregraph(:,:,x))
%    title('Sound pressure in the ukulele through the x axis');
%    axis square
%    grid off
%    shading interp
%    caxis([0 0.25])
%    zlim([-2 2])
%    
%    subplot(132)
%    pcolor(Ypressuregraph(:,:,x))
%    title('Sound pressure in the ukulele through the y axis.');
%    axis square
%    grid off
%    shading interp
%    caxis([0 0.25])
%    zlim([-2 2])
%    
%    subplot(133)
%    pcolor(Zpressuregraph(:,:,x))
%    title('Sound pressure in the ukulele through the z axis');
%    
%    axis square
%    grid off
%    shading interp
%    caxis([0 0.25])
%    zlim([-2 2])
%    
%    pause(0.04);
   
   currentBinFreq  = x / (iterations*delta_t);
   
   if currentBinFreq < helmholtzCutoffFreq * 1.25
        freq_label(x) = currentBinFreq;
   end


x
end

%correct time domain x-axis to read real times
timescale = 0:delta_t:delta_t*(iterations-1);

figure
timeDomainTitle = [iterstring ' timesteps of the resonance of Uke after impulse excitation using Lower Res Model.'];
plot(timescale, measuredPressure);
title(timeDomainTitle);
xlabel('Time(secs)');
ylabel('Wave Pressure');
timeDomainPlotName = ['lores3d Time Domain ' iterstring 'x'];
savefig(timeDomainPlotName);

X = fft(measuredPressure);
X_mag = abs(X);
 

figure
fftTitle = ['FFT using ' iterstring ' timesteps of the resonance of Uke after impulse excitation using Lower Res Model.'];
stem(freq_label(1:length(freq_label)-1), X_mag(2:length(freq_label)));
title(fftTitle);
xlabel('Discrete Frequency components (Hz) up to 1.25x the TLM Cutoff Freq');
ylabel('Magnitude');
fftPlotName = ['lores3d FFT  ' iterstring 'x'];
savefig(fftPlotName);


%DFT at 232.1901hz
figure
subplot(131)
pcolor(abs(XPressureDFT1))
title('DFT of 232.1901Hz in the ukulele through the x axis');
axis square
grid off
shading interp
caxis([0 27])
zlim([-2 2])
   
subplot(132)
pcolor(abs(YPressureDFT1))
title('DFT of 232.1901Hz in the ukulele through the y axis.');
axis square
grid off
shading interp
caxis([0 27])
zlim([-2 2])
   
subplot(133)
pcolor(abs(ZPressureDFT1))
title('DFT of 232.1901Hz in the ukulele through the z axis');
   
axis square
grid off
shading interp
caxis([0 27])
zlim([-2 2])
savefig('lowres uke dft at 232Hz.fig');




%DFT at 751.0346hz
figure
subplot(131)
pcolor(abs(XPressureDFT2))
title('DFT of 751.0346Hz in the ukulele through the x axis');
axis square
grid off
shading interp
caxis([0 15])
zlim([-2 2])
   
subplot(132)
pcolor(abs(YPressureDFT2))
title('DFT of 751.0346Hz in the ukulele through the y axis.');
axis square
grid off
shading interp
caxis([0 15])
zlim([-2 2])
   
subplot(133)
pcolor(abs(ZPressureDFT2))
title('DFT of 751.0346Hz in the ukulele through the z axis');
   
axis square
grid off
shading interp
caxis([0 15])
zlim([-2 2])
savefig('lowres uke dft at 751Hz.fig');





%DFT at 1482.003hz
figure
subplot(131)
pcolor(abs(XPressureDFT3))
title('DFT of 1482.003Hz in the ukulele through the x axis');
axis square
grid off
shading interp
caxis([0 60])
zlim([-2 2])
   
subplot(132)
pcolor(abs(YPressureDFT3))
title('DFT of 1482.003Hz in the ukulele through the y axis.');
axis square
grid off
shading interp
caxis([0 60])
zlim([-2 2])
   
subplot(133)
pcolor(abs(ZPressureDFT3))
title('DFT of 1482.003Hz in the ukulele through the z axis');
   
axis square
grid off
shading interp
caxis([0 60])
zlim([-2 2])
savefig('lowres uke dft at 1482Hz.fig');





%DFT at 2138.4419hz
figure
subplot(131)
pcolor(abs(XPressureDFT4))
title('DFT of 2138.4419Hz in the ukulele through the x axis');
axis square
grid off
shading interp
caxis([0 20])
zlim([-2 2])
   
subplot(132)
pcolor(abs(YPressureDFT4))
title('DFT of 2138.4419Hz in the ukulele through the y axis.');
axis square
grid off
shading interp
caxis([0 20])
zlim([-2 2])
   
subplot(133)
pcolor(abs(ZPressureDFT4))
title('DFT of 2138.4419Hz in the ukulele through the z axis');
   
axis square
grid off
shading interp
caxis([0 20])
zlim([-2 2])
savefig('lowres uke dft at 2138Hz.fig');




%DFT at 2422.2298hz
figure
subplot(131)
pcolor(abs(XPressureDFT5))
title('DFT of 2422.2298Hz in the ukulele through the x axis');
axis square
grid off
shading interp
caxis([0 10])
zlim([-2 2])
   
subplot(132)
pcolor(abs(YPressureDFT5))
title('DFT of 2422.2298Hz in the ukulele through the y axis.');
axis square
grid off
shading interp
caxis([0 10])
zlim([-2 2])
   
subplot(133)
pcolor(abs(ZPressureDFT5))
title('DFT of 2422.2298Hz in the ukulele through the z axis');
   
axis square
grid off
shading interp
caxis([0 10])
zlim([-2 2])
savefig('lowres uke dft at 2422Hz.fig');




%DFT at 2809.2132hz
figure
subplot(131)
pcolor(abs(XPressureDFT6))
title('DFT of 2809.2132Hz in the ukulele through the x axis');
axis square
grid off
shading interp
caxis([0 40])
zlim([-2 2])
   
subplot(132)
pcolor(abs(YPressureDFT6))
title('DFT of 2809.2132Hz in the ukulele through the y axis.');
axis square
grid off
shading interp
caxis([0 40])
zlim([-2 2])
   
subplot(133)
pcolor(abs(ZPressureDFT6))
title('DFT of 2809.2132Hz in the ukulele through the z axis');
   
axis square
grid off
shading interp
caxis([0 40])
zlim([-2 2])
savefig('lowres uke dft at 2809Hz.fig');

save('lowres uke', 'XPressureDFT1', 'XPressureDFT2', 'XPressureDFT3', 'XPressureDFT4', 'XPressureDFT5', 'XPressureDFT6', 'YPressureDFT1', 'YPressureDFT2', 'YPressureDFT3', 'YPressureDFT4', 'YPressureDFT5', 'YPressureDFT6', 'ZPressureDFT1', 'ZPressureDFT2', 'ZPressureDFT3', 'ZPressureDFT4', 'ZPressureDFT5', 'ZPressureDFT6', 'Xpressuregraph', 'Ypressuregraph', 'Zpressuregraph');

% 2809.2132Hz - between F and F sharp
% 2422.2298Hz - between D and D sharp
% 2138.4419Hz - between C and C sharp
% 1482.003Hz - F sharp
% 751.0346Hz - between F sharp and G
% 232.1901Hz - A sharp

