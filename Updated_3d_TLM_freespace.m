clear all;
gridSize = 49;

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
graphcurrentPressure = zeros(gridSize,gridSize);

measuredPressure = (1:200);

reflection = 0;

figure('units','normalized','outerposition',[0 0 1 1])


for x = 1:200
    
    f=0.05;
    
    level=sin(2*pi*f*x);
%     if x < 10
%         level=3;
%     else
%         level=0;
%     end;
    
    incidentBack(gridSize / 2 + 0.5, gridSize / 2 + 0.5, gridSize / 2 + 0.5) = incidentBack(gridSize / 2 + 0.5, gridSize / 2 + 0.5, gridSize / 2 + 0.5) + level;
    incidentFront(gridSize / 2 + 0.5, gridSize / 2 + 0.5, gridSize / 2 + 0.5) = incidentFront(gridSize / 2 + 0.5, gridSize / 2 + 0.5, gridSize / 2 + 0.5)+ level;
    incidentLeft(gridSize / 2 + 0.5, gridSize / 2 + 0.5, gridSize / 2 + 0.5) = incidentLeft(gridSize / 2 + 0.5, gridSize / 2 + 0.5, gridSize / 2 + 0.5)+ level;
    incidentRight(gridSize / 2 + 0.5, gridSize / 2 + 0.5, gridSize / 2 + 0.5) = incidentRight(gridSize / 2 + 0.5, gridSize / 2 + 0.5, gridSize / 2 + 0.5) + level;
    incidentTop(gridSize / 2 + 0.5, gridSize / 2 + 0.5, gridSize / 2 + 0.5) = incidentTop(gridSize / 2 + 0.5, gridSize / 2 + 0.5, gridSize / 2 + 0.5) + level;
    incidentBottom(gridSize / 2 + 0.5, gridSize / 2 + 0.5, gridSize / 2 + 0.5) = incidentBottom(gridSize / 2 + 0.5, gridSize / 2 + 0.5, gridSize / 2 + 0.5) + level;
    
    
    for i = 1:gridSize
       
       for j = 1:gridSize
              
           for z = 1:gridSize
               
               reflectedFront(z, j, i) = (1/3) * (incidentRight(z, j, i) - 2*incidentFront(z, j, i) + incidentLeft(z, j, i) + incidentBack(z, j, i) + incidentTop(z, j, i) + incidentBottom(z, j, i));
               reflectedBack(z, j, i) = (1/3) * (incidentRight(z, j, i) + incidentFront(z, j, i) + incidentLeft(z, j, i) - 2*incidentBack(z, j, i) + incidentTop(z, j, i) + incidentBottom(z, j, i));
               reflectedLeft(z, j, i) = (1/3) * (incidentRight(z, j, i) + incidentFront(z, j, i) - 2*incidentLeft(z, j, i) + incidentBack(z, j, i) + incidentTop(z, j, i) + incidentBottom(z, j, i));
               reflectedRight(z, j, i) = (1/3) * (-2*incidentRight(z, j, i) + incidentFront(z, j, i) + incidentLeft(z, j, i) + incidentBack(z, j, i) + incidentTop(z, j, i) + incidentBottom(z, j, i));
               reflectedTop(z, j, i) = (1/3) * (incidentRight(z, j, i) + incidentFront(z, j, i) + incidentLeft(z, j, i) + incidentBack(z, j, i) -2* incidentTop(z, j, i) + incidentBottom(z, j, i));
               reflectedBottom(z, j, i) = (1/3) * (incidentRight(z, j, i) + incidentFront(z, j, i) + incidentLeft(z, j, i) + incidentBack(z, j, i) + incidentTop(z, j, i) - 2*incidentBottom(z, j, i));
               
           end
           
       end
      
    end 
   
    
    for i = 1:gridSize
       
       for j = 1:gridSize
              
           for z = 1:gridSize
               
               if i == 1 || j == 1 ||z == 1 || i == gridSize || j == gridSize || z == gridSize 
                   
                   
                   
               else
                  %freespace propagation
                   incidentFront(z, j, i) = reflectedBack(z, j+1, i);
                   incidentBack(z, j, i) = reflectedFront(z, j-1, i);
                   incidentLeft(z, j, i) = reflectedRight(z, j, i-1);
                   incidentRight(z, j, i) = reflectedLeft(z, j, i+1);
                   incidentTop(z, j, i) = reflectedBottom(z-1, j, i);
                   incidentBottom(z, j, i) = reflectedTop(z+1, j, i);
               end
               
               currentPressure(z, j, i) = reflectedRight(z, j, i) + reflectedLeft(z, j, i) + reflectedFront(z, j, i) + reflectedBack(z, j, i) + reflectedTop(z, j, i) + reflectedBottom(z, j, i);
       
               if z == 25
                   
                   graphcurrentPressure(j, i) = currentPressure(z, j, i);
               end    
               
           end    
       end  
    end 
    
   measuredPressure(x) = currentPressure(5, 5, 25);
   
   pcolor(graphcurrentPressure)
   axis square
   grid off
   shading interp
   caxis([0 0.15])
   zlim([-2 2])
   pause(0.0001);
   
end


fs = 1;
X = fft(measuredPressure);
X_mag = abs(X);
f = 0 : fs/length(X_mag) : fs - fs/length(X_mag);

figure
plot(f(1:length(f)), X_mag(1:length(measuredPressure)));