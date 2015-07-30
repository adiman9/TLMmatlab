clear;
gridSize = 49;

reflectedLeft = zeros(gridSize,gridSize);
reflectedRight = zeros(gridSize,gridSize);
reflectedFront = zeros(gridSize,gridSize);
reflectedBack = zeros(gridSize,gridSize);

incidentLeft = zeros(gridSize,gridSize);
incidentRight = zeros(gridSize,gridSize);
incidentFront = zeros(gridSize,gridSize);
incidentBack = zeros(gridSize,gridSize);

currentPressure = zeros(gridSize,gridSize);
measuredPressure = (1:500);

reflection = 0;

figure('units','normalized','outerposition',[0 0 1 1])
pause(1);


uke = imread('uke side 2 test.bmp');

    
for x = 1:500
    
    f=0.01;
    
%     level=sin(2*pi*f*x);
    if x < 10
        level=3;
    else
        level=0;
    end;
    
    incidentBack(1201) = incidentBack(1201) + level;
    incidentFront(1201) = incidentFront(1201)+ level;
    incidentLeft(1201) = incidentLeft(1201)+ level;
    incidentRight(1201) = incidentRight(1201) + level;
    
    
   % loop to propogate the wave and place the boundary of the ukulele
   for i = 1:gridSize
       
       for j = 1:gridSize
       
          if uke(j, i) == 0
              
              reflectedFront(j, i) = incidentFront(j, i);
              reflectedBack(j, i) = incidentBack(j, i);
              reflectedLeft(j, i) = incidentLeft(j, i);
              reflectedRight(j, i) = incidentRight(j, i);
              
          else
              
              reflectedFront(j, i) = 0.5 * (incidentRight(j, i) - incidentFront(j, i) + incidentLeft(j, i) + incidentBack(j, i));
              reflectedRight(j, i) = 0.5 * (-incidentRight(j, i) + incidentFront(j, i) + incidentLeft(j, i) + incidentBack(j, i));
              reflectedBack(j, i) = 0.5 * (incidentRight(j, i) + incidentFront(j, i) + incidentLeft(j, i) - incidentBack(j, i));
              reflectedLeft(j, i) = 0.5 * (incidentRight(j, i) + incidentFront(j, i) - incidentLeft(j, i) + incidentBack(j, i));
          
          end
          
       end
      
   end 
   
   
   
   % loop to take care of the outer boundary conditions of the mesh
   for i = 1:gridSize
       
       for j = 1:gridSize
       
          if i == 1 && j == 1
             incidentRight(j, i) = reflectedLeft(j, i+1);
             incidentFront(j, i) = reflectedBack(j+1, i);
             incidentLeft(j, i) = reflection * reflectedLeft(j, i);
             incidentBack(j, i) = reflection * reflectedBack(j, i);
          elseif i == gridSize && j == 1
              incidentRight(j, i) = reflection * reflectedRight(j, i);
              incidentFront(j, i) = reflectedBack(j+1, i);
              incidentLeft(j, i) = reflectedRight(j, i-1);
              incidentBack(j, i) = reflection * reflectedBack(j, i);
          elseif i == 1 && j == gridSize
              incidentRight(j, i) = reflectedLeft(j, i+1);
              incidentFront(j, i) = reflection * reflectedFront(j, i);
              incidentLeft(j, i) = reflection * reflectedLeft(j, i);
              incidentBack(j, i) = reflectedFront(j-1, i);
          elseif i == gridSize && j == gridSize
              incidentRight(j, i) = reflection * reflectedRight(j, i);
              incidentFront(j, i) = reflection * reflectedFront(j, i);
              incidentLeft(j, i) = reflectedRight(j, i-1);
              incidentBack(j, i) = reflectedFront(j-1, i);
          elseif i == 1
              incidentRight(j, i) = reflectedLeft(j, i+1);
              incidentFront(j, i) = reflectedBack(j+1, i);
              incidentLeft(j, i) = reflection * reflectedLeft(j, i);
              incidentBack(j, i) = reflectedFront(j-1, i);
          elseif i == gridSize
              incidentRight(j, i) = reflection * reflectedRight(j, i);
              incidentLeft(j, i) = reflectedRight(j, i-1);
              incidentFront(j, i) = reflectedBack(j+1, i);
              incidentBack(j, i) = reflectedFront(j-1, i);
          elseif j == 1
              incidentRight(j, i) = reflectedLeft(j, i+1);
              incidentFront(j, i) = reflectedBack(j+1, i);
              incidentLeft(j, i) = reflectedRight(j, i-1);
              incidentBack(j, i) = reflection * reflectedBack(j, i);
          elseif j == gridSize
              incidentRight(j, i) = reflectedRight(j, i+1);
              incidentFront(j, i) = reflection * reflectedFront(j, i);
              incidentLeft(j, i) = reflectedRight(j, i-1);
              incidentBack(j, i) = reflectedFront(j-1, i);
          else
              incidentRight(j, i) = reflectedLeft(j, i+1);
              incidentFront(j, i) = reflectedBack(j+1, i);
              incidentLeft(j, i) = reflectedRight(j, i-1);
              incidentBack(j, i) = reflectedFront(j-1, i);
          end

          currentPressure(j, i) = (reflectedRight(j, i) + reflectedLeft(j, i) + reflectedFront(j, i) + reflectedBack(j, i));
       end
   end  
   
   
   measuredPressure(x) = currentPressure(1301);
%    surf(currentPressure)
   pcolor(currentPressure)
   axis square
   grid off
   shading interp
   caxis([0 0.4])
   zlim([-2 2])
   pause(0.1);
end

figure
plot(measuredPressure);

figure
fs = 1;
X = fft(measuredPressure);
X_mag = abs(X);
f = 0 : fs/length(X_mag) : fs - fs/length(X_mag);
plot(f(1:length(f)), X_mag(1:length(measuredPressure)));


