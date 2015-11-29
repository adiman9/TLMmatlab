# TLMmatlab

This code is a Matlab model that uses the Scalar Transmission Line Matrix (TLM) method to simulate the acoustic characteristics of a Ukulele.

Currently the Ukulele is hardcoded into the project, but in the future further development will be done to allow any arbitrary model to be used. 

The methoed used to get the 3D geometry into the model is to quantise the 3D space and create layers of bitmap images where each pixel represents one node in the quantised mesh of space. If the pixel is black, this means that node is part of a hard surface, if it is white then the node is freespace.

##Why

This was done as part of my Final Year Project at Loughborough University where I attempted to use this model to find the Helmholtz resonance of a Ukulele which was then compared against experimental and scholarly data to verify the accuracy of the model within Matlab.
