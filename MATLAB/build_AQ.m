function [ A, Q ] = build_AQ( T, proc_var )
%BUILD_AQ Construct transition matrixes A and Q

d = 0.1;

A = [1 0 T           0          ;
     0 1 0           T          ;
     0 0 exp(-d*T) 0          ;
     0 0 0           exp(-d*T)];
 
Q = proc_var * ...
    [T^3/3  0      T^2/2  0    ;
     0      T^3/3  0      T^2/2;
     T^2/2  0      T      0    ;
     0      T^2/2  0      T    ];

end

