function R = axisangle2rot(omega,theta)
% AXISANGLE2ROT Calculates the 3x3 rotation matrix corresponding to a
% rotation theta about an axis defined by omega.
%
% This function effectively implents Rodrigues' formula.
%
% Inputs: omega - 3D vector representing the axis of rotation
%         theta - scalar representing the rotation angle
%
% Output: R - rotation matrix (SO(3))
%
% RBE 501 - Robot Dynamics - Spring 2021
% Worcester Polytechnic Institute
%
% Instructor: L. Fichera <lfichera@wpi.edu>
% Last modified: 01/29/2020

    omega_bracket = [0,-omega(3),omega(2); omega(3),0,-omega(1); -omega(2), omega(1),0];
    I = eye(3);
    R = I + (sin(theta)*omega_bracket) + ((1-cos(theta))*(omega_bracket*omega_bracket));
end