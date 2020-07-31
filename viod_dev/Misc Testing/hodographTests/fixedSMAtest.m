% from the SMA error relations, we know that position error PERCENTAGE
% varies with respect to the inverse root of the semi-major axis, or the
% inverse of velocity.
%
% However, the ABSOLUTE position error actually varies with the inverse of
% velocity cubed, because the PERCENTAGE is calculated by dividing by the
% semi-major axis. So for example, when the velocity is doubled, the
% PERCENTAGE error becomes 1/2 its original magnitude. However the ABSOLUTE
% error is 1/8 its original magnitude. The discrepancy comes from the fact 
% that the semi-major axis became 1/4 its original length. 
%
% In the eccentricity variation case, we know that the semi-major axis does
% not change in length. Therefore we need to check whether the orbit
% position error actually varies with the inverse of velocity cubed as
% predicted, hence the creation of this test.

close all
clear;clc

