% currently we think it's producing a hyperbolic orbit with negative a and
% therefore imaginary E

% plot velocity vectors to check if collinear using rng(1) on day 22
%   > the error was not due to collinearity. velocity vectors were already
%     imaginary prior to input to IOD3V.m
%   > error still comes from TimeProp.m, therefore it's an error with
%     the hyperbolic orbit
%   > solved by implementing universal solver TimeProp_Universal.m
%     10/11/2019

% check days 12, 15, 16 for eccentricity and arg. of pe. error
% check results of changing seed for rng
%   > it appears that all anomalies still occur at the same days
%   > removing absolute value confirms the mean error trends to zero
%   > FALSE - it was because of a bug in the code causing rng to be
%   re-seeded every time

% check semi-major axis error spike at 55 days
% check Neptune 2000
%   especially true anomaly, arg. of pe, ecc
%   > the spike only occurred because of the small sample size. 10/23/2019


% The mean error of orbital elements for small measurement intervals is
% non-zero, and are biased in the same direction across multiple trials
% (this bias direction is different for different orbital elements)
%
%   - suspect this is caused by the orbit setup (either due to being
%     retrograde, or due to initial true anomaly, or due to other causes)
%   > Not actually - when using the Earth-Mars case the mean approaches
%   zero quickly. The real reason is due to the elliptic vs. hyperbolic
%   sorting, which means most elliptic orbits have eccentricities lower
%   than the true value.
%
%   ALL FOLLOWING ANALYSES ARE FOR ELLIPTIC CASES
%   > Interestingly, the mean semi-major axis error is greater than zero
%   for elliptic cases while the eccentricity error was less than zero.
%       > This is NOT due to the algorithm generating larger and less
%       eccentric orbits
%       > INSTEAD, it is beccause for most cases a smaller and less
%       eccentric orbit is generated, but the SMA change is proportionally
%       smaller than the ECC change when the orbit is smaller, and vice
%       versa. 
%       > This causes the 'mean' SMA error to be greater than zero since a
%       single case of a larger orbit skews the entire dataset; it also
%       causes the 'mean' ECC error to be less than zero because smaller
%       orbit cases have a more significant ECC change. 
%       > As the measurement interval increases, standard deviation
%       decreases until hyperbolic cases are no longer a significant
%       portion of the data. At this point the skewing effect is canceled
%       by the frequency of larger vs. smaller measurements, causing the
%       mean error to go to zero.
%   > ECC Vector rotates against orbital direction for smaller orbits
%       - WHY?
%   > Arg. of Pe is smaller than the mean because most orbit estimates are
%   smaller, which means the s.c. must be closer to p.e. to match velocity
%   > True Anomaly error largely mirrors that of Arg. of Pe, which appears
%   to be due to the consistency in inclination.
%       > since inclination is minimally affected, the angle between the
%       s.c. position vector and the line of nodes (i.e. w+f) must also
%       deviate minimally from the true value. 
%   > no clear interpretation of longitude of ascending node yet.
%
%   - plot each of the 100 cases across all time intervals to get 100
%   strands. the goal is to see trends for each rng sample.
%   - check zero noise gives true orbit
%   - fix orbit direction error in IOD3V_V2
%
%   - find cause for spikes
%   - find max eccentricity / farthest planet where no hyperbolic cases
%   start to happen