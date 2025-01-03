# README for Orbit Integrator
## Creator: Raymann Singh
### Created: 5/7/2023
### Last Edited: 1/3/2025

## About
This code integrates the position and velocity vectors over a set period using Euler Integration. It also finds the classical orbital elements of a given initial position and velocity vector and then uses them to final the final classical orbital elements using Kepler's Equation. These results are then plotted and compared by using each method's specific mechanical energy and specific angular momentum. It is important to note that these are the two elements that are being compared since Kepler's Equation only gives final position and velocity vectors which through Euler Integration will only be reached as dt approaches 0.

The orbit is inetegrated using the Perifocal Coordinate System. The integrator is ran over a default of 5 different timesteps ranging from 0.1 to 300 seconds.

![Perifocal Orbit](https://github.com/RaymannS/Orbit-Integrator/blob/main/orbit_perifocal.png)

The radius and velocity can then be seen over time. The default orbit has an eccentricity of near 0 meaning as the integration steps approach zero the radius and velocity will remain near constant over time, as seen below. 

![Radius over Time](https://github.com/RaymannS/Orbit-Integrator/blob/main/radii_time.png)

![Velocity over Time](https://github.com/RaymannS/Orbit-Integrator/blob/main/vel_time.png)

## How to Use
Download the [orbit_integrator.m](https://github.com/RaymannS/Orbit-Integrator/blob/main/orbit_integrator.m) file and open it in MATLAB.
Change the variables in the _Init Variables_ section at the top, these are the initial position and velocty vectors alongside the integration step sizes and period to integrate over.
Run the code and wait for the integration to be done (step sizes approaching zero are more precise, but increase computation time).

## License
R. Singh, hereby disclaims all copyright interest in the program “Orbit Integrator” written by R. Singh.

signature of R. Singh 3 January 2025
R. Singh, Owner

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
