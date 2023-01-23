# shallow-water-sim
A shallow water simulator in a rectangular tub that uses the Navier-Stokes Equations to create an animation of how the surface of the water will 
move over time when given the initial state of the surface of the water.
## Outline
* [1 - Explanation of Inputs](1---Explanation-of-Inputs)
* [2 - Animations](2---Animations)
  * [(i) Water droplet as an input]((i)-Water-droplet-as-an-input)
  * [(ii) Water wall as an input]((ii)-Water-wall-as-an-input)
* [Improvements to Make](Improvements-to-Make)

## 1 - Explanation of Inputs
The calculations for the simulation as well as the function that does the calculations is contained in the 'shallow_water_simulation.py' file. 
The function that performs the calculation is `shallow_water_simulation`. It takes in the parameters `(Lx, Ly, T, eps, dx, dy, dt, init_cond, H)`, where
* **Lx**: The width of the rectangular tub in meters.
* **Ly**: The length of the rectangular tub in meters.
* **T**: The total duration of the simulation in seconds (i.e. the simulation will show how the surface of the water behaves from t=0 to t=T).
* **eps**: The smoothing parameter. This determines how smoothed out the waves are after every calculation. This makes it so that the waves in 
the animation look smooth rather than jagged. Note that eps is short for epsilon.
* **dx**: Step size for the horizontal spatial dimension. Splits the grid into N = Lx/dx number of cells horizontally. Program calculates surface height at every grid point.
* **dy**: Step size for the vertical spatial dimension. Splits the grid into M = Ly/dy number of cells vertically.
* **dt**: Step size for the time dimension. Program calculates the state of the surface of the water after every time step from t=0 to t=T.
* **init_cond**: A numpy array containing the initial state of the surface of the water. It is an array containing 3 arrays in the form `init_cond = [u_0, v_0, eta_init]`.
  Each of the following arrays has M rows and N columns so that there are N elements horizontally and M elements vertically:
    * **u_0**: Contains the horizontal component of the velocity of the water at each grid point.
    * **v_0**: Contains the vertical component of the velocity of the water at each grid point.
    * **eta_init**: Contains the height of the water at each grid point.
* **H**: The bathymetry function. I.e., this is a function that returns the height (and thus the shape) of the surface of the bottom of the tub.

To input the initial state and create the animation, use the `Animation code and results.py` file.

## 2 - Animations
### (i) Water droplet as an input

The output of the simulation with a still water droplet (i.e., Gaussian peak) at some position in the tub initially:

https://user-images.githubusercontent.com/104711470/213978371-f234efce-2b73-4f2d-b8d4-51f609c59a99.mp4

The same output with a different surface colouring:

https://user-images.githubusercontent.com/104711470/213978173-0954c208-721f-431e-9dc8-689bf9f677be.mp4

### (ii) Water wall as an input

The output of the simulation with a still water "wall" at some position in the tub initially:

https://user-images.githubusercontent.com/104711470/213978383-4bde5cf5-c272-4408-8afa-f72389cc88cd.mp4

The same output with a different surface colouring:

https://user-images.githubusercontent.com/104711470/213978420-0f12a86d-9397-481d-b516-5f7c81d202fc.mp4

## Improvements to Make
* Work on animations (fill in under the plot, reduce height of initial water droplet/wall, and make waves easier to see).
* Combine all plot and animation code into one "create_animation" function.
* Create a user interface where changing the input of the simulation is a lot easier/more accesible to someone not familiar with programming OR 
offer a range of inputs to click on.
