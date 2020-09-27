![](img/cpic.png)

Particle in Cell Method, written in C++

# Installation

You will need a C++ compiler, for example [gcc](https://gcc.gnu.org/), as well as the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library installed on your system. Make sure that the directory `/usr/include/Eigen` exists on your system, if not you can create a symlink to where Eigen was installed, usually it is `/usr/include/eigen3` 
```
ln -s /usr/include/eigen3 /usr/include/Eigen
```
Now you are all set, just run `make` in the root directory and the code should compile. With `make run` you can run a test case from the `tests` directory. Check out the `Makefile` to see how to use `libcpic` in your own code.

# Examples

## Free Electrons Moving Around a Cloud of Oxygen Ions (Collisionless)

![](img/vlasov.gif)

## Two Oxygen Velocity Distributions Merging (DSMC Collisions)

![](img/dsmc.gif)
