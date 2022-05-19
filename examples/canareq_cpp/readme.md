# Equilibrium C++ Instructions

The ```main.cpp``` is a c++ script that will execute all functions included inside the following
code:

```c++
int main(int argc, char** argv) {
    aerodes *aero = new aerodes(argc, argv);    // create new aero object

    double angle=0, d_angle=0.5, angle0=1;
    // canar equillibrium
    aero->canary->readInputParams("canareq.data");
    aero->canary->readInputPolar("canarpolar.dat");

    std::string filename = "results.dat";
    for(int i = 0; i < 5; i++) {
        angle = d_angle*i + angle0; // do not set alpha equal to 0
        aero->canary->setCanardAngle(angle);
        aero->canary->linearModel();
        aero->canary->nonlinearModel();
        aero->canary->outputResults2Dat(filename);
    }

    delete aero;
    return 1;
}
```

**Note:** everything inside the parenthesis is within the scope of code that is executed.

**Note:** the for loop sweeps through various canard setting angles. Afterwards the non-linear solver performs a search for the stability criteria.

## Compiling The Example

In order to compile this code, so that we may run the resulting executable called ```test```, use this command.

On Mac OS X
```g++ -o test -laerolib main.cpp```

On Linux
```g++ -o test main.cpp -laerolib```

Now we may run the executable like so,

```./test```

