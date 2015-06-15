CellProfilerPelkmans -- "refactor_modules" branch
=================================================

The purpose of this branch is to refactor modules in such a way that GUI and input/output handling is separated from the actual image processing.

To this end, each module will be separated into two parts:
- A *main* function (the actual "module") which is concerned with reading input from the **handles** structure and writing output to the handles structure. It mustn't deal with image processing.
- One or several *sub* functions, which receive input and return output in a **handles independent** way. They do the actual image processing.

The subfunctions are
- bundled in a *library* (a location outside of the module files)
- organized as [packages](http://ch.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html)
- [imported](http://ch.mathworks.com/help/matlab/ref/import.html) explicitly by modules

This structure will allow testing of code outside the GUI framework and reuse by other tools, such as [Jterator](https://github.com/pelkmanslab/JtLib).


Modules to start with
---------------------

- IdentifyPrimIterative
- IdentifySecondaryIterative
