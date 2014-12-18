# Notes on using stratimkos

## Goma input deck

The solution algorithm needs to be set to stratimikos

    Solution Algorithm = stratimikos

The matrix storage format needs to be set to epetra

    Matrix storage format = epetra

One new goma input deck card is needed.

    Stratimikos File = stratimikos.xml

Included in this directory is a working stratimkos xml file, stratimkos.xml. These files are rather difficult to setup due to sparsity of documentation. The most helpful reference seems to be the [aztec user guide] (https://trilinos.org/oldsite/packages/aztecoo/AztecOOUserGuide.pdf).

## Switch between AztecOO and Belos Solvers

The value of the "Linear Solver Type" string specifies the use of either Belos or AztecOO. Set the value to the exact name of your Parameterlist name, i.e. stratimkos.xml lines 5 or 15.

## More Documentation

[Belos] (https://trilinos.org/docs/r12.6/packages/belos/doc/html/index.html)
[IFPack] (https://trilinos.org/oldsite/packages/ifpack/IfpackUserGuide.pdf)
[Stratimikos] (https://trilinos.org/docs/r12.2/packages/stratimikos/doc/html/index.html)