# Epetra Matrix format has been added to goma

Available solvers with epetra are Amesos and AztecOO

Manual changes:

# Solver specifications

## Matrix storage format

``Matrix storage format = {msr | vbr | epetra}``

epetra    trilinos native format for matrices

## Solution Algorithm

aztecoo   trilinos updated AztecOO solver (object oriented aztec),
          only available with epetra matrix format. Uses regular aztec
          preconditioner options, adds an additional `AztecOO Solver`
          parameter needed similar to amesos

    AztecOO Solver = {gmres | cgs | cg | tfqmr | bicgstab | y12m}

Example:

    Solution Algorithm = aztecoo
    AztecOO Solver = gmres