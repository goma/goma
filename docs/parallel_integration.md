# Notes on Brk Fix integration in goma

Additions to goma with brk fix parallel integration

## Command Line Arguments

`-brk fn` | Specify a Brk file, this overrides the setting in Problem Description File

## Problem Description File

### General Specifications

#### Brk file

```
Brk file = <file_name>
```

##### Description/Usage

This optional card specifies the name of the Brk file for this problem, if one does not exist goma will attempt to create one. The Brk file is used by the brk utility to break the Exodus II files on parallel runs for each processor.

Brk files can only be created on single processor runs.

`<file_name>` | A Brk file in the brk file syntax with specifications for material blocks

##### Examlpes

Following is a sample card:

```
    Brk file = in.brk
```

#### Fix Frequency

```
Fix Frequency = <integer>
```

##### Description/Usage

This optional card specifies the frequency at which goma should fix, or combine, the parallel pieces of the Output Exodus II file.

Without this card goma will only fix after the problem is solved.

`<integer>` | The timestep interval at which to fix the Output Exodus II file, must be positive.

#### Examples

Following are sample cards:

```
    # Fix every timestep
    Fix Frequency = 1
```

```
    # Fix every 5 timesteps
    Fix Frequency = 5
```

## Technical notes

Currently uses brkfix as a library with a hacky map_names.h which changes function names to not conflict with goma function names.

This was how brkfix was integrated previously in rd_mesh and then brk would be called with a populated *argv and argc to mimic the main executable call.

I've modified the brk.c and fix.c to fit a little better into the current goma structure (still pretty ugly) with the `fix_exo_file.c` and `brk_exo_file.c`
