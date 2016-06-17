# Notes on Brk Fix integration in goma

Additions to goma with brk fix parallel integration

## Command Line Arguments

`-brk fn` | Specify a Brk file, this overrides the setting in Problem Description File

## Problem Description File

### File Specifications

    -- After SOLN file in input file

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

### Time Integration Specifications

    -- After Printing Frequency in input file

#### Fix Frequency

```
Fix Frequency = <integer>
```

##### Description/Usage

This optional card specifies the frequency at which goma should fix, or combine, the parallel pieces of the Output Exodus II file. This frequency is relative to the `Printing Frequency` so if Fix Frequency is set to 10 it will fix the output exodus file on the 10th print.

Without this card goma will only fix after the problem is solved.

Fixing only occurs if goma brk the exodus files (`-brk` or `Brk file`)

`<integer>` | Number of prints needed to fix (0 is disabled)

#### Examples

Following are sample cards:

```
    # Fix every print
    Fix Frequency = 1
```

```
    # Fix every 5 prints
    Fix Frequency = 5
```

### Continuation Specifications

    -- After Continuation Printing Frequency in input file

#### Continuation Fix Frequency

```
Continuation Fix Frequency = <integer>
```

##### Description/Usage

See Fix Frequency
