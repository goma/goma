
def multiname(filename, numproc):
    ''' Returns a list of filenames
    matching the number of processors,
    this follows the goma format of multiname'''
    return [filename + ".%d.%d" % (numproc, i) for i in range(numproc)]
