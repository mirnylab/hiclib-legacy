'''Add the parent directory ../ to one of the bash settings files.
'''
import os
import sys

print
print "To skip checks use python install_linux.py --noCheck"
print

if (len(sys.argv) <= 1) or  (sys.argv[1].lower()[-7:] != "nocheck"):
    print "Checking python version..",
    ver = sys.version_info
    if (ver < (2, 6)) or (ver >= (3, 0)):
        raise RuntimeError("Please use python 2.6+, but not 3.x")
    print 'Correct!'

    print "Checking for mirnylib..",
    try:
        import mirnylib #@UnusedImport @IgnorePep8
    except:
        print "Mirnylib library not installed"
        print "see http://bitbucket.org/mirnylab/mirnylib"
        raise RuntimeError("Mirnylib library not installed")
    print "Found!"

    print "Checking for numpy..",
    try:
        import numpy
        print "Found!"
    except:
        print "Numpy not found"
        raise RuntimeError("Numpy not found")

    print "Checking for numpy version..",
    try:
        nv = numpy.__version__
        nums = tuple([int(i) for i in nv.split('.')[:2]])
        assert nums >= (1, 6)
        print "Correct!"
    except:
        print "numpy version is %s" % nv
        print "Needs at least numpy 1.6"
        print "See manual for numpy installation guide"
        raise RuntimeError("Wrong numpy version")

    print "Checking for mirnylib.h5dict install..",
    from mirnylib.h5dict import h5dict
    a = h5dict()
    b = numpy.empty(1000000, dtype="int16")
    c = "bla bla bla"
    a["numpy"] = b
    a["object"] = c
    assert (a["numpy"] - b).sum() == 0
    print "H5dict test successful!"

    print "Checking for joblib..",
    try:
        import joblib
        print "Found!"
    except:
        print "joblib not found"
        raise RuntimeError("joblib not found")

    print


os.chdir("src")
libpath = os.getcwd()
libpath = os.path.normpath(libpath.replace(os.path.expanduser('~'), '$HOME/'))

export_line = 'export PYTHONPATH="$PYTHONPATH:{0}"'.format(libpath)

profiles = [os.path.expanduser(i) for i in ['~/.bash_profile', '~/.bashrc']]

# Do nothing if the library is already exported.
for profile_path in profiles:
    setflag = 1
    if os.path.isfile(profile_path):
        for line in open(profile_path):
            if export_line in line:
                print 'The PYTHONPATH is already set in {0}'.format(
                    profile_path)
                setflag = 0
                break

    if setflag == 1:

        profile_file = open(profile_path, 'a')
        profile_file.writelines(
            ['\n# Added by the mirnylab install script.\n',
             export_line,
             '\n'])
        print 'PYTHONPATH is added to {0}'.format(profile_path)
