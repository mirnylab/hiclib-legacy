'''Add the parent directory ../ to one of the bash settings files.
'''
import os, sys

libpath = os.getcwd()
if libpath.startswith(os.path.expanduser('~')):
    libpath = os.path.join(*(['$HOME',] + libpath.split(os.sep)[3:-1]))

export_line = 'export PYTHONPATH="$PYTHONPATH:{0}"'.format(libpath)

profiles = [os.path.expanduser(i) 
            for i in ['~/.bash_profile', '~/.bashrc', '~/.profile']]

# Do nothing if the library is already exported.
for profile_path in profiles:
    if os.path.isfile(profile_path):
        for line in open(profile_path):
            if export_line in line:
                print 'The PYTHONPATH is already set in {0}'.format(
                       profile_path)
                sys.exit()
    
# If not, modify the first existing file in the chain of profiles.
for profile_path in profiles:
    if os.path.isfile(profile_path):
        profile_file = open(profile_path, 'a')
        profile_file.writelines(
            ['\n# Added by the mirnylab install script.\n',
             export_line, 
             '\n'])
        print 'PYTHONPATH is added to {0}'.format(profile_path)
        sys.exit()

# Create the first file in the chain if the chain is empty.
profile_path = profiles[0]
profile_file = open(profile_path, 'w')
profile_file.writelines(
    ['\# Added by the mirnylab install script.\n',
     export_line,
     '\n'])
print 'PYTHONPATH is added to {0}'.format(profile_path)
