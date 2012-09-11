#!/bin/bash
make html 
rm -r ../../mirnylab_documentation/mirnylab.bitbucket.org/hiclib
mv  build/html  ../../mirnylab_documentation/mirnylab.bitbucket.org/hiclib
rm -r build
cd ../../mirnylab_documentation/mirnylab.bitbucket.org
hg add * 
hg commit -m "Automatic commit from a new build" 
hg push 

