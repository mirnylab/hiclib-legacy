#!/bin/bash
make html 
rm -r ../../manual_pages/mimakaev.bitbucket.org/iterative_correction_paper_f7f91c861335
mv  build/html  ../../manual_pages/mimakaev.bitbucket.org/iterative_correction_paper_f7f91c861335
rm -r build
cd ../../manual_pages/mimakaev.bitbucket.org
hg add * 
hg commit -m "Automatic commit from a new build" 
hg push 

