#!/bin/bash
sphinx-build . _build
cp -r  _build/*  ../../manual_pages/mimakaev.bitbucket.org/iterative_correction_paper_f7f91c861335
cd ../../manual_pages/mimakaev.bitbucket.org
hg add * 
hg commit -m "Automatic commit from a new build" 
hg push 

