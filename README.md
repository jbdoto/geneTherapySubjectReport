# geneTherapySubjectReport

## Background
This is a Dockerized version of the Gene Therapy Subject Report that can be run as part of the IntSiteCaller pipeline.

## Build Info

The easiest way to build and push this container is to run `./build.sh` you may need to adjust the profile to match
your personal settings.

## Container commands:

Demo reports:

    Rscript geneTherapySubjectReport/report.R --patient "pFR03" --trial "WAS" --outputDir /home/everett/projects/gtVISA_dashboard/gtVISA/AdaptimmuneCART --reportFile report.Rmd
    Rscript geneTherapySubjectReport/report.R --patient "pFR01" --trial "WAS" --outputDir /home/everett/projects/gtVISA_dashboard/gtVISA/AdaptimmuneCART --reportFile report.Rmd
    Rscript geneTherapySubjectReport/report.R --patient "pFR05" --trial "WAS" --outputDir /home/everett/projects/gtVISA_dashboard/gtVISA/AdaptimmuneCART --reportFile report.Rmd

    Rscript geneTherapySubjectReport/report.R --patient "p04511-202" --trial "AdaptimmuneCART" --outputDir /home/everett/projects/gtVISA_dashboard/gtVISA/AdaptimmuneCART --reportFile report.Rmd 


 