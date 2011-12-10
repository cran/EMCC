
.onLoad <-
    function (libname, pkgname)
{
    this.year <- substr(as.character(Sys.Date( )), 1, 4)
    packageStartupMessage('##\n',
                          '## Evolutionary Monte Carlo Clustering Package (EMCC)\n',
                          '##\n',
                          '## Functionality: evolutionary Monte Carlo clustering, temperature\n',
                          '## ladder construction and placement\n',
                          '##\n',
                          '## Use: "help(package = EMCC)" at the R prompt for more info\n',
                          '##\n',
                          '## Copyright (C) 2006-', this.year, ' Gopi Goswami\n',
                          '##\n',
                          '##    Created by: Gopi Goswami <goswami@stat.harvard.edu>\n',
                          '## Maintained by: Gopi Goswami <grgoswami@gmail.com>\n',
                          '##\n')

    library.dynam(pkgname, pkgname, lib.loc=libname)
}






