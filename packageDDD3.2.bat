d:
cd "d:\data\ms"
del DDD_3.2.tar.gz
del DDD_3.2.zip
R CMD build DDD
R CMD INSTALL --build DDD_3.2.tar.gz
pause
R CMD check --timings --as-cran DDD_3.2.tar.gz
pause