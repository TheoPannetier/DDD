d:
cd "d:\data\ms\DDD\"
del DDD_3.3.tar.gz
del DDD_3.3.zip
R CMD build DDD
R CMD INSTALL --build DDD_3.3.tar.gz
pause
R CMD check --timings --as-cran DDD_3.3.tar.gz
pause