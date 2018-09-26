data_gen.exe 2 spec 0 100 100 data
univariate_spline.exe data result 1000 11 4
python ..\vis.py 2 data result result_knots