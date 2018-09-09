data_gen.exe spec 0 100 100 data
spline.exe data result 1000 11 4
python vis.py data result result_knots