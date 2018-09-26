# Approximation of splines
 ## Compilation Univariate spline
```shell
$ make all
```
## Using
**Data generation**
```text
$ data_gen.exe <dimension> <func type> <begin val x> <end val x>
                [<begin val y> <end val y>]** <step count> <out filename>
                [for dimension = 3]"
```
**Create Univariate spline**
```text
$ univariate_spline.exe <data filename> <result filename> <points count> <spline knots count> <spline degree>
```
 **Visualisation**
```text
$ python vis.py <dimension> <data filename> <result filename> <result knots filename>
                * <result knots filename> created auto by slpine.exe
```
 ## Run example 2D-spline
```shell
$ data_gen.exe 2 spec 0 100 100 data
$ univariate_spline.exe data result 1000 11 4
$ python vis.py 2 data result result_knots
```
 ![alt tag](https://github.com/alex-belov/splines/blob/master/resource/example_1.png)
