# Approximation of splines

## Compilation 2D-spline
```shell
$ make all
```
## Using
**Data generation**
```text
$ data_gen.exe <func type*> <begin val> <end val> <step count> <out filename>
               *Allow function types: 'linear', 'exp', 'sin', 'cos', 'spec'
```
**Create spline**
```text
$ spline.exe <data filename> <result filename> <points count> <spline knots count> <spline degree>
```

**Visualisation**
```text
$ python vis.py <data filename> <result filename> <*result knots filename>
                * <result knots filename> created auto by slpine.exe
```

## Run example 2D-spline
```shell
$ data_gen.exe spec 0 100 100 data
$ spline.exe data result 1000 11 4
$ python vis.py data result result_knots
```

![alt tag](https://github.com/alex-belov/splines/blob/master/resource/example.png)