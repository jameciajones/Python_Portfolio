# Python_Portfolio
This is the portfolio for python code.

## Analyzing Patient Data
In this analysis we looked at inflammation data for multiple patients.

```python
import numpy
```


```python
numpy.loadtxt(fname = 'inflammation-01.csv', delimiter = ",")
```




    array([[0., 0., 1., ..., 3., 0., 0.],
           [0., 1., 2., ..., 1., 0., 1.],
           [0., 1., 1., ..., 2., 1., 1.],
           ...,
           [0., 1., 1., ..., 1., 1., 1.],
           [0., 0., 0., ..., 0., 2., 0.],
           [0., 0., 1., ..., 1., 1., 0.]])




```python
data = numpy.loadtxt(fname = 'inflammation-01.csv', delimiter = ",")
```


```python
print(data)
```

    [[0. 0. 1. ... 3. 0. 0.]
     [0. 1. 2. ... 1. 0. 1.]
     [0. 1. 1. ... 2. 1. 1.]
     ...
     [0. 1. 1. ... 1. 1. 1.]
     [0. 0. 0. ... 0. 2. 0.]
     [0. 0. 1. ... 1. 1. 0.]]



```python
print(data.shape)
```

    (60, 40)



```python
print('first value in data:' , data [0,0])
```

    first value in data: 0.0



```python
print('middle value in data:' , data[29,19])
```

    middle value in data: 16.0



```python
print(data[0:4, 0:10])
```

    [[0. 0. 1. 3. 1. 2. 4. 7. 8. 3.]
     [0. 1. 2. 1. 2. 1. 3. 2. 2. 6.]
     [0. 1. 1. 3. 3. 2. 6. 2. 5. 9.]
     [0. 0. 2. 0. 4. 2. 2. 1. 6. 7.]]



```python
print(data[0:5, 0:10])
```

    [[0. 0. 1. 3. 1. 2. 4. 7. 8. 3.]
     [0. 1. 2. 1. 2. 1. 3. 2. 2. 6.]
     [0. 1. 1. 3. 3. 2. 6. 2. 5. 9.]
     [0. 0. 2. 0. 4. 2. 2. 1. 6. 7.]
     [0. 1. 1. 3. 3. 1. 3. 5. 2. 4.]]



```python
small = data[:3, 36:]
```


```python
print('small is:')
```

    small is:



```python
print(small)
```

    [[2. 3. 0. 0.]
     [1. 1. 0. 1.]
     [2. 2. 1. 1.]]



```python
# Let us a numpy function
print(numpy.mean(data))
```

    6.14875



```python
maxval, minval, stdvat = numpy.amax(data), numpy.amin(data), numpy.std(data)


```


```python
print(maxval)
print(minval)
print(stdval)
```

    20.0
    0.0



    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-30-566e6d04f418> in <module>
          1 print(maxval)
          2 print(minval)
    ----> 3 print(stdval)
    

    NameError: name 'stdval' is not defined



```python
maxval = numpy.amax(data)
minval = numpy.amin(data)
stdval = numpy.std(data)
```


```python
print(maxval)
print(minval)
print(stdval)
```

    20.0
    0.0
    4.613833197118566



```python
print('maximum inflammation:' , maxval)
print('minimum inflammation:' , minval)
print('standard deviation:' , stdval)
```


```python
# Sometimes we want to look at variation in statistical values, such as maximum inflammation per patient, 
# or average from day one

patient_0 = data[0, :] # 0 on the first axis (rows), everything on the second (columns)

print('maximum inflammation for patient 0:', numpy.amax(patient_0))
```


```python
print('maximum inflammation for patient:2' , numpy.amax(data[2, :]))
```


```python
print(numpy.mean(data, axis =0))
```


```python
print(numpy.mean(data, axis = 0).shape)
```


```python
print(numpy.mean(data, axis = 1))
```


