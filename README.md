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


## Storing Values in Lists


```python
odds = [1, 3, 5, 7]
print('odds are:', odds)
```

    odds are: [1, 3, 5, 7]



```python
print('first element:', odds [0])
print('last element:' , odds [3])
print('"-1" element:', odds [-1])
```

    first element: 1
    last element: 7
    "-1" element: 7



```python
names = ['Curie' ,'Darwing', 'Turing'] # Typo in Darwin's name

print('names is originally:' , names)

name[1] = 'Darwin' # Correct the name

print('final value of names:' , names)
```

    names is originally: ['Curie', 'Darwing', 'Turing']
    final value of names: ['Curie', 'Darwing', 'Turing']



```python
#name = 'Darwin'
#name[0] = 'd'
```


```python
odds.append(11)
print('odds after adding a value:' , odds)
```

    odds after adding a value: [1, 3, 5, 7, 11]



```python
removed_element = odds.pop(0)
print('odds after removing the first element:' , odds)
print('removed_element:' , removed_element)
```

    odds after removing the first element: [3, 5, 7, 11]
    removed_element: 1



```python
odds.reverse()
print('odds after reversing:' , odds)
```

    odds after reversing: [11, 7, 5, 3]



```python
odds = [3,5,7]
primes = odds
primes.append(2)
print('primes;' , primes)
print('odds:' , odds)
```

    primes; [3, 5, 7, 2]
    odds: [3, 5, 7, 2]



```python
odds = [3,5,7]
primes = list(odds)
primes.append(2)
print('primes;' , primes)
print('odds:' , odds)
```

    primes; [3, 5, 7, 2]
    odds: [3, 5, 7]



```python
binomial_name = "Drosophila melangogaster"
group = binomial_name[0:10]
print('group:' , group)

species = binomial_name[11:23]
print('species:' , species)

chromosomes = ['X' , 'Y' , '2' , '3', '4']
autosomes = chromosomes[2:5]
print('autosomes:' , autosomes)

last = chromosomes [-1]
print('last:', last)
```

    group: Drosophila
    species: melangogaste
    autosomes: ['2', '3', '4']
    last: 4



```python
date = 'Monday 4 January 2023'
day = date[0:6]
print('Using 0 to begin range:', day)
day = date[:6]
print('Omitting beginning index:', day)
```

    Using 0 to begin range: Monday
    Omitting beginning index: Monday



```python
months = ['jan', 'feb' , 'mar' , 'apr' ,'may' , 'jun', 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec']
sond = months

print('With known last postion:' , sond)

sond = months[8:len(months)]
print('Using len() to get last entry:' , sond)

sond = months[8:]
print('Omitting ending index:' , sond)


```

    With known last postion: ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
    Using len() to get last entry: ['sep', 'oct', 'nov', 'dec']
    Omitting ending index: ['sep', 'oct', 'nov', 'dec']



